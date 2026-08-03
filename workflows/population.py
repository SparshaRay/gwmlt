"""
Create Master Parameter Database
"""

import resource
resource.setrlimit(resource.RLIMIT_NOFILE, (65536, 65536))

import sys
import yaml
import shutil
import atexit
import signal
import sqlite3
import hashlib
import datetime
import argparse
import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm

import dask
from sqlmodel import create_engine, text

from gwmlt.config import config
from gwmlt.population import population_sampler
from gwmlt.utils.orchestrator import (
    HTCondorOrchestrator, SpawnerLogger, watchdog
)


# Parse input ------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Generate Compact Binary Population Database"
)

parser.add_argument(
    "-s", "--samples", 
    type=int, required=True, 
    help="Total number of samples to generate"
)
parser.add_argument(
    "-r", "--seed", 
    type=int, required=False, 
    default=0, help="Master RNG seed"
)
parser.add_argument(
    "-e", "--ecc-range", 
    type=float, required=False, nargs=2, metavar=('MIN', 'MAX'),
    default=(0.05, 0.50), help="Minimum and maximum eccentricities"
)

args = parser.parse_args()


# Setup variables & paths ------------------------------------------------

# Number of samples

TOTAL_NUM_SAMPLES  = args.samples
SAMPLES_PER_WORKER = config.orchestrator.populations_samples_per_worker

if TOTAL_NUM_SAMPLES % SAMPLES_PER_WORKER != 0 :
    raise ValueError(
        f"Total samples ({TOTAL_NUM_SAMPLES}) must be divisible by "
        f"samples per worker ({SAMPLES_PER_WORKER})."
    )

# Seed and eccentricities

BASE_SEED = args.seed
ECC_RANGE = args.ecc_range

# Setup paths

args_md5 = hashlib.md5(
    f"{ECC_RANGE}{BASE_SEED}{TOTAL_NUM_SAMPLES}".encode('utf-8')
).hexdigest()

run_id   = f"{config.orchestrator.populations_dir_name}_{args_md5}"

BASE_DIR = config.database_root / f"{config.orchestrator.populations_dir_name}s"
RUN_DIR  = BASE_DIR / run_id
LOG_DIR  = RUN_DIR / "logs"
TMP_DIR  = RUN_DIR / "tmp"
DB_DIR   = RUN_DIR / "dbs"

for d in [LOG_DIR, TMP_DIR, DB_DIR] : d.mkdir(parents=True, exist_ok=False) 

# Dump run config to yaml

with open(RUN_DIR / "run_config.yaml", "x") as f :
    yaml.dump({
        "ID"            : run_id,
        "total_samples" : TOTAL_NUM_SAMPLES,
        "ecc_range"     : list(ECC_RANGE),
        "base_seed"     : BASE_SEED,
    }, f, default_flow_style=False)

# Runtime args

TABLE_NAME = config.orchestrator.populations_table_name
DB_PATH    = DB_DIR / f"{TABLE_NAME}.sqlite3"

MAX_CONCURRENT_JOBS = config.orchestrator.max_concurrent_jobs

# Logger 

sys.stdout = SpawnerLogger(LOG_DIR, sys.stdout)
sys.stderr = SpawnerLogger(LOG_DIR, sys.stderr)

print(f"init @ {datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")}\n")


# Initialize orchestrator & exit directives ------------------------------

# Dask environment

dask_env = HTCondorOrchestrator(
    log_dir  = LOG_DIR, 
    max_jobs = MAX_CONCURRENT_JOBS,
    cores    = 1,
    memory   = "2GB",
    disk     = "128MB"
)

# Exit handlers

atexit.register(lambda: print(
    f"Shutdown complete @ {datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")}"
))
atexit.register(dask_env.close)

def handler(*_) :
    print("\nCleaning up and exiting...")
    sys.exit(1)

signal.signal(signal.SIGINT,  handler)
signal.signal(signal.SIGTERM, handler)


# Dask tasks -------------------------------------------------------------

print("Starting population database generation. Check Dask dashboard for progress.\n")

# Generator

@dask.delayed
@watchdog(timeout_seconds=3600)
def gen_df_and_write(
    num_samples: int,
    write_dir: Path,
    chunk_idx: int, 
    seed: int
) -> Path :

    final_file_path = write_dir / f"frag_{chunk_idx:06d}.h5"
    if final_file_path.exists() : return final_file_path
    
    df = population_sampler(
        num_samples = num_samples, 
        seed        = seed, 
        ecc_range   = ECC_RANGE,
        use_q_pdf_interpolator = True
    )
            
    start_uid = chunk_idx * num_samples
    df['UID'] = range(start_uid, start_uid + num_samples)
    df = df.reindex(sorted(df.columns), axis=1)

    tmp_file_path = final_file_path.with_suffix(".tmp")
    df.to_hdf(tmp_file_path, key='pop', mode='w', format='fixed')
    tmp_file_path.rename(final_file_path)
    
    return final_file_path

# Create delayed tasks for Dask workers

delayed_tasks = [
    gen_df_and_write(SAMPLES_PER_WORKER, TMP_DIR, i, i+BASE_SEED)
    for i in range(TOTAL_NUM_SAMPLES // SAMPLES_PER_WORKER)
]


# Create db schema -------------------------------------------------------

meta_df = population_sampler(num_samples=1, seed=0, use_q_pdf_interpolator=True)
meta_df['UID'] = 0
meta_df = meta_df.reindex(sorted(meta_df.columns), axis=1)

sqlite_engine = create_engine(f"sqlite:///{DB_PATH}", echo=False)

create_table_stmt = pd.io.sql.get_schema(meta_df.iloc[0:0], TABLE_NAME, con=sqlite_engine)
create_table_stmt = create_table_stmt.replace("\n)", ",\n\tPRIMARY KEY (`UID`)\n)")

with sqlite_engine.begin() as conn :
    conn.execute(text(f"DROP TABLE IF EXISTS {TABLE_NAME}"))
    conn.execute(text(create_table_stmt))


# Execution --------------------------------------------------------------

print("Spawning workers to create hdf fragments...")
generated_files = dask.compute(*delayed_tasks, retries=4)

print(
    f"Fragment generation complete. Starting to assemble to database...\n"
    f"{datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")}\n"
)

with sqlite3.connect(DB_PATH) as conn :
    
    conn.execute("PRAGMA journal_mode = OFF;")
    conn.execute("PRAGMA synchronous = 0;")
    conn.execute("PRAGMA cache_size = -1048576;")
    conn.execute("PRAGMA locking_mode = EXCLUSIVE;")
    conn.execute("PRAGMA temp_store = MEMORY;")
    
    for file_path in tqdm(generated_files, desc="Streaming HDF fragments to SQLite") :
        chunk_df = pd.read_hdf(file_path, key='pop')
        chunk_df.to_sql(TABLE_NAME, conn, if_exists='append', index=False)

print(f"Completed database generation [ID:{args_md5}]. Finishing up...\n")
shutil.rmtree(TMP_DIR)