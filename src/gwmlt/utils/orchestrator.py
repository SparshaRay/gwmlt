"""
Scheduler Orchestrator and Utilities
"""

import os
import signal
import socket
import threading
import functools
import subprocess

from pathlib import Path
from threadpoolctl import threadpool_limits

import dask
from dask_jobqueue import HTCondorCluster
from dask.distributed import Client


# Orchestrator for Dask+HTCondor ------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

class HTCondorOrchestrator :
    """
    Orchestrator for Dask+HTCondor.
    """

    def __init__(
        self, 
        log_dir: Path, 
        max_jobs: int, 
        cores: int = 1, 
        memory: str = "2GB", 
        disk: str = "128MB"
    ) -> None :
        """
        Initialize the Dask HTCondor orchestrator.

        Parameters
        ----------
        log_dir : Path
            Directory path where Dask and HTCondor logs will be kept.
        max_jobs : int
            Maximum number of concurrent HTCondor workers to scale up to.
        cores : int, optional
            Number of CPU cores to request per worker (default is 1).
            Note that `threadpoolctl` is used to strictly limit the number of threads per worker to 1.
        memory : str, optional
            Amount of RAM to request per worker (default is "2GB").
        disk : str, optional
            Amount of disk space to request per worker (default is "128MB").
        """
        
        self.log_dir = log_dir

        # Dask setup 

        dask.config.set({
            "distributed.scheduler.worker-saturation" : 1.0,
            "distributed.scheduler.work-stealing" : True,
            "distributed.worker.memory.target"    : 0.75,
            "distributed.worker.memory.spill"     : 0.80,
            "distributed.worker.memory.pause"     : 0.85,
            "distributed.worker.memory.terminate" : 0.90,
        })
        
        self.job_extra_directives = {
            "accounting_group"        : "iucaa.gwml0",
            "accounting_group_user"   : "sparsha.ray",
            "should_transfer_files"   : "Yes",
            "when_to_transfer_output" : "ON_EXIT_OR_EVICT",
            "stream_output"           : "True",
            "stream_error"            : "True",
            "getenv"                  : "True",
            "universe"                : "vanilla",
            "output" : str(log_dir / "worker_$(ClusterId).out"),
            "error"  : str(log_dir / "worker_$(ClusterId).err"),
            "log"    : str(log_dir / "condor_master.log"),
        }

        # HTCondor setup

        print("Initializing Dask+HTCondor cluster...")

        self.cluster = HTCondorCluster(
            cores  = cores,
            memory = memory,
            disk   = disk,

            local_directory = "/tmp", 
            log_directory   = log_dir,

            job_extra_directives = self.job_extra_directives,
            scheduler_options    = {"interface": "bond0"},

            job_script_prologue=[
                "export OMP_NUM_THREADS=1",
                "export OPENBLAS_NUM_THREADS=1",
                "export MKL_NUM_THREADS=1",
                "export NUMEXPR_NUM_THREADS=1",
                "export DASK_DISTRIBUTED__LOGGING__DISTRIBUTED=info",
                "export DASK_DISTRIBUTED__LOGGING__BOKEH=critical",
            ]
        )

        # Initialize cluster

        self.cluster.adapt(minimum=0, maximum=max_jobs)
        self.client = Client(self.cluster)

        print("Successfully initialized Dask+HTCondor cluster.")

        dask_port = int(self.client.scheduler_info()['services']['dashboard'])
        username = os.environ.get('USER') or os.environ.get('USERNAME')

        print(
            "\nYou may forward the Dask Dashboard port with:\n"
            f"`ssh -N -L <available_port>:localhost:{dask_port} {username}@{socket.gethostname()}`\n"
            "and then open up localhost:<available_port> in your web browser.\n"
        )


    def close(self) -> None :
        """Shutdown Dask and HTCondor workers"""

        print("Shutting down Dask+HTCondor cluster...")

        try :
            # Consolidate dask logs and close the client and cluster

            scheduler_logs = self.client.get_scheduler_logs()[::-1]
            with open(self.log_dir / "dask_scheduler.log", "w") as f :
                for log_entry in scheduler_logs : f.write(f"{str(log_entry)}\n\n")

            self.client.close()
            self.cluster.close()

            # Close HTCondor workers and consolidate logs

            if Path(self.job_extra_directives['log']).is_file() :
                subprocess.run(
                    ["condor_wait", "-wait", "15", self.job_extra_directives['log']], 
                    check=True, stdout=subprocess.DEVNULL
                )
            
            subprocess.run( # Force shutdown
                ["condor_rm", "-forcex", "-all"], check=False, 
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )

        except Exception as e : print(f"Error while closing Dask cluster: {e}.")

        print("Dask+HTCondor workers closed. Consolidating HTCondor logs...")

        for log_type in ["out", "err"] :
            filepaths = list(self.log_dir.glob(f"worker_*.{log_type}"))
            if not filepaths : continue

            master_file_path = self.log_dir / f"condor_{log_type}.log"
            with open(master_file_path, "w") as f :
                for filepath in filepaths :
                    filename = filepath.name
                    f.write(f"{filename} {'='*100} \n\n")
                    with open(filepath, "r") as chunk : f.write(chunk.read() + "\n\n")
                    filepath.unlink()

        print("Successfully consolidated logs and closed Dask+HTCondor cluster.\n")


# Helper functions --------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Watchdog decorator

def watchdog(timeout_seconds: int) -> callable :
    """
    Kill process after `timeout_seconds` have elapsed.
    This terminates both the running job and the parent dask worker.
    """

    def decorator(func) :
        @functools.wraps(func)
        def wrapper(*args, **kwargs) :

            def single_threaded_func(*args, **kwargs) :
                with (
                    threadpool_limits(limits=1, user_api='blas'),
                    threadpool_limits(limits=1, user_api='openmp')
                ) : return func(*args, **kwargs)

            def kill() :
                print(
                    f"[WATCHDOG] > Process exceeded timeout of {timeout_seconds}s. "
                    "Terminating both the process and the worker with SIGKILL.\n"
                    f"Input arguments and keyword arguments : \n{args, kwargs}",
                    flush=True
                )
                os.kill(os.getpid(), signal.SIGKILL)

            timer = threading.Timer(timeout_seconds, kill)
            timer.start()

            try : return single_threaded_func(*args, **kwargs)
            finally : timer.cancel()

        return wrapper
    return decorator

# Logger class for spawner

class SpawnerLogger :
    def __init__(self, log_dir: Path, stream) -> None :
        self.log_file = open(log_dir / "spawner.log", 'a')
        self.stream = stream

    def write(self, message) -> None :
        self.stream.write(message)
        self.log_file.write(message)
        self.flush()

    def flush(self) -> None :
        self.stream.flush()
        self.log_file.flush()