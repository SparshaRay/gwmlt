from pathlib import Path

def test_config_override() :
    from gwmlt.config import config, config_override

    old_f_ref = config.waveform.f_ref
    old_database_root = config.database_root

    with config_override({
        "waveform.f_ref": 0.0,
        "database_root": Path("/new/root")
    }) :
        new_f_ref = config.waveform.f_ref
        new_database_root = config.database_root

    final_f_ref = config.waveform.f_ref
    final_database_root = config.database_root
        
    assert old_f_ref == final_f_ref != new_f_ref == 0.0
    assert old_database_root == final_database_root != new_database_root == Path("/new/root")