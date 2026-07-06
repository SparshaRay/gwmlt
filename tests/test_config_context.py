from pathlib import Path

def test_config_override() :
    """
    Test that the config override with context manager is working as expected.
    """

    from gwmlt.config import config, config_override

    old_f_ref = config.waveform.f_ref
    old_project_root = config.project_root

    with config_override({
        "waveform.f_ref": 0.0,
        "project_root": Path("/new/root")
    }) :
        new_f_ref = config.waveform.f_ref
        new_project_root = config.project_root

    final_f_ref = config.waveform.f_ref
    final_project_root = config.project_root
        
    assert old_f_ref == final_f_ref != new_f_ref == 0.0
    assert old_project_root == final_project_root != new_project_root == Path("/new/root")