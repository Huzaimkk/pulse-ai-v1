from pulse_ecg import config


def test_config_has_every_section():
    cfg = config.load()
    for key in ("paths", "ptbxl", "recording", "filter", "beat", "bounds", "qc",
                "cohort", "gate2", "parallel"):
        assert key in cfg, key


def test_rules_encoded_in_config():
    cfg = config.load()
    assert cfg["ptbxl"]["records"] == "records500"          # Rule 1
    assert cfg["ptbxl"]["fs"] == 500
    assert cfg["cohort"]["split_by"] == "patient_id"         # Rule 2
    assert cfg["cohort"]["gate1_min_mean_age"] >= 55         # Gate 1
    assert cfg["cohort"]["sinus_codes"]                      # sinus-rhythm cohort
    assert cfg["gate2"]["primary"] == "ludb"                 # Gate 2


def test_ranges_are_ordered():
    cfg = config.load()
    lo, hi = cfg["filter"]["band"]
    assert 0 < lo < hi < cfg["ptbxl"]["fs"] / 2
    for name, (a, b) in cfg["qc"].items():
        assert a < b, name


def test_paths_resolve_inside_project():
    assert config.path("ptbxl").is_relative_to(config.ROOT)
    assert config.path("refs_bib").name == "refs.bib"
