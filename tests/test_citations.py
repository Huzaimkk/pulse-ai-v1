"""Every threshold carries a citation, and every citation exists in refs.bib."""

import re

import pulse_ecg.features  # noqa: F401  (registers features)
from pulse_ecg import config
from pulse_ecg.registry import REGISTRY


def bib_keys():
    text = config.path("refs_bib").read_text(encoding="utf-8")
    return set(re.findall(r"@\w+\{([^,\s]+),", text))


def test_refs_bib_has_the_required_entries():
    keys = bib_keys()
    for k in ("wagner2020ptbxl", "kalyakulina2020ludb", "laguna1997qtdb",
              "martinezselles2020bayes", "wolf1991framingham", "attia2019af",
              "kamel2024arcadia"):
        assert k in keys, k


def test_every_feature_citation_is_in_refs_bib():
    missing = REGISTRY.cite_keys() - bib_keys()
    assert not missing, f"cite keys not in refs.bib: {sorted(missing)}"
