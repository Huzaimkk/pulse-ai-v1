"""One function per paper table/figure number. Every number in the paper comes from
here or from a script that calls it (Rule 8)."""

import math

import numpy as np
import pandas as pd

from .registry import REGISTRY


def prf(y, yhat) -> tuple[float, float, float]:
    """Precision, recall, F1 for the positive class (1). A ratio with an empty
    denominator is NaN, never 0 or 1. F1 is the harmonic mean of P and R, so it can
    never exceed both (the arithmetic error that reached the previous thesis)."""
    y, yhat = np.asarray(y).astype(bool), np.asarray(yhat).astype(bool)
    tp = int(np.sum(y & yhat))
    fp = int(np.sum(~y & yhat))
    fn = int(np.sum(y & ~yhat))
    p = tp / (tp + fp) if tp + fp else math.nan
    r = tp / (tp + fn) if tp + fn else math.nan
    if math.isnan(p) or math.isnan(r):
        return p, r, math.nan
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)


def table_iii() -> pd.DataFrame:
    """Table III: every feature with a published normal range or threshold, and its
    citation — straight from the registry declarations."""
    import pulse_ecg.features  # noqa: F401
    t = REGISTRY.table()
    return t[(t["normal"] != "") | (t["abnormal"] != "")].reset_index(drop=True)
