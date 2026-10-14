"""Figures 1-3 stay inside their boxes and never go below 7.5 pt."""

import matplotlib.pyplot as plt
import pytest

import figstyle as S
import make_figures as F


@pytest.mark.parametrize("name", list(F.FIGURES))
def test_figure_layout_is_clean(name):
    fig = F.FIGURES[name]()
    try:
        assert S.check_layout(fig) == []
        assert fig.get_size_inches()[0] <= S.COL_2 + 1e-9
    finally:
        plt.close(fig)


def test_no_stroke_prediction_wording():
    # Device_Role §6: the model detects a marker; it does not predict stroke.
    src = open(F.__file__, encoding="utf-8").read().lower()
    assert "pre-stroke risk" not in src
