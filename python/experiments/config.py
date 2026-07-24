"""Config-as-data: one resolved object per figure.

A figure is reproducible from a :class:`Config` plus a seed -- not from an
editing session (the failure mode CODE_REVIEW.md 3 describes). XGBoost
hyperparameters are read from the same R tuning JSONs the R runs used, so the
learners are identical across languages.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List

# repo root: python/experiments/config.py -> parents[2]
REPO_ROOT = Path(__file__).resolve().parents[2]
TUNING_ROOT = REPO_ROOT / "data-raw" / "tuning"


@dataclass
class Config:
    name: str
    n: int
    d: int
    xsetting: str
    ysetting: str
    intsetting: str
    strengths: List[float]
    reps: int = 200
    n_boot: int = 100
    nfolds: int = 5
    normalise: bool = False
    adaptive: List[str] = field(default_factory=list)
    competitors: List[str] = field(default_factory=list)

    def xgb_params(self, setting: str) -> dict:
        path = TUNING_ROOT / f"n{self.n}_numclass{self.d}" / f"tune_{setting}_results.json"
        return json.loads(path.read_text())["xgb"]

    def to_dict(self) -> dict:
        return asdict(self)


def power_config(xsetting: str, ysetting: str, intsetting: str, **overrides) -> Config:
    """The standard power-figure block (mirrors data-raw/power/simulation_power.R).

    Adaptive method set depends on the interaction, exactly as in the R driver:
    ``tree`` for binary_tree, ``ordinal`` for step, plus the non-adaptive
    comparators and the three competitors.
    """
    adaptive = ["max", "euclid", "mGCM"]
    if intsetting == "binary_tree":
        adaptive = ["tree"] + adaptive
    elif intsetting == "step":
        adaptive = ["ordinal"] + adaptive

    cfg = dict(
        name=f"power_{xsetting}_{ysetting}_{intsetting}",
        n=1000, d=8, xsetting=xsetting, ysetting=ysetting, intsetting=intsetting,
        strengths=[round(0.2 * k, 1) for k in range(1, 10)],  # 0.2 .. 1.8
        reps=200, n_boot=100, nfolds=5, normalise=False,
        adaptive=adaptive, competitors=["ankan", "chi_sq", "multinomial"],
    )
    cfg.update(overrides)
    return Config(**cfg)
