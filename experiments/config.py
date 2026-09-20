"""Config-as-data: one resolved object per figure.

A figure is reproducible from a :class:`Config` plus a seed -- not from an
editing session.

``n_boot`` is 1000 rather than the 100 the R runs used. The minP calibration is
exact at any ``n_boot``, but its p-values live on the ``1/(n_boot+1)`` grid and
with ``L = dx + dy - 3`` search levels several paths tie at the floor, so small
``n_boot`` costs power: at ``d = 8`` (``L = 13``) ``n_boot = 100`` recovers about
half the power available at ``n_boot = 1000``, and 1000 is where the curve has
flattened. It is also what the ``_bonf`` comparators need to be able to reject at
all (floor ``L/(n_boot+1)``). XGBoost
hyperparameters live beside this file in ``tuning/`` -- one JSON per marginal
setting, carried over from the tuning runs so the package is self-contained.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List

# repo root: experiments/config.py -> parents[1]
REPO_ROOT = Path(__file__).resolve().parents[1]
TUNING_ROOT = Path(__file__).resolve().parent / "tuning"


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
    n_boot: int = 1000
    normalise: bool = False
    adaptive: List[str] = field(default_factory=list)
    competitors: List[str] = field(default_factory=list)

    def xgb_params(self, setting: str) -> dict:
        path = TUNING_ROOT / f"n{self.n}_numclass{self.d}" / f"tune_{setting}_results.json"
        return json.loads(path.read_text())["xgb"]

    def to_dict(self) -> dict:
        return asdict(self)


def size_config(xsetting: str, ysetting: str, **overrides) -> Config:
    """Null-calibration (size) block: strength 0, measure rejection at the level."""
    cfg = dict(
        name=f"size_{xsetting}_{ysetting}",
        n=1000, d=8, xsetting=xsetting, ysetting=ysetting, intsetting="step",
        strengths=[0.0], reps=1000, n_boot=1000, normalise=False,
        adaptive=["tree", "ordinal", "tree_bonf", "ordinal_bonf", "max", "euclid", "mGCM"],
        competitors=["ankan", "chi_sq"],
    )
    cfg.update(overrides)
    return Config(**cfg)


def power_config(xsetting: str, ysetting: str, intsetting: str, **overrides) -> Config:
    """The standard power-figure block.

    The adaptive method set depends on the interaction: ``tree`` for binary_tree,
    ``ordinal`` for step, plus the non-adaptive comparators and the competitors.
    """
    adaptive = ["max", "euclid", "mGCM"]
    if intsetting == "binary_tree":
        adaptive = ["tree", "tree_bonf"] + adaptive
    elif intsetting == "step":
        adaptive = ["ordinal", "ordinal_bonf"] + adaptive

    cfg = dict(
        name=f"power_{xsetting}_{ysetting}_{intsetting}",
        n=1000, d=8, xsetting=xsetting, ysetting=ysetting, intsetting=intsetting,
        strengths=[round(0.2 * k, 1) for k in range(1, 10)],  # 0.2 .. 1.8
        reps=200, n_boot=1000, normalise=False,
        adaptive=adaptive, competitors=["ankan", "chi_sq", "multinomial"],
    )
    cfg.update(overrides)
    return Config(**cfg)
