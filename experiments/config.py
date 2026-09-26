"""Config-as-data: one resolved object per figure.

A figure is reproducible from a :class:`Config` plus a seed -- not from an
editing session.

``n_boot`` is 1000 rather than the 100 the R runs used. The minP calibration is
exact at any ``n_boot``, but its p-values live on the ``1/(n_boot+1)`` grid and
with ``L = dx + dy - 3`` search levels several paths tie at the floor, so small
``n_boot`` costs power: at ``d = 8`` (``L = 13``) ``n_boot = 100`` recovers about
half the power available at ``n_boot = 1000``, and 1000 is where the curve has
flattened. It is also what the ``_bonf`` comparators need to be able to reject at
all (floor ``L/(n_boot+1)``).

``dx`` and ``dy`` are separate axes: the applications are asymmetric (hundreds of
levels on one side, a handful on the other), and the dense ``Sigma`` is
``(dx*dy)^2``, so symmetric ``d`` in the hundreds is out of reach anyway.

XGBoost hyperparameters live beside this file in ``tuning/`` -- one JSON per
``(n, num_class, setting)``, written by ``tune.py``. ``X`` and ``Y`` are tuned
separately because each one's law given ``Z`` depends only on its own setting and
dimension (the interaction preserves both margins).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List

# repo root: experiments/config.py -> parents[1]
REPO_ROOT = Path(__file__).resolve().parents[1]
TUNING_ROOT = Path(__file__).resolve().parent / "tuning"

# The configuration every figure used before d became an axis.
DEFAULT_N, DEFAULT_D = 1000, 8


def tuning_path(n: int, num_class: int, setting: str, root: Path = TUNING_ROOT) -> Path:
    return Path(root) / f"n{n}_numclass{num_class}" / f"tune_{setting}_results.json"


@dataclass
class Config:
    name: str
    n: int
    dx: int
    dy: int
    xsetting: str
    ysetting: str
    intsetting: str
    strengths: List[float]
    reps: int = 200
    n_boot: int = 1000
    normalise: bool = False
    adaptive: List[str] = field(default_factory=list)
    competitors: List[str] = field(default_factory=list)

    def xgb_params(self, setting: str, num_class: int) -> dict:
        path = tuning_path(self.n, num_class, setting)
        if not path.exists():
            raise FileNotFoundError(
                f"No tuned hyperparameters at {path}. Run: "
                f"python experiments/tune.py --n {self.n} --d {num_class} --settings {setting}"
            )
        return json.loads(path.read_text())["xgb"]

    def to_dict(self) -> dict:
        return asdict(self)


def _dims(overrides: dict) -> tuple[int, int]:
    """Pop ``d`` (sets both) / ``dx`` / ``dy`` from ``overrides``."""
    d = overrides.pop("d", DEFAULT_D)
    return overrides.pop("dx", d), overrides.pop("dy", d)


def _suffix(n: int, dx: int, dy: int) -> str:
    """Name suffix for anything off the historical ``n = 1000, d = 8`` default."""
    if (n, dx, dy) == (DEFAULT_N, DEFAULT_D, DEFAULT_D):
        return ""
    return f"_n{n}_dx{dx}_dy{dy}"


def size_config(xsetting: str, ysetting: str, **overrides) -> Config:
    """Null-calibration (size) block: strength 0, measure rejection at the level."""
    dx, dy = _dims(overrides)
    n = overrides.pop("n", DEFAULT_N)
    cfg = dict(
        name=f"size_{xsetting}_{ysetting}{_suffix(n, dx, dy)}",
        n=n, dx=dx, dy=dy, xsetting=xsetting, ysetting=ysetting, intsetting="step",
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

    dx, dy = _dims(overrides)
    n = overrides.pop("n", DEFAULT_N)
    cfg = dict(
        name=f"power_{xsetting}_{ysetting}_{intsetting}{_suffix(n, dx, dy)}",
        n=n, dx=dx, dy=dy, xsetting=xsetting, ysetting=ysetting, intsetting=intsetting,
        strengths=[round(0.2 * k, 1) for k in range(1, 10)],  # 0.2 .. 1.8
        reps=200, n_boot=1000, normalise=False,
        adaptive=adaptive, competitors=["ankan", "chi_sq", "multinomial"],
    )
    cfg.update(overrides)
    return Config(**cfg)


def d_grid(xsetting: str, ysetting: str, intsetting: str, dxs, dy: int, n: int, **overrides) -> List[Config]:
    """The ``d`` axis: one power block per ``dx``, with ``n`` and ``dy`` held fixed.

    This is the experiment for "does the merging advantage widen as ``d`` grows?"
    -- the per-observation signal of the DGP is roughly constant in ``d`` (flat
    for lin/vee/hat, decaying gently for sin/sig), so a test that pays for every
    cell loses power with ``d`` and one that merges should not.
    """
    return [power_config(xsetting, ysetting, intsetting, n=n, dx=dx, dy=dy, **overrides) for dx in dxs]
