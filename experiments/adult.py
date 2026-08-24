"""The UCI adult-income dataset, preprocessed as Ankan & Textor (AAAI-23) describe it.

The paper uses eleven variables (Fig. 1a) and states two discretisations
verbatim (p. 12186):

    we discretized the variable "Age" into the categories < 21, 21-30, ..., 61-70,
    > 70 and the variable "HoursPerWeek" into the categories <= 20, 21-30, 30-40,
    > 40

Everything else about the coding is left implicit, so this module fixes the
remaining choices explicitly and records why. The choices are checkable: the
degrees of freedom printed in Fig. 1b are a deterministic function of the coding,
``df = (k-1)(r-1) * n_strata``, and all six rows are consistent with a single
one -- see :func:`fig1b_df_check` and ``tests/test_adult.py``.

Data provenance
---------------
``data/adult-train.csv`` is the UCI "adult.data" train split (32,561 rows,
headerless), vendored from the ``jbrownlee/Datasets`` mirror so that a run needs
no network and cannot silently change under us. It is the same split pgmpy
ships as its ``adult`` dataset.

Coding decisions not stated by the paper
----------------------------------------
* **Missing values.** ``workclass``, ``occupation`` and ``native-country`` carry
  a literal ``?``. The paper never mentions them -- but Fig. 1b's df column
  settles it. Those six df require Workclass to show 6 levels and Occupation 13
  in the n = 1000 sample of Fig. 1. Keeping ``?`` as a category cannot produce
  those counts (it can only add a level); dropping the ``?`` rows can, and does:
  listwise deletion leaves 30,162 rows with 7 workclasses and 14 occupations, and
  a 1000-row subsample that misses the rarest levels (``Without-pay``,
  ``Never-worked``, ``Armed-Forces``) lands exactly on 6 and 13. Two of ten seeds
  tried reproduce all six df. Hence ``drop_missing=True`` is the default -- it is
  what the paper did, inferred from the only checkable quantity it reports.
* **Education order.** Taken from the dataset's own ``education-num`` column,
  which is the authoritative 1..16 ordering (Preschool .. Doctorate) rather than
  a hand-written one.
* **HoursPerWeek bin 3.** The paper's own list overlaps at 30 ("21-30, 30-40").
  Read as ``31-40``, the only non-overlapping reading.
* **Income** is binarised at $50K/year, as the paper says.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent / "data"
ADULT_CSV = DATA_DIR / "adult-train.csv"

# UCI adult.names column order (the file is headerless).
RAW_COLUMNS = [
    "age", "workclass", "fnlwgt", "education", "education-num",
    "marital-status", "occupation", "relationship", "race", "sex",
    "capital-gain", "capital-loss", "hours-per-week", "native-country", "income",
]

Kind = Literal["ordinal", "categorical"]


@dataclass(frozen=True)
class Variable:
    """One analysis variable: its short code, its type, and its level labels.

    ``kind`` drives which of the paper's statistics applies to a pair (Q1 for
    ordinal-ordinal, Q2 for categorical-ordinal, Q3 for categorical-categorical)
    and which conditional probability model estimates its propensities.
    """

    name: str
    code: str  # the abbreviation used in Fig. 1
    kind: Kind
    levels: tuple[str, ...]

    @property
    def n_levels(self) -> int:
        return len(self.levels)


AGE_LABELS = ("<21", "21-30", "31-40", "41-50", "51-60", "61-70", ">70")
AGE_EDGES = (-np.inf, 20, 30, 40, 50, 60, 70, np.inf)

HOURS_LABELS = ("<=20", "21-30", "31-40", ">40")
HOURS_EDGES = (-np.inf, 20, 30, 40, np.inf)

EDUCATION_LEVELS = (
    "Preschool", "1st-4th", "5th-6th", "7th-8th", "9th", "10th", "11th", "12th",
    "HS-grad", "Some-college", "Assoc-voc", "Assoc-acdm", "Bachelors", "Masters",
    "Prof-school", "Doctorate",
)

INCOME_LEVELS = ("<=50K", ">50K")


def load_raw() -> pd.DataFrame:
    """The vendored CSV, with column names attached and whitespace stripped."""
    if not ADULT_CSV.exists():
        raise FileNotFoundError(
            f"{ADULT_CSV} is missing. Re-vendor it with:\n"
            "  curl -sSo experiments/data/adult-train.csv "
            "https://raw.githubusercontent.com/jbrownlee/Datasets/master/adult-train.csv"
        )
    df = pd.read_csv(ADULT_CSV, header=None, names=RAW_COLUMNS, skipinitialspace=True)
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip()
    return df


def _bin(series: pd.Series, edges, labels) -> pd.Categorical:
    return pd.cut(series, bins=list(edges), labels=list(labels), right=True, ordered=True)


@dataclass(frozen=True)
class AdultData:
    """Preprocessed adult income data.

    ``codes`` holds one column per analysis variable as **1-based** contiguous
    integer labels -- the convention ``catci`` uses throughout -- and ``variables``
    maps each column name to its :class:`Variable` spec. Levels are those actually
    observed, so ``variables[v].n_levels`` always equals ``codes[v].max()``.
    """

    codes: pd.DataFrame
    variables: dict[str, Variable]

    @property
    def n(self) -> int:
        return len(self.codes)

    def kind(self, name: str) -> Kind:
        return self.variables[name].kind

    def n_levels(self, name: str) -> int:
        return self.variables[name].n_levels

    def subsample(self, n: int, seed: int) -> "AdultData":
        """A row subsample, re-coded so labels stay contiguous and 1-based.

        Rare categories can vanish in a small subsample -- which is exactly what
        Fig. 1b shows, where Workclass has 6 observed levels and Occupation 13
        rather than their full-data counts. Re-coding keeps that faithful.
        """
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(self.codes), size=n, replace=False)
        sub = self.codes.iloc[np.sort(idx)].reset_index(drop=True)
        return _recode(sub, self.variables)


def _recode(codes: pd.DataFrame, variables: dict[str, Variable]) -> AdultData:
    """Drop unobserved levels and renumber the survivors 1..m, preserving order."""
    out = {}
    specs = {}
    for name, var in variables.items():
        observed = np.sort(codes[name].unique())
        remap = {old: new for new, old in enumerate(observed, start=1)}
        out[name] = codes[name].map(remap).astype(np.int64)
        specs[name] = Variable(
            name=var.name,
            code=var.code,
            kind=var.kind,
            levels=tuple(var.levels[o - 1] for o in observed),
        )
    return AdultData(codes=pd.DataFrame(out), variables=specs)


def load(drop_missing: bool = True) -> AdultData:
    """Load and preprocess the eleven analysis variables of Fig. 1a.

    Parameters
    ----------
    drop_missing:
        If True (the default), drop rows where ``workclass``, ``occupation`` or
        ``native-country`` is ``?`` (32,561 -> 30,162 rows). This is what the
        paper did; see the module docstring for the Fig. 1b df evidence. Pass
        False to keep ``?`` as its own category and the raw 32,561 rows.
    """
    raw = load_raw()
    if drop_missing:
        missing = (raw[["workclass", "occupation", "native-country"]] == "?").any(axis=1)
        raw = raw.loc[~missing].reset_index(drop=True)

    frame = pd.DataFrame(index=raw.index)
    variables: dict[str, Variable] = {}

    def add_ordinal_from_labels(name, code, values: pd.Series, levels):
        cat = pd.Categorical(values, categories=list(levels), ordered=True)
        if cat.isna().any():
            bad = sorted(set(values[pd.isna(cat)]))
            raise ValueError(f"{name}: values outside the declared level set: {bad}")
        frame[name] = cat.codes + 1
        variables[name] = Variable(name, code, "ordinal", tuple(levels))

    def add_categorical(name, code, values: pd.Series):
        levels = tuple(sorted(values.unique()))
        cat = pd.Categorical(values, categories=list(levels))
        frame[name] = cat.codes + 1
        variables[name] = Variable(name, code, "categorical", levels)

    # --- ordinal variables ---------------------------------------------------
    add_ordinal_from_labels("Age", "Age", _bin(raw["age"], AGE_EDGES, AGE_LABELS).astype(str), AGE_LABELS)
    add_ordinal_from_labels(
        "HoursPerWeek", "HrPW",
        _bin(raw["hours-per-week"], HOURS_EDGES, HOURS_LABELS).astype(str), HOURS_LABELS,
    )
    # education-num is the dataset's own authoritative 1..16 ordering
    frame["Education"] = raw["education-num"].astype(np.int64)
    variables["Education"] = Variable("Education", "Edct", "ordinal", EDUCATION_LEVELS)
    add_ordinal_from_labels("Income", "Incm", raw["income"], INCOME_LEVELS)

    # --- categorical variables -----------------------------------------------
    add_categorical("Workclass", "Wrkc", raw["workclass"])
    add_categorical("MaritalStatus", "MrSt", raw["marital-status"])
    add_categorical("Occupation", "Occp", raw["occupation"])
    add_categorical("Relationship", "Rltn", raw["relationship"])
    add_categorical("Race", "Race", raw["race"])
    add_categorical("Sex", "Sex", raw["sex"])
    add_categorical("NativeCountry", "NtvC", raw["native-country"])

    return _recode(frame, variables)


# --------------------------------------------------------------------------- #
# Fig. 1b cross-check
# --------------------------------------------------------------------------- #
# The paper's Fig. 1b prints p and df for six pairs given Z = {Age, Sex}, from a
# stratified mutual-information test on 1000 samples. df is deterministic given
# the coding, so it validates preprocessing without implementing that baseline.
FIG1B = (
    # (X, Y, df)
    ("Education", "Workclass", 1050),
    ("Occupation", "Workclass", 840),
    ("Relationship", "HoursPerWeek", 210),
    ("Income", "Occupation", 168),
    ("Income", "Workclass", 70),
    ("Income", "HoursPerWeek", 42),
)


def fig1b_df_check(data: AdultData, z=("Age", "Sex")) -> pd.DataFrame:
    """Recompute Fig. 1b's df column from ``data``'s observed level counts.

    ``df = (k-1)(r-1) * n_strata`` where ``n_strata`` is the product of the
    conditioning variables' level counts. Returns a frame with the paper's df
    beside ours, so a caller (or a test) can compare.
    """
    n_strata = int(np.prod([data.n_levels(v) for v in z]))
    rows = []
    for x, y, paper_df in FIG1B:
        k, r = data.n_levels(x), data.n_levels(y)
        ours = (k - 1) * (r - 1) * n_strata
        rows.append(
            dict(X=x, Y=y, k=k, r=r, n_strata=n_strata,
                 df_ours=ours, df_paper=paper_df, match=ours == paper_df)
        )
    return pd.DataFrame(rows)


if __name__ == "__main__":
    full = load()
    print(f"full data: n = {full.n}")
    print(pd.DataFrame([
        dict(variable=v.name, code=v.code, kind=v.kind, levels=v.n_levels)
        for v in full.variables.values()
    ]).to_string(index=False))
    print("\nFig. 1b df check on a 1000-row subsample (paper's Fig. 1 setting, seed 6):")
    print(fig1b_df_check(full.subsample(1000, seed=6)).to_string(index=False))
