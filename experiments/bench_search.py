"""Runtime of the two search directions, as a function of the number of labels.

One search is one bootstrap draw's worth of work, so the numbers here scale
directly into ``n_boot``. ``Sigma`` is shared across draws, which the divisive
search exploits by building its block-sum table once
(:class:`~catci.blocks.SigmaBlocks`); the table build is reported separately and
amortised over ``n_boot + 1`` paths in the "per draw" columns.

Usage:
    python bench_search.py                     # tree, dx = dy
    python bench_search.py --structure ordinal
    python bench_search.py --dy 10             # the paper's motivating dX x dY
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from catci.blocks import SigmaBlocks
from catci.search import divisive_search, greedy_search
from catci.structure import Ordinal, Tree


def random_sigma(p: int, rng: np.random.Generator) -> np.ndarray:
    A = rng.standard_normal((p + 5, p))
    return A.T @ A / (p + 5)


def timed(fn, repeats: int) -> float:
    """Seconds per call, best of ``repeats`` (best-of resists scheduler noise)."""
    best = np.inf
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def structures(name: str, dx: int, dy: int):
    if name == "tree":
        return Tree.binary(dx), Tree.binary(dy)
    if name == "ordinal":
        return Ordinal(), Ordinal()
    raise ValueError(f"Unknown structure {name!r}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", default="tree", choices=["tree", "ordinal"])
    ap.add_argument("--dx", type=int, nargs="+", default=[4, 8, 16, 32, 64])
    ap.add_argument("--dy", type=int, default=None, help="fixed dy; default dy = dx")
    ap.add_argument("--levels", type=int, nargs="+", default=[2, 4],
                    help="truncation budgets to report alongside the full path")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    budgets = args.levels

    head = f"{'dx':>4} {'dy':>4} {'L':>4} {'merge':>10} {'split':>10} {'ratio':>7} {'table':>9}"
    head += "".join(f"{'split@' + str(b):>10}" for b in budgets)
    print(f"\n### search runtime, structure={args.structure}, "
          f"amortised over n_boot={args.n_boot} draws  (ms per draw)")
    print(head)
    print("-" * len(head))

    for dx in args.dx:
        dy = args.dy if args.dy is not None else dx
        p = dx * dy
        Sigma = random_sigma(p, rng)
        T = rng.multivariate_normal(np.zeros(p), Sigma)
        xs, ys = structures(args.structure, dx, dy)

        t_table = timed(lambda: SigmaBlocks(Sigma, dx, dy), args.repeats)
        blocks = SigmaBlocks(Sigma, dx, dy)
        amortised = t_table / (args.n_boot + 1)

        t_merge = timed(lambda: greedy_search(T, Sigma, dx, dy, xs, ys), args.repeats)
        t_split = timed(
            lambda: divisive_search(T, Sigma, dx, dy, xs, ys, sigma_blocks=blocks),
            args.repeats,
        ) + amortised
        trunc = [
            timed(
                lambda b=b: divisive_search(
                    T, Sigma, dx, dy, xs, ys, sigma_blocks=blocks, max_levels=b
                ),
                args.repeats,
            ) + amortised
            for b in budgets
        ]

        row = (f"{dx:>4} {dy:>4} {dx + dy - 3:>4} {1e3 * t_merge:>10.2f} "
               f"{1e3 * t_split:>10.2f} {t_merge / t_split:>6.1f}x {1e3 * t_table:>9.1f}")
        row += "".join(f"{1e3 * t:>10.2f}" for t in trunc)
        print(row)

    print("\n  merge    = greedy_search, full path (dx + dy - 3 levels)")
    print("  split    = divisive_search, full path, Sigma table amortised")
    print("  ratio    = merge / split, >1 means the divisive direction is faster")
    print("  table    = one-off SigmaBlocks build, ms (shared by all n_boot + 1 draws)")
    print("  split@L  = divisive_search truncated to L splits, table amortised")


if __name__ == "__main__":
    main()
