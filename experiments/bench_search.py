"""Wall time of the greedy search: per-candidate loop vs the vectorised batch.

One search is one bootstrap draw's worth of work. The loop column times the
original implementation (``_greedy_search_loop``) on a few draws and reports ms
per draw; the batch columns time ``greedy_search_paths`` on ``B`` draws at once.

Usage:
    python bench_search.py                          # tree, dx = dy
    python bench_search.py --structure greedy --dx 4 8 16
    python bench_search.py --dx 32 --n-boot 1000 10000 --no-loop
    python bench_search.py --dx 30 --dy 10          # the paper's motivating dX x dY
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from catci.search import _greedy_search_loop, greedy_search_paths
from catci.structure import Ordinal, Saturated, Tree


def random_sigma(p: int, rng: np.random.Generator) -> np.ndarray:
    A = rng.standard_normal((p + 5, p))
    return A.T @ A / (p + 5)


def structures(name: str, dx: int, dy: int):
    if name == "tree":
        return Tree.binary(dx), Tree.binary(dy)
    if name == "ordinal":
        return Ordinal(), Ordinal()
    if name == "greedy":
        return Saturated(), Saturated()
    raise ValueError(f"Unknown structure {name!r}")


CLOCK = time.perf_counter


def timed(fn) -> float:
    t0 = CLOCK()
    fn()
    return CLOCK() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", default="tree", choices=["tree", "ordinal", "greedy"])
    ap.add_argument("--dx", type=int, nargs="+", default=[4, 8, 16, 32])
    ap.add_argument("--dy", type=int, default=None, help="fixed dy; default dy = dx")
    ap.add_argument("--n-boot", type=int, nargs="+", default=[100, 1000])
    ap.add_argument("--loop-draws", type=int, default=5, help="draws timed for the loop")
    ap.add_argument("--no-loop", action="store_true", help="skip the (slow) loop reference")
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cpu", action="store_true",
                    help="time process CPU, not wall clock (for a loaded machine; "
                         "pair with VECLIB_MAXIMUM_THREADS=1 / OMP_NUM_THREADS=1)")
    args = ap.parse_args()
    global CLOCK
    if args.cpu:
        CLOCK = time.process_time

    rng = np.random.default_rng(args.seed)
    head = f"{'dx':>4} {'dy':>4} {'p':>6} {'L':>4} {'loop ms/draw':>13}"
    for B in args.n_boot:
        head += f" {'B=' + str(B) + ' s':>12} {'ms/draw':>9} {'speedup':>8}"
    clock = "CPU" if args.cpu else "wall"
    print(f"\n### greedy search, structure={args.structure}, n_jobs={args.n_jobs}, {clock} time")
    print(head)
    print("-" * len(head))

    for dx in args.dx:
        dy = args.dy if args.dy is not None else dx
        p = dx * dy
        Sigma = random_sigma(p, rng)
        xs, ys = structures(args.structure, dx, dy)
        L = dx + dy - 4

        loop = np.nan
        if not args.no_loop:
            T = rng.standard_normal((p, args.loop_draws))
            loop = timed(lambda: [_greedy_search_loop(T[:, b], Sigma, dx, dy, xs, ys)
                                  for b in range(args.loop_draws)]) / args.loop_draws

        row = f"{dx:>4} {dy:>4} {p:>6} {L:>4} {1e3 * loop:>13.2f}"
        for B in args.n_boot:
            T = rng.standard_normal((p, B))
            t = timed(lambda: greedy_search_paths(T, Sigma, dx, dy, xs, ys, n_jobs=args.n_jobs))
            row += f" {t:>12.3f} {1e3 * t / B:>9.3f} {loop / (t / B):>7.0f}x"
        print(row, flush=True)

    print("\n  loop    = original per-candidate search, one draw at a time")
    print("  B=... s = greedy_search_paths on B draws, total seconds (includes setup)")
    print("  speedup = loop ms/draw over batch ms/draw")


if __name__ == "__main__":
    main()
