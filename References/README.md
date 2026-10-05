# References/ — machine-checkable complexity derivations

**Policy, set by the user 2026-09-25 and binding on all future work:**

> Every channel-estimation scheme and every receiver design in this project
> ships with a Lean file in `References/` encoding its complexity derivation,
> so the derivation can be machine-verified rather than taken on trust.

This is not optional documentation. A new estimator or receiver is not
finished until its `flops_*.m` and its `References/*.lean` both exist and
agree.

## Why this exists

Complexity derivations here are produced by reading MATLAB source and
counting operations by hand. That process is error-prone in a specific and
dangerous way: a wrong operation count produces a *plausible* number, and a
plausible wrong number is worse than no number, because nothing downstream
flags it. This project has already been burned repeatedly by plausible
arithmetic — see `notes/DTMUSIC_PORT_AUDIT.md` for the catalogue.

Lean removes one whole class of that error: the algebra.

## What a Lean file here DOES and DOES NOT verify

Read this before trusting one.

**It DOES verify:**
- that the stated total is the sum of the stated terms;
- that concrete evaluations at the shipped configuration are arithmetically
  correct (so the smoke-test numbers in the `.m` file are machine-checked);
- that claimed scaling laws follow from the term structure (e.g. "cost is
  proportional to `1 + NGS`", "dominant term is linear in `Phat`");
- that claimed dominance ratios hold (e.g. "subdominant terms are under 1%").

**It DOES NOT verify:**
- **that the operation counts match the MATLAB source.** That link is a human
  audit step and Lean cannot see it. If someone miscounts `match_proj` calls,
  Lean will happily prove theorems about the wrong algorithm.
- that the MATLAB implements the published algorithm.
- anything about numerical accuracy, only about operation counts.

So the Lean file makes the *derivation* trustworthy given the counts. Keeping
the *counts* trustworthy still requires reading the code, and each `.lean`
file must cite the exact source file and line range each count came from.

## Required contents of each file

1. A `Cfg` structure whose fields are the real config parameters
   (`L`, `N`, `Qtau`, `Qnu`, `NGS`, cost constants …), named as in
   `oddm_config.m`.
2. One definition per cost term, each carrying a comment naming the **source
   file and lines** it was counted from.
3. A `total` that is literally the sum of those terms.
4. A definition of the shipped configuration (e.g. Figure 7's point).
5. Concrete theorems fixing the evaluated numbers at that configuration,
   proved by `rfl` or `decide`. **These must match the `.m` file's smoke
   test exactly** — that agreement is the point.
6. At least one structural theorem capturing the derivation's headline claim
   (scaling law, dominance ratio, or both).
7. A header block stating what is *not* verified, and any outstanding
   empirical validation.

## Conventions

- **Lean 4, no Mathlib.** Everything is `Nat` arithmetic closable by `rfl`,
  `decide`, or `simp` with core lemmas. This keeps the files checkable with a
  bare `lean` install and no project scaffolding. Do not add a Mathlib
  dependency for convenience.
- **Naturals only.** Operation counts are naturals. Where the MATLAB uses a
  fractional `Phat` (DT-MUSIC's expected order is ~4.4), state theorems at
  integer `Phat` and give the polynomial in `Phat` so any value can be
  evaluated by hand. Do not smuggle in rationals.
- **Ratios as cross-multiplications.** `a / b = k` becomes `a = k * b`, to
  stay in `Nat` and avoid truncating division.
- One file per algorithm, named `<Algorithm>Flops.lean`.

## Status

| algorithm | `flops_*.m` | `References/*.lean` | clean timing validation |
|---|---|---|---|
| DD-RELAX | `flops_dd_relax.m` | `DDRelaxFlops.lean` | outstanding |
| SAGE | `flops_sage.m` | `SAGEFlops.lean` | outstanding |
| OMP | `flops_omp.m` | `OMPFlops.lean` | outstanding |
| DT-MUSIC (win3) | `flops_dtmusic_ord.m` | queued | outstanding |

### First cross-algorithm result

DD-RELAX and SAGE cost almost exactly the same per frame at the shipped
configuration, and for a reason the Lean files make explicit:

| term | DD-RELAX | SAGE |
|---|---|---|
| grid correlation (dominant) | 354,931,200 | 354,931,200 |
| atom generation | 725,760 | 259,200 |
| joint LS / Gram / solve | 1,218,024 | — (none in SAGE) |
| residual recompute | (inside joint LS) | 497,664 |
| **total** | **356,874,984** | **355,688,064** |

The dominant terms are *identical* — both make `Phat*(1+NGS) = 45` calls to
the same `match_proj` against the same dictionary, so both pay
`45 * Qtau*Qnu*L*N * 8`. Everything that distinguishes the two algorithms
lives in the remaining 0.3%.

The practical consequence: at fixed `NGS`, per-frame cost is set almost
entirely by **how many grid correlations an estimator performs**, not by what
it does with the results. Any future optimisation that does not reduce the
correlation count is optimising the 0.3%.

OMP confirms this sharply. It is single-phase — no `NGS`, no cyclic
refinement — so it makes `Phat = 9` correlations where the other two make
`Phat*(1+NGS) = 45`:

| estimator | grid correlations | total flops/frame | vs OMP |
|---|---|---|---|
| OMP | 9 | 71,705,232 | 1.00x |
| SAGE | 45 | 355,688,064 | 4.96x |
| DD-RELAX | 45 | 356,874,984 | 4.98x |

The cost ratio is `45/9 = 5`, and the measured ratios are 4.96 and 4.98 — the
entire difference between these three estimators is the correlation count,
to within half a percent.

Read alongside accuracy, that is the actual complexity story of Figure 7:
**DD-RELAX buys its accuracy with 5x the correlations.** At gamma_p = 20 dB
its channel RMSE is 0.066 against OMP's 0.10, for 4.98x the arithmetic. SAGE
pays DD-RELAX's price (4.96x) for 0.082. Whether that trade is worth it is a
judgement for the paper; the numbers for making it are now machine-checked.

A caution carried over from the OMP file: its subdominant terms are **1.003%**
of total, not under 1%. The agent writing it checked rather than inheriting
DD-RELAX's threshold, and the claim would have been false. Do not copy
dominance bounds between files — recheck each one.

### The complete four-way result

DT-MUSIC sits outside the correlation-count pattern entirely — it has no grid
correlation in its Doppler stage, and its cost is dominated (94.0%) by the
delay-dictionary extraction. With its Lean file done, all four are comparable.

FLOP counts are from the `.lean` files, at `Phat = Pmax = 9` (worst case,
integer per policy). Accuracy is from Profile 5 at Eb/N0 = 16 dB — **one
harness, same rows, same frames** — because the per-estimator numbers quoted
elsewhere in this project come from different diagnostic rigs and are not
comparable with each other.

| estimator | flops/frame | vs cheapest | chan RMSE | BER |
|---|---|---|---|---|
| DT-MUSIC (win3) | 53,440,464 | 1.00x | 0.1546 | 6.49e-3 |
| OMP | 71,705,232 | 1.34x | 0.1066 | 4.95e-4 |
| SAGE | 355,688,064 | 6.66x | 0.0823 | 2.62e-4 |
| DD-RELAX | 356,874,984 | 6.68x | 0.0695 | 2.06e-4 |
| perfect CSI | — | — | — | 8.66e-5 |

**Cost and accuracy are monotonically ordered.** No estimator is dominated on
the pair, so each genuinely buys accuracy with arithmetic, and the figure has
a real cost/accuracy frontier rather than a winner.

**With one exception, and it is the strongest single statement available
here: DD-RELAX dominates SAGE outright.** Same cost to within 0.3% — they
perform the identical 45 grid correlations — but 16% better channel RMSE
(0.0695 vs 0.0823) and 21% better BER. That is not a trade-off; at equal
arithmetic DD-RELAX is simply better, and the Lean files establish the "equal
arithmetic" half rather than asserting it.

Read the other direction, DD-RELAX costs **6.68x** DT-MUSIC for 2.2x the
channel accuracy, and **4.98x** OMP for 1.53x. Whether those are good trades
is a judgement for the paper; the numbers for making it are machine-checked.

Caveat on the BER column: at these levels BER is dominated by rare bad frames
(see Profile 10's header for the worked case), so use the RMSE column for
comparisons. Frame counts also differ between rows (5,000 to 12,850).
| PT-MMSE receiver | not started | not started | — |
| SIC-MMSE receiver | not started | not started | — |
| CMC-MMSE receiver | not started | not started | — |

"Clean timing validation" means a `cputime` measurement on an otherwise idle
machine with realized `Phat` recorded. Both timing datasets gathered so far
were `tic`/`toc` wall clock taken while collection workers saturated the
machine, and are too noisy to adjudicate model differences of under ~10%.
Do not cite them as validation.
