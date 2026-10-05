/-
  DTMUSICOrdFlops.lean — machine-checkable complexity derivation for
  sliding-window DT-MUSIC with Eq. (30) order selection, as implemented in
  this repository.

  Mirrors:  Comm Functions/ODDM Functions/DD-RELAX-paper/flops_dtmusic_ord.m
  Algorithm: Comm Functions/ODDM Functions/DD-RELAX-paper/est_dtmusic_ord.m
  Policy:    References/README.md

  WHAT THIS FILE VERIFIES
    * that `totalBeta`/`totalAtom` are exactly the sum of the stated terms;
    * that the evaluated numbers at the SHIPPED Profile 5 configuration
      (F = 10, nAtom = 6, ordCrit = "beta", Pc = Pmax = 9) are correct --
      these are the same numbers flops_dtmusic_ord.m's smoke test prints,
      cross-checked term-by-term against a fresh independent Python
      re-evaluation of the same integer arithmetic before this file was
      written;
    * DOMINANCE: extract_dict (the per-path, per-delay-grid-point
      dictionary build) is more than 94% of the total at the shipped
      config -- NOT the >99% seen in DDRelaxFlops.lean; this estimator's
      subdominant terms are real, not negligible, and the theorem below
      states the true ~6% bound rather than reusing DD-RELAX's <1% figure;
    * the ordCrit BRANCH RATIO: choosing ordCrit="atom" instead of the
      shipped "beta" costs between 5x and 6x more (measured 5.700x) at this
      configuration -- the quantitative form of why the shipped choice
      (accuracy-motivated, per saved_profiles.m's Profile 5 comment: "atom"
      and "beta" pick the same order ~79/80 trials, so "beta" was kept
      because it is what Eq. 30 specifies, not for cost reasons) also
      happens to be the cheap one.

  WHAT THIS FILE DOES NOT VERIFY
    * THAT THE OPERATION COUNTS MATCH est_dtmusic_ord.m. Lean cannot see the
      MATLAB. Each term below cites the source lines it was counted from
      (both in the algorithm file and in flops_dtmusic_ord.m); that reading
      is a human audit step and is where an error would hide.
    * anything about numerical accuracy -- only operation counts.
    * the model against measured runtime. NO timing data exists for this
      estimator at all (not even noisy data) -- outstanding, same as OMP.
    * the K-FLATNESS claim in flops_dtmusic_ord.m's header. That is an
      explicitly UNVALIDATED STRUCTURAL PREDICTION (total FLOPs grow only
      mildly with nAtom/K because extract_dict, K-independent, dominates)
      -- this file's `extractRefitFlops` term makes the K-dependence
      concrete and machine-checkable in principle, but no theorem below
      claims flatness is validated. The K=1..8 s/frame timing data
      originally offered as support for that prediction was WITHDRAWN by
      desktop after re-inspection: it was a scheduling artifact (one arm's
      own blocks spanned a 2.4x range purely from other concurrent arms
      finishing and freeing cores mid-run), not a property of the
      estimator, and MUST NOT be cited as evidence for or against anything
      here.
    * this file only models ordSel = true (the shipped default). The
      ordSel = false single-trial branch exists in the .m file but is not
      encoded here.

  Lean 4, no Mathlib (see References/README.md).
-/

namespace CWS.DTMUSICOrd

/-- Configuration parameters, named as in `oddm_config.m` where they
overlap with the other estimators, and as in `est_dtmusic_ord.m`'s own
argument list otherwise. -/
structure Cfg where
  /-- delay-domain support of the pilot window, `L = L1 + L2 + 1`. -/
  L       : Nat
  /-- Doppler bins per frame. -/
  N       : Nat
  /-- delay grid points -- the SAME dictionary DD-RELAX/SAGE/OMP search,
  reused here only in the final tau/phi extraction, never for Doppler. -/
  Qtau    : Nat
  /-- max model order (also the ordSel sweep's upper bound before the
  `min` with `mSmooth - 1`). -/
  Pmax    : Nat
  /-- sliding-window length, `num_init_frames` in saved_profiles.m --
  NOT part of `oddm_config.m`; the single biggest cost lever, so it is a
  required field with no silent default. -/
  F       : Nat
  /-- spatial-smoothing subarray length (shipped default: `N`, no
  smoothing). -/
  mSmooth : Nat
  /-- THIS ESTIMATOR'S OWN Doppler search-grid size, from `gridOv` --
  a DIFFERENT grid from `Qnu` (DD-RELAX/SAGE/OMP's dictionary grid, built
  from `PiNu`). Do not conflate the two. -/
  Ngrid   : Nat
  /-- K, delay atoms per retained Doppler (shipped: 6). -/
  nAtom   : Nat
  /-- flops per complex multiply-accumulate (6 mul + 2 add). -/
  cCMAC   : Nat
  /-- calibrated flops per grid point of one `s_atom` evaluation. -/
  cAtom   : Nat
  /-- calibrated flops constant for a dense Hermitian eigendecomposition
  of size `n`, taken as `cEig * n^3` -- a numerical-library detail, not
  derived exactly, same spirit as `cAtom`. -/
  cEig    : Nat

/-- `N - mSmooth + 1`, the number of smoothing sub-windows per lag row per
frame (est_dtmusic_ord.m line 232). At the shipped `mSmooth = N` this is 1. -/
def nSub (c : Cfg) : Nat := c.N - c.mSmooth + 1

/-- `min(Pmax, mSmooth - 1)`, the order-selection sweep's upper bound
(est_dtmusic_ord.m line 302). -/
def PnMax (c : Cfg) : Nat := min c.Pmax (c.mSmooth - 1)

/-! ### Cost terms, in the order est_dtmusic_ord.m executes them -/

/-- WINDOW RESHAPE. `Yall(:,:,f) = ifftshift(reshape(...),2) * conj(F_N)`
for each of `F` window frames (est_dtmusic_ord.m lines 240-243) -- an
`L x N` times `N x N` DFT-matrix product, `L*N^2` cmacs, done `F` times
(flops_dtmusic_ord.m line 237). -/
def windowFlops (c : Cfg) : Nat := c.F * c.L * (c.N * c.N) * c.cCMAC

/-- COVARIANCE ACCUMULATION. Nested loop `f=1..F, d=1..L, s=1..nSub`, each
iteration an `mSmooth`-length outer product (est_dtmusic_ord.m lines
263-276; flops_dtmusic_ord.m line 238). Often co-dominant with the window
term: at the shipped `mSmooth = N`, `nSub = 1`, so this collapses to the
same `F*L*N^2` shape. -/
def covFlops (c : Cfg) : Nat := c.F * c.L * (nSub c) * (c.mSmooth * c.mSmooth) * c.cCMAC

/-- EIGENDECOMPOSITION. One `mSmooth x mSmooth` Hermitian `eig()` per
frame, REGARDLESS of `ordSel` (est_dtmusic_ord.m line 293, before the
order-selection loop; flops_dtmusic_ord.m line 239). -/
def eigFlops (c : Cfg) : Nat := c.cEig * (c.mSmooth * c.mSmooth * c.mSmooth) * c.cCMAC

/-- `Sum_{Pn=1}^{PnMax} (mSmooth - Pn) = PnMax*mSmooth - PnMax*(PnMax+1)/2`,
the ordSel=true sweep's total `Esteer*Vn` work (est_dtmusic_ord.m lines
315-344, the `Vn = Vec(:,Pn+1:mSmooth)` / `Esteer*Vn` step at line 396 of
est_dtmusic_ord.m's helper; flops_dtmusic_ord.m lines 244-249). -/
def sumMSmoothMinusPn (c : Cfg) : Nat :=
  (PnMax c) * c.mSmooth - (PnMax c) * (PnMax c + 1) / 2

/-- ORDER-SELECTION MUSIC-SPECTRUM SWEEP. `Ngrid x mSmooth` times
`mSmooth x (mSmooth-Pn)`, summed over the `PnMax` trial orders
(flops_dtmusic_ord.m line 249). This is where `ordSel` multiplies a COST,
but a cheap one -- see the header note on the branch ratio. -/
def musicFlops (c : Cfg) : Nat := c.Ngrid * c.mSmooth * (sumMSmoothMinusPn c) * c.cCMAC

/-- Per-trial `beta_hat` solve (Gram + regularised solve + `Efull'*Y_last`),
approximated at size `Pc` for every one of the `PnMax` trials
(est_dtmusic_ord.m lines 416-419; flops_dtmusic_ord.m line 254). -/
def betaSolveFlops (c : Cfg) (Pc : Nat) : Nat :=
  (PnMax c) * (Pc*Pc*c.N + Pc*Pc*Pc + Pc*c.N*c.L) * c.cCMAC

/-- DOMINANT. FINAL TAU/PHI EXTRACTION, always runs ONCE at the winning
order's `Pc`, regardless of `ordSel`/`ordCrit`. For each of `Pc` retained
paths, EVERY `Qtau` delay-grid atom is generated (`s_atom`, `cAtom*L*N`)
AND transformed into the delay-time domain via an `L x N` times `N x N`
`conj(F_N)` product (`L*N^2` cmacs) BEFORE scoring -- a genuinely larger
per-atom cost than DD-RELAX/SAGE/OMP's plain `s_atom` calls
(est_dtmusic_ord.m lines 463-466, `dtmusicw_extract_multi`, and the
byte-equivalent `dtmusicw_extract_tau_phi` at lines 500-502;
flops_dtmusic_ord.m line 257). -/
def extractDictFlops (c : Cfg) (Pc : Nat) : Nat :=
  Pc * c.Qtau * (c.cAtom * c.L * c.N + c.L * (c.N * c.N)) * c.cCMAC

/-- `Sum_{k=1}^{K} k`, `Sum k^2`, `Sum k^3` for the per-path OMP-style
growing joint refit inside `dtmusicw_extract_multi` (only when
`nAtom > 1`; est_dtmusic_ord.m lines 470-479; flops_dtmusic_ord.m lines
262-266). Assumes `nAtom <= Qtau` (true at every shipped/tested config:
`nAtom` is at most 8, `Qtau = 79`), matching `min(K,Qtau) = K`. -/
def sumK  (c : Cfg) : Nat := c.nAtom * (c.nAtom + 1) / 2
def sumK2 (c : Cfg) : Nat := c.nAtom * (c.nAtom + 1) * (2 * c.nAtom + 1) / 6
def sumK3 (c : Cfg) : Nat := (sumK c) * (sumK c)

/-- Subdominant refit cost on top of `extractDictFlops`, present only for
`nAtom > 1` (the shipped config has `nAtom = 6`, so this is NOT zero here
-- unlike DD-RELAX/SAGE/OMP's `nAtom`-free structure). Absent entirely at
`nAtom = 1` (`dtmusicw_extract_tau_phi`'s closed-form per-atom phi, already
inside `extractDictFlops`'s per-atom accounting). -/
def extractRefitFlops (c : Cfg) (Pc : Nat) : Nat :=
  Pc * c.nAtom * c.Qtau * c.L * c.cCMAC +
  Pc * ((sumK2 c) * c.L + (sumK3 c) + (sumK c) * c.L) * c.cCMAC

/-- ordCrit == "atom" ONLY: the expensive extraction re-runs INSIDE the
order-selection loop, once per trial, evaluated at that trial's own `Pn`
(est_dtmusic_ord.m lines 325-334; flops_dtmusic_ord.m lines 277-282).
Confirmed by TWO independent readings (this file's author and desktop,
2026-09-25) that ordCrit == "beta" (the shipped default) does NOT run this
term at all -- it is what makes "beta" the cheap branch. -/
def ordCritAtomFlops (c : Cfg) : Nat :=
  c.cCMAC * c.Qtau * (c.cAtom * c.L * c.N + c.L * (c.N * c.N)) *
    ((PnMax c) * (PnMax c + 1) / 2)

/-- Total per-frame cost, SHIPPED branch (`ordCrit = "beta"`). -/
def totalBeta (c : Cfg) (Pc : Nat) : Nat :=
  windowFlops c + covFlops c + eigFlops c + musicFlops c + betaSolveFlops c Pc +
    extractDictFlops c Pc + extractRefitFlops c Pc

/-- Total per-frame cost, `ordCrit = "atom"` branch. -/
def totalAtom (c : Cfg) (Pc : Nat) : Nat :=
  totalBeta c Pc + ordCritAtomFlops c

/-! ### The shipped configuration

Profile 5's DT-MUSIC-WIN config as of `impl = "win3"` (saved_profiles.m,
the struct with `'ordSel', true, 'ordCrit', "beta", 'nAtom', 6`):
`N = 16, Pmax = 9, PiTau = 32, PiNu = 16` give `Qtau = 79, L = 12` via
`oddm_config.m`, same as DDRelaxFlops.lean's `fig7`. `F = 10` is
`num_init_frames` from the SAME struct, confirmed to be the actual sliding-
window length by reading `sim_fun_ODDM_DTMUSICW_PTMMSE.m` line 268
(`yp_win = zeros(P.L*P.N, nif)` where `nif = cs.num_init_frames`) --
`frames_per_trial` (20) is a separate parameter controlling window reuse
between estimator calls, NOT the window length itself, and is not part of
this per-frame cost model. `mSmooth = 16` (unset in the struct, so the
estimator's own default `= N` applies) and `Ngrid = 127` are computed from
`gridOv`'s default (32) via `floor(2*nu_max*gridOv/F0) + 1` at this
project's shipped EVA/v=500 channel profile, read directly from
`flops_dtmusic_ord.m`'s own printed `params.Ngrid` at this configuration
rather than hand-derived. -/
def fig7Ord : Cfg :=
  { L := 12, N := 16, Qtau := 79, Pmax := 9, F := 10, mSmooth := 16,
    Ngrid := 127, nAtom := 6, cCMAC := 8, cAtom := 30, cEig := 9 }

/-- The worst-case order, `min(Pmax, mSmooth-1) = min(9,15) = 9`. DT-MUSIC's
early "expected" order is measured at ~4.40 (desktop, 2026-09-24 audit) --
a FRACTIONAL value with no `Nat` representation, so per policy this file
states only the worst-case bound. The polynomial in `Pc` is given
explicitly in each `def` above so any integer value (e.g. 4, rounding the
measured 4.40 down) can be evaluated by substitution. -/
def PcFig7Ord : Nat := 9

/-! ### Concrete checks

These fix the numbers flops_dtmusic_ord.m's smoke test prints at the
SHIPPED configuration (F=10, nAtom=6, ordCrit="beta", Pc=Pmax=9), and were
cross-checked against an independent Python re-evaluation of the same
integer arithmetic before this file was written -- not merely asserted. -/

example : nSub fig7Ord = 1 := rfl
example : PnMax fig7Ord = 9 := rfl
example : sumMSmoothMinusPn fig7Ord = 99 := rfl
example : sumK fig7Ord = 21 := rfl
example : sumK2 fig7Ord = 91 := rfl
example : sumK3 fig7Ord = 441 := rfl

theorem window_fig7ord        : windowFlops fig7Ord              = 245760    := rfl
theorem cov_fig7ord           : covFlops fig7Ord                 = 245760    := rfl
theorem eig_fig7ord           : eigFlops fig7Ord                 = 294912    := rfl
theorem music_fig7ord         : musicFlops fig7Ord               = 1609344   := rfl
theorem betaSolve_fig7ord     : betaSolveFlops fig7Ord PcFig7Ord = 270216    := rfl
theorem extractDict_fig7ord   : extractDictFlops fig7Ord PcFig7Ord = 50236416 := rfl
theorem extractRefit_fig7ord  : extractRefitFlops fig7Ord PcFig7Ord = 538056 := rfl
theorem ordCritAtom_fig7ord   : ordCritAtomFlops fig7Ord         = 251182080 := rfl

theorem totalBeta_fig7ord : totalBeta fig7Ord PcFig7Ord = 53440464 := rfl
theorem totalAtom_fig7ord : totalAtom fig7Ord PcFig7Ord = 304622544 := rfl

/-- The seven "beta"-branch terms account for the total exactly. -/
theorem beta_terms_sum_to_total :
    windowFlops fig7Ord + covFlops fig7Ord + eigFlops fig7Ord + musicFlops fig7Ord +
      betaSolveFlops fig7Ord PcFig7Ord + extractDictFlops fig7Ord PcFig7Ord +
      extractRefitFlops fig7Ord PcFig7Ord
      = totalBeta fig7Ord PcFig7Ord := rfl

/-- The "atom" branch is exactly the "beta" total plus the one extra term. -/
theorem atom_is_beta_plus_extra :
    totalAtom fig7Ord PcFig7Ord = totalBeta fig7Ord PcFig7Ord + ordCritAtomFlops fig7Ord := rfl

/-! ### Structural claims -/

/-- DOMINANCE. `extract_dict` alone is more than 94% of the total at the
shipped configuration -- equivalently, everything else combined is under
6%. THIS IS A WEAKER BOUND THAN DDRelaxFlops.lean's `<1%`, stated
correctly rather than reused: DT-MUSIC's subdominant terms (the covariance
accumulation and the order-selection MUSIC sweep in particular) are real
costs at this window length and grid size, not negligible the way
DD-RELAX's joint-LS terms are. -/
-- TIGHT: 320,404,800 vs 320,642,784, a 0.07% margin. True and provable, but
-- do not assume slack -- a term edit that adds even a small amount of
-- subdominant cost could flip this without extract_dict's own share moving
-- much at all. Re-check the margin, not just the truth value, after any
-- edit to a subdominant term.
theorem subdominant_under_six_percent :
    100 * (totalBeta fig7Ord PcFig7Ord - extractDictFlops fig7Ord PcFig7Ord)
      < 6 * totalBeta fig7Ord PcFig7Ord := by decide

/-- ORDCRIT BRANCH RATIO. Choosing `ordCrit = "atom"` instead of the
shipped `"beta"` costs BETWEEN 5x and 6x more at this configuration
(measured 5.700x) -- stated as a cross-multiplied bracket to stay in `Nat`
and avoid truncating division. This is the quantitative form of why the
shipped choice, made on accuracy grounds (saved_profiles.m: "atom" and
"beta" pick the same order in ~79/80 trials, so "beta" was kept because it
is what Eq. 30 specifies), also happens to be the cheap one. -/
theorem ordcrit_atom_between_5x_and_6x_beta :
    5 * totalBeta fig7Ord PcFig7Ord < totalAtom fig7Ord PcFig7Ord ∧
      totalAtom fig7Ord PcFig7Ord < 6 * totalBeta fig7Ord PcFig7Ord := by decide

end CWS.DTMUSICOrd
