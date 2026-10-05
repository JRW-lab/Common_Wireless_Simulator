/-
  OMPFlops.lean — machine-checkable complexity derivation for the OMP
  channel estimator, as implemented in this repository.

  Mirrors:  Comm Functions/ODDM Functions/DD-RELAX-paper/flops_omp.m
  Algorithm: Comm Functions/ODDM Functions/DD-RELAX-paper/est_omp.m
  Exemplar:  References/DDRelaxFlops.lean (structure, naming, tactics)
  Policy:    References/README.md

  WHAT THIS FILE VERIFIES
    * that `total` is exactly the sum of the stated terms;
    * that the evaluated numbers at the shipped configuration are correct
      (these are the same numbers `flops_omp.m`'s own smoke test prints,
      confirmed 2026-09-25 by running
      `flops_omp(P,'Phat',9)` at N=16,M=64,fc=4e9,sub=15000,alpha=0.3,Q=4,
      v_kmh=500,Pmax=9,PiTau=32,PiNu=16,frame_layout="guard_end", which
      gives Qtau=79, Qnu=65, L=12);
    * that the closed-form triangular/square/cube sums used by the growing
      joint-LS term equal the actual recursive sums, at that concrete Phat;
    * that the dominant (`match_proj`) term is exactly linear in `Phat`;
    * a dominance ratio for the subdominant terms — see the note below,
      because for OMP (unlike DD-RELAX) it does NOT clear the 1% bar.

  WHAT THIS FILE DOES NOT VERIFY
    * THAT THE OPERATION COUNTS MATCH est_omp.m. Lean cannot see the MATLAB.
      Each count below cites the source lines it was read from; that
      reading is a human audit step and is where an error would hide.
    * that the MATLAB implements the published algorithm (OMP.pdf,
      Tropp & Gilbert 2007, Algorithm 3).
    * anything about numerical accuracy — only operation counts.
    * the model against measured runtime. NO timing validation exists for
      OMP at all: the only OMP data collected so far is a PiTau/PiNu RMSE
      grid sweep, which is not a timing measurement. `flops_omp.m`'s own
      header (2026-09-25) records this as outstanding.

  TWO FACTS ESTABLISHED BY THE 2026-09-25 AUDIT, load-bearing for the
  constants below:
    * measured Phat = 9.00 of Pmax = 9 (early-stop fires in only ~2% of
      frames), so worst case and expected coincide almost exactly here —
      the opposite situation from DT-MUSIC, whose order selection typically
      retains far fewer than Pmax.
    * OMP is strictly on-grid: there is no `parabola_refine` call anywhere
      in `est_omp.m`, unlike a sub-grid-refined estimator.

  Lean 4, no Mathlib (see References/README.md). Only `rfl`, `by decide`,
  and `by simp [...]` with core `Nat` lemmas are used, matching the
  exemplar — no `induction`, no `ring`, no `omega`.
-/

namespace CWS.OMP

/-- Configuration parameters, named as in `oddm_config.m`. Unlike
`CWS.DDRelax.Cfg` there is no `NGS` field: OMP has no Gauss-Seidel /
cyclic-refinement phase at all (`est_omp.m`'s own docstring, lines 13-14:
"no Gauss-Seidel polish -- this is what keeps OMP the coarsest of the three
baselines"). -/
structure Cfg where
  /-- delay-domain support of the pilot window, `L = L1 + L2 + 1`. -/
  L     : Nat
  /-- Doppler bins per frame. -/
  N     : Nat
  /-- delay grid points, `ceil(tau_max/Dtau) + 1`. -/
  Qtau  : Nat
  /-- Doppler grid points, `2*ceil(nu_max/Dnu) + 1`. -/
  Qnu   : Nat
  /-- flops per complex multiply-accumulate (6 mul + 2 add). -/
  cCMAC : Nat
  /-- calibrated flops per grid point of one `s_atom` evaluation. -/
  cAtom : Nat

/-- Pilot-window length, the inner dimension of every matvec below. -/
def LN (c : Cfg) : Nat := c.L * c.N

/-! ### Call counts

`est_omp.m` has exactly ONE phase (lines 39-52): `for t = 1:Pmax` calling
`match_proj` once (line 40) then `s_atom` once (line 44) then
`joint_ls_omp` once (line 48), with an early break that sets `Phat = t`
(line 51); if it never breaks, line 53 sets `Phat = Pmax`. So — same
"always exactly Phat calls" pattern as DD-RELAX's Phase 1 — the loop makes
exactly `Phat` calls to each of `match_proj` and `s_atom`, with NO `(1+NGS)`
multiplier anywhere, because there is no second phase. -/

/-- `match_proj` calls per frame: one phase, no NGS multiplier. -/
def numProj (c : Cfg) (Phat : Nat) : Nat := Phat

/-- `s_atom` calls per frame: one per iteration, no duplicate re-evaluation
(contrast DD-RELAX's `2 + 3*NGS` per path). -/
def numAtom (c : Cfg) (Phat : Nat) : Nat := Phat

/-! ### Cost terms -/

/-- DOMINANT. Each `match_proj` call computes a dense `(Qtau*Qnu) x LN` by
`LN x 1` complex matvec against the precomputed dictionary, identical
per-call cost to DD-RELAX's and SAGE's `match_proj` (`match_proj.m` line
18). -/
def projFlops (c : Cfg) (Phat : Nat) : Nat :=
  numProj c Phat * (c.Qtau * c.Qnu) * LN c * c.cCMAC

/-- Subdominant, approximate: `s_atom` evaluation is `O(L*N)` per call, same
calibrated `cAtom` caveat as DD-RELAX (interpolation cost is a library
detail, not an algorithmic one). -/
def atomFlops (c : Cfg) (Phat : Nat) : Nat :=
  numAtom c Phat * LN c * c.cAtom

/-! ### The growing joint-LS refit

What distinguishes OMP from plain on-grid matching pursuit is the JOINT
least-squares refit over ALL currently-selected atoms after every new atom
(`est_omp.m` line 48, `joint_ls_omp` defined lines 60-66: Gram `S'*S` at
line 63, ridge-regularized solve against `S'*yp` at line 65). Unlike
DD-RELAX/SAGE, whose cyclic phase always solves a FIXED-size `Phat x Phat`
system a fixed number of times, OMP's system has `t` columns at loop
iteration `t`, so the Gram/solve/`S'*yp` costs are SUMMED over `t = 1..Phat`
rather than multiplied by a constant. These sums are given here by genuine
recursion (the actual `Sum_{t=1}^{n} t^k`, not the closed-form shortcut),
so each definition is transparently "the sum `flops_omp.m` describes"; the
closed-form identities `flops_omp.m` uses instead (`sum_t = n(n+1)/2`, etc.)
are checked against these below as concrete facts. -/

/-- `Sum_{t=1}^{n} t`. -/
def sumUpTo (n : Nat) : Nat :=
  match n with
  | 0 => 0
  | k + 1 => (k + 1) + sumUpTo k

/-- `Sum_{t=1}^{n} t^2`. -/
def sumSqUpTo (n : Nat) : Nat :=
  match n with
  | 0 => 0
  | k + 1 => (k + 1) * (k + 1) + sumSqUpTo k

/-- `Sum_{t=1}^{n} t^3`. -/
def sumCubeUpTo (n : Nat) : Nat :=
  match n with
  | 0 => 0
  | k + 1 => (k + 1) * (k + 1) * (k + 1) + sumCubeUpTo k

/-- Gram `S'*S` at iteration `t` costs `t^2 * LN` cmacs; summed over the
whole loop. -/
def gramFlops (c : Cfg) (Phat : Nat) : Nat :=
  c.cCMAC * LN c * sumSqUpTo Phat

/-- The regularized solve at iteration `t` costs `O(t^3)`; summed over the
whole loop. -/
def solveFlops (c : Cfg) (Phat : Nat) : Nat :=
  c.cCMAC * sumCubeUpTo Phat

/-- `S'*yp` at iteration `t` costs `t * LN` cmacs; summed over the whole
loop. -/
def styFlops (c : Cfg) (Phat : Nat) : Nat :=
  c.cCMAC * LN c * sumUpTo Phat

/-- REDUNDANT CALL, kept visible on purpose. `est_omp.m` line 56 makes ONE
EXTRA `joint_ls_omp(S,yp)` call AFTER the loop, at the final `Phat`-column
`S` — this recomputes exactly what the last in-loop iteration (`t = Phat`)
already produced. It is a genuine, if minor, redundancy in the source, and
`flops_omp.m` counts it (lines ~109-111) rather than folding it silently
into the summed terms above: it is a single call at size `Phat`, not part
of the `t = 1..Phat` sum. -/
def finalExtraFlops (c : Cfg) (Phat : Nat) : Nat :=
  c.cCMAC * LN c * (Phat * Phat) + c.cCMAC * (Phat * Phat * Phat)
    + c.cCMAC * LN c * Phat

def jointLsFlops (c : Cfg) (Phat : Nat) : Nat :=
  gramFlops c Phat + solveFlops c Phat + styFlops c Phat + finalExtraFlops c Phat

/-- Total per-frame estimation cost. -/
def total (c : Cfg) (Phat : Nat) : Nat :=
  projFlops c Phat + atomFlops c Phat + jointLsFlops c Phat

/-! ### The shipped configuration

N=16, M=64, fc=4e9, sub=15000, alpha=0.3, Q=4, v_kmh=500, Pmax=9, PiTau=32,
PiNu=16, frame_layout="guard_end". The derived grid sizes `Qtau = 79`,
`Qnu = 65` and `L = 12` were read from `oddm_config.m` at those parameters
by running `flops_omp` directly (2026-09-25); they match the values
independently recorded for Figure 7 elsewhere in this project, since both
configs share the same underlying `N`/`M`/`v_kmh`/grid parameters. -/
def shippedCfg : Cfg :=
  { L := 12, N := 16, Qtau := 79, Qnu := 65, cCMAC := 8, cAtom := 30 }

/-- `Pmax = 9`, and the measured `Phat` is 9.00 of 9 (early stop fires in
only ~2% of frames), so this is both the worst case and the expected case
almost exactly. -/
def PhatShipped : Nat := 9

/-! ### Concrete checks

These fix the numbers `flops_omp.m`'s smoke test prints
(`flops_omp(P,'Phat',9)` at the shipped configuration). If a term in the
MATLAB is edited without updating the count here, these stop being `rfl`. -/

example : LN shippedCfg = 192 := rfl
example : numProj shippedCfg PhatShipped = 9 := rfl
example : numAtom shippedCfg PhatShipped = 9 := rfl

theorem proj_shipped  : projFlops  shippedCfg PhatShipped = 70986240 := rfl
theorem atom_shipped  : atomFlops  shippedCfg PhatShipped = 51840    := rfl
theorem gram_shipped  : gramFlops  shippedCfg PhatShipped = 437760   := rfl
theorem solve_shipped : solveFlops shippedCfg PhatShipped = 16200    := rfl
theorem sty_shipped   : styFlops   shippedCfg PhatShipped = 69120    := rfl
theorem finalExtra_shipped :
    finalExtraFlops shippedCfg PhatShipped = 144072 := rfl
theorem jointLs_shipped : jointLsFlops shippedCfg PhatShipped = 667152 := rfl
theorem total_shipped   : total       shippedCfg PhatShipped = 71705232 := rfl

/-- The three top-level terms account for the total exactly — no hidden
remainder. -/
theorem terms_sum_to_total :
    projFlops shippedCfg PhatShipped + atomFlops shippedCfg PhatShipped
        + jointLsFlops shippedCfg PhatShipped
      = total shippedCfg PhatShipped := rfl

/-- The four joint-LS sub-terms (Gram, solve, `S'*yp`, the redundant final
call) account for `jointLsFlops` exactly. -/
theorem jointLs_terms_sum :
    gramFlops shippedCfg PhatShipped + solveFlops shippedCfg PhatShipped
        + styFlops shippedCfg PhatShipped + finalExtraFlops shippedCfg PhatShipped
      = jointLsFlops shippedCfg PhatShipped := rfl

/-! ### Closed-form sum identities, checked at the shipped `Phat`

`flops_omp.m` computes these sums via the textbook closed forms
(`sum_t = n(n+1)/2`, `sum_t2 = n(n+1)(2n+1)/6`, `sum_t3 = sum_t^2`) rather
than by recursion. A fully general proof that the recursive sums above
equal those closed forms needs induction, which the exemplar
(`DDRelaxFlops.lean`) does not use and which is avoided here per
`References/README.md`'s tactic restriction; instead each identity is
checked as a concrete, `decide`-able fact at `Phat = 9`, stated as a
cross-multiplication so no `Nat` division appears. -/

/-- `2 * Sum_{t=1}^{9} t = 9 * 10`. -/
theorem sumUpTo_closed_9 : 2 * sumUpTo 9 = 9 * 10 := by decide

/-- `6 * Sum_{t=1}^{9} t^2 = 9 * 10 * 19`. -/
theorem sumSqUpTo_closed_9 : 6 * sumSqUpTo 9 = 9 * 10 * 19 := by decide

/-- Nicomachus's identity, `Sum t^3 = (Sum t)^2`, checked at `n = 9`
(`45^2 = 2025`) rather than proved in general, for the same reason as
above. -/
theorem sumCubeUpTo_eq_sumUpTo_sq_9 :
    sumCubeUpTo 9 = sumUpTo 9 * sumUpTo 9 := by decide

/-! ### Structural claims -/

/-- The dominant term is exactly linear in `Phat` — same shape as
DD-RELAX's `proj_linear_in_Phat`, since `numProj` here is just `Phat`
itself (no `(1+NGS)` factor to also vary). -/
theorem proj_linear_in_Phat (c : Cfg) (a b : Nat) :
    projFlops c (a + b) = projFlops c a + projFlops c b := by
  simp [projFlops, numProj, Nat.add_mul, Nat.mul_add]

/-- Same linearity for the atom-generation term. -/
theorem atom_linear_in_Phat (c : Cfg) (a b : Nat) :
    atomFlops c (a + b) = atomFlops c a + atomFlops c b := by
  simp [atomFlops, numAtom, Nat.add_mul, Nat.mul_add]

/-- DOMINANCE, WITH A CAVEAT. At the shipped configuration the joint-LS
growth plus atom generation is under 2% of the total — nowhere near
overtaking `match_proj`. -/
theorem subdominant_under_two_percent :
    50 * (atomFlops shippedCfg PhatShipped + jointLsFlops shippedCfg PhatShipped)
      < total shippedCfg PhatShipped := by decide

/-- HONESTY CHECK: unlike DD-RELAX (whose subdominant terms are under 1% of
its total, see `DDRelaxFlops.subdominant_under_one_percent`), OMP's
subdominant terms here are NOT under 1% — they are about 1.003% of the
total, just over the line, because the growing joint-LS sum (`O(Phat^3)`
peak term) is a larger share of a smaller `Qtau*Qnu`-independent total than
DD-RELAX's flat `(NGS+1)`-multiplied version. This theorem records the
actual direction of the inequality so the file cannot silently overstate
dominance the way a copy-pasted "`< 1%`" claim would. -/
theorem subdominant_not_under_one_percent :
    total shippedCfg PhatShipped
      ≤ 100 * (atomFlops shippedCfg PhatShipped + jointLsFlops shippedCfg PhatShipped) := by
  decide

end CWS.OMP
