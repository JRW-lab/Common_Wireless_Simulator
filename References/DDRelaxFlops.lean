/-
  DDRelaxFlops.lean — machine-checkable complexity derivation for DD-RELAX
  channel estimation, as implemented in this repository.

  Mirrors:  Comm Functions/ODDM Functions/DD-RELAX-paper/flops_dd_relax.m
  Algorithm: Comm Functions/ODDM Functions/DD-RELAX-paper/dd_relax.m
  Policy:    References/README.md

  WHAT THIS FILE VERIFIES
    * that `total` is exactly the sum of the stated terms;
    * that the evaluated numbers at Figure 7's configuration are correct
      (these are the same numbers the .m file's smoke test prints);
    * that the dominant term scales exactly as `Phat * (1 + NGS)`;
    * that the subdominant terms are under 1% of the total at that config.

  WHAT THIS FILE DOES NOT VERIFY
    * THAT THE OPERATION COUNTS MATCH dd_relax.m. Lean cannot see the MATLAB.
      Each count below cites the source lines it was read from; that reading
      is a human audit step and is where an error would hide.
    * anything about numerical accuracy — only operation counts.
    * the model against measured runtime. The only timing data available is
      tic/toc wall clock taken while four collection workers saturated the
      machine; it is consistent with the (1+NGS) prediction to within its own
      noise, which is not the same as validation. A clean cputime measurement
      is outstanding.

  Lean 4, no Mathlib (see References/README.md).
-/

namespace CWS.DDRelax

/-- Configuration parameters, named as in `oddm_config.m`. -/
structure Cfg where
  /-- delay-domain support of the pilot window, `L = L1 + L2 + 1`. -/
  L     : Nat
  /-- Doppler bins per frame. -/
  N     : Nat
  /-- delay grid points, `ceil(tau_max/Dtau) + 1`. -/
  Qtau  : Nat
  /-- Doppler grid points, `2*ceil(nu_max/Dnu) + 1`. -/
  Qnu   : Nat
  /-- Gauss-Seidel sweeps. -/
  NGS   : Nat
  /-- flops per complex multiply-accumulate (6 mul + 2 add). -/
  cCMAC : Nat
  /-- calibrated flops per grid point of one `s_atom` evaluation. -/
  cAtom : Nat

/-- Pilot-window length, the inner dimension of every matvec below. -/
def LN (c : Cfg) : Nat := c.L * c.N

/-! ### Call counts

`dd_relax.m` has two phases. Phase 1 (lines 106-124) is the matching-pursuit
initialisation: `for t = 1:Pmax` calling `match_proj` once per iteration, with
an early break that sets `Phat = t`; if it never breaks, line 125 sets
`Phat = Pmax`. So Phase 1 makes exactly `Phat` calls in either case — NOT
`Pmax`. Phase 2 (lines 131-141) is `for cyc = 1:NGS { for t = 1:Phat }`, each
inner iteration calling `relax_component`, which calls `match_proj` once. -/

/-- `match_proj` calls per frame. -/
def numProj (c : Cfg) (Phat : Nat) : Nat := Phat * (1 + c.NGS)

/-- `s_atom` calls per frame: one per path in Phase 1, three per path per
sweep in Phase 2 (refine-then-reassign), plus one final per path. -/
def numAtom (c : Cfg) (Phat : Nat) : Nat := Phat * (2 + 3 * c.NGS)

/-! ### Cost terms -/

/-- DOMINANT. Each `match_proj` call computes `sd2' * r`, a dense
`(Qtau*Qnu) x LN` by `LN x 1` complex matvec against the precomputed
dictionary (`match_proj.m` line 18). The dictionary itself is built once per
configuration in `build_ctx.m`, not per frame, and is correctly excluded. -/
def projFlops (c : Cfg) (Phat : Nat) : Nat :=
  numProj c Phat * (c.Qtau * c.Qnu) * LN c * c.cCMAC

/-- Subdominant, approximate: `s_atom` evaluation is `O(L*N)` per call via
`apg_interp`'s three-shift Akima interpolation. `cAtom` is calibrated rather
than derived, because the interpolation cost is a library detail, not an
algorithmic one. -/
def atomFlops (c : Cfg) (Phat : Nat) : Nat :=
  numAtom c Phat * LN c * c.cAtom

/-- `joint_ls` runs `NGS + 1` times: Gram `S'*S` (dense, not exploiting
Hermitian symmetry, matching what the code actually computes), `S'*yp`, and
an `O(Phat^3)` regularised solve. Separately, the Phase 2 inner loop
recomputes the residual `S * phi_c(1:Phat)` once per path per sweep — `NGS`
sweeps only, since the final `joint_ls` call does not re-enter the loop. -/
def gramFlops  (c : Cfg) (Phat : Nat) : Nat := (c.NGS + 1) * (Phat * Phat) * LN c * c.cCMAC
def solveFlops (c : Cfg) (Phat : Nat) : Nat := (c.NGS + 1) * (Phat * Phat * Phat) * c.cCMAC
def styFlops   (c : Cfg) (Phat : Nat) : Nat := (c.NGS + 1) * Phat * LN c * c.cCMAC
def residFlops (c : Cfg) (Phat : Nat) : Nat := c.NGS * (Phat * Phat) * LN c * c.cCMAC

def jointLsFlops (c : Cfg) (Phat : Nat) : Nat :=
  gramFlops c Phat + solveFlops c Phat + styFlops c Phat + residFlops c Phat

/-- Total per-frame estimation cost. -/
def total (c : Cfg) (Phat : Nat) : Nat :=
  projFlops c Phat + atomFlops c Phat + jointLsFlops c Phat

/-! ### The shipped configuration

Figure 7 / Profile 5: `N = 16, M = 64, v = 500 km/h, alpha = 0.3, Q = 4,
PiTau = 32, PiNu = 16, Pmax = 9, NGS = 4`. The derived grid sizes
`Qtau = 79`, `Qnu = 65` and `L = 12` come from `oddm_config.m` at those
parameters and match values independently recorded elsewhere in the project. -/
def fig7 : Cfg :=
  { L := 12, N := 16, Qtau := 79, Qnu := 65, NGS := 4, cCMAC := 8, cAtom := 30 }

/-- `Pmax = 9`, and DD-RELAX's early stop rarely fires, so this is both the
worst case and close to the expected case. -/
def PmaxFig7 : Nat := 9

/-! ### Concrete checks

These fix the numbers the `.m` file's smoke test prints. If a term in the
MATLAB is edited without updating the count here, these stop being `rfl`. -/

example : LN fig7 = 192 := rfl
example : numProj fig7 PmaxFig7 = 45 := rfl
example : numAtom fig7 PmaxFig7 = 126 := rfl

theorem proj_fig7    : projFlops    fig7 PmaxFig7 = 354931200 := rfl
theorem atom_fig7    : atomFlops    fig7 PmaxFig7 = 725760    := rfl
theorem jointLs_fig7 : jointLsFlops fig7 PmaxFig7 = 1218024   := rfl
theorem total_fig7   : total        fig7 PmaxFig7 = 356874984 := rfl

/-- The three terms account for the total exactly — no hidden remainder. -/
theorem terms_sum_to_total :
    projFlops fig7 PmaxFig7 + atomFlops fig7 PmaxFig7 + jointLsFlops fig7 PmaxFig7
      = total fig7 PmaxFig7 := rfl

/-! ### Structural claims -/

/-- DOMINANCE. Everything other than the grid correlation is under 1% of the
total, which is what licenses quoting the leading-order term alone. -/
theorem subdominant_under_one_percent :
    100 * (atomFlops fig7 PmaxFig7 + jointLsFlops fig7 PmaxFig7)
      < total fig7 PmaxFig7 := by decide

/-- The dominant term is exactly linear in `Phat`. -/
theorem proj_linear_in_Phat (c : Cfg) (a b : Nat) :
    projFlops c (a + b) = projFlops c a + projFlops c b := by
  simp [projFlops, numProj, Nat.add_mul, Nat.mul_add]

/-- SCALING IN NGS. At fixed `Phat`, the dominant term is proportional to
`1 + NGS`. Stated as a cross-multiplication to stay in `Nat` and avoid
truncating division. This is the claim that was compared against measured
runtime; the measurement was too noisy to confirm it, but the algebra is
exact. -/
theorem proj_ngs_scaling (c : Cfg) (Phat n m : Nat) :
    projFlops { c with NGS := n } Phat * (1 + m)
      = projFlops { c with NGS := m } Phat * (1 + n) := by
  simp [projFlops, numProj, LN, Nat.mul_comm, Nat.mul_left_comm, Nat.mul_assoc]

/-- Concrete instance of the above: going from `NGS = 4` to `NGS = 64`
multiplies the dominant term by exactly 13. -/
theorem proj_ngs_4_to_64 (Phat : Nat) :
    projFlops { fig7 with NGS := 64 } Phat
      = 13 * projFlops { fig7 with NGS := 4 } Phat := by
  simp [projFlops, numProj, LN, Nat.mul_comm, Nat.mul_left_comm, Nat.mul_assoc]

end CWS.DDRelax
