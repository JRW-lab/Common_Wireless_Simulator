/-
  SAGEFlops.lean — machine-checkable complexity derivation for the SAGE
  channel estimator, as implemented in this repository.

  Mirrors:  Comm Functions/ODDM Functions/DD-RELAX-paper/flops_sage.m
  Algorithm: Comm Functions/ODDM Functions/DD-RELAX-paper/est_sage.m
  Policy:    References/README.md
  Sibling:   References/DDRelaxFlops.lean (same shape, different constants —
             read that file's header too; this one only restates what
             differs for SAGE).

  WHAT THIS FILE VERIFIES
    * that `total` is exactly the sum of the stated terms;
    * that the evaluated numbers at Figure 7's configuration are correct
      (these are the same numbers flops_sage.m's own smoke test prints,
      obtained here by actually running it, not by hand arithmetic — see
      the concrete-checks section below for the command used);
    * that the dominant (`match_proj`) and atom-generation terms both scale
      exactly as `Phat * (1 + NGS)` — the SAME pure shape as DD-RELAX's
      dominant term, and, unlike DD-RELAX, that these two call counts are
      literally EQUAL to each other (SAGE has no extra atom re-evaluation);
    * that the subdominant term is under 1% of the total at that config.

  WHAT THIS FILE DOES NOT VERIFY
    * THAT THE OPERATION COUNTS MATCH est_sage.m. Lean cannot see the
      MATLAB. Each count below cites the source lines it was read from;
      that reading is a human audit step and is where an error would hide.
    * that the MATLAB implements the published SAGE algorithm (Fessler &
      Hero 1994) correctly.
    * anything about numerical accuracy — only operation counts.
    * the model against measured runtime. flops_sage.m's own header records
      a comparison against a desktop NGS-ceiling sweep (tic/toc wall clock,
      four collection workers saturating the machine simultaneously) with
      deviations up to -10.2% at NGS=32 — too noisy to adjudicate a pure
      (1+NGS) model to that precision, and explicitly NOT to be read as a
      confirmed discrepancy requiring a new term. No clean cputime
      validation exists yet for SAGE, same as for DD-RELAX.

  TWO THINGS ESTABLISHED BY A PRIOR AUDIT (both already documented in
  flops_sage.m's own header, restated here rather than re-derived):
    * SAGE HAS NO joint-LS / Gram-matrix / linear-solve term. An earlier
      version of est_sage.m had a periodic full joint re-solve; a
      2026-09-14 correction removed it as not part of literal SAGE (Fessler
      & Hero's own pseudocode updates one parameter at a time with every
      other one held exactly fixed — there is no periodic joint step at
      all). So unlike DD-RELAX there is no `O(Phat^3)` solve and no
      `O(Phat^2 * L*N)` Gram term here — do not add one just because
      DD-RELAX has one.
    * Measured `Phat` is about 8.99 of `Pmax = 9` at this project's configs
      (flops_sage.m header, 2026-09-24), so the worst-case (`Phat = Pmax`)
      and expected-value breakdowns nearly coincide. Per References/README.md
      §"Naturals only", theorems below are stated only at the integer
      `Phat = Pmax = 9` (worst case); the fractional 8.99 expected value is
      not encoded as a `rfl` theorem, only noted here for context.

  Lean 4, no Mathlib (see References/README.md).
-/

namespace CWS.SAGE

/-- Configuration parameters, named as in `oddm_config.m`. Identical field
set to `CWS.DDRelax.Cfg` — SAGE needs no extra fields because it has no
Gram/solve step to parameterize. -/
structure Cfg where
  /-- delay-domain support of the pilot window, `L = L1 + L2 + 1`. -/
  L     : Nat
  /-- Doppler bins per frame. -/
  N     : Nat
  /-- delay grid points, `ceil(tau_max/Dtau) + 1`. -/
  Qtau  : Nat
  /-- Doppler grid points, `2*ceil(nu_max/Dnu) + 1`. -/
  Qnu   : Nat
  /-- Gauss-Seidel-style cyclic refinement sweeps. -/
  NGS   : Nat
  /-- flops per complex multiply-accumulate (6 mul + 2 add). -/
  cCMAC : Nat
  /-- calibrated flops per grid point of one `s_atom` evaluation. -/
  cAtom : Nat

/-- Pilot-window length, the inner dimension of every matvec below. -/
def LN (c : Cfg) : Nat := c.L * c.N

/-! ### Call counts

`est_sage.m` has the same two-phase shape as `dd_relax.m`. Phase 1
(lines 57-65) is the greedy Matching-Pursuit initialisation: `for t = 1:Pmax`
calling `match_proj` once per iteration (line 58) and `Omega` (`s_atom`)
once per iteration (line 61), with an early break that sets `Phat = t`; if
it never breaks, line 66 sets `Phat = Pmax`. So Phase 1 makes exactly `Phat`
calls to each, NOT `Pmax` — same reasoning as DD-RELAX's Phase 1. Phase 2
(lines 73-82) is `for cyc = 1:NGS { for t = 1:Phat }`, each inner iteration
calling `match_proj` once (line 76) and `Omega` once (line 79) — ONE atom
evaluation per inner iteration, not two: est_sage.m has no `refine_atom`-
then-reassignment pattern and no parabolic sub-grid refinement anywhere
(the file's own docstring, lines 15-16: "on-grid selection only ... this is
what keeps SAGE coarser than DD-RELAX"). -/

/-- `match_proj` calls per frame. -/
def numProj (c : Cfg) (Phat : Nat) : Nat := Phat * (1 + c.NGS)

/-- `Omega` / `s_atom` calls per frame: exactly one per path per phase-visit
(Phase 1: one per path; Phase 2: one per path per sweep) — literally the
SAME expression as `numProj`, unlike DD-RELAX where the atom count
(`Phat * (2 + 3*NGS)`) exceeds the projection count. -/
def numAtom (c : Cfg) (Phat : Nat) : Nat := Phat * (1 + c.NGS)

/-! ### Cost terms -/

/-- DOMINANT. Each `match_proj` call computes `sd2' * r` (or, in Phase 2,
against the freshly-recomputed `yt`), a dense `(Qtau*Qnu) x LN` by `LN x 1`
complex matvec against the precomputed dictionary (`match_proj.m` line 18,
same routine DD-RELAX calls). The dictionary itself is built once per
configuration in `build_ctx.m`, not per frame, and is correctly excluded. -/
def projFlops (c : Cfg) (Phat : Nat) : Nat :=
  numProj c Phat * (c.Qtau * c.Qnu) * LN c * c.cCMAC

/-- Subdominant, approximate: each `s_atom` evaluation is `O(L*N)` via
`apg_interp`'s interpolation, same calibrated-constant convention as
DD-RELAX's `atomFlops` (`cAtom` is a library detail, not derived exactly). -/
def atomFlops (c : Cfg) (Phat : Nat) : Nat :=
  numAtom c Phat * LN c * c.cAtom

/-- Subdominant. The cyclic phase recomputes the FULL residual
`yt = yp - S*phi_c(:) + S(:,t)*phi_c(t)` fresh for every `t`
(`est_sage.m` line 75) — `S*phi_c(:)` costs `Phat*LN` complex
multiply-accumulates, recomputed `Phat` times per sweep over `NGS` sweeps
(the `NGS`-cycle loop at lines 73-82; there is no extra post-loop call, so
this is `NGS` times, not `NGS + 1`). Costed the same way as DD-RELAX's own
Gauss-Seidel residual recompute. NO joint-LS / Gram / solve term exists
here — see the header note above: the 2026-09-14 correction removed the
periodic joint re-solve as not part of literal SAGE, so unlike DD-RELAX
there is no `O(Phat^3)` and no `O(Phat^2 * L*N)` Gram term to add. -/
def residFlops (c : Cfg) (Phat : Nat) : Nat :=
  c.NGS * (Phat * Phat) * LN c * c.cCMAC

/-- Total per-frame estimation cost. -/
def total (c : Cfg) (Phat : Nat) : Nat :=
  projFlops c Phat + atomFlops c Phat + residFlops c Phat

/-! ### The shipped configuration

Figure 7 / Profile 5 config, same as `CWS.DDRelax.fig7`: `N = 16, M = 64,
v = 500 km/h, alpha = 0.3, Q = 4, PiTau = 32, PiNu = 16, Pmax = 9, NGS = 4`.
The derived grid sizes were obtained by actually running

  matlab -batch "cd('C:/MATLAB Projects/Common Wireless Simulator'); ^
    addpath(genpath('Comm Functions')); ^
    P = oddm_config('N',16,'M',64,'fc',4e9,'sub',15000,'alpha',0.3,'Q',4, ^
      'v_kmh',500,'Pmax',9,'NGS',4,'PiTau',32,'PiNu',16, ^
      'frame_layout','guard_end'); ^
    out = flops_sage(P,'Phat',8.99); disp(out.worst_case); disp(out.params);"

which printed `L = 12, Qtau = 79, Qnu = 65, Nbins = 16` — identical to
DD-RELAX's grid at the same config (as expected: both read the same
`oddm_config.m` fields, and SAGE adds no grid parameters of its own) — and
`worst_case = {proj: 354931200, atom_gen: 259200, resid_recompute: 497664,
total: 355688064}` at `Phat = Pmax = 9`, which is exactly what the concrete
theorems below fix. -/
def sageFig7 : Cfg :=
  { L := 12, N := 16, Qtau := 79, Qnu := 65, NGS := 4, cCMAC := 8, cAtom := 30 }

/-- `Pmax = 9`; measured `Phat` is ~8.99 of 9, so this is both the worst
case and (per the header note) essentially the expected case too. -/
def PmaxFig7 : Nat := 9

/-! ### Concrete checks

These fix the numbers `flops_sage.m`'s `worst_case` breakdown printed for
the run above. If a term in the MATLAB is edited without updating the
count here, these stop being `rfl`. -/

example : LN sageFig7 = 192 := rfl
example : numProj sageFig7 PmaxFig7 = 45 := rfl
example : numAtom sageFig7 PmaxFig7 = 45 := rfl

theorem proj_fig7  : projFlops  sageFig7 PmaxFig7 = 354931200 := rfl
theorem atom_fig7  : atomFlops  sageFig7 PmaxFig7 = 259200    := rfl
theorem resid_fig7 : residFlops sageFig7 PmaxFig7 = 497664    := rfl
theorem total_fig7 : total      sageFig7 PmaxFig7 = 355688064 := rfl

/-- The three terms account for the total exactly — no hidden remainder. -/
theorem terms_sum_to_total :
    projFlops sageFig7 PmaxFig7 + atomFlops sageFig7 PmaxFig7 + residFlops sageFig7 PmaxFig7
      = total sageFig7 PmaxFig7 := rfl

/-! ### Structural claims -/

/-- SAGE-SPECIFIC STRUCTURAL FACT: unlike DD-RELAX (where the atom count
strictly exceeds the projection count), SAGE calls `match_proj` and
`Omega`/`s_atom` EXACTLY the same number of times per frame — both phases
pair one atom evaluation with one projection, with no separate refine step.
This is definitional (both are literally `Phat * (1 + NGS)`), but stating
it as a theorem makes the equality machine-checked rather than merely
"the same by inspection of two definitions". -/
theorem numAtom_eq_numProj (c : Cfg) (Phat : Nat) :
    numAtom c Phat = numProj c Phat := rfl

/-- DOMINANCE. The atom-generation and residual-recompute terms together
are under 1% of the total, which is what licenses quoting the leading-order
(`match_proj`) term alone. -/
theorem subdominant_under_one_percent :
    100 * (atomFlops sageFig7 PmaxFig7 + residFlops sageFig7 PmaxFig7)
      < total sageFig7 PmaxFig7 := by decide

/-- The dominant term is exactly linear in `Phat`. -/
theorem proj_linear_in_Phat (c : Cfg) (a b : Nat) :
    projFlops c (a + b) = projFlops c a + projFlops c b := by
  simp [projFlops, numProj, Nat.add_mul, Nat.mul_add]

/-- SCALING IN NGS. At fixed `Phat`, the dominant (`match_proj`) term is
proportional to `1 + NGS` — the SAME pure shape as DD-RELAX's
`proj_ngs_scaling`. Stated as a cross-multiplication to stay in `Nat` and
avoid truncating division. -/
theorem proj_ngs_scaling (c : Cfg) (Phat n m : Nat) :
    projFlops { c with NGS := n } Phat * (1 + m)
      = projFlops { c with NGS := m } Phat * (1 + n) := by
  simp [projFlops, numProj, LN, Nat.mul_comm, Nat.mul_left_comm, Nat.mul_assoc]

/-- SCALING IN NGS, atom term. Because `numAtom = numProj` (see
`numAtom_eq_numProj`), the atom-generation term has EXACTLY the same
`(1 + NGS)`-proportional shape as the dominant term — a second instance of
the claim above, not merely an analogous one. -/
theorem atom_ngs_scaling (c : Cfg) (Phat n m : Nat) :
    atomFlops { c with NGS := n } Phat * (1 + m)
      = atomFlops { c with NGS := m } Phat * (1 + n) := by
  simp [atomFlops, numAtom, LN, Nat.mul_comm, Nat.mul_left_comm, Nat.mul_assoc]

/-- Concrete instance of `proj_ngs_scaling`: going from `NGS = 4` to
`NGS = 64` multiplies the dominant term by exactly 13, same multiplier as
DD-RELAX's `proj_ngs_4_to_64` (the multiplier depends only on `(1+NGS)`
ratios, which are identical between the two estimators). -/
theorem proj_ngs_4_to_64 (Phat : Nat) :
    projFlops { sageFig7 with NGS := 64 } Phat
      = 13 * projFlops { sageFig7 with NGS := 4 } Phat := by
  simp [projFlops, numProj, LN, Nat.mul_comm, Nat.mul_left_comm, Nat.mul_assoc]

end CWS.SAGE
