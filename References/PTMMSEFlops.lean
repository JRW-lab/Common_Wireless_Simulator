/-
  PTMMSEFlops.lean - machine-checkable complexity derivation for the
  partitioned time-domain MMSE receiver (PT-MMSE).

  Mirrors:   Comm Functions/TX RX Functions/flops_ptmmse.m
  Algorithm: Comm Functions/TX RX Functions/equalizer_ptmmse.m,
             which implements Section V / Algorithm 2 of the project's
             own DD-RELAX paper (Eq. 48-55).
  Policy:    References/README.md

  WHAT THIS FILE VERIFIES
    * that `total` is exactly the sum of the six stated terms;
    * the evaluated numbers at the shipped configuration, matching the
      .m file digit for digit;
    * that the Gram term alone exceeds 77% of one call;
    * that the O^3/3 in the LU term divides exactly at this configuration,
      so the Nat division below truncates nothing.

  WHAT THIS FILE DOES NOT VERIFY
    * THAT THE OPERATION COUNTS MATCH equalizer_ptmmse.m. Lean cannot see
      the MATLAB. Each count cites the lines it was read from; that
      reading is a human audit step and is where an error would hide.
    * anything about numerical accuracy - only operation counts.
    * that FLOPs predict runtime. They demonstrably do not here: the Gram
      is a dense BLAS3 product that runs near peak, while the estimator's
      grid correlations are BLAS2. A term being larger in FLOPs does NOT
      make it larger in seconds. See README.md's timing caveat.

  SCOPE. |J(m)| is taken at its Corollary-1 bound of 2L-1 blocks. With no
  guard (DD-RELAX-SI) every candidate block carries data, so the bound is
  EXACT; with a guard the true count is lower near the guard and this is
  an upper bound, matching the worst-case convention used elsewhere.
  Constellation slicing inside decide_block and the LDPC decoder are
  excluded - the decoder is not part of the receiver.

  Lean 4, no Mathlib (see References/README.md).
-/

namespace CWS.PTMMSE

/-- Configuration parameters, named as in oddm_config.m. -/
structure Cfg where
  /-- subcarriers / delay blocks per frame. -/
  M       : Nat
  /-- Doppler bins per frame, the size of one delay block. -/
  N       : Nat
  /-- delay-support width, L = L1 + L2 + 1. -/
  L       : Nat
  /-- data delay blocks. M with no guard; 52 for the 1L guard and 41 for
      the 2L guard at this project's geometry. -/
  Md      : Nat
  /-- candidate blocks coupling into one observation, 2L-1 by Corollary 1. -/
  Jblocks : Nat
  /-- detection iterations (the profile's N_iters). -/
  NI      : Nat
  /-- flops per complex multiply-accumulate. -/
  cCMAC   : Nat
  /-- flops for a squared magnitude plus an accumulate. -/
  cABS2   : Nat

/-- Observation rows per block, |O(m)| = L*N (Eq. 49). -/
def O (c : Cfg) : Nat := c.L * c.N

/-- Columns in J(m) (Eq. 50). -/
def J (c : Cfg) : Nat := c.Jblocks * c.N

/-! ### Terms -/

/-- Delay-time transform, Pm = kron(I_M, F_N'), applied on each side of H
    plus once to y. equalizer_ptmmse.m lines 226-228. -/
def transform (c : Cfg) : Nat :=
  2 * (c.M * c.N * c.N * c.M * c.N) * c.cCMAC + (c.M * c.N * c.N) * c.cCMAC

/-- Gram per data block, Gm = Hm*(colVar.*Hm'): an (O x J)(J x O) product
    plus the per-column scaling. Line 329. -/
def gram (c : Cfg) : Nat :=
  c.Md * (O c * J c * O c) * c.cCMAC + c.Md * (J c * O c) * c.cCMAC

/-- LU of the O x O system in W1 = Es*(Hown'/(Gm + N0 I)), line 333.
    Nat division; see lu_exact below. -/
def solve (c : Cfg) : Nat :=
  c.Md * (O c * O c * O c / 3 + c.N * O c * O c) * c.cCMAC

/-- h2 = sum(|Hown|^2,1), line 330. -/
def hnorm (c : Cfg) : Nat := c.Md * (O c * c.N) * c.cABS2

/-- Iteration 1: the LMMSE apply, the decision transform, the soft-output
    transform. Lines 358-368. -/
def iter1 (c : Cfg) : Nat :=
  c.Md * (c.N * O c + 2 * c.N * c.N) * c.cCMAC

/-- Iterations 2..NI: cancellation Hm*xJ, matched filter Hown'*z, and the
    two F_N transforms. Lines 375-397. -/
def itersN (c : Cfg) : Nat :=
  c.Md * (c.NI - 1) * (O c * J c + c.N * O c + 2 * c.N * c.N) * c.cCMAC

/-- Cost of ONE equalizer call. -/
def total (c : Cfg) : Nat :=
  transform c + gram c + solve c + hnorm c + iter1 c + itersN c

theorem total_is_sum (c : Cfg) :
    total c = transform c + gram c + solve c + hnorm c + iter1 c + itersN c := rfl

/-- Per frame, for a scheme that calls the receiver `calls` times. -/
def perFrame (c : Cfg) (calls : Nat) : Nat := calls * total c

/-! ### The shipped configuration

M=64, N=16, Q=4 give L1=4, L2=7, L=12, so O = 192 and J = 23*16 = 368.
Md = 64 for the guard-free superimposed frame; N_iters = 8. -/

def shipped : Cfg :=
  { M := 64, N := 16, L := 12, Md := 64, Jblocks := 23, NI := 8,
    cCMAC := 8, cABS2 := 4 }

/-- The 1L-guard frame reserves L = 12 delay blocks. -/
def shippedL : Cfg := { shipped with Md := 52 }
/-- The 2L-guard frame reserves 2L-1 = 23 delay blocks. -/
def shipped2L : Cfg := { shipped with Md := 41 }

/-! ### Concrete evaluations - these must match flops_ptmmse.m exactly -/

theorem O_shipped : O shipped = 192 := rfl
theorem J_shipped : J shipped = 368 := rfl

theorem transform_shipped : transform shipped = 268566528 := rfl
theorem gram_shipped      : gram shipped      = 6981943296 := rfl
theorem solve_shipped     : solve shipped     = 1509949440 := rfl
theorem hnorm_shipped     : hnorm shipped     = 786432 := rfl
theorem iter1_shipped     : iter1 shipped     = 1835008 := rfl
theorem itersN_shipped    : itersN shipped    = 266076160 := rfl

theorem total_shipped : total shipped = 9029156864 := rfl

/-- DD-RELAX-SI calls the receiver once per SIC round and once per turbo
    round: 8 at the measured Eb/N0 = 6 dB round count, 17 worst case. -/
theorem perFrame_si_expected : perFrame shipped 8 = 72233254912 := rfl
theorem perFrame_si_worst    : perFrame shipped 17 = 153495666688 := rfl

/-- The guard-isolated baselines call it once. -/
theorem perFrame_1L : perFrame shippedL 1 = 7386546176 := rfl
theorem perFrame_2L : perFrame shipped2L 1 = 5880819712 := rfl

/-! ### Structural results -/

/-- The LU term's division is EXACT at this configuration, so the Nat
    division in `solve` truncates nothing. -/
theorem lu_exact : 3 * (O shipped * O shipped * O shipped / 3)
    = O shipped * O shipped * O shipped := by decide

/-- The Gram alone is over 77% of one call. The receiver's cost is
    essentially the cost of forming one (L*N) x (L*N) Gram per data
    block. -/
theorem gram_dominant :
    gram shipped * 100 > total shipped * 77 := by decide

/-- Everything other than the Gram and the LU solve is under 7%. -/
theorem rest_small :
    (transform shipped + hnorm shipped + iter1 shipped + itersN shipped) * 100
      < total shipped * 7 := by decide

/-- Cost is exactly linear in the number of calls, which is why calling the
    receiver once per SIC round dominates DD-RELAX-SI's frame budget. -/
theorem perFrame_linear (c : Cfg) (a b : Nat) :
    perFrame c (a + b) = perFrame c a + perFrame c b := by
  simp [perFrame, Nat.add_mul]

end CWS.PTMMSE
