import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Analysis.SpecialFunctions.Log.Base
import Mathlib.Tactic

/-!
The final algebraic balance in equation (total) of the manuscript.
This file proves bounds on a specified real-valued expression. It does not
assert that an executable has this cost; that refinement remains separate.
-/

namespace BBP

noncomputable def balanceParameter (L : ℝ) : ℝ := (⌈Real.sqrt L⌉₊ : ℝ)

noncomputable def costExpression (N L u : ℝ) : ℝ :=
  N * (u * L ^ 2 + L ^ 3 / u + L * u ^ 2 + L ^ 2 * (Real.log L) ^ 2)

theorem balanceParameter_bounds {L : ℝ} (hL : 1 ≤ L) :
    Real.sqrt L ≤ balanceParameter L ∧ balanceParameter L ≤ 2 * Real.sqrt L := by
  have hs : 1 ≤ Real.sqrt L := Real.one_le_sqrt.mpr hL
  constructor
  · exact Nat.le_ceil _
  · have h := Nat.ceil_lt_add_one (Real.sqrt_nonneg L)
    dsimp [balanceParameter]
    linarith

theorem costExpression_bound {N L u : ℝ} (hN : 0 ≤ N) (hL : 1 ≤ L)
    (hu₁ : Real.sqrt L ≤ u) (hu₂ : u ≤ 2 * Real.sqrt L)
    (hlog : (Real.log L) ^ 2 ≤ Real.sqrt L) :
    costExpression N L u ≤ 8 * N * L ^ 2 * Real.sqrt L := by
  have hLp : 0 < L := lt_of_lt_of_le zero_lt_one hL
  have hs : 1 ≤ Real.sqrt L := Real.one_le_sqrt.mpr hL
  have hs0 : 0 ≤ Real.sqrt L := Real.sqrt_nonneg L
  have hsq : (Real.sqrt L) ^ 2 = L := Real.sq_sqrt hLp.le
  have hup : 0 < u := lt_of_lt_of_le (by linarith : 0 < Real.sqrt L) hu₁
  have hfirst : u * L ^ 2 ≤ 2 * (L ^ 2 * Real.sqrt L) := by
    nlinarith [mul_le_mul_of_nonneg_right hu₂ (sq_nonneg L)]
  have hsecond : L ^ 3 / u ≤ L ^ 2 * Real.sqrt L := by
    apply (div_le_iff₀ hup).mpr
    have h := mul_le_mul_of_nonneg_left hu₁ (show 0 ≤ L ^ 2 * Real.sqrt L by positivity)
    nlinarith [h]
  have husq : u ^ 2 ≤ 4 * L := by
    have h := mul_self_le_mul_self hup.le hu₂
    nlinarith
  have hthird : L * u ^ 2 ≤ 4 * (L ^ 2 * Real.sqrt L) := by
    have h := mul_le_mul_of_nonneg_left husq hLp.le
    have h' := mul_le_mul_of_nonneg_left hs (show 0 ≤ 4 * L ^ 2 by positivity)
    nlinarith
  have hfourth : L ^ 2 * (Real.log L) ^ 2 ≤ L ^ 2 * Real.sqrt L :=
    mul_le_mul_of_nonneg_left hlog (sq_nonneg L)
  have hsum : u * L ^ 2 + L ^ 3 / u + L * u ^ 2 + L ^ 2 * (Real.log L) ^ 2 ≤
      8 * (L ^ 2 * Real.sqrt L) := by linarith
  have h := mul_le_mul_of_nonneg_left hsum hN
  dsimp [costExpression]
  nlinarith [h]

open Filter Asymptotics

theorem eventually_log_sq_le_sqrt :
    ∀ᶠ L : ℝ in atTop, (Real.log L) ^ 2 ≤ Real.sqrt L := by
  have h := (isLittleO_log_rpow_rpow_atTop (2 : ℝ)
    (s := (1 / 2 : ℝ)) (by norm_num)).eventuallyLE
  filter_upwards [h] with L hL
  simpa only [Real.rpow_two, ← Real.sqrt_eq_rpow,
    Real.norm_of_nonneg (sq_nonneg (Real.log L)),
    Real.norm_of_nonneg (Real.sqrt_nonneg L)] using hL

theorem pow_two_mul_sqrt {L : ℝ} (hL : 0 < L) :
    L ^ 2 * Real.sqrt L = L ^ (5 / 2 : ℝ) := by
  rw [Real.sqrt_eq_rpow, ← Real.rpow_two,
    ← Real.rpow_add hL]
  norm_num

/-- This is an asymptotic theorem about the manuscript's cost expression,
not a machine-cost theorem about the C++ or Lean implementation. -/
theorem balanced_cost_isBigO {α : Type*} {l : Filter α} {N L : α → ℝ}
    (hN : ∀ᶠ x in l, 0 ≤ N x) (hL : Tendsto L l atTop) :
    (fun x => costExpression (N x) (L x) (balanceParameter (L x))) =O[l]
      (fun x => N x * (L x) ^ (5 / 2 : ℝ)) := by
  apply Asymptotics.IsBigO.of_bound 8
  filter_upwards [hN, hL.eventually (eventually_ge_atTop (1 : ℝ)),
    hL.eventually eventually_log_sq_le_sqrt] with x hn hl hlog
  obtain ⟨hu₁, hu₂⟩ := balanceParameter_bounds hl
  have hup : 0 < balanceParameter (L x) :=
    lt_of_lt_of_le (Real.sqrt_pos.mpr (lt_of_lt_of_le zero_lt_one hl)) hu₁
  have hcost : 0 ≤ costExpression (N x) (L x) (balanceParameter (L x)) := by
    dsimp [costExpression]
    positivity
  have htarget : 0 ≤ N x * (L x) ^ (5 / 2 : ℝ) := by positivity
  rw [Real.norm_of_nonneg hcost, Real.norm_of_nonneg htarget]
  have h := costExpression_bound hn hl hu₁ hu₂ hlog
  have heq := pow_two_mul_sqrt (lt_of_lt_of_le zero_lt_one hl)
  calc
    costExpression (N x) (L x) (balanceParameter (L x)) ≤
        8 * N x * L x ^ 2 * Real.sqrt (L x) := h
    _ = 8 * (N x * (L x) ^ (5 / 2 : ℝ)) := by rw [← heq]; ring

theorem logarithmic_cost_isBigO :
    (fun N : ℝ => costExpression N (Real.log N) (balanceParameter (Real.log N)))
      =O[atTop] (fun N : ℝ => N * (Real.log N) ^ (5 / 2 : ℝ)) :=
  balanced_cost_isBigO (eventually_ge_atTop 0) Real.tendsto_log_atTop

/-- The manuscript's discrete base-two size parameter, represented in reals
only for the asymptotic estimates. -/
noncomputable def bitLog (N : ℕ) : ℝ := max 1 (⌈Real.logb 2 (N : ℝ)⌉₊ : ℝ)

theorem bitLog_tendsto : Tendsto bitLog atTop atTop := by
  have hlog : Tendsto (fun N : ℕ => Real.logb 2 (N : ℝ)) atTop atTop :=
    (Real.tendsto_logb_atTop (by norm_num : (1 : ℝ) < 2)).comp
      tendsto_natCast_atTop_atTop
  apply tendsto_atTop_mono (f := fun N : ℕ => Real.logb 2 (N : ℝ)) _ hlog
  intro N
  exact (Nat.le_ceil _).trans (le_max_right _ _)

/-- The same bound for natural positions and the actual integer-rounded
logarithm and square-root parameter in the manuscript. -/
theorem discrete_balanced_cost_isBigO :
    (fun N : ℕ => costExpression (N : ℝ) (bitLog N) (balanceParameter (bitLog N)))
      =O[atTop] (fun N : ℕ => (N : ℝ) * (bitLog N) ^ (5 / 2 : ℝ)) :=
  balanced_cost_isBigO (Eventually.of_forall (fun N : ℕ => Nat.cast_nonneg N)) bitLog_tendsto

end BBP

