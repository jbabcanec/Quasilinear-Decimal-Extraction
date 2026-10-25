import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Exact splitting identities and the geometric BBP tail estimate

The splitting results are integer algorithms with rational specifications.
The BBP tail theorem keeps `HasSum bbpTerm Real.pi` as an explicit premise:
this file does not prove the BBP identity or any running-time bound.
-/

namespace BBP
namespace Series

/-- Unreduced integer data for a product and a sum of prefix products. -/
structure SplitData where
  P : ℤ
  Q : ℤ
  U : ℤ
  deriving DecidableEq, Repr

def merge (left right : SplitData) : SplitData :=
  ⟨left.P * right.P, left.Q * right.Q,
    left.U * right.Q + left.P * right.U⟩

def leaf (pq : ℤ × ℤ) : SplitData := ⟨pq.1, pq.2, pq.1⟩

def empty : SplitData := ⟨1, 1, 0⟩

def productValue (d : SplitData) : ℚ := (d.P : ℚ) / d.Q

def sumValue (d : SplitData) : ℚ := (d.U : ℚ) / d.Q

theorem merge_assoc (x y z : SplitData) :
    merge (merge x y) z = merge x (merge y z) := by
  cases x; cases y; cases z
  simp only [merge, SplitData.mk.injEq]
  constructor
  · ring
  constructor <;> ring

theorem merge_product (left right : SplitData) :
    productValue (merge left right) = productValue left * productValue right := by
  simp [productValue, merge, mul_div_mul_comm]

theorem merge_sum (left right : SplitData)
    (hl : left.Q ≠ 0) (hr : right.Q ≠ 0) :
    sumValue (merge left right) =
      sumValue left + productValue left * sumValue right := by
  have hl' : (left.Q : ℚ) ≠ 0 := by exact_mod_cast hl
  have hr' : (right.Q : ℚ) ≠ 0 := by exact_mod_cast hr
  simp only [sumValue, productValue, merge, Int.cast_add, Int.cast_mul]
  field_simp

/-- Product of all consecutive ratios in a segment. -/
def ratioProduct : List (ℤ × ℤ) → ℚ
  | [] => 1
  | pq :: rest => ((pq.1 : ℚ) / pq.2) * ratioProduct rest

/-- Sum of the nonempty prefix products of a segment. -/
def prefixSum : List (ℤ × ℤ) → ℚ
  | [] => 0
  | pq :: rest => ((pq.1 : ℚ) / pq.2) * (1 + prefixSum rest)

theorem ratioProduct_append (xs ys : List (ℤ × ℤ)) :
    ratioProduct (xs ++ ys) = ratioProduct xs * ratioProduct ys := by
  induction xs with
  | nil => simp [ratioProduct]
  | cons x xs ih => simp [ratioProduct, ih, mul_assoc]

theorem prefixSum_append (xs ys : List (ℤ × ℤ)) :
    prefixSum (xs ++ ys) = prefixSum xs + ratioProduct xs * prefixSum ys := by
  induction xs with
  | nil => simp [prefixSum, ratioProduct]
  | cons x xs ih =>
    simp only [List.cons_append, prefixSum, ratioProduct, ih]
    ring

/-- A splitting tree may use any parenthesization. -/
inductive SplitTree where
  | nil
  | singleton (pq : ℤ × ℤ)
  | node (left right : SplitTree)

def SplitTree.flatten : SplitTree → List (ℤ × ℤ)
  | .nil => []
  | .singleton pq => [pq]
  | .node left right => left.flatten ++ right.flatten

def SplitTree.eval : SplitTree → SplitData
  | .nil => empty
  | .singleton pq => leaf pq
  | .node left right => merge left.eval right.eval

/-- Every splitting tree computes the product and prefix-product sum of its
leaves exactly. Its denominator remains nonzero. -/
theorem splitting_correct (tree : SplitTree)
    (hq : ∀ pq ∈ tree.flatten, pq.2 ≠ 0) :
    tree.eval.Q ≠ 0 ∧
    productValue tree.eval = ratioProduct tree.flatten ∧
    sumValue tree.eval = prefixSum tree.flatten := by
  induction tree with
  | nil => simp [SplitTree.eval, SplitTree.flatten, empty, productValue, sumValue,
      ratioProduct, prefixSum]
  | singleton pq =>
    have hpq := hq pq (by simp [SplitTree.flatten])
    simpa [SplitTree.eval, SplitTree.flatten, leaf, productValue, sumValue,
      ratioProduct, prefixSum] using hpq
  | node left right ihl ihr =>
    have hl := ihl (fun pq hpq => hq pq (by simp [SplitTree.flatten, hpq]))
    have hr := ihr (fun pq hpq => hq pq (by simp [SplitTree.flatten, hpq]))
    refine ⟨?_, ?_, ?_⟩
    · exact mul_ne_zero hl.1 hr.1
    · simpa only [SplitTree.eval, SplitTree.flatten, ratioProduct_append,
        merge_product] using congrArg₂ (· * ·) hl.2.1 hr.2.1
    · simp only [SplitTree.eval, SplitTree.flatten, prefixSum_append,
        merge_sum _ _ hl.1 hr.1, hl.2.1, hl.2.2, hr.2.2]

/-- The integer factors for the real correction series: the initial factor
is `1/1`; later factors are `-a*r / (15*(A+a*r))`. -/
def hypergeometricFactor (A a r : ℕ) : ℤ × ℤ :=
  if r = 0 then (1, 1) else (-((a * r : ℕ) : ℤ), 15 * ((A + a * r : ℕ) : ℤ))

theorem hypergeometric_denominator_pos (A a r : ℕ) (hA : 0 < A) :
    0 < (hypergeometricFactor A a r).2 := by
  unfold hypergeometricFactor
  split_ifs
  · norm_num
  · simp only
    positivity

/-- Specialization to a finite truncation of the real correction series.
The premise specifies the leaves and permits any splitting-tree shape. -/
theorem hypergeometric_splitting_correct (A a R : ℕ) (hA : 0 < A)
    (tree : SplitTree)
    (hleaves : tree.flatten = (List.range R).map (hypergeometricFactor A a)) :
    tree.eval.Q ≠ 0 ∧
    productValue tree.eval = ratioProduct ((List.range R).map (hypergeometricFactor A a)) ∧
    sumValue tree.eval = prefixSum ((List.range R).map (hypergeometricFactor A a)) := by
  rw [← hleaves]
  apply splitting_correct
  intro pq hpq
  rw [hleaves] at hpq
  obtain ⟨r, _, rfl⟩ := List.mem_map.mp hpq
  exact ne_of_gt (hypergeometric_denominator_pos A a r hA)

/-- A noninitial ratio of the correction series, indexed from zero. -/
noncomputable def correctionRatio (A a : ℝ) (r : ℕ) : ℝ :=
  -(a * (r + 1)) / (15 * (A + a * (r + 1)))

/-- The finite product defining the coefficient `h_R`, including `h_0=1`. -/
noncomputable def correctionCoefficient (A a : ℝ) (R : ℕ) : ℝ :=
  ∏ r ∈ Finset.range R, correctionRatio A a r

theorem correctionRatio_abs (A a : ℝ) (r : ℕ) (hA : 0 < A) (ha : 0 ≤ a) :
    |correctionRatio A a r| = a * (r + 1) / (15 * (A + a * (r + 1))) := by
  unfold correctionRatio
  rw [abs_div, abs_neg, abs_of_nonneg (by positivity), abs_of_pos (by positivity)]

/-- The uniform strict contraction used to terminate the block's real
remainder search, including inputs outside the asymptotic regime. -/
theorem correctionRatio_lt (A a : ℝ) (r : ℕ) (hA : 0 < A) (ha : 0 ≤ a) :
    |correctionRatio A a r| < 1 / 15 := by
  rw [correctionRatio_abs A a r hA ha]
  apply (div_lt_iff₀ (by positivity)).mpr
  nlinarith

/-- The sharper finite coefficient bound used for the `ell/u` truncation
count. This proves the bound without invoking the infinite-series identity. -/
theorem correctionCoefficient_bound (A a : ℝ) (R : ℕ) (hA : 0 < A) (ha : 0 ≤ a) :
    |correctionCoefficient A a R| ≤ (a * R / (15 * A)) ^ R := by
  unfold correctionCoefficient
  rw [Finset.abs_prod]
  calc
    (∏ r ∈ Finset.range R, |correctionRatio A a r|) ≤
        ∏ _r ∈ Finset.range R, a * R / (15 * A) := by
      apply Finset.prod_le_prod
      · intro r _; exact abs_nonneg _
      · intro r hr
        have hr' : (r : ℝ) + 1 ≤ (R : ℝ) := by
          exact_mod_cast Nat.succ_le_of_lt (Finset.mem_range.mp hr)
        rw [correctionRatio_abs A a r hA ha]
        apply div_le_div₀ (by positivity) (mul_le_mul_of_nonneg_left hr' ha)
          (by positivity)
        nlinarith [show 0 ≤ a * ((r : ℝ) + 1) by positivity]
    _ = (a * R / (15 * A)) ^ R := by simp [div_pow]

/-- The unshifted term of the binary BBP series. -/
noncomputable def bbpTerm (k : ℕ) : ℝ :=
  (1 / 16 : ℝ) ^ k *
    (4 / (8 * (k : ℝ) + 1) - 2 / (8 * (k : ℝ) + 4) -
      1 / (8 * (k : ℝ) + 5) - 1 / (8 * (k : ℝ) + 6))

theorem bbpTerm_bound (k : ℕ) :
    ‖bbpTerm k‖ ≤ 8 * (1 / 16 : ℝ) ^ k := by
  have hk : 0 ≤ (k : ℝ) := Nat.cast_nonneg k
  have h1 : 0 ≤ 4 / (8 * (k : ℝ) + 1) := by positivity
  have h2 : 0 ≤ 2 / (8 * (k : ℝ) + 4) := by positivity
  have h3 : 0 ≤ 1 / (8 * (k : ℝ) + 5) := by positivity
  have h4 : 0 ≤ 1 / (8 * (k : ℝ) + 6) := by positivity
  have u1 : 4 / (8 * (k : ℝ) + 1) ≤ 4 := by
    apply (div_le_iff₀ (by positivity)).mpr; nlinarith
  have u2 : 2 / (8 * (k : ℝ) + 4) ≤ 2 := by
    apply (div_le_iff₀ (by positivity)).mpr; nlinarith
  have u3 : 1 / (8 * (k : ℝ) + 5) ≤ 1 := by
    apply (div_le_iff₀ (by positivity)).mpr; nlinarith
  have u4 : 1 / (8 * (k : ℝ) + 6) ≤ 1 := by
    apply (div_le_iff₀ (by positivity)).mpr; nlinarith
  have hb : |4 / (8 * (k : ℝ) + 1) - 2 / (8 * (k : ℝ) + 4) -
      1 / (8 * (k : ℝ) + 5) - 1 / (8 * (k : ℝ) + 6)| ≤ 8 := by
    apply abs_le.mpr
    constructor <;> linarith
  rw [bbpTerm, Real.norm_eq_abs, abs_mul, abs_of_nonneg (by positivity)]
  nlinarith [mul_le_mul_of_nonneg_left hb (show 0 ≤ (1 / 16 : ℝ) ^ k by positivity)]

/-- Geometric tail estimate for BBP, conditional only on the series identity.
The finite sum contains terms `0` through `K`, inclusive. -/
theorem bbp_tail (hbbp : HasSum bbpTerm Real.pi) (K : ℕ) :
    |Real.pi - ∑ k ∈ Finset.range (K + 1), bbpTerm k| ≤
      (128 / 15 : ℝ) * (1 / 16 : ℝ) ^ (K + 1) := by
  have h := norm_sub_le_of_geometric_bound_of_hasSum
    (by norm_num : (1 / 16 : ℝ) < 1) bbpTerm_bound hbbp (K + 1)
  rw [Real.norm_eq_abs, abs_sub_comm] at h
  convert h using 1
  ring

/-- The tail bound after multiplying by the decimal shift `10^n`. -/
theorem shifted_bbp_tail (hbbp : HasSum bbpTerm Real.pi) (n K : ℕ) :
    |(10 : ℝ) ^ n * Real.pi -
      (10 : ℝ) ^ n * ∑ k ∈ Finset.range (K + 1), bbpTerm k| ≤
      (128 / 15 : ℝ) * (10 : ℝ) ^ n * (1 / 16 : ℝ) ^ (K + 1) := by
  rw [← mul_sub, abs_mul, abs_of_nonneg (by positivity)]
  have h := mul_le_mul_of_nonneg_left (bbp_tail hbbp K)
    (show 0 ≤ (10 : ℝ) ^ n by positivity)
  nlinarith

/-- The integer cutoff from the manuscript, using `(G+3)/4` for `ceil(G/4)`. -/
def cutoff (n G : ℕ) : ℕ := n + (G + 3) / 4 + 2

/-- The cutoff makes the purely geometric upper bound smaller than `2^(-G)`.
There is no floating-point estimate in this inequality. -/
theorem cutoff_bound (n G : ℕ) :
    (128 / 15 : ℝ) * (10 : ℝ) ^ n * (1 / 16 : ℝ) ^ (cutoff n G + 1) <
      1 / (2 : ℝ) ^ G := by
  let q := (G + 3) / 4
  have hG : G ≤ 4 * q := by dsimp [q]; omega
  have hpow : (2 : ℝ) ^ G ≤ (16 : ℝ) ^ q := by
    calc
      (2 : ℝ) ^ G ≤ (2 : ℝ) ^ (4 * q) := pow_le_pow_right₀ (by norm_num) hG
      _ = (16 : ℝ) ^ q := by rw [pow_mul]; norm_num
  have hinv : (1 / 16 : ℝ) ^ q ≤ 1 / (2 : ℝ) ^ G := by
    simpa only [one_div_pow] using one_div_le_one_div_of_le (by positivity) hpow
  have hsmall : (10 : ℝ) ^ n * (1 / 16 : ℝ) ^ n ≤ 1 := by
    rw [← mul_pow]
    exact pow_le_one₀ (by norm_num) (by norm_num)
  have hcut : cutoff n G + 1 = n + q + 3 := by simp [cutoff, q]
  rw [hcut, pow_add, pow_add]
  have hq : 0 < (1 / 16 : ℝ) ^ q := by positivity
  calc
    (128 / 15 : ℝ) * 10 ^ n * (((1 / 16 : ℝ) ^ n * (1 / 16 : ℝ) ^ q) *
        (1 / 16 : ℝ) ^ 3) =
        (1 / 480 : ℝ) * (10 ^ n * (1 / 16 : ℝ) ^ n) * (1 / 16 : ℝ) ^ q := by ring
    _ ≤ (1 / 480 : ℝ) * (1 / 16 : ℝ) ^ q := by nlinarith
    _ < (1 / 16 : ℝ) ^ q := by nlinarith
    _ ≤ 1 / (2 : ℝ) ^ G := hinv

/-- Certified absolute error for the manuscript's integer cutoff, with the
BBP identity remaining an explicit mathematical dependency. -/
theorem shifted_bbp_cutoff (hbbp : HasSum bbpTerm Real.pi) (n G : ℕ) :
    |(10 : ℝ) ^ n * Real.pi -
      (10 : ℝ) ^ n * ∑ k ∈ Finset.range (cutoff n G + 1), bbpTerm k| <
      1 / (2 : ℝ) ^ G :=
  lt_of_le_of_lt (shifted_bbp_tail hbbp n (cutoff n G)) (cutoff_bound n G)

end Series
end BBP
