import Mathlib.RingTheory.Coprime.Lemmas
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Data.Rat.Cast.CharZero
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination

/-!
# Exact arithmetic for weighted-block reconstruction

These theorems prove the rational merge operations and the CRT integrality
step. The local congruences are explicit hypotheses: no claim is made here
that the manuscript's prime-power coefficient algorithm computes them.
-/

namespace BBP

/-- The uncancelled merge used in the prime-fraction tree. -/
theorem rational_merge (b₁ d₁ b₂ d₂ : ℚ) (h₁ : d₁ ≠ 0) (h₂ : d₂ ≠ 0) :
    (b₁ * d₂ + b₂ * d₁) / (d₁ * d₂) = b₁ / d₁ + b₂ / d₂ := by
  field_simp

/-- The merge used for a series of successive products. -/
theorem hypergeometric_merge (p₁ q₁ u₁ p₂ q₂ u₂ : ℚ)
    (h₁ : q₁ ≠ 0) (h₂ : q₂ ≠ 0) :
    (p₁ * p₂ / (q₁ * q₂) = (p₁ / q₁) * (p₂ / q₂)) ∧
    ((u₁ * q₂ + p₁ * u₂) / (q₁ * q₂) = u₁ / q₁ + (p₁ / q₁) * (u₂ / q₂)) := by
  constructor <;> field_simp

/-- An exact weighted merge, including a binary shift on the left child. -/
theorem shifted_rational_merge (b₁ d₁ b₂ d₂ : ℚ) (k : ℕ)
    (h₁ : d₁ ≠ 0) (h₂ : d₂ ≠ 0) :
    ((2 : ℚ)^k * b₁ * d₂ + b₂ * d₁) / (d₁ * d₂) =
      (2 : ℚ)^k * (b₁ / d₁) + b₂ / d₂ := by
  field_simp

section CRT

variable {ι : Type*} [DecidableEq ι]

/-- Product of all moduli except the specified modulus. -/
def cofactor (s : Finset ι) (q : ι → ℤ) (i : ι) : ℤ :=
  ∏ j ∈ s.erase i, q j

/-- The numerator obtained by adding all prime-power fractions over the
uncancelled product denominator. -/
def reconstructionNumerator (s : Finset ι) (q c : ι → ℤ) : ℤ :=
  ∑ i ∈ s, c i * cofactor s q i

theorem modulus_mul_cofactor (s : Finset ι) (q : ι → ℤ) {i : ι} (hi : i ∈ s) :
    q i * cofactor s q i = ∏ j ∈ s, q j := by
  exact Finset.mul_prod_erase s q hi

/-- Every other prime-fraction contribution vanishes modulo a fixed modulus. -/
theorem modulus_dvd_other_contribution (s : Finset ι) (q c : ι → ℤ)
    {i j : ι} (hi : i ∈ s) (hij : i ≠ j) :
    q i ∣ c j * cofactor s q j := by
  apply dvd_mul_of_dvd_right
  exact Finset.dvd_prod_of_mem q (Finset.mem_erase.mpr ⟨hij, hi⟩)

/-- CRT turns the individual local numerator congruences into a single
divisibility statement by the whole product. -/
theorem product_dvd_reconstruction_error (s : Finset ι) (q c : ι → ℤ) (A : ℤ)
    (hcop : (s : Set ι).Pairwise (fun i j => IsCoprime (q i) (q j)))
    (hlocal : ∀ i ∈ s, q i ∣ A - c i * cofactor s q i) :
    (∏ i ∈ s, q i) ∣ A - reconstructionNumerator s q c := by
  apply Finset.prod_dvd_of_coprime hcop
  intro i hi
  have hother : q i ∣ ∑ j ∈ s.erase i, c j * cofactor s q j := by
    apply Finset.dvd_sum
    intro j hj
    exact modulus_dvd_other_contribution s q c hi (Finset.mem_erase.mp hj).1.symm
  have hs : reconstructionNumerator s q c =
      c i * cofactor s q i + ∑ j ∈ s.erase i, c j * cofactor s q j := by
    exact (Finset.add_sum_erase s (fun j => c j * cofactor s q j) hi).symm
  rw [hs]
  simpa only [sub_add_eq_sub_sub] using dvd_sub (hlocal i hi) hother

/-- Rational interpretation of the numerator assembled with integer cofactors. -/
theorem reconstructionNumerator_div_product (s : Finset ι) (q c : ι → ℤ)
    (hnz : ∀ i ∈ s, q i ≠ 0) :
    (reconstructionNumerator s q c : ℚ) / (∏ i ∈ s, q i : ℤ) =
      ∑ i ∈ s, (c i : ℚ) / q i := by
  have hD : (∏ i ∈ s, q i : ℤ) ≠ 0 := Finset.prod_ne_zero_iff.mpr hnz
  have hDq : ((∏ i ∈ s, q i : ℤ) : ℚ) ≠ 0 := by exact_mod_cast hD
  simp only [reconstructionNumerator, Int.cast_sum, Int.cast_mul, Finset.sum_div]
  apply Finset.sum_congr rfl
  intro i hi
  have hqi : (q i : ℚ) ≠ 0 := by exact_mod_cast hnz i hi
  have hc : (q i : ℚ) * cofactor s q i = (∏ j ∈ s, q j : ℤ) := by
    exact_mod_cast modulus_mul_cofactor s q hi
  field_simp
  linear_combination (c i : ℚ) * hc

/-- The difference between a rational with the common product denominator
and its CRT prime-fraction reconstruction is an integer. Local congruences
and pairwise coprimality are checked separately from this theorem. -/
theorem crt_integer_correction (s : Finset ι) (q c : ι → ℤ) (A : ℤ)
    (hnz : ∀ i ∈ s, q i ≠ 0)
    (hcop : (s : Set ι).Pairwise (fun i j => IsCoprime (q i) (q j)))
    (hlocal : ∀ i ∈ s, q i ∣ A - c i * cofactor s q i) :
    ∃ z : ℤ, (A : ℚ) / (∏ i ∈ s, q i : ℤ) -
      (∑ i ∈ s, (c i : ℚ) / q i) = z := by
  obtain ⟨z, hz⟩ := product_dvd_reconstruction_error s q c A hcop hlocal
  refine ⟨z, ?_⟩
  rw [← reconstructionNumerator_div_product s q c hnz, ← sub_div]
  have hD : (∏ i ∈ s, q i : ℤ) ≠ 0 := Finset.prod_ne_zero_iff.mpr hnz
  have hDq : ((∏ i ∈ s, q i : ℤ) : ℚ) ≠ 0 := by exact_mod_cast hD
  have hzq : (A : ℚ) - reconstructionNumerator s q c =
      (∏ i ∈ s, q i : ℤ) * (z : ℚ) := by exact_mod_cast hz
  rw [hzq]
  exact mul_div_cancel_left₀ (z : ℚ) hDq

end CRT

/-- The local inverse calculation behind a prime-power coefficient.
Here a source denominator is `r * v`, the prime-power modulus is `r * s`,
and their respective quotients in the common denominator are `k` and `e`.
The inverse is supplied by its exact modular equation. -/
theorem local_inverse_congruence (r v s e k inv : ℤ)
    (hr : r ≠ 0) (hquot : (r * v) * k = (r * s) * e)
    (hinv : r * s ∣ v * inv - 1) :
    r * s ∣ k - e * s * inv := by
  have hcancel : v * k = s * e := by
    apply mul_left_cancel₀ hr
    linear_combination hquot
  have heq : k - e * s * inv = (-k) * (v * inv - 1) := by
    linear_combination inv * hcancel
  rw [heq]
  exact dvd_mul_of_dvd_right hinv (-k)

section SourceSum

variable {ι κ : Type*} [DecidableEq ι]

/-- Exact common-denominator conversion for a finite weighted source sum. -/
theorem source_sum_common_denominator (t : Finset κ) (w μ : κ → ℤ) (D : ℤ)
    (hD : D ≠ 0) (hμ : ∀ j ∈ t, μ j ≠ 0)
    (hdiv : ∀ j ∈ t, μ j ∣ D) :
    (∑ j ∈ t, (w j : ℚ) / μ j) =
      ((∑ j ∈ t, w j * (D / μ j) : ℤ) : ℚ) / D := by
  simp only [Int.cast_sum, Int.cast_mul, Finset.sum_div]
  apply Finset.sum_congr rfl
  intro j hj
  have hμq : (μ j : ℚ) ≠ 0 := by exact_mod_cast hμ j hj
  have hDq : (D : ℚ) ≠ 0 := by exact_mod_cast hD
  have hquot : ((D / μ j : ℤ) : ℚ) * μ j = D := by
    exact_mod_cast Int.ediv_mul_cancel (hdiv j hj)
  field_simp
  linear_combination -(w j : ℚ) * hquot

/-- The partial-fractions-modulo-integers conclusion for an arbitrary
finite source sum. The local residues are hypotheses expressed solely
using integer arithmetic and the actual weighted common numerator. -/
theorem weighted_sum_crt_integer_correction
    (t : Finset κ) (w μ : κ → ℤ) (s : Finset ι) (q c : ι → ℤ)
    (hnz : ∀ i ∈ s, q i ≠ 0)
    (hcop : (s : Set ι).Pairwise (fun i j => IsCoprime (q i) (q j)))
    (hμ : ∀ j ∈ t, μ j ≠ 0)
    (hdiv : ∀ j ∈ t, μ j ∣ ∏ i ∈ s, q i)
    (hlocal : ∀ i ∈ s, q i ∣
      (∑ j ∈ t, w j * ((∏ k ∈ s, q k) / μ j)) - c i * cofactor s q i) :
    ∃ z : ℤ, (∑ j ∈ t, (w j : ℚ) / μ j) -
      (∑ i ∈ s, (c i : ℚ) / q i) = z := by
  rw [source_sum_common_denominator t w μ (∏ i ∈ s, q i)
    (Finset.prod_ne_zero_iff.mpr hnz) hμ hdiv]
  exact crt_integer_correction s q c _ hnz hcop hlocal

/-- A version that derives the local congruences from modular inverses.
For each modulus and source denominator, `r` is the stripped factor,
`v` the remaining unit, and `a = q / r`. The theorem includes all source
terms; when `r = 1`, their coefficient is automatically zero modulo `q`.
Prime-power valuation facts that supply these factorizations remain outside
this statement, but the complete finite modular algebra is proved here. -/
theorem partial_fractions_of_local_inverses
    (t : Finset κ) (w μ : κ → ℤ) (s : Finset ι) (q c : ι → ℤ)
    (r v a inv : ι → κ → ℤ)
    (hnz : ∀ i ∈ s, q i ≠ 0)
    (hcop : (s : Set ι).Pairwise (fun i j => IsCoprime (q i) (q j)))
    (hμ : ∀ j ∈ t, μ j ≠ 0)
    (hdiv : ∀ j ∈ t, μ j ∣ ∏ i ∈ s, q i)
    (hr : ∀ i ∈ s, ∀ j ∈ t, r i j ≠ 0)
    (hμfactor : ∀ i ∈ s, ∀ j ∈ t, μ j = r i j * v i j)
    (hqfactor : ∀ i ∈ s, ∀ j ∈ t, q i = r i j * a i j)
    (hinv : ∀ i ∈ s, ∀ j ∈ t, q i ∣ v i j * inv i j - 1)
    (hcoeff : ∀ i ∈ s, q i ∣ (∑ j ∈ t, w j * a i j * inv i j) - c i) :
    ∃ z : ℤ, (∑ j ∈ t, (w j : ℚ) / μ j) -
      (∑ i ∈ s, (c i : ℚ) / q i) = z := by
  apply weighted_sum_crt_integer_correction t w μ s q c hnz hcop hμ hdiv
  intro i hi
  have hterm : ∀ j ∈ t, q i ∣
      w j * ((∏ k ∈ s, q k) / μ j) -
        cofactor s q i * (w j * a i j * inv i j) := by
    intro j hj
    have hquot : (r i j * v i j) * ((∏ k ∈ s, q k) / μ j) =
        (r i j * a i j) * cofactor s q i := by
      rw [← hμfactor i hi j hj, ← hqfactor i hi j hj]
      rw [mul_comm (μ j), Int.ediv_mul_cancel (hdiv j hj)]
      exact (modulus_mul_cofactor s q hi).symm
    have hlocal := local_inverse_congruence
      (r i j) (v i j) (a i j) (cofactor s q i)
      ((∏ k ∈ s, q k) / μ j) (inv i j) (hr i hi j hj) hquot
      (by simpa only [← hqfactor i hi j hj] using hinv i hi j hj)
    rw [← hqfactor i hi j hj] at hlocal
    have heq : w j * ((∏ k ∈ s, q k) / μ j) -
        cofactor s q i * (w j * a i j * inv i j) =
        w j * (((∏ k ∈ s, q k) / μ j) - cofactor s q i * a i j * inv i j) := by
      ring
    rw [heq]
    exact dvd_mul_of_dvd_right hlocal (w j)
  have hsum : q i ∣ ∑ j ∈ t, (w j * ((∏ k ∈ s, q k) / μ j) -
      cofactor s q i * (w j * a i j * inv i j)) := Finset.dvd_sum hterm
  have hc := dvd_mul_of_dvd_right (hcoeff i hi) (cofactor s q i)
  have hall := dvd_add hsum hc
  have heq : (∑ j ∈ t, w j * ((∏ k ∈ s, q k) / μ j)) - c i * cofactor s q i =
      (∑ j ∈ t, (w j * ((∏ k ∈ s, q k) / μ j) -
        cofactor s q i * (w j * a i j * inv i j))) +
        cofactor s q i * ((∑ j ∈ t, w j * a i j * inv i j) - c i) := by
    simp only [Finset.sum_sub_distrib, ← Finset.mul_sum]
    ring
  rw [heq]
  exact hall

end SourceSum

end BBP

#print axioms BBP.rational_merge
#print axioms BBP.hypergeometric_merge
#print axioms BBP.shifted_rational_merge
#print axioms BBP.crt_integer_correction
#print axioms BBP.local_inverse_congruence
#print axioms BBP.weighted_sum_crt_integer_correction
#print axioms BBP.partial_fractions_of_local_inverses
