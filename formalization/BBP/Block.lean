import BBP.Reconstruction
import BBP.Certificate

/-!
# Assembling exact weighted-block reconstruction

This file joins the finite CRT decomposition to nearest-integer recovery.
The approximation error is an explicit premise. Thus it verifies the returned
integer and rational once that error bound and the modular data are available;
it does not claim that an executable has generated or certified those inputs.
-/

namespace BBP

/-- Once the prime-fraction discrepancy is an integer and the real estimate
is accurate to less than one half, the rounded correction and the returned
integer numerator over the original denominator are both exact. -/
theorem exact_block_from_integer_correction
    (S : ℚ) (b d : ℤ) (approx : ℝ) (hd : d ≠ 0)
    (hinteger : ∃ z : ℤ, S - (b : ℚ) / d = z)
    (happrox : |approx - ((S - (b : ℚ) / d : ℚ) : ℝ)| < 1 / 2) :
    ((Int.floor (approx + 1 / 2) : ℚ) = S - (b : ℚ) / d) ∧
    (((b + Int.floor (approx + 1 / 2) * d : ℤ) : ℚ) / d = S) := by
  obtain ⟨z, hz⟩ := hinteger
  have hround : Int.floor (approx + 1 / 2) = z :=
    nearest_integer_correction approx z (by simpa only [hz, Rat.cast_intCast] using happrox)
  have hdq : (d : ℚ) ≠ 0 := by exact_mod_cast hd
  constructor
  · simpa only [hround] using hz.symm
  · rw [hround]
    push_cast
    calc
      ((b : ℚ) + (z : ℚ) * d) / d = (b : ℚ) / d + z := by field_simp
      _ = S := by rw [← hz]; ring

/-- Integrated block theorem. The finite weights may in particular be the
manuscript's `16^(t - 1 - j)`, and the denominators may be `A + a*j`.

The factorization and modular inverse premises supply the CRT decomposition;
the four-term paper error budget supplies nearest-integer recovery. The
conclusion identifies both the recovered integer and the exact returned
numerator over the uncancelled product denominator.
-/
theorem weighted_block_crt_reconstruction
    {ι κ : Type*} [DecidableEq ι]
    (t : Finset κ) (w μ : κ → ℤ) (s : Finset ι) (q c : ι → ℤ)
    (r v a inv : ι → κ → ℤ) (approx : ℝ)
    (hnz : ∀ i ∈ s, q i ≠ 0)
    (hcop : (s : Set ι).Pairwise (fun i j => IsCoprime (q i) (q j)))
    (hμ : ∀ j ∈ t, μ j ≠ 0)
    (hdiv : ∀ j ∈ t, μ j ∣ ∏ i ∈ s, q i)
    (hr : ∀ i ∈ s, ∀ j ∈ t, r i j ≠ 0)
    (hμfactor : ∀ i ∈ s, ∀ j ∈ t, μ j = r i j * v i j)
    (hqfactor : ∀ i ∈ s, ∀ j ∈ t, q i = r i j * a i j)
    (hinv : ∀ i ∈ s, ∀ j ∈ t, q i ∣ v i j * inv i j - 1)
    (hcoeff : ∀ i ∈ s, q i ∣ (∑ j ∈ t, w j * a i j * inv i j) - c i)
    (happrox : |approx - (((∑ j ∈ t, (w j : ℚ) / μ j) -
      (∑ i ∈ s, (c i : ℚ) / q i) : ℚ) : ℝ)| <
        1 / 15 + 1 / 64 + 1 / 256 + 1 / 256) :
    ((Int.floor (approx + 1 / 2) : ℚ) =
      (∑ j ∈ t, (w j : ℚ) / μ j) - (∑ i ∈ s, (c i : ℚ) / q i)) ∧
    (((reconstructionNumerator s q c +
        Int.floor (approx + 1 / 2) * (∏ i ∈ s, q i) : ℤ) : ℚ) /
      (∏ i ∈ s, q i : ℤ) = ∑ j ∈ t, (w j : ℚ) / μ j) := by
  have hint := partial_fractions_of_local_inverses t w μ s q c r v a inv
    hnz hcop hμ hdiv hr hμfactor hqfactor hinv hcoeff
  have hT := reconstructionNumerator_div_product s q c hnz
  have h := exact_block_from_integer_correction
    (∑ j ∈ t, (w j : ℚ) / μ j)
    (reconstructionNumerator s q c) (∏ i ∈ s, q i) approx
    (Finset.prod_ne_zero_iff.mpr hnz)
    (by simpa only [hT] using hint)
    (by simpa only [hT] using happrox.trans correction_error_budget)
  simpa only [hT] using h

#print axioms exact_block_from_integer_correction
#print axioms weighted_block_crt_reconstruction

end BBP
