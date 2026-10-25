import BBP.Certificate
import BBP.Series
import BBP.Guard

/-!
Connect the finite BBP approximation, its proved tail bound, and the
integer stopping test. The BBP identity and the exact finite representation
are explicit inputs. This is not a refinement theorem for the C++ program.
-/

namespace BBP

/-- If the computed rational really represents the shifted finite BBP sum,
the integer stopping test proves the requested decimal window is correct.
No assumption about random carries or the distribution of pi's digits is used. -/
theorem certified_bbp_window
    (hbbp : HasSum Series.bbpTerm Real.pi)
    (n w G : ℕ) (Q U d r a : ℤ)
    (hQ : 0 < Q) (hU0 : 0 ≤ U) (hUQ : U < Q)
    (hrepresentation :
      (10 : ℝ) ^ n * ∑ k ∈ Finset.range (Series.cutoff n G + 1), Series.bbpTerm k =
        (a : ℝ) + (U : ℝ) / Q)
    (hdiv : (10 : ℤ) ^ w * U = d * Q + r) (hr0 : 0 ≤ r) (hrQ : r < Q)
    (hleft : (10 : ℤ) ^ w * Q < (2 : ℤ) ^ G * r)
    (hright : (10 : ℤ) ^ w * Q < (2 : ℤ) ^ G * (Q - r)) :
    Int.floor ((10 : ℝ) ^ w * Int.fract ((10 : ℝ) ^ n * Real.pi)) = d := by
  apply decimal_output_certificate w G Q U d r a _
    hQ hU0 hUQ hdiv hr0 hrQ hleft hright
  rw [← hrepresentation]
  exact Series.shifted_bbp_cutoff hbbp n G

/-- The stated initial guard passes, provided the published exponent-eight
rational-approximation bound holds at this denominator and the finite sum
is represented exactly. The quantitative hypothesis is deliberately exposed. -/
theorem initial_pass_certified_bbp_window
    (hbbp : HasSum Series.bbpTerm Real.pi)
    (n w : ℕ) (Q U d r a : ℤ)
    (hquant : ∀ p : ℤ, 1 / ((10 : ℝ) ^ (n + w)) ^ 8 <
      |Real.pi - (p : ℝ) / (10 : ℝ) ^ (n + w)|)
    (hQ : 0 < Q) (hU0 : 0 ≤ U) (hUQ : U < Q)
    (hrepresentation :
      (10 : ℝ) ^ n * ∑ k ∈ Finset.range (Series.cutoff n (initialGuard n w) + 1),
        Series.bbpTerm k = (a : ℝ) + (U : ℝ) / Q)
    (hdiv : (10 : ℤ) ^ w * U = d * Q + r) (hr0 : 0 ≤ r) (hrQ : r < Q) :
    ((10 : ℤ) ^ w * Q < (2 : ℤ) ^ initialGuard n w * r ∧
      (10 : ℤ) ^ w * Q < (2 : ℤ) ^ initialGuard n w * (Q - r)) ∧
      Int.floor ((10 : ℝ) ^ w * Int.fract ((10 : ℝ) ^ n * Real.pi)) = d := by
  have hsep := initial_guard_separates_boundaries Real.pi n w hquant
  have happrox : |(10 : ℝ) ^ n * Real.pi - ((a : ℝ) + (U : ℝ) / Q)| <
      1 / (2 : ℝ) ^ initialGuard n w := by
    rw [← hrepresentation]
    exact Series.shifted_bbp_cutoff hbbp n (initialGuard n w)
  have hc := certificate_of_grid_separation Q U ((10 : ℤ) ^ w)
    ((2 : ℤ) ^ initialGuard n w) d r a ((10 : ℝ) ^ n * Real.pi)
    hQ (by positivity) (by positivity) hdiv hr0 hrQ
    (by intro k; simpa using (hsep k).le) (by simpa using happrox)
  exact ⟨hc, certified_bbp_window hbbp n w (initialGuard n w) Q U d r a
    hQ hU0 hUQ hrepresentation hdiv hr0 hrQ hc.1 hc.2⟩

end BBP
