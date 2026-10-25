import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

/-!
Explicit arithmetic used to turn the published exponent-eight rational
approximation bound into a linear guard. The rational approximation bound
is a hypothesis here: this does not formalize the published theorem for pi.
-/

namespace BBP

def initialGuard (n w : ℕ) : ℕ := 32 * (n + w + 2)

theorem guard_exponent_margin (n w : ℕ) :
    1 + 4 * (7 * n + 8 * w) < initialGuard n w := by
  simp only [initialGuard]
  omega

theorem guard_power_margin (n w : ℕ) :
    2 * (10 : ℕ) ^ (7 * n + 8 * w) < 2 ^ initialGuard n w := by
  calc
    2 * (10 : ℕ) ^ (7 * n + 8 * w) ≤
        2 * 16 ^ (7 * n + 8 * w) := by gcongr; norm_num
    _ = 2 ^ (1 + 4 * (7 * n + 8 * w)) := by
      rw [pow_add (2 : ℕ) 1 (4 * (7 * n + 8 * w)), pow_one,
        pow_mul (2 : ℕ) 4 (7 * n + 8 * w)]
      rfl
    _ < 2 ^ initialGuard n w := by
      exact pow_lt_pow_right₀ (by norm_num) (guard_exponent_margin n w)

theorem guard_error_margin (n w : ℕ) :
    2 / (2 : ℝ) ^ initialGuard n w < 1 / (10 : ℝ) ^ (7 * n + 8 * w) := by
  have h : 2 * (10 : ℝ) ^ (7 * n + 8 * w) < 2 ^ initialGuard n w := by
    exact_mod_cast guard_power_margin n w
  apply (div_lt_div_iff₀ (by positivity) (by positivity)).mpr
  simpa using h

/-- Applying an exponent-eight rational approximation hypothesis at q=10^(n+w).
The premise is deliberately visible and is not supplied by an axiom. -/
theorem scaled_boundary_separation (x : ℝ) (n w : ℕ) (p : ℤ)
    (h : 1 / ((10 : ℝ) ^ (n + w)) ^ 8 <
      |x - (p : ℝ) / (10 : ℝ) ^ (n + w)|) :
    1 / (10 : ℝ) ^ (7 * n + 8 * w) <
      |(10 : ℝ) ^ n * x - (p : ℝ) / (10 : ℝ) ^ w| := by
  have hm := mul_lt_mul_of_pos_left h (show 0 < (10 : ℝ) ^ n by positivity)
  have hscale : (10 : ℝ) ^ n *
      |x - (p : ℝ) / (10 : ℝ) ^ (n + w)| =
      |(10 : ℝ) ^ n * x - (p : ℝ) / (10 : ℝ) ^ w| := by
    calc
      (10 : ℝ) ^ n * |x - (p : ℝ) / (10 : ℝ) ^ (n + w)| =
          |(10 : ℝ) ^ n| * |x - (p : ℝ) / (10 : ℝ) ^ (n + w)| := by
            rw [abs_of_pos (show 0 < (10 : ℝ) ^ n by positivity)]
      _ = |(10 : ℝ) ^ n * (x - (p : ℝ) / (10 : ℝ) ^ (n + w))| :=
        (abs_mul _ _).symm
      _ = _ := by
        congr 1
        rw [pow_add]
        field_simp
  have hpow : ((10 : ℝ) ^ (n + w)) ^ 8 =
      (10 : ℝ) ^ n * (10 : ℝ) ^ (7 * n + 8 * w) := by
    rw [← pow_mul, ← pow_add]
    congr 1
    omega
  rw [hscale, hpow] at hm
  have hcancel : (10 : ℝ) ^ n *
      (1 / ((10 : ℝ) ^ n * (10 : ℝ) ^ (7 * n + 8 * w))) =
      1 / (10 : ℝ) ^ (7 * n + 8 * w) := by field_simp
  rwa [hcancel] at hm

theorem initial_guard_separates_boundaries (x : ℝ) (n w : ℕ)
    (h : ∀ p : ℤ, 1 / ((10 : ℝ) ^ (n + w)) ^ 8 <
      |x - (p : ℝ) / (10 : ℝ) ^ (n + w)|) :
    ∀ p : ℤ, 2 / (2 : ℝ) ^ initialGuard n w <
      |(10 : ℝ) ^ n * x - (p : ℝ) / (10 : ℝ) ^ w| := by
  intro p
  exact (guard_error_margin n w).trans (scaled_boundary_separation x n w p (h p))

end BBP
