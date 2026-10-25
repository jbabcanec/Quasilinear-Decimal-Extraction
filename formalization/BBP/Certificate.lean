import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Algebra.Order.Floor.Ring
import Mathlib.Analysis.Real.Pi.Irrational
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Positivity

/-!
# Decimal output certificates and nearest-integer correction

This file proves the soundness part of Proposition `prop:certificate` and
the rounding step in equation `eq:z` in the manuscript.
The real approximation estimate is an explicit hypothesis. In particular,
this file does not assert a BBP identity, a tail estimate, or termination of
an executable program. Its eventual-passing theorem assumes valid real
approximations at arbitrarily high precision.
-/

namespace BBP

/-- An approximation within one half of an integer recovers that integer
by adding one half and taking the floor. -/
theorem nearest_integer_correction (approx : ℝ) (z : ℤ)
    (h : |approx - (z : ℝ)| < 1 / 2) :
    Int.floor (approx + 1 / 2) = z := by
  apply Int.floor_eq_iff.mpr
  obtain ⟨hlo, hhi⟩ := abs_lt.mp h
  constructor <;> linarith

/-- The four error allowances used for the exact CRT correction fit strictly
inside the nearest-integer recovery radius. -/
theorem correction_error_budget :
    (1 / 15 + 1 / 64 + 1 / 256 + 1 / 256 : ℝ) < 1 / 2 := by
  norm_num

/-- The paper's explicit total error allowance implies exact recovery. -/
theorem nearest_integer_correction_of_budget (approx : ℝ) (z : ℤ)
    (h : |approx - (z : ℝ)| < 1 / 15 + 1 / 64 + 1 / 256 + 1 / 256) :
    Int.floor (approx + 1 / 2) = z :=
  nearest_integer_correction approx z (h.trans correction_error_budget)

/-- The integer stopping test guarantees the indicated digit block.

Here `b` is the output radix multiplier (the manuscript uses `10^w`), and
`E` is the reciprocal error scale (the manuscript uses `2^G`). The integers
`d,r` are the quotient and remainder of `b*U` on division by `Q`.
No assumption about the floor or fractional part of `x` is needed: the
proof also excludes wraparound at either endpoint of the unit interval.
-/
theorem decimal_certificate
    (Q U b E d r a : ℤ) (x : ℝ)
    (hQ : 0 < Q) (hb : 0 < b) (hE : 0 < E)
    (hU0 : 0 ≤ U) (hUQ : U < Q)
    (hdiv : b * U = d * Q + r) (hr0 : 0 ≤ r) (hrQ : r < Q)
    (hleft : b * Q < E * r) (hright : b * Q < E * (Q - r))
    (happrox : |x - ((a : ℝ) + (U : ℝ) / Q)| < 1 / (E : ℝ)) :
    Int.floor ((b : ℝ) * Int.fract x) = d := by
  have hd0 : 0 ≤ d := by
    by_contra h
    have hd : d ≤ -1 := by omega
    nlinarith [mul_nonneg (le_of_lt hb) hU0,
      mul_le_mul_of_nonneg_right hd (le_of_lt hQ)]
  have hdb : d + 1 ≤ b := by
    by_contra h
    have hd : b ≤ d := by omega
    nlinarith [mul_pos hb (sub_pos.mpr hUQ),
      mul_le_mul_of_nonneg_right hd (le_of_lt hQ)]
  have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
  have hbr : (0 : ℝ) < b := by exact_mod_cast hb
  have hEr : (0 : ℝ) < E := by exact_mod_cast hE
  have hdr : (0 : ℝ) ≤ d := by exact_mod_cast hd0
  have hdbr : (d : ℝ) + 1 ≤ b := by exact_mod_cast hdb
  have hdivr : (b : ℝ) * U = (d : ℝ) * Q + r := by exact_mod_cast hdiv
  have hlr : (b : ℝ) * Q < (E : ℝ) * r := by exact_mod_cast hleft
  have hrr : (b : ℝ) * Q < (E : ℝ) * ((Q : ℝ) - r) := by
    exact_mod_cast hright
  have hfrac : (b : ℝ) * ((U : ℝ) / Q) = (d : ℝ) + (r : ℝ) / Q := by
    calc
      (b : ℝ) * ((U : ℝ) / Q) = ((b : ℝ) * U) / Q := by ring
      _ = (d : ℝ) + (r : ℝ) / Q := by
        field_simp
        nlinarith only [hdivr]
  have hgapL : (b : ℝ) / E < (r : ℝ) / Q := by
    apply (div_lt_div_iff₀ hEr hQr).mpr
    nlinarith
  have hgapR : (b : ℝ) / E < ((Q : ℝ) - r) / Q := by
    apply (div_lt_div_iff₀ hEr hQr).mpr
    nlinarith
  have hrem : ((Q : ℝ) - r) / Q = 1 - (r : ℝ) / Q := by
    field_simp
  obtain ⟨heL, heR⟩ := abs_lt.mp happrox
  have hscaledL := mul_lt_mul_of_pos_left heL hbr
  have hscaledR := mul_lt_mul_of_pos_left heR hbr
  rw [hrem] at hgapR
  simp only [div_eq_mul_inv] at hfrac hgapL hgapR hscaledL hscaledR
  have hcellL : (d : ℝ) < (b : ℝ) * (x - a) := by
    nlinarith only [hfrac, hgapL, hscaledL]
  have hcellR : (b : ℝ) * (x - a) < (d : ℝ) + 1 := by
    nlinarith only [hfrac, hgapR, hscaledR]
  have hx0 : 0 ≤ x - (a : ℝ) := by nlinarith
  have hx1 : x - (a : ℝ) < 1 := by nlinarith
  have hfloor : Int.floor x = a := Int.floor_eq_iff.mpr ⟨by linarith, by linarith⟩
  have hfract : Int.fract x = x - (a : ℝ) := by
    rw [Int.fract, hfloor]
  rw [hfract]
  exact Int.floor_eq_iff.mpr ⟨le_of_lt hcellL, hcellR⟩

/-- Specialization of `decimal_certificate` to the exact powers appearing
in the manuscript's decimal stopping test. -/
theorem decimal_output_certificate
    (w G : ℕ) (Q U d r a : ℤ) (x : ℝ)
    (hQ : 0 < Q) (hU0 : 0 ≤ U) (hUQ : U < Q)
    (hdiv : (10 : ℤ) ^ w * U = d * Q + r) (hr0 : 0 ≤ r) (hrQ : r < Q)
    (hleft : (10 : ℤ) ^ w * Q < (2 : ℤ) ^ G * r)
    (hright : (10 : ℤ) ^ w * Q < (2 : ℤ) ^ G * (Q - r))
    (happrox : |x - ((a : ℝ) + (U : ℝ) / Q)| < 1 / (2 : ℝ) ^ G) :
    Int.floor ((10 : ℝ) ^ w * Int.fract x) = d := by
  have h := decimal_certificate Q U ((10 : ℤ) ^ w) ((2 : ℤ) ^ G) d r a x
    hQ (pow_pos (by norm_num) w) (pow_pos (by norm_num) G)
    hU0 hUQ hdiv hr0 hrQ hleft hright (by simpa using happrox)
  simpa using h

/-- A nonintegral real number has a positive distance from the integer grid.
This lower bound holds for every integer, including integers outside the
two neighboring grid points. -/
theorem integer_grid_margin (y : ℝ) (hy : y ≠ (Int.floor y : ℝ)) :
    ∃ δ : ℝ, 0 < δ ∧ ∀ k : ℤ, δ ≤ |y - (k : ℝ)| := by
  refine ⟨min (Int.fract y) (1 - Int.fract y), ?_, ?_⟩
  · exact lt_min (Int.fract_pos.mpr hy) (sub_pos.mpr (Int.fract_lt_one y))
  · intro k
    by_cases hk : k ≤ Int.floor y
    · have hkr : (k : ℝ) ≤ Int.floor y := by exact_mod_cast hk
      have hf := Int.floor_add_fract y
      have hm := min_le_left (Int.fract y) (1 - Int.fract y)
      have ha := le_abs_self (y - (k : ℝ))
      linarith
    · have hk' : Int.floor y + 1 ≤ k := by omega
      have hkr : (Int.floor y : ℝ) + 1 ≤ k := by exact_mod_cast hk'
      have hf := Int.floor_add_fract y
      have hm := min_le_right (Int.fract y) (1 - Int.fract y)
      have ha := neg_le_abs (y - (k : ℝ))
      linarith

/-- Every irrational real has a strictly positive distance from a fixed
radix grid. Irrationality alone suffices for this qualitative statement;
it supplies no quantitative complexity bound. -/
theorem irrational_grid_margin (x : ℝ) (hx : Irrational x)
    (b : ℤ) (hb : 0 < b) :
    ∃ δ : ℝ, 0 < δ ∧ ∀ k : ℤ, δ < |x - (k : ℝ) / b| := by
  have hbr : (0 : ℝ) < b := by exact_mod_cast hb
  have hbx : Irrational ((b : ℝ) * x) := hx.intCast_mul (ne_of_gt hb)
  obtain ⟨μ, hμ, hmargin⟩ := integer_grid_margin ((b : ℝ) * x)
    (hbx.ne_int (Int.floor ((b : ℝ) * x)))
  refine ⟨μ / (2 * (b : ℝ)), div_pos hμ (by positivity), ?_⟩
  intro k
  have hmul : (b : ℝ) * (x - (k : ℝ) / b) = (b : ℝ) * x - k := by
    field_simp
  have habs : |(b : ℝ) * x - k| = (b : ℝ) * |x - (k : ℝ) / b| := by
    rw [← hmul, abs_mul, abs_of_pos hbr]
  have hm := hmargin k
  rw [habs] at hm
  apply (div_lt_iff₀ (by positivity : (0 : ℝ) < 2 * b)).mpr
  nlinarith

/-- Actual irrationality of pi, from mathlib, provides the nonzero decimal
grid distance used in the manuscript's qualitative termination argument. -/
theorem shifted_pi_decimal_margin (n w : ℕ) :
    ∃ δ : ℝ, 0 < δ ∧ ∀ k : ℤ,
      δ < |(10 : ℝ) ^ n * Real.pi - (k : ℝ) / (10 : ℝ) ^ w| := by
  have hi : Irrational ((10 : ℝ) ^ n * Real.pi) := by
    have h := irrational_pi.intCast_mul
      (by positivity : (10 : ℤ) ^ n ≠ 0)
    simpa using h
  have h := irrational_grid_margin ((10 : ℝ) ^ n * Real.pi) hi
    ((10 : ℤ) ^ w) (by positivity)
  simpa using h

/-- If the exact target is at least twice the approximation error from every
grid point, both integer inequalities in the stopping test must hold. -/
theorem certificate_of_grid_separation
    (Q U b E d r a : ℤ) (x : ℝ)
    (hQ : 0 < Q) (hb : 0 < b) (hE : 0 < E)
    (hdiv : b * U = d * Q + r) (hr0 : 0 ≤ r) (hrQ : r < Q)
    (hsep : ∀ k : ℤ, 2 / (E : ℝ) ≤ |x - (k : ℝ) / b|)
    (happrox : |x - ((a : ℝ) + (U : ℝ) / Q)| < 1 / (E : ℝ)) :
    b * Q < E * r ∧ b * Q < E * (Q - r) := by
  have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
  have hbr : (0 : ℝ) < b := by exact_mod_cast hb
  have hEr : (0 : ℝ) < E := by exact_mod_cast hE
  have hr0r : (0 : ℝ) ≤ r := by exact_mod_cast hr0
  have hrQr : (r : ℝ) < Q := by exact_mod_cast hrQ
  have hdivr : (b : ℝ) * U = (d : ℝ) * Q + r := by exact_mod_cast hdiv
  have hL : ((a : ℝ) + (U : ℝ) / Q) - ((a * b + d : ℤ) : ℝ) / b =
      (r : ℝ) / ((b : ℝ) * Q) := by
    push_cast
    field_simp
    nlinarith only [hdivr]
  have hR : ((a * b + d + 1 : ℤ) : ℝ) / b - ((a : ℝ) + (U : ℝ) / Q) =
      ((Q : ℝ) - r) / ((b : ℝ) * Q) := by
    push_cast
    field_simp
    nlinarith only [hdivr]
  have hLAbs : |((a : ℝ) + (U : ℝ) / Q) - ((a * b + d : ℤ) : ℝ) / b| =
      (r : ℝ) / ((b : ℝ) * Q) := by
    rw [hL, abs_of_nonneg (div_nonneg hr0r (le_of_lt (mul_pos hbr hQr)))]
  have hRAbs : |((a : ℝ) + (U : ℝ) / Q) - ((a * b + d + 1 : ℤ) : ℝ) / b| =
      ((Q : ℝ) - r) / ((b : ℝ) * Q) := by
    rw [abs_sub_comm, hR, abs_of_pos (div_pos (sub_pos.mpr hrQr) (mul_pos hbr hQr))]
  have hdistL : 1 / (E : ℝ) < (r : ℝ) / ((b : ℝ) * Q) := by
    have ht := abs_sub_le x ((a : ℝ) + (U : ℝ) / Q) (((a * b + d : ℤ) : ℝ) / b)
    rw [hLAbs] at ht
    have hs := hsep (a * b + d)
    rw [show (2 : ℝ) / E = 2 * (1 / (E : ℝ)) by ring] at hs
    linarith only [ht, hs, happrox]
  have hdistR : 1 / (E : ℝ) < ((Q : ℝ) - r) / ((b : ℝ) * Q) := by
    have ht := abs_sub_le x ((a : ℝ) + (U : ℝ) / Q) (((a * b + d + 1 : ℤ) : ℝ) / b)
    rw [hRAbs] at ht
    have hs := hsep (a * b + d + 1)
    rw [show (2 : ℝ) / E = 2 * (1 / (E : ℝ)) by ring] at hs
    linarith only [ht, hs, happrox]
  have hl := (div_lt_div_iff₀ hEr (mul_pos hbr hQr)).mp hdistL
  have hr := (div_lt_div_iff₀ hEr (mul_pos hbr hQr)).mp hdistR
  constructor
  · have h : (b : ℝ) * Q < (E : ℝ) * r := by nlinarith only [hl]
    exact_mod_cast h
  · have h : (b : ℝ) * Q < (E : ℝ) * ((Q : ℝ) - r) := by nlinarith only [hr]
    exact_mod_cast h

/-- Qualitative eventual passing for every sequence of approximations to an
irrational target with error less than `2^(-G)`.

This is a mathematical statement about arbitrary valid approximations. It
does not assert that an implementation computes those approximations or
prove an execution-time bound for a program's adaptive loop.
-/
theorem irrational_eventual_dyadic_certificate
    (x : ℝ) (hx : Irrational x) (b : ℤ) (hb : 0 < b) :
    ∃ G₀ : ℕ, ∀ G : ℕ, G₀ ≤ G → ∀ Q U d r a : ℤ,
      0 < Q → b * U = d * Q + r → 0 ≤ r → r < Q →
      |x - ((a : ℝ) + (U : ℝ) / Q)| < 1 / (2 : ℝ) ^ G →
      b * Q < (2 : ℤ) ^ G * r ∧ b * Q < (2 : ℤ) ^ G * (Q - r) := by
  obtain ⟨δ, hδ, hmargin⟩ := irrational_grid_margin x hx b hb
  obtain ⟨G₀, hG₀⟩ := pow_unbounded_of_one_lt (2 / δ) (by norm_num : (1 : ℝ) < 2)
  refine ⟨G₀, ?_⟩
  intro G hG Q U d r a hQ hdiv hr0 hrQ happrox
  have hpow : (0 : ℝ) < 2 ^ G := by positivity
  have hlarge : 2 / δ < (2 : ℝ) ^ G :=
    hG₀.trans_le (pow_le_pow_right₀ (by norm_num : (1 : ℝ) ≤ 2) hG)
  have hsmall : 2 / (2 : ℝ) ^ G < δ := by
    apply (div_lt_iff₀ hpow).mpr
    have h := (div_lt_iff₀ hδ).mp hlarge
    nlinarith only [h]
  apply certificate_of_grid_separation Q U b ((2 : ℤ) ^ G) d r a x
    hQ hb (by positivity) hdiv hr0 hrQ
  · intro k
    simpa using (hsmall.trans (hmargin k)).le
  · simpa using happrox

/-- Specialization of qualitative eventual passing to the shifted pi target
and decimal grid in the manuscript. -/
theorem shifted_pi_eventual_decimal_certificate (n w : ℕ) :
    ∃ G₀ : ℕ, ∀ G : ℕ, G₀ ≤ G → ∀ Q U d r a : ℤ,
      0 < Q → (10 : ℤ) ^ w * U = d * Q + r → 0 ≤ r → r < Q →
      |(10 : ℝ) ^ n * Real.pi - ((a : ℝ) + (U : ℝ) / Q)| < 1 / (2 : ℝ) ^ G →
      (10 : ℤ) ^ w * Q < (2 : ℤ) ^ G * r ∧
        (10 : ℤ) ^ w * Q < (2 : ℤ) ^ G * (Q - r) := by
  have hi : Irrational ((10 : ℝ) ^ n * Real.pi) := by
    have h := irrational_pi.intCast_mul
      (by positivity : (10 : ℤ) ^ n ≠ 0)
    simpa using h
  exact irrational_eventual_dyadic_certificate _ hi ((10 : ℤ) ^ w) (by positivity)

#print axioms nearest_integer_correction
#print axioms correction_error_budget
#print axioms nearest_integer_correction_of_budget
#print axioms decimal_certificate
#print axioms decimal_output_certificate
#print axioms integer_grid_margin
#print axioms irrational_grid_margin
#print axioms shifted_pi_decimal_margin
#print axioms certificate_of_grid_separation
#print axioms irrational_eventual_dyadic_certificate
#print axioms shifted_pi_eventual_decimal_certificate

end BBP
