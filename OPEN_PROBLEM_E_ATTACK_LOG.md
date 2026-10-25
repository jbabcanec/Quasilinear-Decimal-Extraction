# Open Problem E — Complete Attack Log

> Historical research log. For the current interpretation of proof drafts,
> implementation limits, and measurements, read [RESEARCH_STATUS.md](RESEARCH_STATUS.md).
> Earlier speed and completion claims are superseded there. The entries below
> are retained as a chronological record, not a current readiness assessment.

**2026-09-06 update:** §19 gives a constructive O(N log^(5/2) N) bound for
the exact weighted core, using short-block LCM compression and certified
real-series integer reconstruction. §18 resolves the unweighted atom in
O(N log²N) and identifies a gap in the unrestricted interface theorem.
The weighted O(N log²N) target remains open. These are mathematical bounds;
the timing experiments are not evidence for their logarithmic exponents.
Section 20 records a further analytic audit and a certified decimal-output
construction in the same bound, using exact word aggregation.

**Goal.** Decimal digit extraction of π from binary BBP in \(O(N\log^2 N)\) (AGM parity).
**Original manuscript/default engine.** \(O(N\log^3 N)\) via the v3 division-free pipeline.
**New construction.** \(O(N\log^{5/2}N)\) for V; full proof in
`notes/weighted_block_crt_2026-09-06.md`, independent review/priority pending.
**Barrier.** Open Problem E: the additive core
\[
V=\sum_t\pm 2^{J_{\max}-J_t}(\mu_t^{-1}\bmod 2^{J_{\max}})\bmod 2^{J_{\max}}
\]
in \(O(\mathsf{M}(N)\log N)\). The unweighted surrogate
\(H=\sum_{k=1}^{K}(8k+1)^{-1}\bmod 2^{\Theta(N)}\) is solved in that budget;
the old claim of cost-preserving equivalence to V was not established.

**Philosophy.** Don’t like the game, change the rules — but every rule-change gets a kill law or a lemma.

**Status.** The log³ upper bound has been improved to log^(5/2) by a new
construction. The log² target is not solved. Earlier sections retain the
historical claims and must be read with the corrections in §§18–19.

---

## 0. What was already cracked before this campaign (Theorem crack)

| Object | Status | Cost |
|---|---|---|
| \(P_0=\prod_{k\le K}(8k+1)\bmod 2^m\) | CRACKED (Borwein factorization / sieve mass \(\Theta(N)\)) | \(O(N\log^2 N)\) |
| \(E=\sum\operatorname{Log}(8k+1)=\operatorname{Log}P_0/\operatorname{Log}9\) | CRACKED | same |
| Dyadic power sums \(\sum(1+8k)^{2^r}\bmod 2^{2r}\) | CRACKED | same |
| Diamond log-gamma at args \(8/(1+8K)\) (via EM reverse of \(E\)) | CRACKED | same |
| Full extraction pipeline consuming \(V\) | DONE | \(O(N\log^3 N)\) |

**Razor.** Differentiation destroys factorization: \(e_T=\prod\mu_i\) smooth; \(e_{T-1}=\sum_i\prod_{j\neq i}\mu_j\) generic. Leave-one-out is the barrier.

---

## 1. Campaign map (this conversation, chronological)

### Phase A — Fresh doors / rule-change hunt

**Scripts:** `rule_change_hunt.py`, then `door_beta_gamma.py`

| Attempt | Idea | Verdict |
|---|---|---|
| Approx \(V\) (limb support) | Digit only needs \(O(\log N)\) bits of fractional part | **KILL.** Top \(G\) bits of \(XV\) flip unless \(V\) correct through limb \(J-G\). Still \(\Theta(N)\) bits of \(V\). |
| Diff quotient of cracked \(\Pi\) | \(H\approx\varepsilon^{-1}\mathrm{Log}(\Pi_\varepsilon/\Pi_0)\) | **KILL (rate law).** Peak yield \(\sim a\) bits at sieve-cost \(\sim 1\), then cost explodes. |
| LCM-denominator binary splitting | Compress merge mass via lcm of denoms | **KILL.** Operand mass still \(\sim N\log^3 N\). |
| Fixed-dimension Landen / digamma duplication | Keep Landen dimension bounded | **KILL.** Dimension grows; no low-height collapse. |
| Archimedean short-window | \(O(\log N)\) real terms for digit | **KILL.** Every term matters mod 1 in the 2-adic pipeline. |
| **Door α** (Landen) | Classical AGM/Landen on digamma | **KILL.** Branches; no collapse. |
| **Door β** (admissible deformation) | Recover \(H\) from \(\mathrm{Log}\Pi(u)\) at sieve-mass-\(\Theta(N)\) points | **KILL.** Finite-difference stencils rate-limited to \(O(\log N)\) bits of \(H\). LLL: \(H\) not in span of admissible Logs, poly-weighted Logs, or \(O(\log N)\) Bernoulli/Mazur moments. |
| **Door γ** (merge model lower bound) | Prove third log in restricted SLP / merge model | **Partial theorem.** Merge depth / mass \(\Omega(N\log^3 N)\) in the model; \(e_T\) vs \(e_{T-1}\) razor formalized. |

**Winning shapes named but not constructed:** B1 (digamma via bit-burst — restates E), B2/B3 (Mazur pairing / unknown).

### Phase B — Exact identities for the atom \(H\)

**Scripts:** `diff_quotient_H.py`, `diff_approach.py`, `faulhaber_H.py`, `em_2adic.py`, `h_probe.py`

| Identity | Statement | Complexity |
|---|---|---|
| Difference quotient | \(H\equiv(\mathrm{Log}\,P(2^r)-\mathrm{Log}\,P(0))/2^r\pmod{2^m}\) for \(r\ge m+2-v_2(S_2)\) | Large-\(r\): not soft-linear; admissible \(r\): only \(O(\log N)\) bits |
| **Dual / high-word** | \(\prod_k(R+8k+1)\equiv P_0(1+RH)\pmod{R^2}\), \(R=2^m\) | \(H=\) high word of \(P_R\cdot P_0^{-1}\). \(P_0\) cheap; \(P_R\) is the barrier |
| Faulhaber-geometric | \(H=\sum_j(-8)^j\sum_{k\le K}k^j\) (2-adic) | Bit-exact; needs \(\Theta(m)\) power sums → Bernoulli mass / \(O(N^2)\) |
| 2-adic Euler–Maclaurin | \(E\) ↔ Diamond log-gamma (cracked side); digamma derivative is \(H\) | EM+Bernoulli matches \(H\) (sanity); no exp/log stand-in for Bernoulli correction |

**Special-function hunt (`diff_approach.py`):** no low-height exp/log combination replaces Bernoulli correction for \(H\).

### Phase C — Dual product pinned as the barrier

**Scripts:** `attack_dual_product.py`, `odd_prod_lift.py`, `egf_power_sums.py`

**Pinned form.** Compute
\[
P_R=B(m,K)=\prod_{k=1}^{K}(2^m+8k+1)\bmod 2^{2m}
\]
in \(O(N\log^2 N)\) (digit regime \(K,m=\Theta(N)\)).

| Attack | Result |
|---|---|
| Half-word dual multiply | Correct for recovering \(H\); cost tracks \(\Theta(K\,\mathsf{M}(m))\) |
| AP-halving / multifactorial | Restates \(H\) or wrong regime |
| LLL vs \(P_R\) | No cheap algebraic relation found |
| OEIS A395199 \((2^n-1)!!\) lift by squaring | **Works for FULL odd double-factorial** in poly(\(m\)); **wrong regime** for short window length \(K=\Theta(m)\) at height \(2^m\) |
| **new_block lemma** | \(B_m=\prod_{t=0}^{2^{m-3}-1}(2^m+8t+1)\) satisfies \(v_2(B_{m+1}/B_m^2-1)=2m-3\) — squaring for **dyadic-length** 1-mod-8 blocks (exponential length); not digit-\(K\) |
| Dual factorial \(m!=(m_0!)(1+m_1 H_{m_0}\varepsilon)\) | Restates dual picture |
| EGF batch power sums \(\sum e^{kt}\) | **Blocked:** \(1/n!\notin\mathbb{Z}_2\). \(\exp(8t)\) OK valuation-wise, but \(\sum_k\exp(kt)\) needs \(\exp(t)\). OGF collapses to digamma/\(H\) |
| Caruso bit-burst / digamma | Digamma not holonomic; Bernoulli wall remains |
| Product-tree claims of \(O(\mathsf{M}(n)\log n)\) for \(n\) full-size residues mod \(2^{O(n)}\) | Wrong model; empirical cost \(\sim K\cdot\mathsf{M}(m)\) |

### Phase D — BSGS, Stirling, Q-structure, Borwein dual

**Scripts:** `bsgs_dual.py`, `batch_Sj_attack.py`, `q_structure.py`, `dual_borwein_hunt.py`

| Attack | Result |
|---|---|
| BSGS block poly \(Q(x)=\prod_{t<s}(A+8t+8x)\), \(s\sim\sqrt{K}\) | **Correct** dual product / \(H\) recovery. Complexity rung \(\widetilde{O}(N^{3/2})\) — below naive \(N^2\), **above** target \(O(N\log^2 N)\) and above existing \(O(N\log^3 N)\) extract |
| Stirling-2nd batch \(S_j\) | Bit-exact via hockey-stick falling sums. Still \(\Theta(J^2)=\Theta(N^2)\). Soft-linear Stirling triangle not known; EGF blocked in \(\mathbb{Z}_2\) |
| Soft-linear \(S_j\) ↔ soft-linear Bernoulli (Harvey \(O(m^2\log^{2+\varepsilon}m)\)) ↔ \(H\) | Equivalence pinned |
| \(Q(x)\) coeff \(v_2\) profile | Unit coeff present; Stirling-1-like; no sparse 2-adic gift |
| Truncated \(x^2\) product tree | Ops \(\Theta(K)\); **equals leave-one-out / Open Problem E** |
| Dual Borwein doubling monomials | Algebraic pairing \((R+16j+1)(R+16j+9)=R^2+\cdots\) only. Sieve mass of dual factors is \(\Theta(N^2/\log)\) — **not smooth**, Borwein does not transfer |
| \(n!\bmod 2^{2m}\) / Codeforces \(p^k\) log-expansion | Factorial vanishes for \(n\ge 2m\); log-expansion of AP product collapses to Faulhaber/power sums |

### Phase E — Dual squaring lemma (THE NEW LEMMA)

**Scripts:** `dual_lift.py`, `dual_lift2.py`, `dual_jump.py`
**Paper:** `decimal_extraction_paper.tex` — Prop.~E-equiv, Lemma dual-square (after Open Problem E)

---

## 2. THE LEMMA (write this on the wall)

### Lemma (dual squaring)

Let \(R=2^m\), \(\mu_k=8k+1\), \(B(m,K)=\prod_{k=1}^{K}(R+\mu_k)\), \(P_0(K)=\prod_{k=1}^{K}\mu_k\). Then **exactly**
\[
\frac{B(m,K)^{2}}{P_0(K)\,B(m+1,K)}
=\prod_{k=1}^{K}\left(1+\frac{R^{2}}{\mu_k(\mu_k+2R)}\right),
\]
and consequently
\[
v_2\bigl(B(m,K)^{2}-P_0(K)\,B(m+1,K)\bigr)\ge 2m.
\]

**One-line proof of the ratio.** Per factor,
\[
\frac{(1+R/\mu)^{2}}{1+2R/\mu}=\frac{(\mu+R)^{2}}{\mu(\mu+2R)},
\]
and
\[
\frac{(\mu+R)^{2}}{\mu(\mu+2R)}-1=\frac{R^{2}}{\mu(\mu+2R)}.
\]
Product over \(k\) gives the claim. Since each \(\mu\) is odd and \(\mu+2R\) is odd, \(v_2(R^{2}/(\mu(\mu+2R)))=2m\), so the product is \(1+O(2^{2m})\).

**Machine check.** \(6\le m\le 21\), \(K\in\{m,2m,8,16\}\); ratio identity exact; \(v_2(\mathrm{err})\ge 2m\) (often \(=2m\) when \(K\sim m\)).

### What the lemma does **not** give

- \(m\mapsto m+1\) lift at dual-barrier width \(2(m+1)\) needs **surplus** valuation. Surplus is \(O(1)\).
- Correcting with the expansion
  \[
  \mathrm{ratio}=1+R^{2}\sum\mu^{-2}+O(R^{4})
  \]
  reintroduces inverse-power sums \(\sum\mu^{-2},\sum\mu^{-4},\ldots\) (the **additive tower**).
- \(\Theta(N)\) unit steps × \(\mathsf{M}(N)\) = \(O(N^{2}\log N)\) — same **rate-law family** as door β / Postnikov / admissible \(t=2^r\).
- Geometric jump \(m\mapsto 2m\): chained \(B(m)^{2^s}\approx P_0^{2^s-1}B(m+s)\) loses margin linearly in \(s\); fails full width for \(s\ge 3\).
- \(B(2m)\equiv P_0\pmod{2^{2m}}\) is **trivial** (\(2^{2m}+\mu\equiv\mu\)) and loses \(H\).

---

## 3. Equivalence dictionary (Prop. E-equiv in the paper)

In digit regime \(K,m=\Theta(N)\), soft-linear for **any one** ⇒ soft-linear for **all** (via \(O(\mathsf{M}(N)\log N)\) post-processing + Theorem crack for \(P_0\)):

| ID | Atom |
|---|---|
| E0 | \(H=\sum_{k=1}^{K}(8k+1)^{-1}\bmod 2^m\) |
| E1 | Dual product \(B(m,K)\bmod 2^{2m}\) |
| E2 | Batch power sums \(S_j(K)\) for \(j\le m/3\) |
| E3 | Bernoulli batch \(B_0,\ldots,B_m\) to precision \(2^{O(m)}\) (Harvey: \(O(m^2\log^{2+\varepsilon}m)\)) |
| E4 | \(2\)-adic digamma \(\psi_2(K+1+1/8)\) to \(m\) bits (not holonomic ⇒ no Caruso DE) |
| E5 | Short rising factorial \((a)_K\bmod 2^{2m}\), \(a=2^m+1\), \(K=\Theta(m)\) |

---

## 4. Complexity rungs (honest ladder)

| Method | Cost for atom / dual | Notes |
|---|---|---|
| Direct inverse sum | \(\Theta(K\,\mathsf{M}(m)\log m)\) | Baseline |
| Product tree / merge (full extract) | \(O(N\log^3 N)\) | **Current champion** for digits |
| BSGS block dual product | \(\widetilde{O}(N^{3/2})\) | Below \(N^2\), above extract |
| Dual \(m\)-lift (+ corrections) | \(O(N^{2}\log N)\) | Rate law |
| Faulhaber + Bernoulli | \(O(N^{2}\mathrm{polylog})\) | Harvey |
| **Target** | \(O(N\log^2 N)\) | Open Problem E |
| AGM all-digits | \(O(N\log^2 N)\) | Different problem (not extraction) |

---

## 5. Literature / external pulls (what we borrowed)

| Source | What we took | Outcome |
|---|---|---|
| Borwein 1985 factorial method | Crack for smooth \(P_0\) | Already Theorem crack |
| OEIS A395199 | \((2^k-1)!!\) lift by squaring | Full odds only; wrong regime for short window |
| Caruso–Mezzarobba–Takayama–Vaccon bit-burst | Soft-linear \(p\)-adic special functions | Holonomic DEs only; digamma excluded |
| Harvey Bernoulli multimodular | \(O(k^2\log^{2+\varepsilon}k)\) | Confirms E3 wall |
| Codeforces / \(n!\bmod p^k\) | Log-expansion via power sums | Same Faulhaber wall; \(p=2\) factorial vanishes |
| Diamond \(G_p\) / Morita \(\Gamma_p\) | Multiplication / difference equations | Log-gamma cracked side; derivative is \(H\) |
| Johansson / binary splitting harmonics | Best known for isolated \(H_n\bmod p^k\) | Still not soft-linear |
| Postnikov / Hiary composition | Rate laws for deformations | Door β family |
| Costa–Kedlaya–Roe remainder trees | Amortized \(p\)-adic \(\Gamma\) over many \(p\) | Not single-prime precision at \(p=2\) |

---

## 6. Files touched / added this campaign

### Validation scripts (new or heavily used for E)

```
validation/rule_change_hunt.py
validation/door_beta_gamma.py
validation/diff_quotient_H.py
validation/diff_approach.py
validation/faulhaber_H.py
validation/attack_dual_product.py
validation/odd_prod_lift.py
validation/egf_power_sums.py
validation/bsgs_dual.py
validation/batch_Sj_attack.py
validation/q_structure.py
validation/dual_borwein_hunt.py
validation/dual_lift.py
validation/dual_lift2.py
validation/dual_jump.py
validation/bet1_rising_factorial.py
validation/bet2_inv_power_batch.py
validation/bet3_leave_one_out.py
validation/bet4_iwasawa_mazur.py
validation/bet5_rulechange_approx.py
OPEN_PROBLEM_E_ATTACK_LOG.md          # this catalog
```

### Paper updates

- `decimal_extraction_paper.tex`: Prop. E-equiv (five equivalent atoms), Lemma dual-square, rate-law paragraph after Open Problem E.

### Related earlier validators still in the picture

`merge_kary.py`, `postnikov_compose.py`, `em_2adic.py`, `h_probe.py`, `borwein_crack.py`, `question_e.py`, `atom_check.py`, fences / Dwork / orbit scripts (mostly killed faces from earlier campaigns).

---

## 7. Live bets — results (continued after the log was first written)

### Bet 1 — short rising factorial at height \(2^m\)
**Script:** `validation/bet1_rising_factorial.py`

| Sub-attack | Result |
|---|---|
| Nilpotent Log expansion | Mod \(R^2\) only \(Q_1=H\) survives — tautology with E0/E1 |
| Stirling-1 row evaluation | Row costs \(\Theta(K^2)\); no compression |
| Downward dual-square lift | Agree \(v_2\) typically \(2m-1\) or \(2m-2\); **1–2 bits short** of need \(2m\) |
| Factorial / odd-part ratio | \(v_2(n!)\gg 2m\) at height \(2^m\) — vanishes mod modulus |
| Exact rational CRT | Numerator \(\Theta(N^2)\) bits — \(\Omega(N^2)\) |
| Dyadic split \(K=K_1+K_2\) | Shift \(8K_1\) not a high power of 2; door-β rate law |

**Verdict:** no soft-linear rising-factorial algorithm. Bet 1 closed as a *separate* door (collapses to E0–E5).

### Bet 2 — inverse-power batch / additive tower
**Script:** `validation/bet2_inv_power_batch.py`

| Sub-attack | Result |
|---|---|
| Reverse Newton / reciprocal roots | \(e_{K-1}/e_K=H\); \(e_{K-1}\) is leave-one-out = E |
| Freshman dream on \(Q_j\) | Squaring produces cross terms (partition products of inverses) |
| Higher dual nilpotency \(B\bmod R^{J+1}\) | Encodes \(Q_1..Q_J\) but is **harder**, not easier |
| \(Q_2\) vs \(Q_1\) timing | Same \(\Theta(K)\) class |
| Faulhaber for all \(Q_j\) | One \(S_j\) batch yields all \(Q_j\) — **tower ≡ Bernoulli wall ≡ E** |

**Verdict:** soft-linear inverse-power *batch* is soft-linear-equivalent to Open Problem E. No separation of \(Q_2\) from \(Q_1\). Bet 2 closed as independent door.

### Bet 3 — leave-one-out Borwein-style compression
**Script:** `validation/bet3_leave_one_out.py`

| Sub-attack | Result |
|---|---|
| Smoothness of \(e_{T-1}\) | **Not smooth.** After stripping \(p\le 8K+1\), cofactor has large primes (e.g. \(K=10\): primes \(7309\times 4546631\); \(K=25+\): cofactor \(>100\) bits) |
| Duplication \(e_{T-1}(2K)\) | Exact split into first-half + shifted-half leave-one-out — **same hardness** |
| Charpoly / \(e_{K-1}\) | Product tree; no soft-linear for single elem-sym without \(H\cdot P_0\) |
| \(v_p(e_{T-1})\) from AP counts | Cancellation possible; residues needed — exponent calculus fails for **sums** |
| Fantasy: Borwein for sums | **KILL.** \(v_p(\sum a_i)\) not determined by \(\{v_p(a_i)\}\) |

**Verdict:** structural block — product compression (exponent vectors) does not transfer to sums. Leave-one-out compression remains the name of E, with evidence it is structurally blocked for Borwein-type methods.

### Still open after bets 1–3

4. Iwasawa / Mazur measure side (beyond \(O(\log N)\) moments already killed in door β).
5. Deliberate rule-change: drop exactness for something still digit-correct (limb/approx already killed for \(O(\log N)\) — needs a *new* approximate invariant).
6. Any representation of \(e_{T-1}\) outside the bases named in Conjecture mass (the paper’s guidance).

### Bet 4 — Iwasawa / Mazur deeper
**Script:** `validation/bet4_iwasawa_mazur.py`

| Sub-attack | Result |
|---|---|
| Volkenborn of \(1/(1+8(x\bmod K))\) | Not equal to \(H\); \(H\) is a **partial** AP sum, not an integral over \(\Z_2\) |
| Extending to \(\Z_2\) | Tail = Diamond log-gamma (**cracked**) + partial = \(H\) (**barrier**) |
| Mazur/Faulhaber moments | Need \(R=\Theta(m)\) for full \(H\); mass \(\Theta(N^2)\). Extends door-β \(O(\log N)\) kill |
| Kubota–Leopoldt at fixed \(s\) | Misses moving argument \(K=\Theta(N)\) |
| Gamma-transform | Same razor as \(e_T\) vs \(e_{T-1}\) |

**Verdict:** deeper Iwasawa does not soft-solve \(H\). Bet 4 closed.

### Bet 5 — approximate digit-correct rule-change
**Script:** `validation/bet5_rulechange_approx.py`

| Sub-attack | Result |
|---|---|
| Bit-sensitivity of top-\(G\) window of \(XH\bmod 2^m\) | Bits of \(H\) down to \(\sim m-G\) flip the window; \(G=O(\log N)\) ⇒ still \(\Theta(N)\) bits of \(H\) |
| Sparse low-prec + boundary corrections | Equidistribution ⇒ \(\Theta(N/\mathrm{polylog})\) dangerous limbs without new structure |
| Soft-linear MSB predicate | Equivalent to solving E for \(G=1\) |
| New invariant hunt | No candidate satisfying (cheap)+(determines digit)+(not \(\equiv\Theta(N)\) bits of \(V\)) |

**Verdict:** approximate rule-change still blocked. Limb-support kill reconfirmed with sensitivity law. Bet 5 closed.

### Still open after bets 1–5

The attack surface is now extremely narrow:

- A representation of \(e_{T-1}\) / \(H\) **outside** Conjecture-mass bases (moments, Bernoulli/Iwasawa, inverse-power sums, symmetric functions, per-prime residue systems) that is still soft-linear-evaluable.
- A **geometric** dual jump with surplus valuation \(\Theta(m)\) (dual squaring only has \(O(1)\)).
- Something that changes what “digit extraction” means while remaining BBP-extraction (not AGM+radix) — no such invariant found.

---

## 8. Stress-test campaign (validation/stress_dual_claims.py)

Every claim above was re-tested at larger sizes and against sharper predictions.
No verdict overturned. Three claims sharpened.

### 8a. Sharp valuation law (upgrade of the dual squaring lemma)

\[
v_2\bigl(B(m,K)^2 - P_0(K)\,B(m+1,K)\bigr) \;=\; 2m + v_2(K)\quad\text{exactly.}
\]

Checked: 4480 pairs \(6\le m\le 40\), \(1\le K\le 128\) — zero exceptions;
spot checks to \(m=1000\), \(K=2m\); \(K=2^j\) to \(j=12\). Reason: the
correction sum \(C=\sum 1/(\mu(\mu+2R))\equiv Q_2=\sum\mu^{-2}\pmod{2^{m+1}}\),
and \(\mu^{-2}\equiv 1\pmod{16}\) gives \(Q_2\equiv K\pmod{16}\) with all
Faulhaber corrections of strictly larger valuation, so \(v_2(C)=v_2(K)\).

### 8b. Jump margin law (upgrade of the geometric-jump kill)

\[
v_2\bigl(B(m)^{2^s} - P_0^{2^s-1}B(m+s)\bigr) \;=\; 2m + (s-1) + v_2(K),
\]
so a full-width jump of \(s\) levels works **iff** \(s\le v_2(K)-1\)
(54 cases, zero exceptions). The old statement "fails for \(s\ge 3\)" was
K-specific. Since \(K\) is paddable to high 2-valuation (correct \(H\) with
\(O(2^j)\) extra terms), arbitrarily long full-width jumps are available —
**but they don't help**: a \(2^s\)-th power costs \(s\) squarings, so the
multiplication count from \(m_0\) to \(M\) stays \(\Theta(M)\) regardless of
jump size. Rate law intact; the only escape would be \(B(m)\to B(2m)\) in
\(O(1)\) ring ops, and the widened relation hunt (8d) found none.

### 8c. The one-step lift is free to repair; the chain is the wall

Bet1-C's "1–2 bits short" is exactly \(2-v_2(K)\) bits short, and repairable
for free: the needed correction is \(C\equiv K\pmod{16}\) — no computation.
Each further Faulhaber term of \(Q_2=\sum_j(j+1)(-8)^jS_j(K)\) buys ~3 more
bits of surplus (measured: depth 1/2/4/8 → surplus ≈ +5/+7/+13/+25 over
\(2m\)). What is blocked is the **chain**: truncation error at level \(m'\)
gets squared \(M-1-m'\) times, so cheap corrections give final precision
\(\approx M + \text{const}\) with slope exactly 1 in \(M\) (measured
\(M=20/26/32/40\) → 35/41/47/55 at depth 1), never \(2M\). With first-order
exact corrections the chain achieves exactly \(M+3m_0-1\) (measured 49 at
\(m_0=8,M=26\)) — the residue is the second-order cross terms
\(\sum_{i<j}x_ix_j\), i.e. more of the same tower. Full width needs
correction depth \(\Theta(M-m')\) at each level = Bernoulli/Faulhaber mass =
**E3 wall**. Also noted: \(B(m_0,K)=P_0(K+2^{m_0-3})/P_0(2^{m_0-3})\)
exactly, so the ladder can start Borwein-cheap at \(m_0\approx\log_2 K\)
(equivalently: \(B(m,K)\) is a short segment of the 1-mod-8 product at
offset \(2^{m-3}\) — exponential sieve position, which is why Borwein dies).

### 8d. Widened relation hunt: no geometric jump in monomials

Meet-in-the-middle over 7 atoms \(\{P_0, B(m), B(m{+}1), B(m{+}2), B(2m),
B(m,2K), B(m,K/2)\}\), exponents \(\pm4\), two parameter pairs
\((m,K)=(24,24),(28,20)\) — versus the old hunt's 5 atoms, \(\pm2\), one pair.
At width \(2m-4\): 158 hits, all in the lemma family / trivial
\(B(2m)\equiv P_0\) class (sanity pass). At width \(3m\): **zero** nontrivial
relations. No monomial geometric jump exists in this atom set.

### 8e. Non-smoothness of \(e_{T-1}\) is decisive

Extended bet3-A from \(K\le 30\) (strip \(8K{+}1\)) to \(K\le 64\) (strip
\(10^6\), Miller–Rabin): at \(K=40\) the cofactor is a **161-bit prime**; at
\(K=64\) a 290-bit cofactor with no prime factor below \(10^6\). Borwein
exponent assembly for \(e_{T-1}\) is dead at any poly smoothness bound, not
just \(8K\).

### 8f. Paper items verified at scale

Prop. E-equiv (ii) \(B\equiv P_0(1+2^mH)\pmod{2^{2m}}\) at
\((m,K)=(64,64),(128,101),(256,300)\): exact. (iii) geometric expansion
\(H=\sum_j(-8)^jS_j\) converges at \(j\approx m/3\) as claimed. Suggested
paper edits: state the lemma with the sharp equality \(2m+v_2(K)\); replace
"surplus is \(O(1)\)" with the depth-accounting of 8c; make the
geometric-jump sentence the margin law of 8b.

---

## 9. Bets 6–7 — the flank and the ledger

### Bet 6 — tree/FFT flank: ALREADY CLOSED (cross-reference)

Proposed: shave the third log inside the multiplication layer of the V merge
tree (transform caching, scaled/truncated trees, Borwein-assisted P-sides).
Found on inspection: **rounds 1–2 of the transform-reuse attack
(2026-08-01, `merge_levels.py` / `merge_reuse.py`) already closed this.**
The decomposition is definitive:

\[
N\log^3 = (\text{leaf mass } T\log N)\times(\text{depth }\log T)\times(\text{FFT log per multiply}),
\]

within-merge caching is 9→7 transforms, wrapped (mod \(2^{\text{cap}}\))
multiplication halves lengths where the cap binds — both constants. Free
(linear) multiplication would land exactly at \(N\log^2\) = AGM parity, so
the prize lives in leaf mass or depth: conjecture-mass territory, not FFT
engineering. Additional checks this campaign: truncated trees and graded
Montgomery batch inversion re-derived — constants and \(N^2\) respectively;
Borwein-compressed P-sides don't cheapen the S×P mass.

### Bet 7 — the doubling ledger (NEW: closes the 2-adic-AGM door)

**Script:** `validation/bet7_doubling_ledger.py`

The archimedean world crosses this exact one-log gap (Machin splitting
\(N\log^3\) vs AGM \(N\log^2\)) with a quadratically convergent iteration
carrying \(O(1)\) state. Bet 7 asks whether any 2-adic analogue can exist
for the extraction atom, and answers by exact accounting:

| Fact | Measured |
|---|---|
| Naive doubling weak | \(v_2(B(w)^2-P_0B(2w)) = w+1+v_2(K)\) exactly (12/12) — longer jumps DECAY (unit step has \(2w\)) |
| Exp/Log doubling relation | exact to full width \(4w\) with the full \(Q\)-tower |
| H-depth threshold, jump \(w\to 2w\) | \(d^\* = 3w-1\) exactly (4/4 cases) |
| General jump \(w\to cw\) | \(d^\* = (2c-1)w - O(1)\) for \(c=2,3,4\) |
| Subsidy law | holding \(w\) bits of \(H\), a \(c\)-fold jump emits \(cw\) but consumes \((2c-1)w\): **every new certified bit of \(H\) costs 2 external bits of \(H\)** (measured 1.90–1.93, uniform in \(c\)) |

Even given \(B(w)\) **exactly** and the whole higher tower \(Q_2,Q_3,\ldots\)
**exactly**, the doubled state cannot be reached without externally supplying
twice the H-information gained. No state-bounded quadratically-convergent
iteration (AGM, Landen, Newton-on-\(B\), bootstrapping of any arity) can
exist: the door α family is closed by conservation, not by failure to find
the right iteration. This is the 2-adic explanation of why the gap is
exactly one log on both sides of the completion.

### Doors re-derived and collapsed during bets 6–7 (for the record)

Independent re-derivations that each landed on an existing kill — evidence
the kill-map is complete over the standard idea space:

- Partial 2-adic Hurwitz interpolation \(H=F(-1)\), \(F(n)=\sum\mu^n\):
  circular (powering \(\mu^{2^m-1}\) IS the inversion).
- Character/coset expansion of the cutoff: interval \([1,8K{+}1]\) is not a
  union of few cosets; complete coset sums are trivially closed
  (inversion permutes cosets ⇒ arithmetic series), but the window is
  exponentially short of completeness. One phenomenon: **window ≪ modulus**.
- Weighted Faulhaber/Eulerian tower for \(V\) (weights \(2^{-ct}\)): same
  3-bits-per-term convergence, same tower, same tree.
- p=2 block self-similarity (A395199 face): the shifted-block correction
  **is** the dual product — the squaring lift works precisely where the dual
  regime is trivial.
- 2-adic convergence acceleration (Euler/binomial transform): valuations are
  discrete; \(x/(1-x)\) preserves \(v_2\) — no rate gain available.

### Live tickets after bets 1–7 (the honest shortlist)

1. **Bet 8 (positive kill):** beat BGS \(\widetilde O(\sqrt K)\) for the
   AP-factorial segment \(B(m,K)\bmod 2^{2m}\) (= E5) using p=2-specific
   structure. Note: \(B(m,K)=P_0(K+2^{m-3})/P_0(2^{m-3})\) exactly — a
   short window of the 1-mod-8 product at exponential offset; completeness
   (Wilson closures) dies exactly at window ≪ modulus. Recognized open
   frontier; any advance is a standalone result.
2. **Bet 9 (negative kill):** prove Conjecture mass unrestricted.
   Proof-lead: the **1-unit Log isomorphism** — mod \(2^m\), the 1-unit
   group is additively isomorphic via \(\mathrm{Log}\), squaring = doubling
   of logs (dual-squaring lemma = the Frobenius relation, and its ledger is
   2:1 negative). Transport any ring-ops algorithm toward an additive-only
   model where mass counting could become a theorem. Combined with the
   merge decomposition and the subsidy law, all three walls look like faces
   of one conservation phenomenon.

---

## 10. Bet 9 — EXECUTED: the conservation theorems

**Script:** `validation/bet9_formal_ledger.py`
**Paper:** new block after Conjecture mass — Lemma unit-log flattening
(`lem:flatten`), Theorem scale conservation (`thm:conservation`),
Corollary no-lift/no-2-adic-AGM (`cor:no-lift`); Lemma dual-square upgraded
to the exact law \(2m+v_2(K)\); jump law added; Prop E-equiv (iv) fixed
(batch → single Bernoulli shadow; the batch has \(\Theta(m^2)\) output bits
and could never be soft-linear); pointer added in the structural-exclusions
remark. Compiles clean.

### What is now THEOREM (with proofs in the paper)

1. **Flattening (T1).** \(\mathrm{Log}:(1+8\Z_2,\times)\to(8\Z_2,+)\) is a
   valuation-isometric isomorphism; any SLP with 1-unit intermediates is
   step-for-step an additive SLP at identical widths. The isomorphism that
   cracked \(E\) cuts both ways: multiplicative tricks cannot beat additive
   accounting mod \(2^W\); the sum-of-products core is exactly what Log
   does not flatten (the razor, in valuation form).
2. **Scale conservation (T2).** Holding atoms \(B(w_i)\) at precision
   \(p_i\), plus \(P_0\) and the whole tower \(Q_{\ge2}\) exact and FREE,
   every ring combination certifies \(H\) to at most
   \(\max_i(p_i-w_i)\) bits. Five-line valuation proof. No combination of
   scales beats the single best atom.
3. **Corollary (no lift, no 2-adic AGM).** Scale-side ring ops are
   precision-neutral for \(H\); jump \(w\to cw\) emits \(cw\), consumes
   \((2c-1)w-O(1)\): two external bits per certified bit, uniformly in
   \(c\). All precision is bought on the range axis, where merge-mass
   accounting rules. **Conjecture mass is a theorem in the two-axis
   (range × scale) model.** Open: existence of a third axis.

### Machine checks (all pass)

- T1: 200 random pairs at W=256, zero violations of isometry/homomorphism.
- T2 adversarial: 2000 random Z-combinations × 3 precision profiles,
  worst-case staggered error draws: bound attained, never exceeded.
- Methodological note that matters: the first run of the T2 check
  "violated" the bound by measuring the LUCKIEST error draw instead of the
  guaranteed (worst-case-certified) precision, and a second subtlety hides
  in parity locks between equal-valuation noise terms (the \(O(1)\) slack
  of the ledger). Both are now handled; both belong in any referee-proof
  statement: the theorem bounds what an algorithm can CERTIFY from its
  precision promises, not what a lucky error can coincide with.

### Honest scope

Covered: all algorithms whose intermediates stay in the 1-unit closure
(products/quotients/powers of atoms — every lift, AGM, Landen, Newton,
bootstrapping shape ever proposed for E), combined with any
range-decomposition merge. Not covered: essential use of non-unit
intermediates, bit-level (non-ring) operations, or a genuinely new
decomposition axis. That residue is where Open Problem E now lives.

---

## 11. arXiv-prep audit (2026-08-05)

Independent re-verification battery: `validation/arxiv_recheck.py` — written
from the paper's statements, not from the campaign scripts. All 21 checks
pass, including: dual-square sharp law on 1964 pairs (now covering m=5 and
m=1..4, where it also holds — the stated m≥5 is conservative), jump law on
245 triples, v2(Q2−K)>v2(K) to K=1024, E-equiv (ii) to (m,K)=(256,300),
Log isometry/homomorphism, conservation bound on 600 adversarial combos,
Euler starters for all odd primes < 20000, 2-adic EM at K=12 with exact
von-Staudt handling, and end-to-end digits vs mpmath at N=10^3/10^4 (π) and
10^3 (ln 2).

Paper fixes made for submission:

1. **Prop E-equiv restructured.** Old item (iii) (power-sum batch) had the
   same Θ(m²)-output-bits defect the log's §10 caught for (iv) — a batch can
   never be soft-linear, so the claimed two-way reduction was false as
   stated. Now: (i) H, (ii) dual product, (iii) Bernoulli shadow (single
   value, absorbing the geometric expansion), (iv) digamma **difference**
   ψ₂(K+1+1/8) − ψ₂(1/8) (a lone ψ₂ value is only equivalent given ψ₂(1/8)
   to m bits, itself Bernoulli-mass). Preamble now gives the ε-coefficient
   form of V's numerator (weighted leave-one-out) and scopes the unweighted
   distillation honestly.
2. **Dual-square lemma upgraded to a full proof** via von Staudt–Clausen:
   (j+1)S_j = Σ binom(j+1,i)B⁺_i K^{j+1−i}, v2(B⁺_i) ≥ −1, every monomial
   carries K, so the j-th correction has v2 ≥ 3j−1+v2(K) ≥ v2(K)+2.
3. **Jump law now proved** (induction/telescoping, i=0 factor dominates
   alone), not just machine-checked.
4. **Cor no-lift scoped**: "1-unit closure" qualifier added to the iteration
   exclusion; final claim now "mass law holds for the atom family of Prop
   E-equiv" (T2 is proven for H's tower, not for weighted V), with the
   third-axis residue spelled out.
5. Abstract sentence added for the conservation theorem; \pdfoutput=1;
   date → August 2026; overfull boxes fixed. Compiles clean, 22 pp.

Citations verified against sources: Zudilin v1 (2024-09-16, titled "…in
base 10") → v2 (2024-09-17, retitled base 5) exactly as footnoted;
Kristensen–Mathiasen 2208.07696; Cloitre 2508.03743; CKR II published Res.
Number Theory (Springer, 10.1007/s40993-024-00593-8).

Study companion: `notes/walkthrough.tex/pdf` (18 pp) — from-scratch
rederivation notes, private, not part of the arXiv package.

---

## 12. Door δ (Bézout / elimination) — 2026-08-07

**Script:** `validation/bezout_degeneracy.py` (all checks pass)

Entered from the question "what is the topology of these digits". Set
\(f(x)=\prod_{k\le K}(x+\mu_k)\). Then \(f(0)=P_0\) is the cracked side and
\(f'(0)=e_{K-1}=P_0H\) is the barrier, so **\(H\) is the logarithmic
derivative of \(f\) at \(0\)** and the \(e_T\) vs \(e_{T-1}\) razor is
"the curve is cheap, its polar is not". Elimination theory offers exactly
the missing step: \(f\) is squarefree, so \(uf+vf'=\operatorname{Res}(f,f')
=\operatorname{disc}f\), and evaluating at \(x=0\) solves for \(H\) from
\(f\)-side data. Both normalizations are 2-adically vacuous.

| Sub-attack | Result |
|---|---|
| Integral certificate \(uf+vf'=\operatorname{disc}f\) | **KILL.** \(v_2(\operatorname{disc}f)=3K(K-1)+2\sum_{d<K}(K-d)v_2(d)\approx 5K^2\) (verified directly, \(K=2..25\)). At \(m=\Theta(K)\) the RHS is \(\equiv 0\): the certificate carries no information about \(H\) |
| Normalized certificate \(uf+vf'=1\) | **KILL.** Degeneracy moves into the cofactors: \(v_2(v(0))\approx-7.45K\), a 2-adic **pole** ~6x deeper than the whole target width. Reducing mod \(2^m\) needs more guard bits than \(H\) has bits, and those guard bits are \(H\) |
| Cofactor output size | \(\Theta(Km)=\Theta(N^2)\) bits (measured). Same wall as the \(e_j\) batch |
| Bezoutian (Hermite) | \(\det\operatorname{Bez}(f,f')=\pm\operatorname{disc}=\det\operatorname{Hankel}(p_{i+j})\), verified \(K=2..8\). **Bézout's own matrix construction is the batch power sums \(S_j\)** — the repaired route is E2. Superfast Hankel solvers are \(O(\mathsf M(n)\log n)\) but need the entries, and the entries are the barrier |
| Extension rigidity | **KILL.** Unramified \(W(\mathbb F_{2^s})\): the \(\mu_k\) are rational, Frobenius fixes them, all still reduce to one point (verified in \(W(\mathbb F_4)\), distance spectra identical). Ramified: \(v\) rescales by \(e\) uniformly. The configuration's tree is invariant under **every** finite extension of \(\mathbb Q_2\) |

**Root cause (the topological statement).** \(|\mathbb F_2^\times|=1\). Every
odd number reduces to the same point of the special fibre, so *any* set of
units mod \(2^m\) is a maximally degenerate configuration: all \(\binom K2\)
intersections of \(f\) with its polar pile up at one point. Bézout's theorem
still counts them, it just counts them all in the same place. At \(p=2\)
there is no other kind of instance — the AP is not a hard case of the
inverse-sum problem, it is the only case. Combined with the rigidity row,
this is why "go to a bigger residue field" is unavailable.

**Verdict:** door δ closed. Another independent re-derivation landing on the
razor (§9 pattern), now with a geometric reason for the razor rather than
just a valuation count.

---

## 13. Bet 10 — window rigidity: the index-decomposition axis is CLOSED

**Script:** `validation/bet10_window_rigidity.py` (all checks pass)

**The question.** Every merge scheme in this project decomposes \(\{1..K\}\)
into a *nested* family: a partition refinement, i.e. a tree, i.e. a cover
whose nerve is 0-dimensional. The merge-mass accounting (leaf mass × depth)
assumes that nesting. Topology does not require it. Can a non-nested
(overlapping / signed) index decomposition beat a partition?

**First answer — the naive form is empty.** \(H(S)=\sum_{k\in S}\mu_k^{-1}\)
is finitely additive, so the Čech obstruction vanishes identically:
inclusion–exclusion over any cover is exact, and therefore never obstructed
and never cheaper. Measured (check F): partition mass \(=K=40\); a 5-block
overlapping cover costs 188; heavy overlap costs 2322. All exact, all worse.

**Second answer — the real form, and the theorem.** The saving would have to
come from the *blocks*: use blocks whose local sum has a closed form (cost
\(O(1)\) rather than proportional to size) and write the short window as a
signed combination of those. This is the "virtual cover", strictly more
general than the union-of-cosets question §9 already noted as dead.

The available blocks are pinned exactly. \(1+8\Z_2 \bmod 2^m\) is **cyclic**
of order \(M=2^{m-3}\) (generator 9, verified \(m=6..16\)), so its only
subgroups are the filtration \(G_j=1+2^j\Z_2\); inversion permutes any
subgroup, so a complete coset sums in closed form (verified); and in index
coordinates a coset of \(G_j\) is an AP of step \(2^{j-3}\) and length
\(2^{m-j}\) — **index span exactly \(M\), whatever the step** (verified).

> **Theorem (window rigidity).** Let \(M=2^{m-3}\) and \(0<K<M/2\). Suppose
> \[
> \mathbf 1_{[1,K]}=\sum_b c_b\mathbf 1_{B_b}+\sum_{k\in S}d_k\mathbf 1_{\{k\}}
> \]
> with \(B_b\) any inversion-closed blocks (any steps, any offsets, any
> rational \(c_b\), any number of terms) and \(S\) finite. Then \(|S|\ge K\).
>
> *Proof.* Pass to generating polynomials in \(\Z[z]\) and reduce mod
> \(\Phi_M(z)=z^{M/2}+1\). A closed block of step \(s\) has length
> \(M/s\ge2\), so \(s\mid M/2\), and its elements pair off as
> \(x\leftrightarrow x+M/2\) inside the block: its image is \(0\). The window
> has degree \(K<M/2\), so its image is itself, with \(K\) nonzero
> coefficients. Folding never increases support, so \(|S|\ge K\). ∎

Machine-verified on **all** \(2^{m-3}-1\) closed blocks for \(m=7..14\).

**Equivalent form.** The \(\Z\)-span of all closed blocks is exactly the
\(M/2\)-periodic functions, dimension \(M/2\); the step-\((M/2)\) blocks
\(\{a,a+M/2\}\) are already a basis and every other closed block is a sum of
those (verified). Closed blocks see only the \(M/2\)-periodic part of an
index set, and a window of length \(K<M/2\) has none.

**Consequence.** Complete cosets contribute **exactly nothing** to a short
window: any exact representation still needs \(K\) singletons, i.e. \(K\)
separate inversions, i.e. the naive algorithm. No index decomposition —
nested or overlapping, partition or signed cover — beats term-by-term.

This upgrades the informal "window ≪ modulus" observation of §9 to a theorem
with a two-line proof, and it narrows §10's open residue: **the third axis,
if it exists, is not a way of cutting up the range.** It is also the
structural reason behind 8c's remark that \(B(m,K)\) is a short segment of
the 1-mod-8 product at exponential offset — the closed objects all live at
length \(M\), and \(K=\Theta(m)\) is logarithmically short of that.

**Files:** `validation/bezout_degeneracy.py`, `validation/bet10_window_rigidity.py`

---

## 14. What “done” would look like

An algorithm that, on input \(N\), outputs the \(N\)th decimal digit of \(\pi\) (BBP extraction, not AGM+radix), with bit complexity \(O(N\log^2 N)\), by computing \(V\) (or any E0–E5 atom) in \(O(\mathsf{M}(N)\log N)\).

Until then: keep the lemma, keep the kills, keep attacking the live bets.

---

## 15. External sweep + doors ε, ζ, η, θ — 2026-08-23

Four parallel reconnaissance vectors (arXiv computer algebra, exponential
sums / p-adic analysis, complexity-theory models, community sources: MO,
cstheory, OEIS, Tao/Johansson/Harvey, ANTS/ISSAC programs), everything
triaged against the kill map. Outcome: four new doors named and closed,
one new lemma, one control theorem, a stack of verified negative intel,
and two genuinely live tickets. Nothing in the outside world evades the
kill map; the problem is externally confirmed open AND unclaimed.

**New scripts:** `validation/door_eps_zeta.py`, `validation/carlitz_control.py`
(all checks pass).

### 15a. Doors ε (completion) and ζ (geometric blocks): bet 10's two loopholes, measured shut

Bet 10's rigidity theorem covers inversion-closed blocks with rational
coefficients. Two families slip its hypotheses; both are now closed.

**Door ζ.** In dlog coordinates \(u=\mathrm{dlog}_9(1+8x)\) on
\(\Z/M\), \(M=2^{m-3}\): the value is a pure geometric
\((8x+1)^{-1}=9^{-u}\), the map \(u\mapsto x\) is a 2-adic **isometry**,
and any u-AP (any odd step \(\lambda\)) is a closed-form block
(geometric series) that is NOT a coset. Measured laws:

| Fact | Measured |
|---|---|
| Exact run law | \(\mathrm{runs}_\lambda(U)=K-\lvert W\cap 9^{r}W\rvert\), \(r=\lambda^{-1}\), \(W=1+8[1,K]\) — all odd \(\lambda\), three sizes, zero exceptions |
| Best structured base | \(\lambda=1\) (base 9): \(K-\lfloor (K-1)/9\rfloor\) exactly while \(9K+1\le M\) |
| Min over ALL bases | \(\Theta(K)\), tracking the extreme-value floor of a random \(K\)-set over \(M/2\) dilations (the dip below \(8K/9\) is birthday noise) |
| Longest block, any base | \(3\)–\(4\) at \(K\le100\) = \(\log_9(8K)+O(1)\) (geometric escape: a u-AP is a \(9^{\lambda^{-1}}\)-geometric progression trapped in the value window) |

Consequence: any mixed signed cover by geometric blocks needs
\(\ge K/O(\log K)\) blocks. Best case door ζ saves one log off naive
\(\Theta(K)\) mass — \(\Theta(N^2/\log)\), far above even BGS. Shut.

**Door ε.** Fourier completion with \(2^{m-3}\)-th roots of unity
(irrational coefficients, outside bet 10's scope). Three independent kills:
(i) the folded transforms \(\hat g(a)\) are **dense at every odd
frequency**, indistinguishable from the transform of random values —
sparse only at 2-power frequencies, which are exactly the coset part bet
10 already governs (coherent picture: the closed-form-accessible spectrum
IS the \(M/2\)-periodic part, and the window has none); (ii) exact inverse
DFT divides by \(M=2^{m-3}\) — \(m-3\) guard bits, the dual-width doubling
yet again; (iii) external: complete-sum p-adic stationary phase is a
solved exact theory including \(p=2\) (Cochrane–Zheng; Gurak Salié forms),
and its collapse mechanism provably requires the complete inner range —
the incomplete literature (Postnikov → Korolev, Bourgain–Garaev) is
bounds-only and empty below length \(\exp(c(\log q)^{2/3})\); our window
\(K=\Theta(\log q)\) sits exponentially below every regime ever touched.
Bonus symmetry (not a door): \(\hat g_{8k+1}(a)=a\cdot\hat g_{8k+a}(1)\)
— frequency trades into AP residue class; permutes atoms, reduces nothing.

### 15b. Door η — exact covers / the Umans–Wang program, and the vanishing-slack lemma

The one live external program on the BGS \(\sqrt{}\)-barrier:
Umans–Wang divisor covers (arXiv:2511.10851, Nov 2025) — replace interval
products by products over structured sets \(S-T\), aiming below exponent
\(1/2\) for polynomial/integer factorization. He–Sahai (2608.06681, Aug
2026) refutes the rank-1 (AP) version — length \(\ge n^{3/4-o(1)}\) — but
rank \(\ge2\) survives. Transfer analysis to our atom:

**Lemma (vanishing slack).** For odd \(u\):
\(v_2(u^{2^t}-1)=t-1+v_2(u^2-1)\ge t+2\). Hence junk indices with
multiplicity \(\equiv0\pmod{2^{2m-2}}\) **vanish** from the dual product
mod \(2^{2m}\): a cover need only be multiplicity-exact mod
\(2^{2m-2}\), with arbitrarily huge integer multiplicities at
\(O(\log)\)-squaring cost. (Machine-checked; the \(p=2\) analogue of
UW's divisibility slack. New freedom — the exact-cover question never had it.)

**Touch accounting (the pinch).** Evaluating a rank-\(r\) sumset product
\(\prod_{s_1+\cdots+s_r}(R+8d)\) at one point pinches at the balanced
level: the binary-sumset cover \([1,2^j)=\{0,1\}+\{0,2\}+\cdots\)
re-derives BGS exactly (\(\sum_i\max(2^i,2^{j-i})=\Theta(2^{j/2})\));
generally cost \(\ge\sqrt{\text{touches}}\) per atom. Every window index
needs a unit-multiplicity touch even with slack, so
\(\sum t_i\ge K\) and total cost \(\ge\sqrt K\) by concavity. **In the
cover model — any rank, signs, congruence slack, plus coset/geometric
side-blocks — BGS is optimal.** UW's sub-\(\sqrt{}\) escape uses
divisibility slack (junk divisors don't hurt a gcd) that exact-value
windows don't have. Reading the rank-2 mechanism in full stays a ticket
(bet 8's external anchor); the pinch analysis predicts a kill.

Community re-derivations landing on existing kills (§9 pattern, three
more): dyadic-\(K\) rounding (the difference window is the same atom);
mod-8 class-product doubling from OEIS A395199 (Jianing Song, Apr 2026:
\((2^{k+1}-1)!!\equiv((2^k-1)!!)^2\bmod 2^{3k-1}\), exact — complete
windows enjoy a *growing* surplus, the cleanest external display of the
complete/window contrast; class products are coset products = cracked
side); Catalan/sparse-binomial reformulation (ratios of complete windows).

### 15c. Door θ — coefficient extraction / Graeffe (the non-unit linear axis)

\(H(K)=[y^K]\,G(y)/(1-y)\) with \((8yD+1)G=y/(1-y)\): \(G\) is holonomic
of order 2, so the atom is "\(K\)-th coefficient of an explicit holonomic
series mod \(2^{2m}\)" — P-recursive term extraction, i.e. BGS again. The
2023–26 near-linear composition wave (Kinoshita–Li FOCS 2024 via
Bostan–Mori Graeffe; NSSV JACM 2023, \(n^{1.343}\) Jan 2026) is
constant-coefficient/rational or amortized-over-output machinery;
P-recursive halving doubles coefficient-degree per level and re-pinches
at \(\sqrt K\) (balanced level, same as 15b). \(\Z/2^{2m}\) is not
characteristic 2, so Graeffe runs verbatim — and changes nothing. Verdict:
the non-unit *linear* axis (residue (i) of Cor no-lift) collapses onto the
pinch. Formalizing "linear-algebraic reorganizations with non-unit
intermediates obey the same mass bound" as an extension of T2 is a paper
ticket, not an algorithm.

### 15d. The function-field control: E lives on the archimedean/2-adic seam

**Script:** `validation/carlitz_control.py`. Move the atom to
\(\F_2[T]\) (place \(2\to T\), \(8\to T^3\), window \(\to\) the ball
\(V=\{\deg a<d\}\), \(|V|=2^d=K\)):

\[
H_{FF}=\sum_{\deg a<d}(1+T^3a)^{-1}\bmod T^{mm}
=T^{-3}\,c_0/\ell_V(T^{-3}),
\]

because \(\ell_V(z)=\prod_{a\in V}(z+a)\) is a Carlitz **additive
polynomial** — \(d+1\) monomials \(z^{2^i}\), built by \(d\) Ore steps —
and in char 2 the derivative of an additive polynomial is the constant
\(c_0\), so the log-derivative identity closes the whole sum. Measured:
exact match mod \(T^{64}\), **223 ring ops vs 2048 brute inversions at
\(d=11\)** — polylog in \(K\). The function-field E is trivial: there the
window IS a ball, and there is no other kind of bounded window. Over
\(\Z\) the interval is not a ball and Frobenius is not \(\Z\)-linear;
E cannot even be stated over \(\F_q[T]\). Independent echoes: Igusa /
cell-decomposition closed forms stop exactly at the definable = cell =
coset boundary (Cluckers; Zuniga-Galindo's poly-time univariate Igusa
zeta), and the interval is not definable in Macintyre's language; the
exact Lerch/Glaisher incomplete-sum evaluations
(\(H_{[p/2]}\equiv-2q_p(2)\) etc.) exist only at windows that are fixed
*fractions* of the modulus — completion-compatible ones. Model theory and
function fields independently reproduce window rigidity.

### 15e. Verified negative intel (state of the art, checked 2026-08)

| Front | Status |
|---|---|
| P-recursive \(N\)-th term | \(\sqrt N\) (Chudnovsky '88 / BGS '07) unmoved; Bostan–Mori is constant-coefficient only; ANTS 2024/2026 + ISSAC 2025 programs: zero relevant titles |
| Isolated modular values | Johansson (MO 444893, 2023): nothing below \(n^{1/2+o(1)}\) for \(p(n)\bmod\) small primes — the expert-acknowledged frontier |
| p-adic elementary functions | Even \(p\)-adic log sits at \(n\log^3 n\) with "faster?" open (MO 423035, Johansson/Kulkov). Our \(\log^3\) pipeline matches the class-wide wall |
| Lower bounds | No unconditional superlinear bound even for multiplication (network-coding-conditional only, AFKL '19); online mult \(\Omega(n\log n)\) is cell-probe/online only; Farach-Colton–Tsai: \(\Omega(\M(n\log^{4/7-\epsilon}n))\) for \(n!\) in a mult-only model — restricted-model precedent for T2 |
| Catalytic / TreeEval / Williams | All space, never time (CL ⊆ ZPP); Cook–Mertz register updates are RING ops — inside T2's accounting; Williams STOC '25 gives our pipeline \(O(\sqrt{N}\log^2 N)\) *space* with no time guarantee (Gourdon '97 already had polylog space at \(\sim n^2/\log^2 n\) time); the arXiv √t-without-log claim was withdrawn Jan 2026 |
| Quantum | Gauss/Kloosterman sums: *estimation* at cost \(1/\epsilon\) (van Dam–Seroussi) — we need \(2^{-m}\), exponential; no quantum P-recursive/mod-factorial algorithm exists (searched); the Jul 2026 "quantum modular factorials" paper states the classical sub-\(\sqrt{}\) problem open |
| Factoring shield | Strassen/Lipton/Cheng/Shub–Smale reductions (fast factorial-type products ⇒ fast factoring) are all **vacuous at fixed modulus \(2^{2m}\)** — nobody has studied the fixed-smooth-modulus case: E is unprotected by crypto hardness and uncharted |
| The bet-7 question in public | MO 476452 (joro, Aug 2024): doubling formulas for P-recursive sequences over finite rings — **zero answers in two years**. The territory is unclaimed |
| Base-10 formula shopping | Zudilin 2409.10097: v1 "base 10" (Sep 16 2024) self-retracted to "base 5 … we joke about" in v2 within 24h; §3 documents the failure (the mod-\(5^d\) reduction identity dies once \(2n+1>d-k\)); zero citations, no rescue in 23 months. Cloitre 2025: base 16 / base φ, wrong constant or base, no cost analysis |

### 15f. Live tickets after the sweep

1. **Bet 8** (unchanged, now externally anchored): beat BGS for E5 with
   \(p=2\) structure. The UW rank-2 cover mechanism is the one active
   outside program; read it in full and run the transfer against the 15b
   pinch. Expected: a named kill, upgradeable to a clean model theorem
   ("BGS is optimal in the slack-cover model" — the vanishing-slack lemma
   makes the model statement honest at \(p=2\)).
2. **Bet 11 (NEW): the p-adic Dedekind cascade.** The only
   polylog-computable modular-inverse sums in the literature are the
   Dedekind/cotangent/Vasyunin family — exact reciprocity, Euclidean
   cascade, \(O(\log c)\) algorithms (word-rewriting, arXiv:2210.01172),
   with Bettin–Conrey's computable period-function correction. p-adically
   VALUED Dedekind–Apostol sums exist (Rosen–Snyder; Snyder; Kudo — built
   on Washington's p-adic Hurwitz zeta) and satisfy exact p-adic
   reciprocity — and **nobody has ever asked whether the reciprocity
   supports an evaluation cascade for the p-adic values.** First
   experiment: implement at \(p=2\), small \(m\); verify reciprocity
   numerically; execute one Euclid step as an evaluation move and measure
   the 2-adic precision cost of the correction term. Known obstruction
   shape to respect: every archimedean Dedekind algorithm shrinks modulus
   and window TOGETHER through \(SL_2(\Z)\) acting on \(h/k\), while E
   fixes the modulus — if the p-adic correction is precision-conserving,
   it is the first exact 2-adic reciprocity cascade on record, and the
   specialization question toward \(H\) opens; if precision halves per
   step, it is the ledger again and door-closes the last known family.
3. **Paper remarks pack** (cite-ready, from the sweep): (i) Williams
   √-space corollary + the honest frontier statement (simultaneous
   near-linear time and sublinear space is the open profile); (ii) lemma:
   catalytic register programs are ring-op programs, hence inside T2 —
   closes residue axis (ii) against all known bit-level machinery;
   (iii) the fixed-smooth-modulus hardness vacuum (own the question);
   (iv) the function-field control remark (15d); (v) Sawin/Mahler
   irrationality-measure carry bound (MO 433675) for the rigor section;
   (vi) Harvey–van der Hoeven FOCS 2025 (multiplication ≥ transposition)
   for the \(\M(N)\) discussion.
4. **Community moves** (user's call): post the atom to MO/cstheory
   ("\(\sum(8k+1)^{-1}\bmod 2^m\) below \(\sqrt K\)?") — Johansson and
   Jeřábek answer in exactly this zone and the 476452 silence says it's
   unclaimed; watch Bostan–Yurkevich for matrix-factorial successors to
   "beating binary powering" (2302.04299); pull Harvey's ICM 2026
   accumulating-remainder-tree survey when public; pull Jianing Song's
   A069954 proof PDF and check whether the \(3k-1\) argument
   character-twists to the \(\langle 5^2\rangle\) quarter-group.

**Status after the sweep.** Not cracked. The barrier now has four more
named faces (ε, ζ, η, θ), a control theorem pinning it to the
archimedean/ultrametric seam, and external confirmation on every front
that no known technique evades the kill map — including the fact that
the nearest open questions sit unanswered in public. The two live
directions are bet 8 (UW rank-2 transfer) and bet 11 (p-adic Dedekind
cascade); everything else the outside world offers is already inside the
fence.

---

## 16. Bets 11 and 8-external driven to verdict; door ι — 2026-08-23 (solve-it session)

Same-day continuation: both live tickets of §15f executed to closure,
plus one more classic door checked and fenced. **New scripts:**
`validation/bet11_dedekind_bridge.py`, `validation/door_iota_adjoint.py`
(all checks pass).

### 16a. Bet 11 EXECUTED — the Dedekind/reciprocity family is closed

**Script:** `validation/bet11_dedekind_bridge.py`

Three-step structural closure, each step now backed by measurement:

1. **The cascade transports only complete sums.** Every state of the
   reciprocity recursion is a complete weighted sum
   \(\sum_{a\bmod c} w(a/c)\); reciprocity shrinks modulus and window
   TOGETHER. Verified exact (cascade == brute on coprime pairs;
   reciprocity exact; \(O(\log)\) steps including \(s(h,2^t)\) at
   \(t=512\)).
2. **2-power-cusp values are cracked-side.** \(6\cdot2^s\,s(h,2^s)\bmod
   2^s\) reduces to small rational constants (quadratic polynomial in
   \(2^s\); mod \(2^s\) only the constant survives) — weighted complete
   sums, window-blind by rigidity + the §9 weighted-tower kill.
3. **Odd-cusp values are information-thin.** Each \(s_p(h,q)\) is an
   \(O(\log q)\)-bit rational; no fixed number spans \(\Theta(m)\) bits
   of \(H\); a \(\Theta(K)\)-family is per-\(k\) residue data feeding
   the same merge tree.

**The bridge hunt** (lattice reduction, pool of 10 structured couplings
\(\{s(2^s\bmod q^*,q^*), s(q^*,2^s), s((q^*)^{-1}\bmod 2^s,2^s),\)
Apostol-3, neighbors\(\}\) + \(H\) + \(P_0\) + 1; validated by two
witnesses the reducer had to find, and did: the true identity
\(H+8\sum k q_k^{-1}-K=0\) at height 8, and a planted fake):

| Measurement | Result |
|---|---|
| First run | Rediscovered the classical theorem \(s(\bar h,k)=s(h,k)\) from raw residues at height 1 (machinery validation; duplicate atoms then excluded) |
| Grid \((K,s)\in\{(12,512),(16,512),(16,768)\}\) | **27/27** small null vectors have \(c_H=0\) — all are internal small-rational identities among the information-thin atoms |
| \(H\) | Untouched at height bound \(2^{20}\), random floor \(2^{44}\)–\(2^{64}\) |

**p-adic versions closed by the same argument.** Rosen–Snyder (J. reine
angew. Math. 361, 1985), Snyder (BAMS 37, 1988), Kudo (Nagoya 144,
1996); construction recapped in Kozuka 2023 (open access): Koblitz
measures \(M_u\), Gamma-transforms \(\int_{\Z_p^\times}\langle
x\rangle^s dM_u\), interpolation \(S_{p,m}(n{-}1)=\tilde n
k^n(s_{m,n}(h,k)-\psi^{-1}(p)p^n s_{m,n}(p'h,k))\), \(p'=p^{-1}\bmod
k\). This is weight-direction interpolation (Iwasawa/Mazur machinery =
bet 4's dead zone) of the SAME complete Apostol sums; the argument
transform \(h\to p'h\) stays complete. Nothing incomplete exists
anywhere in the family. **Bet 11 closed.**

### 16b. Door η FINAL — the Umans–Wang engine does not transfer (mechanism in hand)

Full mechanism extracted from arXiv:2511.10851: the framework needs
(P1) *coverage* (every squarefree degree-\(n\) polynomial divides
\(\prod g_i\); integer side: every \(i\le n\) divides some \(s-t\)),
(P2) cheap interval products **mod h** (Frobenius powers
\(X^{q^u}\) along APs via modular-composition chains + the same
BSGS/multipoint core), (P3) a gcd/DDF base case. Factors are extracted
by **gcd-pruning against the target**; each \(s-t\) may carry
arbitrarily many junk divisors and huge extraneous primes;
multiplicities are irrelevant. Cost \(n^{1+\max(\alpha,\beta)}\);
at \((\alpha,\beta)=(1/3,1/3)\): exponents 4/3 and 1/6. He–Sahai
(2608.06681): rank-1 APs need length \(\ge n^{3/4-o(1)}\) (prime-band
blocks pairwise covering with \(\le1\) intersection + friendship
bound); the \(\le1\)-intersection step dies at rank 2 — their words:
"the corresponding difference contains two independent coefficients,
and the at-most-one-intersection argument does not follow."

**Transfer verdict.** Open Problem E has no gcd-pruning semantics: the
output IS the exact value of the prescribed window product. The UW
slack (over-coverage, junk divisors, multiplicity-blindness) is a
factoring resource, unavailable for exact evaluation. The slacks E does
have — quotient-junk (itself a window-type product: circular) and the
§15b vanishing-slack multiplicities — are floored at \(\sqrt K\) by the
touch accounting. **Door η closed for E**, precisely scoped: this
closes the *transfer*; it does not close bet 8 (non-cover \(p=2\)
mechanisms), nor UW's own rank-\(\ge2\) conjecture for factoring.
Echo for the collection: UW's Lemma 6.2 obstruction ("a small affine
transformation of a large quantity hitting \(0 \bmod U\)") is the
window-vs-modulus phenomenon wearing their clothes.

### 16c. Door ι — the gradient trick (Baur–Strassen), checked and shut

**Script:** `validation/door_iota_adjoint.py`

The temptation: \(e_{K-1}=\sum_k\partial(\prod x)/\partial
x_k\big|_{x=\mu}=P_0H\), and Baur–Strassen says all partials cost
\(\le5\times\) the circuit — so a soft-linear \(\mu\)-VARIABLE circuit
for \(P_0\) would crack E outright. Two independent reasons it does
not, both machine-checked:

1. **Value/polynomial gap.** The Borwein crack is a sieve on the
   specific integers \(8k+1\) (shared primes across \(k\)), taking
   \(k\) as input — not a straight-line program in the variables
   \(\mu_k\). The cheapest \(\mu\)-variable circuit for \(\prod x_k\)
   has size \(K\) (it must read its inputs), and Baur–Strassen on it
   gives the naive \(\Theta(K\M(m))\). The razor ("differentiation
   destroys factorization") is exactly this gap in circuit language.
2. **Adjoint width law.** Reverse-mode through the value tree mod
   \(2^{2m}\): forward values shrink toward the leaves (measured 10
   bits/leaf at \(K=256,m=256\)); adjoints are cofactor products —
   full width at EVERY level (measured 511 bits/node throughout;
   adj/fwd ratio 50.9 at the leaves). Reverse-mode mass =
   \(\Theta(K\,m)\). Propagating the adjoint SUM instead gives, node
   for node, \(e(T)=P_R e_L+P_L e_R\) — verified identical to the
   S×P merge, i.e. the existing \(O(N\log^3N)\) pipeline.

**Consequence.** Automatic differentiation, adjoint methods, and
transposition (Tellegen) on any value-DAG for \(P_0\) land inside the
merge-mass accounting. Door ι shut.

### 16d. Where Open Problem E stands after this session

Closed with proofs, measurements, or both: doors α–ι, bets 1–11,
completion (rational AND cyclotomic coefficients), all index
decompositions (nested, signed, coset, geometric, convolution covers of
any rank with congruence slack), all iteration/lift shapes (conservation
theorem), holonomic and coefficient-extraction routes, per-prime and
amortized residue systems, the entire Dedekind/reciprocity family
including its p-adic interpolations, gradient/adjoint/transposition
methods, catalytic register programs, quantum estimation, and the
function-field mirror (where the problem does not exist).

The live surface is now exactly two options:

1. **Bet 8, by non-cover means:** a \(p=2\)-specific evaluation
   mechanism for \(B(m,K)\bmod 2^{2m}\) below \(\widetilde O(\sqrt K)\)
   ring operations that is not an index-decomposition — nothing in the
   2026 literature contains a candidate, and the factoring reductions
   that would forbid it are vacuous at fixed smooth modulus.
2. **Prove Conjecture mass beyond the two-axis model.** Concrete next
   theorem target, now well-motivated by θ+ι: extend T2 to the
   linear/non-unit axis — two lemmas that look provable: (i)
   transposition-closure (Tellegen: any linear reorganization of a
   bilinear evaluation DAG preserves its multiplication mass), (ii) the
   adjoint width law of 16c as a formal statement. That would make the
   conjecture a theorem over ring-linear algorithms with any
   intermediates, leaving only genuinely bit-level arithmetic (carries
   as boolean functions) outside.

**Files this session:** `validation/bet11_dedekind_bridge.py`,
`validation/door_iota_adjoint.py`; §16 of this log.

---

## 17. THE PROOF SESSION — the linear-closure theorems (2026-08-23)

**Note:** `notes/linear_closure.tex` (compiles clean, 6 pp).
**Checks:** `validation/linear_closure_ledger.py` (all pass).

Option 2 of §16d executed: the mass conjecture is now a THEOREM over the
linear-closed two-axis model. Three theorems, fully proved in the note,
each machine-checked.

### 17a. Theorem 1 — the write floor (unconditional, in-model)

**Model** (merge-forest): union forests over \([K]\) carrying
\((P_S, e_S)\) jet records, any arity, exact-or-truncated grade-1
slots, arbitrary \(R\)-linear glue with coefficients from the free
multiplicative oracle (ALL products/cofactors/Borwein objects free),
faithful representation. Cost axiom: pay the bits you read/write.

**Theorem.** Any in-model computation of \(e_{K-1}\bmod 2^W\) costs
\(\ge \tfrac{1}{16} bK\log_2 T\), \(T=\min(W/b, K)\) — i.e.
\(\Omega(N\log^2 N)\) in the digit regime.

Proof skeleton (all steps in the note): (i) **unit transfer** — every
leaf-to-root transfer coefficient is forced to be a unit because the
target's cofactors are odd and \(v_2\) never decreases under ring ops
(2 is not invertible: the 2-adic side is load-bearing for the lower
bound); (ii) **truncation floor** — a slot truncated below \(W\) on a
live path breaks certification (two-world argument through the unit
transfer); (iii) **node cost** \(\ge (d_v-1)\min(\beta(S_v), W)\);
(iv) **binary optimality** \((d-1)/\log_2 d \ge 1\); (v) the
**small/big counting**: maximal small subtrees pay Kraft
(\(\sum L_i\log_2 L_i\)), big nodes pay \(W\) each; minimizing over
partitions gives \(K\log_2 T + 0.91K - T\).

Machine checks: the charge depends only on block sizes, so the exact
minimum over ALL binary trees obeys a size DP — verified
\(\ge 0.5K\log_2 T\) on a \((K,T)\) grid to \(K=4096\), with random,
caterpillar, and 3-ary schedules all above the DP optimum. **Measured
exact law:** the DP optimum EQUALS the balanced-tree cost in every
tested case (\(=K(\log_2 T+1)-T\) at powers of two) — balanced is
exactly optimal, not just asymptotically.

**Consequences.** The floor \(N\log^2 N\) IS the stated target of Open
Problem E (AGM parity = the write floor of the natural model). The
floor is achieved at grade 0: Borwein's crack multiplies at write cost
via sieve/exponent bookkeeping. The champion \(N\log^3 N\) is the
balanced forest paying one multiplication-overhead log at grade 1.
**Open Problem E, restricted to the model, is exactly: can grade-1
merges run at write cost?** Bet 3 (exponent vectors cannot see sums) is
why the grade-0 mechanism does not transfer. Conditional corollary:
if per-merge multiplication costs \(\Omega(n\log n)\) (network-coding)
and transforms don't amortize across merges (bet 6's measurement), the
model floor rises to \(\Omega(N\log^3 N)\) — the pipeline is
conditionally optimal in-model.

### 17b. Theorem 2 — Tellegen closure (linear reorganizations preserve the floor)

Any linear SLP over \(R\) computing \(x\mapsto Mx\) transposes to one
computing \(M^{\mathsf T}\) with the IDENTICAL multiset of
constant-multiplications at identical widths (adjoint edge reversal;
additions shift by \(\le n+p\)). Hence the write floor is
transposition-invariant, and reverse-mode/adjoint algorithms (door ι)
obey Theorem 1 verbatim. Machine check: 5 random 40-op programs over
\(\Z/2^{64}\) with mixed unit/non-unit constants —
\(\langle Mx,y\rangle=\langle x,M^{\mathsf T}y\rangle\) and multiset
equality, all pass.

### 17c. Theorem 3 — interface ledger, information form (T2 for ARBITRARY computation)

If an algorithm's only access to the index structure is (jet values
\(J(2^{w_i})\bmod 2^{p_i}\); the tower \(Q_{\ge2}\) exact; the
multiplicative oracle), then — with NO restriction on its operations:
ring, bit-level, lookup, anything — it certifies at most
\(t^*=\max_i(p_i-w_i)\) bits of \(H\). Proof: two worlds
\((Q_1, Q_{\ge2})\) and \((Q_1+2^{t^*}, Q_{\ge2})\), regenerating
atoms through the exact identity
\(J(t)=\exp(\sum_j (-1)^{j-1}t^jQ_j/j)\): every promised value
unchanged mod \(2^{p_i}\) since \(w_i+t^*\ge p_i\). This SUBSUMES T2 +
Cor no-lift and closes residue axes (i) non-unit intermediates and
(ii) bit-level ops for the entire scale axis. Machine check: identity
exact against direct products; two worlds agree at \(t^*=32\) across a
three-atom profile; perturbation at \(t^*-1\) is detected by the
tightest atom (ledger sharp).

### 17d. The conjecture's residue after the proofs

**Now theorem:** every algorithm mixing merge-forest aggregation (any
linear glue, any transpose/adjoint form) with scale-interface
consumption costs \(\Omega(N\log^2 N)\); conditionally
\(\Omega(N\log^3 N)\).

**The honest residue, exactly two questions:**
1. *Inside the model:* multiplication at write cost for grade-1 merges
   — would land exactly on \(N\log^2 N\) = the floor = the target.
   Grade 0 does it by smoothness; sums provably resist exponent
   calculus (bet 3). This is where a POSITIVE resolution of E must
   live if it stays in-model.
2. *Outside the model:* an index-access channel that is neither
   block-aggregation of grade-1 data, nor jet evaluation at 2-adic
   scales, nor in the multiplicative span. Doors α–ι and bets 1–11
   exhaust every channel named in the literature; completeness of that
   list is the remaining conjectural content (and faithfulness /
   no-redundant-encoding is the technical condition of the model to
   keep honest in any referee exchange).

**Paper integration (pending):** the note is written as a standalone;
its Theorems 1–3 + the grade-0-tightness remark slot in after
Conjecture mass, replacing the "two-axis model" corollary with the
linear-closed version, and the AGM-parity remark gains the "floor"
interpretation. Not yet merged into `decimal_extraction_paper.tex`.

**Files this session:** `notes/linear_closure.tex` (+pdf),
`validation/linear_closure_ledger.py`; §17 of this log.

---

## 18. Prime-power compression and weighted Euler acceleration — 2026-09-06

**Status:** O(N log²N) for the unweighted harmonic atom; the full weighted V
still has the O(N log³N) bound. No claim of external novelty or a world record.

**Derivations:** `notes/prime_crt_and_euler_2026-09-06.md`.

### 18a. The harmonic atom has a prime-indexed algorithm

Set B=8K+1 and q_p=p^floor(log_p B), for odd primes p<=B. For mu_k=8k+1,
compute the word-size residues

    c_p = sum_{k:p|mu_k} (q_p/p^v_p(mu_k))
          * (mu_k/p^v_p(mu_k))^(-1) mod q_p.

Then H - sum_p c_p/q_p is an INTEGER. This follows by multiplying by the
odd part D of lcm(1,...,B) and reducing modulo each q_p. The integer is
recoverable by certified fixed-point rounding using only O(log K)-bit
operands. The prime incidence sieve has O(K log log K) entries, and the
sum tree now has O(K) total denominator bits instead of O(K log K).

With fast word inversions the bit cost is
O(M(K) log K + K M(log K)(log log K)^2), hence O(K log²K).
The exact proof and rounding bound are in the note.

`validation/harmonic_prime_crt.py` passes 250 exact Fraction cases, product-tree
comparisons through a million terms, and tests of the higher-moment congruence
below. Using the existing fast P_0 also gives the unweighted dual product
B(m,K)=P_0(1+2^m H) mod 2^(2m) within the target bound.

The old LCM merge-tree experiment does not cover this reorganization over
primes. Bet 1's claim that K log K is Theta(N²), when K=Theta(N), was also
incorrect. Large prime factors in an unreduced sum numerator do not prevent
this algorithm.

### 18b. The weights leave an integer correction

The same construction with weights w_k gives

    sum_k w_k/mu_k = Z_w + sum_p c_p(w)/q_p,  Z_w in Z.

For word-size weights Z_w can still be recovered cheaply. For exponential
dyadic weights it has Theta(N) bits in ordinary examples, and computing it
to the precision needed for V is unresolved. No cost-preserving reduction
from weighted V to one unweighted H was established by stripping weights.

`validation/weighted_crt_gap.py` verifies 100 random weighted identities.
For w_k=16^(K-k), K=10,20,40,80,160, the integer correction has
33,73,153,313,633 bits. These observations are not lower bounds.

### 18c. A working acceleration for the actual weighted core

For odd A and a=2,4,8 the exact 2-adic identity is

    sum_{j>=0}16^j/(A-a*j)
      = -1/(15A) sum_{r>=0} r!(-16a)^r /
                       (15^r product_{h=1}^r(A-a*h)).

It follows from Newton finite differences and the binomial generating
function. Term r has valuation (5+v2(a))*r-popcount(r): 6–8 bits per term.
Across the four BBP families, the term-count ratio tends to 47/84.
Thus the old statement that Euler/binomial transforms cannot increase the
convergence rate was too broad: it counted the transformed argument but
missed the factorial valuation of the finite differences.

`decimal_extract_v5.py` and `decimal_extract_v5.cpp` implement this identity,
combine the four rational families before one Newton inversion, and generate
the band from endpoints. They preserve the original word phase. The extra
logarithm has **not** been removed from their complexity.

Validation: 200 random series identities, 64 exact weighted-V comparisons,
64 Python and 64 C++ digit comparisons with integer Chudnovsky, and 96
additional C++ core cases. Large sequential C++ comparisons and independent
MPFR digit checks are recorded in
`validation/cpp_euler_benchmark_2026-09-06.json`.

The final C++ build also truncates the denominator at each Newton step.
Its source hashes and final paired measurements are in
`validation/cpp_euler_final_benchmark_2026-09-06.json`. At N=10^7:
87.55 s -> 70.04 s total wall time, 56.41 s -> 39.83 s in the core,
814 MB -> 79 MB peak process working set. Both engines return
`725915133612`, confirmed independently with MPFR at two guard precisions.
These final-build figures are one paired run; the earlier file has two
alternating-order pairs before the last Newton optimization.

### 18d. The unrestricted interface ledger needs repair

Theorem 3 in `notes/linear_closure.tex` grants arbitrary operations and a free
oracle containing every individual mu_k^(-1). Reading and summing those values
already contradicts its claimed operation-independent precision ceiling.
That observation alone is not a fast algorithm.

Independently, the exact complete higher tower cannot be fixed while H varies:
for m>=5, Q_(1+2^(m-3))=Q_1 mod 2^m. The formal second world in that proof
need not correspond to any admissible reciprocal configuration. The numerical
series check does not test this input-consistency requirement.

The earlier valuation inequality for independently supplied truncated log
atoms can be retained with its narrower scope. It does not justify closing
every algorithmic route listed in §§10–17.

**Next mathematical target:** compute the exponential-weight integer correction
Z_w mod 2^m in O(M(N) log N), or compress the weighted sum by a different
proved construction. The unweighted atom is no longer the missing step.

---

## 19. A half-logarithm saving for the actual weighted core — 2026-09-06

**Result:** deterministic exact V in O(N log^(5/2) N), using the same
M(b)=O(b log b) convention as the manuscript. A pass of the decimal
pipeline inherits this bound; decimal carry certification retains its
separate qualification. No external priority or runtime record is claimed.

**Full proof:** `notes/weighted_block_crt_2026-09-06.md`.
**Constructive implementation:** `validation/weighted_block_crt.py`.

Set L=log₂N, u=ceil(sqrt(L)), ell=ceil(N/2^u). Within each BBP family,
normalize a consecutive block of t<=ell terms as

    S = sum_{j=0}^{t-1} 16^(t-1-j)/(A+a*j),   A=Theta(N).

For its denominator LCM D, the prime-valuation argument proves

    D divides [product_j(A+a*j)/oddpart(t!)] * lcm_odd(1,...,t),
    log₂D = O(t(1+log(N/t))).

This is a uniform bound, not observed cancellation. Reconstruct the
denominator-prime fractions T as in §18, now with local maximum prime
powers and the exponential weights. The exact identity S-T in Z holds.
The total local denominator bit mass for full blocks is O(ell*u).

Recover the integer by a real, alternating Euler series:

    S_infinity = 16^t/(15A) * sum_{r>=0}
                 (-a)^r*r! / [15^r product_{h=1}^r(A+a*h)].

The finite/infinite difference is below 1/15. Since A/ell grows like
2^u, R=O(ell/u+1) terms approximate S_infinity to absolute error below
1/64. The first omitted term supplies an exact integer certificate.
Approximate T with O(log N) fractional bits; the combined error in S-T
is strictly below 1/2. Because S-T is already known to be an integer,
nearest-integer rounding recovers it exactly, without a carry heuristic.

The two collections of binary-splitting trees cost, respectively,

    (N/ell) * O(M(ell*u)*L)        = O(N*u*L²),
    (N/ell) * O(M(ell*L/u)*L)      = O(N*L³/u).

Every tree level is charged. Choosing u=sqrt(L) makes both
O(N L^(5/2)). The small-prime sieve, small-modulus powers/inverses,
fixed-point divisions, global assembly, and final dyadic inverse fit
below that bound. Consecutive forests keep the stored data O(N) bits,
plus the workspace of the arithmetic primitives.

This changes the algorithm used for the expensive weighted aggregation.
It is not an inference from fast subroutines or a fitted timing exponent.
The estimate NuL²+NL³/u exposes a remaining half-log gap; optimizing this
estimate is not a lower-bound theorem about other algorithms.

Validation is recorded separately from the proof: 300 exact Fraction
blocks, 1,680 LCM-divisibility cases, 44 full-core and independent digit
positions, and full-residue agreement at N=100,000 and N=1,000,000.
The million-position full extractor also matches MPFR's independent pi
window, saved in `validation/weighted_block_crt_million_2026-09-06.json`.

---

## 20. Analytic audit and a certified decimal-output construction

**Derivation:** `notes/weighted_block_crt_audit.md`. No new numerical
experiments are used as premises in this section.

The core proof was rederived with an explicit uniform truncation count
R=ceil((8*ell+16)/u), complete per-level numerator/denominator size bounds,
and a sequential record-sorting implementation of factorization and data
routing. Sorting costs O(N L² log L), with L=log₂N, and small-modulus
arithmetic costs O(N L²(log L)²). Both are below N L^(5/2). No failed
identity or hidden log³ aggregation cost was found in this audit.

The previous decimal-certification qualification can be removed by a
further construction. Aggregate all word fractions sigma*r/mu as ONE
exact rational R_0, using global denominator-prime partial fractions.
Because |r|<mu=O(N), the integer correction for this sum is recoverable
with O(log N)-bit arithmetic. Its full rational denominator has O(N) bits,
and the entire aggregation costs O(N L²(log L)²), independently of the
later guard precision G. This avoids N separate G-bit divisions.

For fixed output width w, put n=N-1 and initially choose

    G = 32(N+w+1),
    K = n + ceil(G/4) + 2.

The absolute tail of the decimal-shifted BBP sum after K is below 2^-G.
Compute its fractional part f=U/Q EXACTLY from R_0 and the weighted V.
Writing 10^w U=dQ+r with 0<=r<Q, return d only if

    2^G r > 10^w Q,    and    2^G(Q-r) > 10^w Q.

If either test fails, double G and repeat. These inequalities are an
integer certificate that the full pi value lies in the same decimal cell.

The published irrationality bound mu(pi)<8 implies, for all sufficiently
large N, that 10^n*pi is more than 10^(-7n-8w) from every w-place decimal
cell boundary. The initial G above makes twice the tail bound smaller
than that distance. Thus the first pass succeeds for all sufficiently
large N. The finitely many remaining positions terminate by irrationality
and geometric tail decay; they do not change the asymptotic worst-case
bound. The stopping test does not require knowing the irrationality
theorem's threshold.

All of G,K,m remain Theta(N) on the asymptotic first pass, so the complete
certified extraction algorithm has the proposed deterministic
O(N log^(5/2) N) time bound. It never constructs earlier decimal digits.
The current fixed-guard Python/C++ engines do not implement this extension.
External novelty and independent proof review remain unestablished.

---

## 21. Shared word residues and a compressed Euler denominator

**Derivation:** `notes/shared_residues_and_euler_denominators_2026-09-07.md`.

For any odd BBP denominator mu, set q=(mu+1)/2, T=10^(N-1) mod mu,
and Y=2^(-q) mod mu. The first/third, fourth, and second families respectively
use the residues 8TY, 8TY^2, and 8TY^4 modulo mu. Sharing these quantities
and using a simultaneous bit scan for the powers of 5 and 2 preserves every
individual floor in the word accumulator. This is an exact algebraic
reorganization, valid for composite moduli including multiples of 5.

The native implementation is `decimal_extract_shared_residues.cpp`.
Three-run medians at the millionth position are 3.677 seconds for the
previous Euler implementation and 2.739 seconds for the combined method,
a 25.5% reduction in wall time. Both were rebuilt with the same compiler
flags and run in rotating order. Separate component measurements and exact
accumulator checks are recorded. An independent directed-rounding MPFR
reference verifies the decimal windows outside extraction timing. Larger
trials encountered concurrent CPU work and are marked exploratory.

For the real Euler series, let P_r=product_(1<=h<=r)(A+ah) and
f_r=oddpart(r!), with A odd and a a power of two. A prime-valuation argument
proves that P_r/f_r is an integer dividing lcm(A+a,...,A+a(R-1)) for r<R.
The truncated sum consequently has a common denominator of
O(R(1+log(N/R))) bits, which is O(ell) at the manuscript's parameters.
The former unreduced product size is O(R log N).

This size improvement does not yet supply a faster numerator-construction
algorithm. Straight cancellation in the same balanced tree still leaves
an O(R log^3 N) upper-bound calculation. The complete proved bound remains
O(N log^(5/2) N); the new timings concern a constant-factor improvement to
the separate fixed-guard Euler implementation. No world-record claim is
established by this experiment.

---

## 22. Analytic recursion for the Euler correction — 2026-09-07

**Proof note:** `notes/recursive_euler_reconstruction_2026-09-07.md`.
This investigation uses algebraic identities, error inequalities, and
bit-cost estimates; it adds no benchmark or extractor optimization.

Define H_r as the first r terms of the real Euler series and let
B_r(x)=(15^r-(x-1)^r)/(16-x)=sum_k B_(r,k)x^k. The quotient is an integer
polynomial. Partial fractions at A=-ak give the exact identity

    15^(r-1) H_r(A) = 15^(r-1) - a sum_(1<=k<r) k B_(r,k)/(A+ak).

The required coefficients can therefore be supplied as residues modulo the
local denominator prime powers. Prime reconstruction leaves an integer;
that integer is recovered exactly from a smaller Euler sum H_s, where

    d=floor(log_2(15A/(ar))),  s=ceil((4r+16)/d).

For d>=16 and r>=8, s<=r/2 and the scaled approximation error is below
2^-16. The integer-rounding certificate is explicit in the note. The
recursion's numerator, denominator, storage, and tail costs are accounted for.

The remaining batch problem is to obtain B_(r,k) modulo each required
prime power without constructing the full coefficient list. That list has
Theta(r^2) bits, as proved in the note. A standard integer matrix-product
analysis does not establish a faster batch, and factorial divisions can
be noninvertible modulo the requested prime powers.

**Conditional transfer only:** a batch bound O(r L^2 log L), L=ceil(log_2 N),
with the stated workspace condition would give full certified extraction
in O(N L^2 (log L)^2). More generally, a batch bound O(r L^(2+delta)) for
fixed 0<delta<1 would give O(N L^(2+delta/2)). The batch algorithms are
not proved. The unconditional O(N log^(5/2) N) theorem is unchanged.

The coefficient polynomial is independent of A. Common truncation lengths
allow its value-tree precomputation to be shared across all blocks at each
recursive depth. This weakens what an eventual amortized batch algorithm
must accomplish, but leaves the prime-power query costs unresolved.

---

## 23. Lucas reduction and sparse prime-power lifting — 2026-09-07

**New proof draft:**
`notes/lucas_sparse_prime_power_batch_2026-09-07.md`.
This extends Section 22 analytically. It introduces no benchmark premise.

At an occurrence of valuation v, the coefficient is needed only modulo
p^v, because its prime-fraction contribution is multiplied by q_p/p^v.
Therefore precision above p is needed only where p^2 divides a denominator.

For a reduced recurrence transfer P_(l,t)=U/D, the proof gives

    D divides 16^(t-l) binomial(r,l),
    v_p(D) <= floor(log_p r) for every odd prime p.

A remainder tree with explicitly inflated precision can consequently
perform every required nonunit division by exact gcd removal followed
by inversion of a unit. Its query count is small after excluding primes
at most (L+2)^8, whose complete modular rows are shared across blocks.

The remaining prime queries use a weighted Lucas identity. All occurrences
at the same prime in one block share a single low-index query. Across
blocks and recursive lengths the large-prime request mass is O(N) bits.
Their common integer value trees are constructed once; their full cost
and repeated forest traversals are charged in the proof. Smaller-prime
rows and block states are processed by sequential scans, avoiding an
unaccounted random-access assumption.

Taking u=ceil(32 log_2(L+2)), L=ceil(log_2 N), the construction first
proves O(N L^2 (log L)^2). A further exact bookkeeping refinement computes
each word residue once per denominator and merges prime contributions
through an O(N)-bit sequential registry in O(L^2) chunks. It also reuses
each original weighted-block residue across its prime divisors. This gives
the final proof-draft bound

    O(N log^2 N log log N),
    with O(N) stored bits plus arithmetic workspace.

The new note provides the algebra, the division invariant, the prime-power
count, the recursion transfer, and the time and storage accounting. It
does not establish the stronger per-block hypothesis from Section 22,
does not establish O(N log^2 N), and makes no external priority or measured
speed claim. Independent review is outstanding. The main manuscript and
arXiv source archive have not been changed to assert this new theorem.

---

## 24. Native implementation and actual measurements — 2026-09-07

`decimal_extract_lucas.cpp` implements the shared coefficient value tree,
Lucas remainder forests, sparse prime-power lifting with checked exact
division, recursive Euler reconstruction, and sequential exact word
aggregation. It uses the integer BBP cutoff and an adaptive integer
certificate for its 12-digit no-prefix decimal output.

The algebra self-test passes 3,560 coefficient checks, 80 reconstructed
Euler sums, and five exact modular-core and rational-word-sum comparisons.
It executes 876 nonunit divisions. Independent MPFR intervals verify 38
smoke positions, four forced low-guard certificate retries, and all 21
timed outputs. Every timed new recursive run at 100,000 or greater
executes recursive, ordinary-prime, higher-power, and nonunit-division paths.

Three-run wall medians at one million are 16.387 seconds for the new
recursive construction, 12.310 seconds for the matched certified direct
control, and 2.840 seconds for the existing fixed-guard Euler executable.
At ten million, single runs take 281.553, 189.082, and 38.195 seconds,
respectively. The new run peaks at 1.740 GiB process working set and
returns the certified and independently verified window 725915133612.

These measurements do not show a practical speed advantage for the new
method. The certified control is matched in all work other than recursive
Euler reconstruction; the existing Euler executable performs less output
certification work. Both comparisons are reported explicitly. No billion,
trillion, or quadrillion-position run has been performed.

Details: `notes/lucas_implementation_and_measurements_2026-09-07.md`.
Raw data and source hashes: `validation/lucas_benchmark_2026-09-07.json`.

## 25. Proof-to-implementation audit — 2026-09-07

The slowdown prompted an audit. The original measured implementation used
serial products and independent word reductions within shared query leaves,
where the proof requires balanced product and remainder trees. The current
source corrects these operations, with 4,100 additional repeated-leaf
checks. The old source is preserved at
`validation/decimal_extract_lucas_measured_2026-09-07.cpp`; original timings
and their hashes remain intact.

The abstract proof's literal parameters give ell=R=1 at ten million,
while the experimental build used ell=78,125 and R=39,065. Fixed
rescalings and correct outputs do not establish a complete correspondence
between the executable and the bit-complexity proof. Neither a practical
speed record nor an eventual runtime crossover has been established.

Full findings: `notes/lucas_proof_implementation_audit_2026-09-07.md`.
Separate corrected-run data: `validation/lucas_leaf_audit_2026-09-07.json`.

## 26. Two prime-power cores and smaller Lucas parameters — 2026-09-07

The native Machin extractor now provides the fastest local measured
implementation. It strips the 5- and 239-adic valuations of each odd
denominator, aggregates each prime-power contribution in one reversed
geometric sum, and uses a rigorous interval certificate for every
returned decimal window. Its analytic core length is approximately
0.425564 N, versus 1.299174 N in the BBP Euler implementation.

Three-run median wall times are 1.256 versus 2.756 seconds at one million,
and 19.086 versus 38.195 seconds at ten million, against the matched
certified Euler control. At ten million the new code uses 138.1 MiB
peak working set versus 75.1 MiB. All 27 timed outputs, 224 smoke windows
and 10 forced guard-retry windows passed independent checks.

These term-count and performance gains do not lower the binary-splitting
exponent. Its low-guard bound remains O(N log^3 N); unrestricted retries
also expose the known wide-accumulator cost unless an exact-word fallback
is supplied. This executable does not inherit the Lucas theorem.

Separately, the Lucas draft admits B=Theta(L^2/log L) and
Y=Theta(L^3/(log L)^2) with the same candidate O(N L^2 log L) bound.
The derivation pays for amortized row generation, large-operand reductions,
bounded sparse forests, shared value-tree storage and the sharper AP
incidence bound needed for local factor records. It still depends on the
earlier reconstruction lemmas and awaits independent review.

Trying to remove the last log log factor by globally merging shifted
prime fractions loses an integer carry needed by the 2-adic result.
That shortcut does not establish O(N log^2 N). The practical gain and
the smaller asymptotic bound have not yet been realized together in one
fully analyzed native implementation.

Derivations: `notes/machin_prime_core_extraction_2026-09-07.md` and
`notes/lucas_smaller_parameters_2026-09-07.md`.
Code: `decimal_extract_machin.cpp`.
Measurements: `validation/machin_benchmark_2026-09-07.json`.
The original measurement harness and unmodified raw capture are preserved.
The main manuscript and submission archive are unchanged.

## 27. Current-status consolidation — 2026-09-07

The practical success criterion is now recorded in RESEARCH_STATUS.md:
a smaller analytic bound must lead to a verified advantage in a corresponding
implementation at relevant large positions to meet this project's objective.
The measured weighted and Lucas prototypes have not demonstrated that win.
Earlier implications that their bounds established the fastest program,
a measured crossover, or a world record are withdrawn. A failed timing
comparison alone is not a proof of numerical error or a false big-O upper
bound; program-to-proof correspondence and tightness must be assessed
separately. An unobserved crossover is not treated as a demonstrated benefit.

The manuscript retains the weighted theorem and proof unchanged. It now
includes the Machin construction, its integer certificate, and the latest
matched comparison as a separate third performance table. Its current
20-page PDF, source archive, metadata, and verification report are synchronized.
The main README is concise; source mappings and commands are in the
reproduction guide. This log and the old paper audit are explicitly historical.
The earlier manuscript, README, audit, and package were backed up before editing.
No new asymptotic result, external speed record, or arXiv submission is claimed
by this editorial consolidation.
