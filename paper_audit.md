# Audit — `decimal_extraction_paper.tex` → arxiv-ready

## EXECUTION STATUS — ALL RESOLVED (paper rebuilt, compiles clean, 17 pp)
- **BLOCKER-1** ✓ forest Proposition added (`prop:forest`); Theorem 1 + abstract reconciled to Θ(N) space at quasi-linear time, with the O(log²N)-writable/quadratic-time variant stated separately.
- **BLOCKER-2** ✓ false p-adic-Gamma-library claim removed; replaced with the accurate PARI/GP `lngamma`/Morita statement.
- **MAJOR-1** ✓ CKR cited, novelty narrowed to single-prime sub-product-tree. **MAJOR-2** ✓ Gourdon build labelled non-micro-optimized, claim rests on asymptotics. **MAJOR-3** ✓ "subquadratic" → N^{1+o(1)}. **MAJOR-4** ✓ ZZ15 venue fixed. **MAJOR-5** ✓ BBG04 author order fixed.
- **P2** ✓ Brent/Salamin split. **P3** ✓ Bernstein dated. **P6** ✓ HAB02 wired (TC⁰), KL24 removed. **P7** ✓ "n! faster" reworded. (P1 abstract-paragraphing and the optional Shamir/unit-cost parenthetical deliberately left — non-essential.)
- **ENHANCE** ✓ §6 "The additive core, undressed" added: `prop:atom` (harmonic/digamma/Stirling equivalence, validated `atom_check.py`), external witnesses (CGH14/Harvey14/BZ10/Johansson09), hardness neighborhood (CCLL10/Bruin18/KS16/MZ23), Dwork closure in exclusion (a) (BV23 + our `dwork_amplify.py` result). All 3 over-reaches avoided (CKR digamma = "as we observe"; completeness = "suggests"; Johansson = "posed as open").
- **Citations** ✓ 10 new bibitems added, all cited; 37 refs, 0 orphans, 0 missing, 0 dangling. **AI disclosure** ✓ added.
- **Validators**: `atom_check.py` (new, backs prop:atom + the 961-bit figure), `dwork_amplify.py`, `dwork_crossformula.py` — all pass.

---

**Purpose.** Every claim, proof, bound, citation, and framing choice stress-tested before we touch the `.tex`. Line numbers refer to the *pre-edit* file. Fix each item, then edit.

**Original verdict (pre-edit — now executed; see EXECUTION STATUS above).** The core was sound and the results real (re-verified live: correctness 60/60 + spots, v3 = reference, N^1.065 scaling, N=10⁶ digit `130927562832` matches). The pre-edit draft carried **one main-theorem overclaim (space/time)**, **one factual error (p-adic Gamma libraries)**, a **novelty claim missing its nearest prior art (CKR)**, and a **comparison-fairness exposure (Gourdon)**; the barrier section was correct but under-anchored. **All of these are now fixed** — the four defects resolved and the barrier fortified with the external literature — the paper compiles clean (17 pp) and is arxiv-ready. The itemized entries below are retained as the record of what was changed and why.

**Severity legend.** `BLOCKER` = a referee rejects or the author is embarrassed · `MAJOR` = defensibility/correctness gap · `MINOR` = polish · `ENHANCE` = strengthens, not fixing an error.

---

## BLOCKER-1 — Theorem 1 claims Θ(N) space AND quasi-linear time together; no variant delivers both
**Where:** Abstract (l.42–50), Theorem 1 (l.127–139: "using Θ(N) bits of auxiliary space of which all but O(log²N) bits are a read-only string"), proof accounting §3.2 (l.356–393, which itself states workspace "O(L+Tw)" = Θ(N log N)).

**The defect.** Two incompatible regimes are being claimed as one:
- **v1** (per-term band reduction): Θ(N) total space, O(log²N) writable — but **O(N²) time** (per-term prefix read is O(J)=O(N)).
- **v3 / batched tree** (Thm 4): **O(N log³N) time** — but the moduli product tree's own root is Θ(N log N) bits, so workspace is **Θ(N log N)**, not Θ(N). (Checked numerically: at N=10⁷ the tree workspace ≈ 54 MB vs the 2.9 MB string — an honest log factor.)

The proof of Thm 4 (l.382–384) *itself* says workspace is Θ(L+Tw)=Θ(N log N). So the theorem's own proof contradicts its space claim. A referee catches this on first read (time × space product).

**Fix (choose one, (a) preferred).**
- **(a) Prove the simultaneous claim via a remainder *forest*.** Split the T moduli into Θ(log N) blocks of Θ(N/log N); process blocks sequentially, never forming the full product. Each block product is Θ(N) bits ⇒ **Θ(N) space**, at O(N log³N)–O(N log⁴N) time. This is exactly Harvey–Sutherland's remainder-forest device (already cited as [HS14]) and Bernstein scaled trees [Bernstein04]. Currently the code comment calls this "block/forest variant … future work" — so either implement/spell it out, or:
- **(b) Restate Theorem 1 as an explicit time–space tradeoff.** Quasi-linear time O(N log³N) at Θ(N log N) space; or Θ(N) space (O(log²N) writable) at O(N²) time; with the Θ(N)-space-quasi-linear point stated as achievable-via-forests (Prop., with (a)'s argument) rather than asserted in the headline.

**Recommendation:** do (a) as a short Proposition (it's true and citable), so the clean "Θ(N) space, quasi-linear time" headline survives — but *earned*, not asserted. Also reconcile the abstract wording.

---

## BLOCKER-2 — False statement: "no p-adic Gamma routine exists in the standard libraries at all"  [CONFIRMED, and worse than suspected]
**Where:** §6.3 crack remark, l.782–786.

**The defect (confirmed against PARI/GP primary docs).** Morita's Γ_p is implemented in **PARI/GP**, **SageMath**, and **Magma**. Sharper and more damaging: **PARI/GP's `lngamma` computes *Diamond's* log-gamma specifically for |x|>1** — and that is *exactly* the function and domain the paper invokes. The barrier argument is Diamond's log-Γ at argument ≈(1+8K)/8, which has v₂=−3, i.e. |x|₂=8>1 — precisely the |x|>1 branch PARI evaluates. So the parenthetical is false for the very function, at the very argument type, the paper is talking about. A p-adic referee knows this instantly.

**Fix.** Delete the parenthetical. Re-anchor the novelty on *complexity at a specific argument family*, not existence — see MAJOR-1. Suggested replacement clause: "(PARI/GP, Magma and SageMath compute Morita's Γ_p, and PARI/GP's `lngamma` computes Diamond's log-Γ for |x|>1, via Mahler/Dwork expansions whose cost grows with the target precision; none provides a *sub-product-tree* evaluation specialized to the arguments 8/(1+8K), K=Θ(m))."

---

## MAJOR-1 — "First sub-product-tree evaluation of Diamond's log-gamma" needs the CKR caveat
**Where:** §6.3, l.780–786 ("the first recorded sub-product-tree evaluation of Diamond's log-gamma at such arguments").

**The issue.** Costa–Kedlaya–Roe (*Hypergeometric L-functions in average polynomial time* I/II, arXiv:2005.13640, 2310.06971) compute the p-adic (Morita) Gamma function in **average polynomial time**, and their block-lower-triangular matrix already carries *product + sum-of-partial-products* (= a Γ_p value and its derivative). That is the nearest prior art and it is currently uncited.

**Fix.** Reword to the precise, unimpeachable claim: *a single-argument, single-modulus sub-product-tree evaluation of Diamond's 2-adic log-Γ at height-Θ(N) arguments*, and explicitly contrast CKR's **amortized-over-primes** regime (quasi-linear in the prime bound X, fixed small precision) with ours (one fixed prime p=2, precision Θ(N)). Add CKR to the bibliography. This turns a vulnerable claim into a sharp one.

---

## MAJOR-2 — The Gourdon 232× comparison is a fairness exposure
**Where:** Abstract (l.65–69), §7 measurements (l.925–957, table l.931–942).

**The issue.** "232× at N=10⁶" and Gourdon "≈N^1.8" rest on *our* implementation of Gourdon's method ("verification-grade"). A referee will ask whether it was optimized; an unfair constant invites rejection of the whole experimental section.

**Fix.** (i) Lead with the **scaling slopes** (N^1.077 vs ≈N^1.8), which are implementation-robust and follow from the *analyses*, not the constants. (ii) Explicitly label the Gourdon build "a straightforward, non-micro-optimized implementation faithful to the published algorithm; the constant-factor gap is indicative, the asymptotic gap is the claim." (iii) State Gourdon's own complexity O(N² loglog N/log²N) so that *even a perfect* implementation is quadratic-class — making the point independent of our build quality.

---

## MAJOR-3 — Sharpen "subquadratic" so Gourdon can't be read as a counterexample
**Where:** contribution (ii), l.170–171; abstract l.56.

**The issue.** Gourdon's O(N² loglog N/log²N) is *technically* o(N²). Calling ours "the first subquadratic decimal extractor" invites the quibble "so is Gourdon."

**Fix.** Say **"first with complexity N^{1+o(1)} (softly linear)"** and classify Gourdon as **"quadratic-class (N²/polylog)."** Keep the "to our knowledge / implemented and verified" hedge (Gourdon Thm 2 is asserted-not-implemented; Plouffe22 is quadratic-class and table-capped — both confirmed by the literature sweep).

---

## MAJOR-4 — ZZ15 cited to the wrong venue (citation error a moderator flags)
**Where:** bibliography, l.1067–1070 (`ZZ15`), used in the certification remark (l.329–332).

**The defect.** The paper cites "Exp. Math. 24 (2015), 419–423." Wrong: the Zeilberger–Zudilin irrationality-measure paper is **Mosc. J. Comb. Number Theory 9 (2020), no. 4, 407–419** (arXiv:1912.06345). An arXiv preprint dated Dec 2019 cannot appear in a 2015 issue; there is no Experimental Mathematics version. The math content (μ(π) ≤ 7.1032… ⇒ the paper's "≤ 7.11") is **correct** — only the venue is wrong.

**Fix.** Replace the bibitem: `D. Zeilberger, W. Zudilin, "The irrationality measure of π is at most 7.103205334137…," Mosc. J. Comb. Number Theory 9 (2020), no. 4, 407–419; arXiv:1912.06345.`

---

## MAJOR-5 — BBG04 author order is wrong (on the keystone impossibility citation)
**Where:** bibliography `BBG04`, l.~986–989.

**The defect.** The published byline is **Jonathan M. Borwein, David Borwein, William F. Galway** (J. M. Borwein first — confirmed on title page, running headers, and the CUP record). The bibitem lists "D. Borwein, J. M. Borwein, W. F. Galway," swapping the first two authors. This is the citation that motivates the entire paper (no non-binary Machin-BBP), so getting the byline wrong is a bad look; a citation-checking referee flags it. (The *usage* is correctly scoped — "no Machin-type BBP arctangent formula for b≠2^m," not "no BBP formula of any kind" — so only the author order is wrong.)

**Fix.** `J.~M.~Borwein, D.~Borwein, W.~F.~Galway`.

---

## MINOR / POLISH
- **P1 — Abstract is a single 50-line paragraph** (l.34–83). Split into 2–3 for readability; lead with the result, then the method, then the barrier.
- **P2 — Split the bundled Brent/Salamin bibitem** (l.1028–1035): prose says "Brent and Salamin," key is `Brent76` only. Make `\cite{Brent76,Salamin76}`. (From citation agent B.)
- **P3 — Optionally date the Bernstein draft** (2004-08-20). (Agent B.)
- **P4 — Terminology register.** "faces," "the war," "atom" read as campaign-slang. "Faces"/"equivalent formulations" is fine if defined once; scrub "war"-type phrasing for arxiv tone.
- **P5 — Las Vegas certification remark** (l.321–333): correct, but state plainly that binary BBP makes the same heuristic tradeoff, so we're not weaker than the accepted standard.
- **P6 — Orphan bibliography entries.** `HAB02` (Hesse–Allender–Mix Barrington) and `KL24` (Kinoshita–Li) appear in the bibliography but are **never `\cite`d** in the body — arXiv/referees flag dangling references. Two options each: (i) wire HAB02 into the l.745 claim "the mass is depth-polylogarithmic parallelizable" (it actually supports *more* — iterated multiplication/division is in DLOGTIME-uniform TC⁰, i.e. *constant* depth — so tighten that sentence and cite it); resolve KL24 similarly or (ii) delete both entries. The l.745 sentence is currently an *unsupported* claim regardless.
- **P7 — L649 "no algorithm is known to compute n! faster" is loose.** Borwein85 *itself* beats the product tree (by ~log n/log log n). Reword to "no algorithm beats the M(n log n) scale," so it's not self-contradicting the sentence above it.

### Citation-agent refinements to bank (barrier-core batch)
- **Cheng, not Borwein, for the factorial lower bound.** If the barrier section ever cites Ω(M(n log^{4/7−ε} n)), it is **Qi Cheng, "On the Ultimate Complexity of Factorials," STACS 2003 / TCS 326 (2004) 269–283**, in a *restricted* (straight-line) model — never Borwein85 (an upper-bound-only paper). Currently not in the manuscript, so no error today; guard against it if adding. Add Cheng to the bib if used.
- **Soften "non-holonomic ⇒ hard" around Γ_p.** A special-purpose fast Morita Γ_p *does* exist (Rodriguez-Villegas 2007, cited within CMTV21). Where §6 leans on "outside the holonomic class ⇒ no fast evaluation," add a clause: Villegas gives a fast Γ_p, but *not* Diamond's log-Γ at the arguments 8/(1+8K) in quasi-linear time. Keeps the barrier honest.
- **Shamir precision (optional):** add "(in the unit-cost arithmetic model)" — Shamir's theorem is an arithmetic-step statement, not a Turing-machine factoring algorithm. Doesn't affect the vacuous-mod-2^N conclusion.

---

## ENHANCE — Fortify §6 (the barrier) into a fortress. This is the biggest upside.
The barrier is currently self-contained ("we tried fifteen routes"). We now hold the external literature that turns it into "this is the wall the whole field is stuck at, mapped." Add:

1. **The equivalence, stated as a Proposition (elementary, provable):** stripping the geometric weights, the additive core's numerator is Π′(0) for Π(t)=∏_{k≤K}(8k+1+t); hence the core ≡ **an isolated harmonic number / 2-adic digamma mod 2^N** ≡ the x¹-coefficient of a product of N linear forms ≡ (via n!H_n=|s(n+1,2)|) an **unsigned Stirling number of the first kind**. Validated in `validation/atom_check.py`.
2. **External witness that it's open** (this is what "bulletproof" needs) — stated at the right strength: the single-prime factorial barrier is **citable verbatim** to Costa–Gerbicz–Harvey (already in our bib): "(p−1)! mod p² in time p^{1/2+ε}," the "best known for a single w_p" — so beating √p is open. Binary-splitting cost of a rational-term sum → cite **Brent–Zimmermann, MCA (2010)**. The "differentiation destroys the factorization, so the multiplicative shave doesn't cross to the sum" point is **our own observation** (validated, item 3), corroborated by Johansson's blog posing exactly this as open — cite the blog *as a blog*, do not source a theorem to it. (Drop the Kurepa/2×2-accumulator flourish — only weakly citable, an MSc thesis.)
3. **Why the multiplicative crack cannot cross to the additive side, as a Remark:** Π(0) is smooth (Borwein's Legendre-exponent factorization applies); its t-derivative Π′(0) is a *generic* integer with a large prime factor (validated: 961-bit factor at K=200, `atom_check.py`) — differentiation destroys the unique factorization the crack runs on.
4. **The hardness neighborhood** (evidence, explicitly not proof), stated at defensible strength: Cai–Chen–Lipton–Lu dichotomy — quadratic exponential sums in P (via Gauss sums), and **specific degree-≥3 families** #P-hard even at fixed prime power (NOT "every degree-≥3 sum"); Bruin — **no *known* efficient classical** single-Kloosterman-sum algorithm (state-of-the-art, not a lower bound); Kowalski–Sawin / Milićević–Zhang — partial Kloosterman sums converge in law to random processes (pseudorandom) **even to high prime-power moduli**. The inference "so completeness, not field-vs-ring, is what keeps complete sums tame" is **our gloss** — present it as ours, attribute only the convergence-in-law results to them.
5. **Close the Dwork escape in print + our new validated result:** excellent Frobenius lifts reach mod p^{κs} only for **odd p** (1≤k<p) — **vacuous at p=2**; Beukers–Vlasenko could not push the truncated-sum congruence past p^s even conjecturally. Add our machine-validated finding (`validation/dwork_amplify.py`, `dwork_crossformula.py`): the Frobenius-descent precision is **exactly s+v₂(d)−1**, slope 1, reproducing the cross-formula tower law from the Dwork side — so basic Frobenius descent yields only Θ(log N) digits at depth log N. This upgrades exclusion (a)/(b) of the existing Remark.
6. **CKR as the named "closest machine."** Cite them (verbatim) for: a block lower-triangular matrix that records *a product AND the sum of the partial products* simultaneously, via accumulating remainder tree/forest, computing Morita Γ_p in average polynomial time — **amortized over all p ≤ X**. Then state as **our observation** that this product+partial-sum object is exactly a Pochhammer and its logarithmic derivative (the digamma), which CKR do not remark on. Forward pointer: can the amortization-over-primes axis be traded for single-argument precision at fixed p=2? — our stated open direction. (Do NOT attribute the digamma reading to CKR.)

**Net effect:** Open Problem V becomes "equivalent to a cluster of recognized open/hard problems (isolated harmonic numbers, single-prime factorials, incomplete Kloosterman), one log-factor above the product side, with the tame/wild dichotomy proven for the neighboring polynomial case." That is a genuine, citable contribution — far stronger than the current self-referential framing.

---

## Citation status (verification agents; ~33 refs)
- **Algorithm machinery (7): ALL CLEAR** — HvdH21, AFKL19 (genuinely conditional ✓), CGH14, HS14, BorodinMoenck74, Bernstein04, Brent76+Salamin (benchmark ✓). Only P2/P3 cosmetics.
- **Prior-art/impossibility (7): DONE.** BBP97 ✓ (SC open problem confirmed), BBG04 (**author order wrong** — MAJOR-5; usage correctly scoped), Plouffe96 ✓, Bellard97 ✓, **Gourdon03 verified verbatim** (Thm 1/2 bounds, memory, "details soon" promise, asserted-not-proved — the novelty claim's foundation is solid), Plouffe22 ✓ (10⁸ cap verbatim), Zudilin24 ✓ (base-10→5 title change + flaw + denominator all verbatim). One error (BBG04), else publication-ready. Minor: "in its final section" → "final theorem (§4)."
- **Barrier-core (6): DONE — all clean** (bib + usage correct): Borwein85, Strassen77, BGS07, Shamir79 (vacuous mod 2^N ✓), CMTV21 (Γ_p/digamma confirmed outside their scope ✓), BGP17. Three refinements banked below.
- **Special functions (5): DONE.** Diamond77 (bib ✓; BLOCKER-2 library claim **confirmed false**), Washington (✓), HAB02 (bib ✓ but **orphan** — P6), ZZ15 (**wrong venue** — MAJOR-4), FK09 (✓, optionally add AJM 2012 pub). Net: 2 clean, 1 confirmed-error, 1 venue-error, 1 orphan.
- **New to add: DONE — verified, with 3 over-reaches to avoid** (see below). Bibitem-ready:
  - Cai–Chen–Lipton–Lu, *On Tractable Exponential Sums*, FAW 2010, LNCS **6213**, 148–159; arXiv:1005.2632. ✓ (phrase hardness as "specific degree-≥3 families," not "all")
  - Bruin, *On quantum computation of Kloosterman sums*, arXiv:1807.03600 (2018) [preprint]. ✓ (frame as "no *known* efficient classical algorithm," not hardness)
  - Kowalski–Sawin, *Kloosterman paths and the shape of exponential sums*, Compos. Math. **152** (2016) 1489–1516; arXiv:1410.7892. ✓
  - Milićević–Zhang, *Distribution of Kloosterman paths to high prime power moduli*, Trans. AMS Ser. B **10** (2023) 636–669; arXiv:2005.08865. ✓ (the high-prime-power one — load-bearing)
  - Costa–Kedlaya–Roe I, ANTS-XIV, Open Book Ser. **4** (2020) 143–159; arXiv:2005.13640. Costa–Kedlaya–Roe **II**, Res. Number Theory **11** (2025) art. 32; arXiv:2310.06971. ✓ (block-triangular matrix = product + sum-of-partial-products **verbatim**; Morita Γ_p in APT; amortized over p≤X)
  - **Harvey factorial barrier → re-point to CGH14** (already in our bib!): *A search for Wilson primes* gives, verbatim, "(p−1)! (mod p²) in time p^{1/2+ε}" and "best known … for a single w_p" ⇒ beating √p is open. For the *point-counting* framing add Harvey, *Counting points on hyperelliptic curves in APT*, Ann. of Math. **179** (2014) 783–803; arXiv:1210.8239 (never mentions factorials — cite ONLY for the APT framing).
  - **Johansson blog = weakest link.** It does NOT assert binary-splitting-optimal-for-H_n (that bound is stated for n!), and the transfer question is an open *musing* ("not obvious to me how it would be done"), not a claim. Cite **Brent–Zimmermann, *Modern Computer Arithmetic*, CUP 2010** for the binary-splitting complexity; treat "the trick doesn't transfer" as our own observation / folklore (optionally cite the blog *as a blog*). Do NOT source a theorem to Johansson.

### THREE interpretive over-reaches to NOT put in the paper (adversarial-verify catches)
1. **CKR do NOT compute a derivative/digamma.** Their "sum of partial products" feeds the Beukers–Cohen–Mellit trace formula — no Γ′/digamma anywhere in either paper. The reading "this is the Pochhammer *and its log-derivative*" is **our observation**, and should be presented as our contribution, not attributed to CKR.
2. **Johansson's "factorial trick doesn't transfer" is an open musing, not an assertion** — see above.
3. **"Completeness (not field-vs-ring) is what makes complete sums tame" is our gloss**, not a statement in Kowalski–Sawin or Milićević–Zhang. Attribute only the convergence-in-law (pseudorandomness) results to them; keep the completeness contrast as our framing.

---

## AI-assistance disclosure (author chose: include). Draft candidates — pick one:
- **(concise)** "Computational experiments, the literature survey, and manuscript preparation were carried out with assistance from an AI system (Anthropic Claude). All mathematical identities were independently machine-validated by the exact-arithmetic scripts in the accompanying repository (`validation/`); the author is solely responsible for the results and their correctness."
- **(fuller)** adds: "In particular, the barrier characterization of §6 was developed against a systematic search of the p-adic, exponential-sum, and holonomic-summation literature; each cited external result was verified against its primary source."

Place in Acknowledgments (l.975–977), replacing the current placeholder comment.

---

## Prioritized edit checklist
1. **BLOCKER-1** — remainder-forest Proposition (or restate as tradeoff); reconcile abstract + Thm 1 + §3.2. *(math + writing)*
2. **BLOCKER-2** — delete the false library parenthetical. *(one line)*
3. **MAJOR-1** — reword the log-Γ novelty; cite CKR. **MAJOR-2/3** — Gourdon fairness + "N^{1+o(1)}". *(writing + 1 bibitem)*
4. **ENHANCE** — rebuild §6 with items 1–6; add `atom_check.py`, `dwork_amplify.py`, `dwork_crossformula.py` to the validation manifest. *(the big value-add)*
5. **Citations** — fix BBG04 author order (MAJOR-5), ZZ15 venue (MAJOR-4), resolve HAB02/KL24 orphans (P6); add the verified new bibitems (CCLL, Bruin, Kowalski–Sawin, Milićević–Zhang, CKR I/II, Brent–Zimmermann; reuse CGH14 for the factorial barrier); write the fortification with the 3 over-reaches avoided; P2/P3 cosmetics.
6. **AI disclosure** + P1/P4/P5 polish.
7. Recompile (pdflatex ×2), re-check all `\ref`/`\cite` resolve, no overfull boxes.
