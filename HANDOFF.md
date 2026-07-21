# HANDOFF — BBP decimal digit extraction campaign
**Written 2026-07-21 (end of the great closing day). Read this first, then `campaign_graph.md` (the war board), then skim `research_log.md` Part 16 for detail.**

## Protocols (standing, non-negotiable)
- **BINGO protocol**: say BINGO **only** when the main prize lands — V computed in O(M(N) log N), machine-validated, which gives digit extraction at O(N log²N) = parity with full computation. Never for lesser results. No celebratory language outrunning results (this failed once; it's why the protocol exists).
- **Machine-validate every identity before enshrining it** in log/board/paper. Validators live in `validation/` (18 scripts, all exact-arithmetic, all passing as of today).
- Update `campaign_graph.md` and `research_log.md` on every result, positive or negative, with kill-reasons for negatives.

## What is proven and shipped
- **Main theorem**: N-th decimal digit of π in O(N log³N) bit ops, Θ(N)-bit read-only space (the shared string X = 5^(N−1)). First implemented+verified subquadratic decimal extractor. Verified against public corpora to N = 10⁷.
- **Measured**: v3 implementation runs N^1.077 (10⁷ in 111.6 s single-thread Python/gmpy2); 232× faster than compiled Gourdon at 10⁶.
- **Sieve-batched exponentiation**: all modexp data for the family at amortized O(loglog N)/modulus — also speeds classical binary BBP.
- **The Borwein crack**: entire multiplicative cluster — Π(8k+1) mod 2^m, exponent sum E, power sums, Diamond log-gamma at small height — down to O(M(N) log N) = N log²N. Validated to K = 3·10⁶ (`validation/borwein_crack.py`).
- **Paper**: `decimal_extraction_paper.tex`, 15 pp, compiles clean (pdflatex, two passes). Contains main theorem, refined pipeline, generalizations, the barrier section with Open Problem (V), mass-conservation conjecture, and the new **Remark (structural exclusions)**. ~25 bibitems incl. Brent76 (AGM benchmark) and FK09 (Fresnel remainder).

## The benchmark (recalibrated today — critical framing)
Brent–Salamin AGM computes **all** N digits in O(N log²N) since 1976, same Θ(N) memory class. So the prize is **extraction at parity with computation**: match N log²N while keeping extraction semantics (no earlier digits, read-only shared string, word-parallelism, spot-verification). We sit exactly **one log factor** above parity. That log = the additive core V, the sole holdout (every other pipeline phase is at N log²N or below — see the ledger in the conversation log / paper §4).

## The atom and its faces
"Face" = an exactly-equivalent formulation of computing V (cheap conversions both ways, machine-validated). One object, many costumes; crack any face and all fall (this actually happened once: the multiplicative faces fell together under the Borwein crack). Current faces: V's native sum; the truncated 2-adic logarithm trunc_M[Log₂(1−z)], v(z)=½ over ℚ₂(ζ₈) (8-section identity validated 880 bits — complete logs are cheap, truncation is the war); digamma/H-tail; tapered Bernoulli/Iwasawa tables (Plouffe22's dependency); Dwork-tower amplification gap; 1/A-Eulerian (Apostol–Bernoulli) expansion; **the carry face** (newest, best for lower bounds): sliced into G-bit windows, every slice-sum is word-cheap modexp data; ALL hardness = the Θ(N/G) carry bits between slices.

## What is closed (with kill-reasons — do not re-attack without new structure)
- **Formula shopping**: cross-formula tower law min v₂(R_s) = s + v₂(d) − 1 — slope universal for every BBP family d·t+1 (`validation/bellard_tower.py`).
- **AGM/Landen depth-doubling**: Lerch duplication is exact but parameter-branching (1 param → N params in log N doublings = the product tree again). Logs have an AGM because they have ONE parameter.
- **Incomplete reciprocity (L2)**: real-analytic version requires a transcendental Fresnel remainder (Fedotov–Klopp); calibrated LLL over ζ₈-duals at two parameter points: empty (`validation/l2_swap_hunt.py`).
- **Algebraic cross-scale relations generally**: calibrated all-integer LLL (controls recovered at height O(1); noise floors 2^17–2^38): nothing in between — linear, quadratic, exp-side (`validation/agm_law_hunt.py`).
- **Salié twist**: quadratic-character twists do NOT collapse V (`validation/twist_hunt.py`) — V is Kloosterman-like, not Salié-like.
- **Dwork tower amplification**: the tower is real and proven-in-sketch (collapse lemma + complete-block inversion, `validation/tower_engine.py`) but its O(s) precision ceiling IS the atom — cannot amplify.
- Also dead: randomization (unit weights, no concentration), matrix-factorial BSGS (N^1.5 at our aspect ratio — contiguity recurrence V(A+8) = 16V(A) + cheap exists, see log), Gosper (non-summable), Mahler transform (mass), inverse-periodicity Mersenne folding (mass), holonomy (x/(e^x−1) non-holonomic), composition/chirp (no roots of unity in ℤ/2^m), tables (don't fit), information-theoretic bounds (impossible — output small).

## Live leads (the attack queue, in order)
1. **Dwork-tower standalone note** — tower + collapse lemma (1/((di+1)(dk+1)) = [2/(di+1)+1/(dk+1)]/(dj+3) for i+2k=j) + engine + cross-formula law = self-contained arXiv note on a Dwork-type congruence for a 2-ramified ₂F₁. Ready to write; all math validated.
2. **Prove mass conservation in a restricted model, via the carry face** — carry/communication complexity has real lower-bound techniques (unlike algebraic circuits). The structural-exclusions remark is the case-analysis scaffolding.
3. **Kloosterman ↔ V reduction hunt (tame/wild mirror)** — finite-field exponential sums mirror our dichotomy exactly: tame (Gauss/Jacobi = Γ_p products) = our cracked cluster; wild (Kloosterman, no polylog algorithm, decades open) = V. A reduction either way imports their hardness evidence wholesale.
4. **Las-Vegas carry refinement** — compute slices at low precision, refine only near-boundary carries. Naive version blocked (terms span many slices → N²/G); needs a new idea for per-term O(1) slice participation.
5. **Literature watch** — single-Bernoulli/Apostol–Bernoulli (Harvey's frontier), p-adic L-values (KL24 line), excellent lifts (BV-III), exponential-sum computation, multiplication itself.

## File map
- `decimal_extraction_paper.tex/.pdf` — THE paper (unified, results-only, 15 pp).
- `campaign_graph.md` — war board: status line, result stack, war-map, attack queue, source index. **The primary orientation document.**
- `research_log.md` — complete campaign record (Part 16 = this campaign, ~20 dated entries).
- `cracks_dossier.md` — ranked leads, ~25 closed doors with reasons, cost tables, sources.
- `decimal_bbp_extract.py` — reference implementation (v1/v2/v3; --v3-check, --bench etc.). `three_phase_extract.cpp` — independent C++ check. `gourdon_digit_extract.cpp` — the race opponent.
- `validation/` — all 18 validator scripts (copied from session scratchpad today; scratchpad dies with the session).
- Memory (auto-loads): user profile, project state (updated today through the closing-day results), BINGO protocol.

## Where the last session ended
All three greenlit attacks closed in one day (AGM: three kills; L2: empty at two parameter points; conjecture fortified into the paper as the exclusions remark). Then the wide branch sweep: five more kills, the tame/wild mirror discovered, the carry face found (from Joseph's "crack it into pieces" push). Strategic posture shifted from assault to siege: build permanent artifacts (queue #1), prosecute the provable front (queue #2), keep watch (#3, #5). The natural next action is **queue #1: write the Dwork-tower note** — say the word or just go.
