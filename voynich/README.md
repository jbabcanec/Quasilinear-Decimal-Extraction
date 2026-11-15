# Poured or fitted? Physical line-break forensics of the Voynich manuscript

**Voynich lines are fitted to the space, not wrapped at the first word that fails to fit. So are the lines of real 15th-century manuscripts.**

This study measures the manuscript's text in pixel space, using a bounding box for every word, and asks when a line ends: because the next word of an existing text did not fit, or because the scribe filled the space?

The first version of this report read the answer as evidence that the text was not copied from a finished text. A control on real manuscripts (section 4) shows that reading was wrong: scribes copying known French prose leave the same fingerprint. What is left is a description of how Voynich lines end, and one difference from the real manuscripts.

## Findings

1. **Voynich lines are fitted, not mechanically wrapped.** Across 3,183 line breaks, the blank at a line's end is unrelated to the width of the word that starts the next line (coupling 0.05, 95% interval −0.03 to 0.14). Mechanically wrapping the same words on the same pages gives 0.48.
2. **Real scribes do the same.** Six French prose manuscripts of the 15th century, measured the same way, give couplings from −0.44 to +0.29, five of them at or below +0.09. Mechanical wrapping of their own words gives +0.18 to +0.51. The fitted line end is ordinary scribal practice.
3. **One difference survives.** Voynich lines end on shorter-than-average words (−9%). All six real manuscripts end lines on longer-than-average words (+13% to +26%).
4. **Words are not copied from the word physically above them.** At most 0.7% of words can be such copies. The Timm–Schinner generator's own output scores about the same on this test, so this limits the naive picture of copying, not their algorithm.

## Data

Word boxes come from voynichese.com (one box per word of Takahashi's transliteration, on page images about 600 px wide), joined by character alignment to the Zandbergen–Landini transliteration (ZL3b) for line, paragraph, section and hand labels. Running paragraph text only: 33,528 words on 3,968 lines.

All lengths are in **word widths**: the median word-box width on that page (about 45 px).

## 1. Lines are fitted, not wrapped

A scribe writing out a text that already exists breaks the line when the next word will not fit. That leaves two fingerprints:

- the blank at the line end is anywhere from nothing to a whole word, and
- the blank is smaller than the next line's first word, because that word is the one that did not fit. Lines followed by long words therefore end earlier.

A scribe who chooses what to write leaves neither.

To know what each fingerprint looks like on these pages, three simulated scribes re-lay the manuscript's own word stream, with the real box widths, on the real page geometry:

- the **wrapper** puts the next word on the line if it fits, otherwise starts a new line;
- the **squeezer** does the same but may compress the last word by up to 30%;
- the **fitter**, near the line end, takes the first of the next twelve words that fits and leaves the rest in order.

![Blank left at the line end](figures/fig1_slack.png)

![Blank against the next line's first word](figures/fig2_coupling.png)

| | Median blank | Middle-half spread | Coupling to next word |
|---|---|---|---|
| Manuscript | 0.13 | 0.33 | 0.05 (−0.03 to 0.14) |
| Wrapper, 20 runs | 0.32 | 0.68 | 0.48 (0.42 to 0.51) |
| Squeezer, 20 runs | 0.17 | 0.48 | 0.33 (0.29 to 0.36) |
| Fitter, 20 runs | 0.15 | 0.36 | 0.01 (−0.03 to 0.07) |

*Coupling* is the regression slope of the blank on the next line's first-word width, within page. Manuscript interval: bootstrap over pages. Simulations: mean and range.

The result holds in every section and for every hand with enough lines to measure. The pharmaceutical section (166 breaks) is too small to say anything: 0.31, interval −0.15 to 1.16.

![Coupling by section and hand](figures/fig3_sections.png)

### Loopholes tested

| Alternative explanation | Test | Result |
|---|---|---|
| The scribe stretched or squeezed the writing to hit the margin | Repeat with widths predicted from each word's glyphs instead of measured widths | Coupling −0.07 (−0.18 to 0.04); wrapper 0.45 |
| | Is the last word narrower than its glyphs predict? | By 0.05 word widths, and the gap before it shrinks by 0.07. Closing the distance to a wrapped text needs about 0.35 |
| Words are split across the line break | Is "last word + next line's first word" a known word? | 1.0% of the time, against 4.3% for neighbours inside a line |
| The last word is padding added after the real next word failed to fit | Couple the room *before* the last word to the next line's first word | 0.09; wrapper 0.59, fitter 0.05 |
| The last word is just any word that happens to fit | Replace it with a random mid-line word that fits the same room | Spread of the blank rises from 0.33 to 0.57 |
| Short or broken lines distort it | Keep only lines of six or more words | Manuscript 0.01, wrapper 0.47 |
| The inferred right edge is wrong | Use no edge at all: within a paragraph, does a line stop further left when the next line starts with a wider word? | Manuscript 0.00; wrapper −0.45, fitter −0.02 (wrapping predicts about −0.5) |

## 2. The last word is cut to measure

![Last word against the room it had](figures/fig4_lastword.png)

The last word's length tracks the room left for it. This is not distinctive: the real manuscripts in section 4 show the same pattern (about 2.2 characters when room is tight, 6.5 when it is wide), and part of it is built in, because the room includes the word that was written into it.

The last word's length tracks the room left for it. Its ending does not: 18% of line-final words end in the glyph `m`, against 1% elsewhere, and that share is about the same whether the room was tight or generous. The line-end mark follows the count (this is the last word); the length follows the pixels.

Two long-known oddities come out of fitting without further assumptions. In the manuscript, line-first words are 9% wider than mid-line words and line-last words 10% narrower. The fitter gives +9% and −16%; the wrapper gives +10% and +1%.

A line is also closed off statistically. The last glyph of a word predicts the first glyph of the next word (0.17 bits between neighbours inside a line, 0.01 for words one apart). The link is intact between the second-last and last word of a line (0.17 bits) and gone across the break, whether measured from the last word (0.004 bits) or the second-last (−0.001). So the last word belongs to its line, and nothing carries over to the next one.

## 3. No copying from the word overhead

The self-citation hypothesis (Timm & Schinner) holds that the text was produced by copying and modifying words already on the page. A simple version says a word is often a variant of the word above it.

For every word, each earlier word on its page was scored as similar (within one edit) or not, and expressed relative to how often two random words of that page are similar. Words overhead do score higher, but they score equally higher on lines 6 to 20 above, where copying from overhead cannot be the cause. That part is a position-in-line effect.

![Similarity by horizontal offset](figures/fig5_kernel.png)

The *overhead contrast* divides the overhead advantage on the two lines above by the same advantage on lines 7 to 20 above. It is 1.0 when nothing is copied from overhead.

| | Overhead contrast |
|---|---|
| Manuscript | 1.04 (0.98 to 1.09) |
| Synthetic pages, no copying, 12 runs | 1.00 ± 0.04 |
| Synthetic pages, 1% of words copied from overhead | 1.09 |
| Synthetic pages, 5% copied from overhead | 1.70 |
| Timm–Schinner generator output, its own lines, 3 seeds | 1.06, 1.07, 1.11 |
| Naibbe cipher of Latin, wrapped | 1.04 (0.97 to 1.11) |
| Latin (Pliny), wrapped | 0.95 (0.83 to 1.08) |

![Calibration](figures/fig6_calibration.png)

Reading the manuscript's upper limit off the calibration gives at most 0.7% of words copied from overhead, with a best estimate of 0.4%. What the manuscript does show is a diffuse effect: a word is 23% more likely than chance to resemble a word on the line above, 15% for three lines up, 6% for ten, at any horizontal offset.

## 4. Control: real manuscripts

Six manuscripts from `Gallicorpora/HTR-MSS-15e-Siecle` (15th-century French, line-level ground truth, commit `707a106`), main text zones only. Each word gets a synthetic box whose width is proportional to its character count, and the result goes through the same code as the Voynich pages, including the simulated wrapper and fitter on each manuscript's own words.

| Manuscript (Gallica id) | Line breaks | Median blank | Coupling (95%) | Its own wrapper | Its own fitter | Last word vs mid-line length |
|---|---|---|---|---|---|---|
| btv1b100737746 | 462 | 0.27 | −0.08 (−0.23 to 0.12) | 0.37 | 0.12 | +13% |
| btv1b55008562q | 1,046 | 0.28 | 0.07 (0.03 to 0.10) | 0.18 | 0.08 | +15% |
| btv1b90069505 | 1,124 | 0.31 | −0.19 (−0.30 to −0.09) | 0.39 | 0.17 | +14% |
| btv1b9009670t | 438 | 0.56 | 0.29 (0.10 to 0.55) | 0.45 | 0.11 | +26% |
| btv1b90595325 | 394 | 0.37 | −0.44 (−0.75 to −0.15) | 0.42 | 0.19 | +21% |
| btv1b90616967 | 416 | 0.33 | 0.09 (−0.03 to 0.18) | 0.51 | 0.23 | +18% |
| **Voynich** | 3,183 | 0.13 | 0.05 (−0.03 to 0.14) | 0.48 | 0.01 | −9% |

Real scribes copying real text do not break the line at the first word that fails to fit. Their line ends carry little or no trace of the next word, as in the Voynich. The coupling test separates mechanical wrapping from scribal line-filling. It does not separate a copied text from one composed on the page.

Reproduce with `python3 control_real/control_real.py <path to the cloned dataset>`.

## What this means

- The Voynich line end looks like the work of a scribe filling lines, which is what scribes of the period did with ordinary texts. It is not evidence that the text was made up on the page, and the earlier version of this report was wrong to read it that way.
- The remaining difference is in how the line is filled. Real scribes end lines on long words, which is what happens when the word straddling the margin is squeezed in. Voynich scribes end lines on short words and start the next line with long ones (+9%).
- The side-tests in sections 1 and 2 (no squeezing of the last word, no split words, no link between words across the line break) still describe the Voynich. They have not been run on a second real corpus with true word boxes, so how unusual they are is open.
- None of this bears on whether the text is meaningful.

## Limits

- The control uses synthetic word boxes (width from character count) and line ends from transcription baselines, so its blanks are noisier than the Voynich's. It covers one dataset, one language, six manuscripts, taken as prose from their transcriptions.
- A line's available right edge is inferred from neighbouring lines of the same paragraph, not observed. The simulations are measured by the same estimator.
- Herbal pages have irregular right edges where drawings intrude, so their intervals are wide.
- Voynich boxes are quantised to 3 px and follow one transliteration.

## Prior work

A background search (limited: voynich.ninja, arXiv and several blogs were unreachable from this environment) found no earlier test of whether the next word would have fitted the measured blank. The closest work found:

- Currier, and Vogt (2012), on the line as a functional unit; Vogt already modelled wrapping in character counts.
- voynich.ninja threads (2026) re-breaking the text greedily in character counts.
- `d2kx/voynich-research`: last-word length against line fullness in glyph counts, and the flat `-m` rate.
- `eo777/vms-spatial` and `vms-dynamic-citation` (2026): word boxes used to test self-citation with a visibility cone; its same-position test was inconclusive.
- Rozanova & Temerev (2026), whose archive supplies the data used here, used the boxes to validate word gaps.

The overhead-copy result therefore confirms earlier null results with a pixel measurement and adds a calibrated limit. The line-break test and its control on real manuscripts were not found elsewhere. Treat that as "not found", not "proven absent".

## Reproduce

```sh
./fetch_data.sh                 # third-party data, pinned commit
pip install numpy scipy pandas matplotlib rapidfuzz
python3 src/exp2_wrap.py        # line-break test            -> results/exp2_wrap.json
python3 src/exp3_lineend.py     # last-word follow-ups       -> results/exp3_lineend.json
python3 src/exp1_kernel.py      # overhead-copy test (~40 min) -> results/exp1_kernel.json
python3 src/figures.py
```

| File | Purpose |
|---|---|
| `src/vdata.py` | Word boxes joined to ZL3b; line and paragraph structure |
| `src/wrap.py` | Line-end measurements and the three simulated scribes |
| `src/pairs.py`, `src/kernel.py` | Word-pair similarity by position; overhead contrast |
| `src/layout.py` | Synthetic pages on the real geometry |
| `control_real/control_real.py` | The same line-end measurements on real manuscripts |

## Data credits

Not redistributed here; `fetch_data.sh` clones them from the archive of L. Rozanova and A. Temerev (`lrozanova/voynich-units`), which keeps the original terms.

- ZL3b transliteration: R. Zandbergen and G. Landini.
- Word boxes: voynichese.com, over T. Takahashi's transliteration.
- Self-citation generator output: T. Timm and A. Schinner's algorithm, as archived there.
- Naibbe cipher text: M. Greshko.
- Latin text: Pliny, *Naturalis Historia*, Tesserae edition.
- Control manuscripts: Gallicorpora, *Données HTR manuscrits du 15e siècle* (transcribed by N. Leroy), cloned separately.
