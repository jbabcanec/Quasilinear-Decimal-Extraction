# Submission package

The manuscript and source archive are prepared locally. Nothing has been
submitted to arXiv.

Use `decimal_extraction_source.zip` as the source archive and
`decimal_extraction_paper.tex` as the main file, with PDFLaTeX. It contains
one top-level article and an `anc/` directory of reproduction material.
The inline bibliography needs no BibTeX step. Compilation needs no shell
escape, external conversion, or files outside the archive.

The package omits manuscript PDFs, build logs, auxiliary files, private
research notes, binaries, and old drafts. The source format follows
[arXiv's TeX submission guidance](https://info.arxiv.org/help/submit_tex.html).
Review arXiv's compiled PDF during submission; the local build does not
establish compatibility with every TeX installation.

`metadata.txt` contains the title, author information, abstract, and current
page and table counts. The existing affiliation and email are retained.
The proposed classification is math.NT, with cs.DS as a possible cross-list;
the submitting author makes the final category and license selections.
No license or endorsement has been selected or asserted here.

The main claims and limits are:

- The weighted-block theorem gives O(N log^(5/2) N) bit time with O(N)
  stored bits plus arithmetic workspace. Its native prototype implements
  the exact core and a fixed-guard output path, not the full certified
  procedure analyzed in the theorem.
- The numerical section includes three performance tables and the paper
  retains the published-method complexity comparison. Earlier batches
  and the subsequent certified comparison are identified separately.
- The Machin program is fastest in the latest local batch through ten
  million, with higher memory use than the matched certified Euler control.
  It certifies every returned window. Its O(N log^3 N) attempt bound assumes
  G = O(log N); unrestricted guard retries need further cost accounting.
- Neither a practical crossover for the smaller-bound prototypes nor an
  external speed record is established. Full-prefix MPFR computations are
  controls outside the no-prefix decimal-extraction ranking. Gourdon's
  earlier linear-memory discussion remains part of the prior work.
- The AI-assistance acknowledgment is retained. No independent human
  verification of the mathematical proof is represented.

The 33-file archive includes exact measured native sources, raw datasets,
baseline provenance, and a supplementary reproduction guide. The latest
comparison contributes 27 timed runs with independent interval checks.
Bellard's original source is not redistributed; the historical benchmark
driver obtains it with a recorded hash check. The article build and native
Machin and weighted self-tests require no such download. Optional older
reference programs are retained for reproducibility.

`source_manifest.json` records every archived file's hash.
`verification.json` records the isolated build, table and source checks,
PDF checks, and the separately recorded visual inspection. These establish
the prepared files and their rendering, not acceptance or proof correctness.
