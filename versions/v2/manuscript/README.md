# The tea-tiling coordination sequence

The note proves Colin Barker's February 11, 2018 generating-function and recurrence conjectures in the FORMULA fields of [OEIS A299285](https://oeis.org/A299285) and [A299286](https://oeis.org/A299286). It also establishes Stefano Spezia's conditional closed form in A299285, proves the companion tcd formulas in [A299287](https://oeis.org/A299287) and [A299288](https://oeis.org/A299288), and gives least eventual recurrence orders and least quasipolynomial periods. For arbitrary orientations of the alternating special and prism slabs, the tcd and tea counts are sharp lower and upper bounds. Their sphere-count gap is exactly `floor(n^2/8)`; Corollary 6 gives the exact orientation correction and equality conditions.

## Contents

- `note.tex`: self-contained LaTeX source, with its bibliography included.
- `metadata.json`: author, license, problem and mathematical-coverage data.
- `coverage.json`: source locators and exact mathematical coverage.
- `verify.py`: standard-library Python integer checks.
- `verification.json`: the corresponding finite-check output.

The source archive contains these files and this README. The PDF is supplied separately.

## Build

Use Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively, run `pdflatex -interaction=nonstopmode note.tex` twice. The bibliography is in the LaTeX source; no external bibliography database is required.

## Finite checks

Python 3.10 or later requires no additional packages:

```sh
python3 verify.py --radius 32 > verification.json
```

The JSON records exact finite ranges. Breadth-first search checks the graph distances and sphere counts through radius 32 for both tea and tcd independently of the closed-form metrics, visiting 90,153 tea vertices. Further checks compare the 43 A299285, 39 A299286, 44 A299287 and 40 A299288 tabulated values, the layer sums and generating-function coefficients through index 1000, small integer hexagons, and tea graph edges and vertex-transitivity maps.

The arbitrary-orientation checks cover all 256 sign assignments on the eight even boundaries from -8 through 6. For each graph they compare distances, exact sphere corrections, sharp sphere and ball bounds and sphere equality conditions through radius eight, totaling 413,696 vertices. The tea–tcd sphere and ball gaps are also checked through radius 1000. The all-radius conclusions follow from the proofs in the note.

## Mathematical coverage

The pinned [TeaCoordinationSequence.lean](https://github.com/the-omega-institute/trureturing/blob/e7f7c25c7de35e89418be730f50ec8fb94949133/D5/S3/StatisticalMechanics/TeaCoordinationSequence.lean) defines the tea coordinate graph and its spheres. The declaration `D5.S3.StatisticalMechanics.TeaCoordinationSequence.result` proves its four positive-radius sphere-count identities. The geometric identification, generating functions, ball count, least-order and least-period results, all tcd conclusions and the arbitrary-orientation bounds have written proofs in the note. Uniform tilings other than tea and tcd in the common OEIS cross-reference list are outside its scope.

## Authors and license

1. Haobo Ma — ChronoAI Pte Ltd, Singapore.
2. Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang. [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/) (`cc-by-4.0`).
