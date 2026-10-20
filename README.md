# GB-BB Codes Classification

Classification of generalized bicycle (GB) quantum codes with weight-3 and weight-4 polynomials, GB(3,3) and GB(4,4), and search results for bivariate bicycle (BB) codes BB(4,4). The repository contains the Magma algorithms and the resulting lists of codes.

## About GB and BB codes

A generalized bicycle code GB(w, w) is a CSS quantum code `[[2n, k, d]]` built from two binary circulant matrices, defined by two polynomials `A(X)` and `B(X)` of Hamming weight `w` over `GF(2)[X]/(X^n - 1)`. A bivariate bicycle code BB(4,4) is a CSS code of length `N = 2lm` built from two polynomials `A(X, Y)` and `B(X, Y)` with four terms each over `GF(2)[X, Y]/(X^l - 1, Y^m - 1)`.

## Classification methodology

- **GB(3,3) and GB(4,4).** For each length `n`, pairs `(A, B)` are reduced to one representative per equivalence class (cyclic shifts, multipliers `X -> X^m`, swap `A <-> B`), using canonical integer bitmasks. The dimension `k` and the minimum distance `d` are then computed with Magma.
- **BB(4,4).** Pairs `(A, B)` are drawn at random for a given `(l, m)`; the minimum distance of each code is computed with Magma, and the codes are grouped by parameters `[[N, k, d]]`. This search is not exhaustive.

## Repository structure

```
README.md
LICENSE                          MIT license (code)
LICENSE-data                     CC-BY-4.0 notice (results)
GB33/
├── README.md
├── algorithms/gb33_classification.m
└── classification/gb33_classification.md
GB44/
├── README.md
├── algorithms/
│   ├── gb44_all_representatives.m
│   └── gb44_one_representative.m
└── classification/
    ├── gb44_complete_classification_n05-n19.md
    ├── gb44_complete_classification_n20-n21.md
    ├── gb44_complete_classification_n22-n23.md
    ├── gb44_complete_classification_n24_part1.md
    ├── gb44_complete_classification_n24_part2.md
    ├── gb44_complete_classification_n25.md
    └── gb44_one_representative.md
BB44/
├── README.md
├── algorithms/bb44_random_search_dmin.m
└── parameters/bb44_best_codes.md
```

| Folder | Algorithms | Results |
|---|---|---|
| [GB33](GB33) | [gb33_classification.m](GB33/algorithms/gb33_classification.m) | [gb33_classification.md](GB33/classification/gb33_classification.md) |
| [GB44](GB44) | [gb44_all_representatives.m](GB44/algorithms/gb44_all_representatives.m), [gb44_one_representative.m](GB44/algorithms/gb44_one_representative.m) | [GB44/classification](GB44/classification) |
| [BB44](BB44) | [bb44_random_search_dmin.m](BB44/algorithms/bb44_random_search_dmin.m) | [bb44_best_codes.md](BB44/parameters/bb44_best_codes.md) |

## How to use / reproduce

[Magma](http://magma.maths.usyd.edu.au/) is required. Each script documents its usage and parameters in its header. For example, for GB(3,3):

```
magma -b n:=10 gb33_classification.m
```

The scripts of the GB families also have a parallel batch mode (see the header of each script). The BB(4,4) script is configured by editing `l`, `m`, `dimension_threshold` and `num_trials` at the top of the file. See the README of each folder for details.

## Results

The result files are Markdown files with collapsible lists: one list per length (`n` for GB, `N` for BB), then one list per parameter `[[N, k, d]]`, containing the pairs `(A, B)` in a code block. For GB codes, the summary of each length gives its best parameters: for each dimension `k`, the parameters of greatest distance `d`.

## Licenses

- Code (`.m` files): MIT, see [LICENSE](LICENSE).
- Results (`.md` files): CC-BY-4.0, see [LICENSE-data](LICENSE-data).

## Citation

If you use this repository, please cite the associated paper (manuscript in preparation):

```bibtex
@unpublished{AGS2O26,
  title  = {{GB(3,3)} and {GB(4,4)} Codes: Exhaustive Classification and Comparison with Bivariate Bicycle Codes},
  author = {Arnault, Fran{\c{c}}ois and Gaborit, Philippe and Saussay, Nicolas},
  note   = {Manuscript in preparation},
  year   = {2026}
}
```

Arnault, Gaborit and Saussay are with XLIM, UMR 7252, Université de Limoges. The full citation will be added once the paper is published.
