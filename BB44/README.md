# BB(4,4) search

Bivariate bicycle codes `[[N, k, d]]` of length `N = l*m`, built from two polynomials `A(X, Y)`, `B(X, Y)` of four terms each over `GF(2)[X, Y]/(X^l - 1, Y^m - 1)`, with `H_X = [A | B]` and `H_Z = [B^T | A^T]`.

## Contents

| File | Description |
|---|---|
| [algorithms/bb44_random_search_dmin.m](algorithms/bb44_random_search_dmin.m) | Magma script: random search of codes and computation of their minimum distance. |
| [parameters/bb44_best_codes.md](parameters/bb44_best_codes.md) | Best codes found: 33 codes for 11 lengths `N` (24 to 144), with `(l, m)`, `A` and `B`. |

## Usage

Edit `l`, `m`, `dimension_threshold` and `num_trials` at the top of the script, then run:

```
magma bb44_random_search_dmin.m
```

Codes of dimension at most `dimension_threshold` are discarded; `num_trials` is the number of valid draws. The output lists, for each `k` then each `d`, the parameters `[[N, k, d]]` with an example `(A, B)`.

## Results format

One collapsible list per length `N`; inside, one list per parameter `[[N, k, d]]` containing the codes found, each with its `l`, `m`, `A` and `B`. The search is random, so the lists are not exhaustive.
