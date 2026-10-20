# GB(3,3) classification

Classification of the generalized bicycle codes `[[2n, k, d]]` built from two polynomials `A(X)`, `B(X)` of Hamming weight 3 over `GF(2)[X]/(X^n - 1)`.

## Contents

| File | Description |
|---|---|
| [algorithms/gb33_classification.m](algorithms/gb33_classification.m) | Magma script: one representative per equivalence class of pairs `(A, B)`, with `k` and `d`. |
| [classification/gb33_classification.md](classification/gb33_classification.md) | Results: 38 lengths `n`, 380 parameters `[[2n, k, d]]`, 9 525 pairs `(A, B)`. |

## Usage

```
magma -b n:=10 gb33_classification.m                      # interactive mode
magma -b n:=10 slice:=0 nb:=4 gb33_classification.m       # batch mode, slice 0 of 4
```

`n` is the length of each circulant block (code length `2n`); `slice` and `nb` split the list of representatives for parallel runs.

## Results format

One collapsible list per length `n`, whose summary gives the best parameters (for each `k`, the greatest `d`); inside, one list per parameter `[[2n, k, d]]` containing its pairs `(A, B)`.
