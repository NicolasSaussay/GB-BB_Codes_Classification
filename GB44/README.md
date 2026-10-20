# GB(4,4) classification

Classification of the generalized bicycle codes `[[2n, k, d]]` built from two polynomials `A(X)`, `B(X)` of Hamming weight 4 over `GF(2)[X]/(X^n - 1)`.

## Contents

| File | Description |
|---|---|
| [algorithms/gb44_all_representatives.m](algorithms/gb44_all_representatives.m) | Magma script: all non-equivalent representatives per parameter `[[2n, k, d]]`. |
| [algorithms/gb44_one_representative.m](algorithms/gb44_one_representative.m) | Magma script: one representative per parameter `[[2n, k, d]]`. |
| classification/gb44_complete_classification_*.md | Complete classification for `n` = 5 to 25 (all representatives), split into six files (below). |
| [classification/gb44_one_representative.md](classification/gb44_one_representative.md) | One representative per parameter for `n` = 26 to 47, plus the codes of length `n` = 60 and 63. |

Complete classification files (`n` = 5 to 25: 308 parameters, 42 613 pairs `(A, B)`):

- [n = 5–19](classification/gb44_complete_classification_n05-n19.md)
- [n = 20–21](classification/gb44_complete_classification_n20-n21.md)
- [n = 22–23](classification/gb44_complete_classification_n22-n23.md)
- [n = 24, part 1](classification/gb44_complete_classification_n24_part1.md) and [part 2](classification/gb44_complete_classification_n24_part2.md) (split by parameter)
- [n = 25](classification/gb44_complete_classification_n25.md)

For `n` = 28, 29 and 30, only one representative per parameter is listed.

## Usage

```
magma -b n:=10 gb44_all_representatives.m                    # interactive mode
magma -b n:=10 slice:=0 nb:=4 gb44_all_representatives.m     # batch mode, slice 0 of 4
magma -b n:=26 gb44_one_representative.m                     # interactive mode
magma -b n:=26 tranche:=0 nb:=4 gb44_one_representative.m    # batch mode, slice 0 of 4
```

`n` is the length of each circulant block (code length `2n`). In `gb44_one_representative.m` the slice parameter is called `tranche`.

## Results format

One collapsible list per length `n`, whose summary gives the best parameters (for each `k`, the greatest `d`); inside, one list per parameter `[[2n, k, d]]` containing its pairs `(A, B)`.
