## Cloneman
### a package to parse and manage information from Copy Number calling and Clone Inference tools.

There exist many programs that are used to infer the copy number aberrations in samples
They typically infer the copy number segments (either as integer numbers, e.g. 1, 3, 4 copies etc)
or as log ratios to the "normal" coverage.
Each program contains produces its own specific output format, in one or multiple files.

This package was conceived in an attempt to normalize these outputs in a standard R format
with a set of predetermined columns and objects that can be used for downstream analysis.


## Installation

```
devtools::install_github(c('kmavrommatis/cloneman'),ref='main')
```
