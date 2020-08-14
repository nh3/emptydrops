# emptydrops
Python implementation of emptydrops-like cell calling as in CellRanger v3.0.2

## Disclaimer:
All code originally comes from https://github.com/10XGenomics/cellranger with
minimal modifications for packaging and running under python3.

## Usage:
```
from emptydrops import find_nonambient_barcodes
from emptydrops.matrix import CountMatrix

matrix = CountMatrix.from_legacy_mtx(mtx_dir)

find_nonambient_barcodes(
    matrix,          # Full expression matrix
    orig_cell_bcs,   # (iterable of str): Strings of initially-called cell barcodes
    min_umi_frac_of_median=0.01,
    min_umis_nonambient=500,
    max_adj_pvalue=0.01
)
```

Returns:
```
[
    'eval_bcs',      # Candidate barcode indices in addition to those in `orig_cell_bcs` (n)
    'log_likelihood',# Ambient log likelihoods (n)
    'pvalues',       # pvalues (n)
    'pvalues_adj',   # B-H adjusted pvalues (n)
    'is_nonambient', # Boolean nonambient calls (n)
]
```
