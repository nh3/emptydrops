# emptydrops

## Disclaimer:
All code originally comes from https://github.com/10XGenomics/cellranger with
minimal modifications for packaging and running under python3.

## Usage:
```
import emptydrops

emptydrops.find_nonambient_barcodes(matrix, orig_cell_bcs, min_umi_frac_of_median=0.01, min_umis_nonambient=500, max_adj_pvalue=0.01)
```
