# dseqr 0.20.0
* switch to `picker` for grid plots

# dseqr 0.19.0
* switch to `picker` for scatter plots
* add support for .h5 uploads

# dseqr 0.16.0
* updated `scDblFinder` to 1.7.4
* improved upload modal
* support for all ensembl species
* switch to `HDF5Array::TENxMatrix` for fast gene exploration without loading time
* see [convert script](data-raw/convert/0.16.0.R) to update data generated with `[v0.15.0-0.16.0)`
*

# dseqr 0.15.0
* Added `human_lung` reference from azimuth
* can zoom plots with brush and double-click
* gene table has filter for each column
* see [convert script](data-raw/convert/0.15.0.R) to update data generated with `[v0.14.1-0.15.0)`


# dseqr 0.14.18
* pathway analysis is faster and GO terms are clustered


# dseqr 0.14.1
* can delete single-cell datasets
* improved upload interface
* added cluster-free differential expression plot
* uses `dgRMatrix` for faster access to logcounts
* see [convert script](data-raw/convert/0.14.1.R) to update data generated with `v[0.14.0-0.14.1)`


# dseqr 0.14.0
* `logcounts` and `SingleCellExperiment` are saved separately for initial load speed
* see [convert script](data-raw/convert/0.14.0.R) to update data generated with `v[0.13.3-0.14.0)`


# dseqr 0.13.3
* control and test sample can be changed through UI
* added cluster-free differential abundance plot
* see [convert script](data-raw/convert/0.13.3.R) to update data generated with `v[0.13.1-0.13.3)`


# dseqr 0.13.1
* Restored 1v1 cluster comparisons now using AUC and percent differences.


# dseqr 0.13.0
* Added a `NEWS.md` file to track changes to the package.
* Removed `LoomExperiment` in favor of `qs::qsave(..., preset = 'fast')`
* Uses `presto` to calculate markers for speed.
* Implements flexible integrated group specification/comparison.
* Uses `scDblFinder` for detecting doublets.
* Uses UMAP for larger datasets (>5000 cells).
* see [convert script](data-raw/convert/0.13.0.R) to update data generated with `v[0.10.0-0.13.0)`
