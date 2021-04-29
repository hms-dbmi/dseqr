# dseqr 0.13.0

* Added a `NEWS.md` file to track changes to the package.
* Removed `LoomExperiment` in favor of `qs::qsave(..., preset = 'fast')`
* Uses `presto` to calculate markers for speed.
* Implements flexible integrated group specification/comparison.
* Uses `scDblFinder` for detecting doublets.
* Uses UMAP for larger datasets (>5000 cells).
