# studies with both RNA-Seq and Microarray assays of CMAP02 treated drugs

# RNA-Seq    | Microarray    | cell line/tissue  | drugs
# -----------|---------------|-------------------|-----
# GSE115609  |  GSE115564    | MCF7              | fulvestrant
# GSE55347   |  GSE47875     | liver             | estradiol, bezafibrate, clofibrate, clotrimazole, econazole, gemfibrozil, ifosfamide, leflunomide, lovastatin, miconazole, pirinixic acid, rosiglitazone, simvastatin
# GSE41586   |  GSE41364     | HT29              | azacitidine
# GSE43526   |  GSE43899     | cortical neurons  | topotecan (structurally similar to irinotecan and camptothecin)


setwd('data-raw/benchmarks/rnaseq_microarray_cmap')
library(crossmeta)

crossmeta::get_raw('GSE47875')
eset <- crossmeta::load_raw('GSE47875')
