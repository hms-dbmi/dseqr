mock_scseq <- function() {

    sce <- scuttle::mockSCE()
    sce <- scuttle::logNormCounts(sce)

    sce@metadata$mrna <- row.names(sce)[1:10]
    sce@metadata$rrna <- row.names(sce)[11:20]

    sce <- add_scseq_qc_metrics(sce, for_qcplots = TRUE)

    return(sce)
}


test_that("evaluating a custom metric based on colData works", {
    scseq <- mock_scseq()
    cutoff <- stats::median(scseq$log10_sum)
    metric <- paste0('log10_sum>', cutoff)

    res <- evaluate_custom_metric(metric, scseq)

    # used to check that wasn't error
    expect_s4_class(res, 'DFrame')

    # same cell order
    expect_equal(row.names(res), colnames(scseq))

    # column name is metric
    expect_equal(colnames(res), metric)

    # value is correct
    expect_equal(res[[1]], scseq$log10_sum > cutoff)
})

test_that("evaluating a custom metric based on a single gene works", {
    scseq <- mock_scseq()
    expr <- SingleCellExperiment::logcounts(scseq)

    gene <- row.names(expr)[1]
    gene_cutoff <- stats::median(expr[gene, ])
    log10_sum_cutoff <- stats::median(scseq$log10_sum)
    metric <- paste0('log10_sum>', log10_sum_cutoff, '&', gene, '<', gene_cutoff)

    res <- evaluate_custom_metric(metric, scseq)

    # used to check that wasn't error
    expect_s4_class(res, 'DFrame')

    # same cell order
    expect_equal(row.names(res), colnames(scseq))

    # column name is metric
    expect_equal(colnames(res), metric)

    # value is correct
    expect_equal(res[[1]], unname(expr[gene, ] < gene_cutoff & scseq$log10_sum > log10_sum_cutoff))
})


test_that("evaluating a custom metric works with duplicate cell ids", {
    scseq <- mock_scseq()
    colnames(scseq)[1:2] <- 'Cell_001'
    expect_equal(colnames(scseq)[1], colnames(scseq)[2])

    expr <- SingleCellExperiment::logcounts(scseq)

    gene <- row.names(expr)[1]
    gene_cutoff <- stats::median(expr[gene, ])
    log10_sum_cutoff <- stats::median(scseq$log10_sum)
    metric <- paste0('log10_sum>', log10_sum_cutoff, '&', gene, '<', gene_cutoff)

    res <- evaluate_custom_metric(metric, scseq)

    # used to check that wasn't error
    expect_s4_class(res, 'DFrame')

    # same cell order
    expect_equal(row.names(res), colnames(scseq))

    # column name is metric
    expect_equal(colnames(res), metric)

    # value is correct
    expect_equal(res[[1]], unname(expr[gene, ] < gene_cutoff & scseq$log10_sum > log10_sum_cutoff))
})


test_that("evaluating a nonsensical metric throws an error", {
    scseq <- mock_scseq()
    cutoff <- stats::median(scseq$log10_sum)
    metric <- paste0('log10_sum>')

    expect_error(evaluate_custom_metric(metric, scseq))
})
