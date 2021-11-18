mapQuery <- function (exp_query, metadata_query, ref_obj, vars = NULL, verbose = TRUE,
                      do_normalize = TRUE, do_umap = TRUE, sigma = 0.1, return_type = c('symphony', 'Seurat')) {
    if (return_type == 'Seurat') {
        que <- Seurat::CreateSeuratObject(
            counts=exp_query,
            meta.data=metadata_query,
            assay='SymphonyQuery'
        )
    }

    if (do_normalize) {
        if (verbose)
            message("Normalizing")
        exp_query = normalizeData(exp_query, 10000, "log")
    }
    if (verbose)
        message("Scaling and synchronizing query gene expression")
    idx_shared_genes = which(ref_obj$vargenes$symbol %in% rownames(exp_query))
    shared_genes = ref_obj$vargenes$symbol[idx_shared_genes]
    if (verbose)
        message("Found ", length(shared_genes), " reference variable genes in query dataset")
    exp_query_scaled = scaleDataWithStats(exp_query[shared_genes,
    ], ref_obj$vargenes$mean[idx_shared_genes], ref_obj$vargenes$stddev[idx_shared_genes],
    1)
    exp_query_scaled_sync = matrix(0, nrow = length(ref_obj$vargenes$symbol),
                                   ncol = ncol(exp_query))
    exp_query_scaled_sync[idx_shared_genes, ] = exp_query_scaled
    rownames(exp_query_scaled_sync) = ref_obj$vargenes$symbol
    colnames(exp_query_scaled_sync) = colnames(exp_query)
    if (verbose)
        message("Project query cells using reference gene loadings")
    Z_pca_query = t(ref_obj$loadings) %*% exp_query_scaled_sync
    if (verbose)
        message("Clustering query cells to reference centroids")
    Z_pca_query_cos = cosine_normalize_cpp(Z_pca_query, 2)
    R_query = soft_cluster(ref_obj$centroids, Z_pca_query_cos,
                           sigma)
    if (verbose)
        message("Correcting query batch effects")
    if (!is.null(vars)) {
        design = droplevels(metadata_query)[, vars] %>% as.data.frame()
        onehot = design %>% purrr::map(function(.x) {
            if (length(unique(.x)) == 1) {
                rep(1, length(.x))
            }
            else {
                stats::model.matrix(~0 + .x)
            }
        }) %>% purrr::reduce(cbind)
        Xq = cbind(1, intercept = onehot) %>% t()
    }
    else {
        Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))),
                    sparse = TRUE)
    }
    Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), as.matrix(Xq),
                              as.matrix(R_query), as.matrix(ref_obj$cache[[1]]), as.matrix(ref_obj$cache[[2]]))
    colnames(Z_pca_query) = row.names(metadata_query)
    rownames(Z_pca_query) = paste0("PC_", seq_len(nrow(Zq_corr)))
    colnames(Zq_corr) = row.names(metadata_query)
    rownames(Zq_corr) = paste0("harmony_", seq_len(nrow(Zq_corr)))
    umap_query = NULL
    if (do_umap & !is.null(ref_obj$save_uwot_path)) {
        if (verbose)
            message("UMAP")
        ref_umap_model = uwot::load_uwot(ref_obj$save_uwot_path,
                                         verbose = FALSE)

        ## UMAP may have been learned on subset of columns
        umap_query = uwot::umap_transform(t(Zq_corr)[, 1:ref_umap_model$norig_col], ref_umap_model)
        #         umap_query = uwot::umap_transform(t(Zq_corr), ref_umap_model)
        colnames(umap_query) = c("UMAP1", "UMAP2")
        rownames(umap_query) <- row.names(metadata_query)
    }
    if (verbose)
        message("All done!")

    if (return_type == 'Seurat') {
        que@assays$SymphonyQuery@data <- exp_query
        que@assays$SymphonyQuery@scale.data <- exp_query_scaled_sync
        que[['pca']] <- Seurat::CreateDimReducObject(
            embeddings = t(Z_pca_query),
            loadings = ref_obj$loadings,
            stdev = as.numeric(apply(Z_pca_query, 1, stats::sd)),
            assay = 'SymphonyQuery',
            key = 'pca_'
        )
        que[['harmony']] <- Seurat::CreateDimReducObject(
            embeddings = t(Zq_corr),
            stdev = as.numeric(apply(Zq_corr, 1, stats::sd)),
            assay = 'SymphonyQuery',
            key = 'harmony_',
            misc=list(R=R_query)
        )
        que <- Seurat::ProjectDim(que, reduction = 'harmony', overwrite = TRUE, verbose = FALSE)
        if (do_umap) {
            que[['umap']] <- Seurat::CreateDimReducObject(
                embeddings = umap_query,
                assay = 'SymphonyQuery',
                key = 'umap_'
            )
        }
        return(que)
    } else if (return_type == 'symphony') {
        return(list(Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query,
                    Xq = Xq, umap = umap_query, meta_data = metadata_query))
    } else {
        stop(glue('The return type = \"{return_type}\" is not available.'))
    }

}

environment(mapQuery) <- environment(symphony::mapQuery)
