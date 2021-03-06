#' UI for docs page (add new sections here)
#'
#' @inheritParams scPageUI
#'
#' @return shiny.tag with html for docs tab
#'
#' @export
docsPageUI <- function(id, tab, active) {

  # Single Cell section ----
  singleCellSection <- docsSection(
    id = 'sc-docs', name = 'Single Cell',
    content = list(
      # Workflow ----
      docsSubsection(id = 'sc-workflow',
                     name = 'Overview',
                     content = tagList(
                       HTML(
                         "<p><a href='https://pachterlab.github.io/kallisto/about'>kallisto</a> is used for
                         quantification and workflow recommendations from <a href='https://osca.bioconductor.org/'>osca.bioconductor.org</a> are followed.
                         Choices are described below.</p>"
                       )
                     )

      ),
      # Feature selection and clustering ----
      docsSubsection(
        id = 'sc-preprocessing',
        name = 'Feature selection and clustering',
        content = tagList(
          HTML("
        <p>Highly variable genes are taken as those within the top 10% of biological variance. A custom set of genes can also be specified
        when subsetting by clicking the <span class='bs-docs-btn'><i class='fa fa-upload fa-fw'></i></span> button. The uploaded file
        should have one HGNC symbol per line. The rationale and approach for doing so is described
        <a href='http://bioconductor.org/books/release/OSCA/feature-selection.html#apriori-hvgs'>here</a>.</p>
        <p><code>igraph::cluster_walktrap</code> is used for datasets with less than 10,000 cells. For larger datasets, <code>igraph::cluster_louvain</code> is used for speed.
        <a href='https://osca.bioconductor.org/dimensionality-reduction.html#based-on-population-structure'>getClusteredPCs</a> is used to pick the number of principle components to use for clustering. This is quite slow but generally provides satisfactory clustering.
        </p>")
        )
      ),
      # Automatic cluster labeling ----
      docsSubsection(
        id = 'sc-label-auto',
        name = 'Automatic cluster labeling',
        content = tagList(
          HTML("
        <p>To automate cell cluster labeling, labels are transfered from a previously labeled dataset to the current dataset using <a href='https://bioconductor.org/packages/release/bioc/html/SingleR.html'>SingleR</a>.
        If the two datasets share an ancestor through subsetting, cluster labels are transfered by using the most frequent label in the reference dataset. Because this is substantially faster, these datasets are the first choices shown.</p>

        <p>To transfer labels, first <b>Select a dataset</b> to transfer labels to and click the label transfer toggle.
        Next, select another single cell dataset to <b>Transfer labels from</b>. If you are satisfied, go ahead and overwrite the previous labels.</p>
             "),
          div(class = 'bs-callout bs-callout-info',
              h4('Mix automatic and manual labeling'),
              p('Automatic cell-type labeling is often a good starting point for further refinements through manual labeling.')
          )
        )
      ),
      # Integrate datasets ----
      docsSubsection(
        id = 'sc-integrate-samples',
        name = 'Dataset integration',
        content = tagList(
          HTML("
             <p>Integrating single cell samples allows for cross-sample differential expression analyses, differential abundance analyses, pathway analyses, and drug queries.
             To integrate single cell samples, click the dataset integation toggle and fill in the inputs. Available integration types include <a href='https://github.com/immunogenomics/harmony'>harmony</a>
             and <a href='https://osca.bioconductor.org/integrating-datasets.html#performing-mnn-correction'>fastMNN</a>. You can select both in order to compare results.</p>
             <p>Integrated datasets with at least three samples will use the same <code>limma</code> pipeline used for bulk analyses with <a href='https://osca.bioconductor.org/multi-sample-comparisons.html#differential-expression-between-conditions'>pseudobulk</a> expression profiles.
             Aggregative methods <a href='https://www.biorxiv.org/content/biorxiv/early/2019/07/26/713412.full.pdf'>outperform</a> methods that treat each cell as an independant replicate.</p>

             <p>You can also download, fill out, and upload a sample pairs csv in order to model for an additional variable (e.g. individual) for pseudobulk analyses. The additional
             variable is also added as a covariate for the <code>vars_use</code> argument to <code>harmony::HarmonyMatrix</code> during integration.</p>

             "),
          div(class = 'bs-callout bs-callout-danger',
              h4('Identify incorrectly integrated clusters'),
              HTML('<p>After integrating datasets, you should look out for clusters that were incorrectly grouped (e.g. RBCs from one sample
              cluster with B-cells from another sample). One way to do this is to select the newly integrated dataset and inspect the source clusters for each new integrated cluster by selecting the <b>labels</b> toggle
                   from the the <b>Perform comparisons between:</b> choices.</p><p>Another useful exercise is to select the <code>is_test</code>, <code>is_ctrl</code>, and individual sample
                   names from the <b>Feature to plot</b> input in order to inspect the extent of mixing. Uneven mixing may reflect changes in cell-type abundances or issues with integration. Trying
                   multiple integration methods can be helpful for deciding between the two.'),
          ),
          div(class = 'bs-callout bs-callout-danger',
              h4('Exclude incorrectly integrated clusters'),
              HTML('<p>If cell labels have changed with integration, you will have to decide if the new labels
                 are correct by inspecting markers. In my experience, a cell type that exists in one sample but not another can cause issues. If you discover clusters that failed to integrate correctly, re-run integration and exclude any offending clusters.</p>')
          )
        )
      ),
      # Excluding ambient genes ----
      docsSubsection(
        id = 'sc-differential-expression',
        name = 'Excluding ambient genes',
        content = tagList(

          HTML("<p>Single cell samples are often plagued by contamination from ruptured cells (seen for example as high baseline expression of
                 RBC specific genes in non RBCs). We originally used <code>soupx</code> to correct for this ambient expression but found that correction
                 would often result in misleading differences between samples. Ambient genes are now identified as those that are outliers in
                 empty droplets with less than ten counts and can be excluded by clicking <span class='bs-docs-btn'><i class='fa fa-ban fa-fw'></i></span>.
                 This simple approach doesn't alter the expression profile and so avoids misleading corrections. This approach is also suggested
               by the authors of the <a href='https://osca.bioconductor.org/multi-sample-comparisons.html#ambient-problems'>osca</a> handbook.</p>
               <p>We additionally only call a gene as ambient if the sample it is in would change the direction of differential expression. For example,
               if a gene is up-regulated in test samples and called as ambient in one or more control samples, then it won't be excluded. The high background
               expression in the control samples results in a smaller difference in expression but does not change the direction of change.</p>")
        )
      )
    ))


  # Bulk Data section ----
  datasetsSection <- docsSection(
    id = 'ds-docs', name = 'Bulk Data',
    content = list(
      # Workflow ----
      docsSubsection(id = 'bulk-overview',
                     name = 'Overview',
                     content = tagList(
                       HTML(
                         "<p><a href='https://pachterlab.github.io/kallisto/about'>kallisto</a> is used for
                         quantification, <a href='https://bioconductor.org/packages/release/bioc/html/sva.html'>surrogate variable analysis</a> is used to discover and account for sources of variation,
                         <code>limma</code> is used for fixed and mixed effect differential expression analysis as well as pathway analysis
                         with <code>cameraPR</code>. Standardized effect sizes (moderated t-statistics in units of standard deviation) are reported
                         and are amenable to meta-analyses."
                       )
                     )

      ),
      # Quantifying RNA-seq datasets ----
      docsSubsection(
        id = 'ds-quantification',
        name = 'Quantifying RNA-seq datasets',
        content = tagList(
          HTML("<p>To add a new dataset, first make sure the <b>Dataset name</b> input is empty. Then type
                 a name for your new bulk dataset and press enter. You will be prompted to select a subfolder containing fastq.gz files.</p>"
          ),
          div(class = 'bs-callout bs-callout-danger',
              h4('Expand subfolders on left only'),
              p('Subfolders can only be expanded on the left side of the popup when selecting a folder with fastq.gz files.')
          ),
          p("If a folder with bulk fastq.gz files is selected, you will be prompted to confirm the auto-detected end type.
              If the end type is pair-ended, you will also need to indicate pairs. To do this, click the table rows corresponding to the
              paired files and click Paired. If present, you must also indicate any files that are replicates of the same sample
              before you run quantification.")
        )),
      # Sample-level gene expression ----
      docsSubsection(
        id = 'ds-gene-plots',
        name = 'Sample-level gene expression',
        content = tagList(p("Sample level gene expression plots allow you to inspect normalized expression values for genes of
                              interest split by groups of interest. These plots allow you to see expression values for all
                              samples in order to avoid being misled by summary statistics."),
                          div(class = 'bs-callout bs-callout-info',
                              h4('Expression normalization algorithm for visualizations'),
                              HTML("<p>Gene expression is normalized with the function <code>DESeq2::rlog</code> with <code>blind = FALSE</code> for RNA-Seq data.
                                     For further details see
                                     <a href='http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization' target='_blank'>Data transformations and visualization</a>
                                     in the DESeq2 vignette. For microarray data, <a href='https://github.com/alexvpickering/crossmeta'>crossmeta</a> uses `affy::rma` or `oligo::rma` for Affymetrix platforms and
                                   `limma::neqc` for Agilent and Illumina platforms. These normalized gene expression values are also used for MDS plots and cell-type deconvolution.</p>")
                          )
        )
      ),
      # Cell-type deconvolution ----
      docsSubsection(
        id = 'ds-cell-plots',
        name = 'Cell-type deconvolution',
        content = tagList(p("Cell type deconvolution uses gene expression from a single-cell dataset to estimate the proportion
                               of different cell types in bulk RNA-Seq samples. Large differences in cell-type proportions may
                               confound differential expression but also point to important biology."),
                          div(class = 'bs-callout bs-callout-danger',
                              h4('Requires comparable single-cell dataset'),
                              HTML("<p>Cell-type deconvolution requires a single-cell dataset that contains at least the cell types that you expect
                                     to observe in the bulk dataset that you are deconvoluting.</p>")
                          ),
                          HTML("<p>To run cell-type deconvolution, select a bulk dataset and click the <span class='bs-docs-btn'><i class='fa fa-object-ungroup fa-fw'></i></span></span> button to show inputs for cell
                               type deconvolution.</p>
                               <p>Next, select a reference single cell dataset that contains the cell types that you expect to observe in the bulk
                                 dataset that you are deconvoluting. Also select the clusters to use for deconvolution. You should only select cluster types
                                 that you expect within your bulk dataset.</p>"),
                          div(class = 'bs-callout bs-callout-info',
                              h4('Cell-type deconvolution algorithm'),
                              HTML("<p>Cell-type deconvolution is performed using <code>dtangle</code> as described in the
                                     <a href='https://github.com/gjhunt/dtangle/blob/master/vign/sc_vignette.md' target='_blank'>single-cell vignette</a>
                                     except using <code>DESeq2::rlog</code> normalized and surrogate variable/pair adjusted bulk counts. For the reference single cell dataset,
                                     log normalized counts are used. See the
                                     <a href='https://academic.oup.com/bioinformatics/article/35/12/2093/5165376?guestAccessKey=ac40b15d-bec0-48c1-be94-fbef567f63ec' target='_blank'>dtangle manuscript</a>
                                     for further details.</p>")
                          )

        ))
    ))

  # Drugs section ----
  drugsSection <- docsSection(
    id = 'drugs-docs', name = 'Drugs',
    content = list(
      # Overview ----
      docsSubsection(
        id = 'drugs-overview',
        name = 'Overview',
        content = tagList(
          HTML("
        <p>CMAP02/L1000 drug and genetic perturbations are matched to a reference differential expression signature based on pearson correlation with the top 200 query signature genes.
        The most negatively correlated drug signatures are predicted to oppose the disease query signature and may thus act as potential therapeutics. The most strongly correlated (positive
        or negative) genetic perturbations may reflect an importance of the perturbed gene to the disease.
        </p>
             "),
          div(class = 'bs-callout bs-callout-info',
              h4('Similarity metrics and data preparations'),
              HTML("In our benchmarks where we use external data that assayed drugs also measured by CMAP02/L1000, pearson correlation and cosine similarity on <code>limma</code> differential
            expression values outperform other approaches including Kolmogorov-Smirnov, eXtreme Sum, and cosine similarity performed on characteristic direction differential expression
            values.")
          )
        )

      ),
      # Custom queries -----
      docsSubsection(
        id = 'drugs-custom-queries',
        name = 'Custom queries',
        content = tagList(
          HTML("<p>Custom queries allow you to specify a custom query signature.
          To run a custom query, toggle the custom signature button, provide a name to save the custom query as, then click the <span class='bs-docs-btn'><i class='fa fa-upload fa-fw'></i></span>
             button to upload a csv with your query signature. The first column of the uploaded csv must be unique HGNC symbols. The uploaded csv must also include either 1) a column named
               one of dprime or logFC, or 2) a column of effect size values as it's second column. This query signature can then be selected to view the results as normal.</p>"),
          div(class = 'bs-callout bs-callout-info',
              h4('Custom query differences'),
              HTML("For custom queries, all uploaded genes that were also measured by L1000/CMAP02 are used. In contrast, query signatures from the Single-Cell or Bulk tabs
                   use only the top 200 most differentially expressed genes.")
          ),

        )
      )
    ))



  # list of sections ----
  docsSections <- list(
    # Single cell documentation
    singleCellSection,

    # Datasets documentation
    datasetsSection,

    # Drugs documentation
    drugsSection
  )


  section_info <- extract_section_info(docsSections)


  withTags({
    tabPane(tab, active,
            div(class = "bs-docs-container",
                div(class = 'row',
                    div(class = 'col-md-6 col-md-offset-1', role="main", docsSections),
                    div(class = "col-md-4 col-md-offset-1 bd-toc hidden-xs hidden-sm", role="complementary", docsSideNav(section_info))
                )
            )
    )
  })
}

#' A primary section in the Docs page
#'
#' @keywords internal
#' @noRd
docsSection <- function(id, name, content) {

  withTags({
    div(class = 'bs-docs-section',
        h1(id = id, class='page-header',
           tags$a(href = paste0('#', id)),
           name
        ),
        content
    )
  })

}

#' A subsections of a docsSection
#'
#' @keywords internal
#' @noRd
docsSubsection <- function(id, name, content) {
  withTags({
    tagList(
      h2(id = id,
         tags$a(href = paste0('#', id)),
         name

      ),
      content
    )
  })
}

#' Navigation on right side of Docs page
#'
#' @keywords internal
#' @noRd
docsSideNav <- function(section_info) {

  section_id_names <- section_info$section_id_names
  subsection_id_names <- section_info$subsection_id_names

  nsections <- length(section_id_names)

  withTags({
    tags$nav(class = "bs-docs-sidebar well-form well-bg",
        tags$ul(class="nav bs-docs-sidenav",
           lapply(seq_len(nsections), function(i) {

             # get id and name for current section
             section_id_name <- section_id_names[[i]]

             # get id and names list of subsections for current section
             subsectioni_id_names <- subsection_id_names[[i]]

             tags$li(
               # link for current section
               a(href=paste0("#", section_id_name[1]), section_id_name[2]),
               tags$ul(class="nav",
                  # link for each subsection of current section
                  lapply(subsectioni_id_names, function(subsectioni_id_name) {
                    tags$li(
                      a(href = paste0("#", subsectioni_id_name[1]), subsectioni_id_name[2])
                    )
                  }
                  )
               )
             )
           })
        ))
  })
}


#' Extracts info needed to construct docsSideNav
#'
#' @keywords internal
#' @noRd
extract_section_info <- function(docsSections) {
  section_id_names <- lapply(docsSections, function(docsSection) {

    # title is in first child
    # get id and name
    id <- docsSection$children[[1]]$attribs$id
    name <- docsSection$children[[1]]$children[[2]]

    return(c(id, name))
  })

  subsection_id_names <- lapply(docsSections, function(docsSection) {

    # subsections are in the second child
    subsections <- docsSection$children[[2]]

    lapply(subsections, function(subsection) {
      # subsection titles are in the first child
      id <- subsection[[1]]$attribs$id
      name <- subsection[[1]]$children[[2]]

      return(c(id, name))
    })

  })

  return(list(
    section_id_names = section_id_names,
    subsection_id_names = subsection_id_names
  ))
}
