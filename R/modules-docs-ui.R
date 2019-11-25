# Datasets section ----
datasetsSection <- docsSection(
  id = 'ds-docs', name = 'Datasets',
  content = list(
    # Quantifying RNA-seq datasets ----
    docsSubsection(
      id = 'ds-quantification',
      name = 'Quantifying RNA-seq datasets',
      content = tagList(
        HTML("<p>To add a new dataset, first make sure the dataset name input is empty. Then type
                 a name for your new bulk or single-cell dataset and hit enter. You will be prompted to select a subfolder of the <b>bulk</b>
                 folder containing fastq.gz files or a subfolder of the <b>single-cell</b> folder containing either fastq.gz or Cell Ranger files.</p>"
        ),
        img(src = 'docs-ds-dataset_name.png', class = 'bs-docs-img'),
        div(class = 'bs-callout bs-callout-danger',
            h4('Folder selection'),
            p('Subfolders can only be selected on the left side of the popup when selecting a folder with fastq.gz or cellranger files.')
        ),
        p("If a folder with bulk fastq.gz files is selected, you will be prompted to confirm the auto-detected end type.
              If the end type is pair-ended, you will also need to indicate pairs. To do this, click the table rows corresponding to the
              paired files and click Paired. If present, you should also indicate any files that are replicates of the same sample
              before you run quantification."),
        p("If a folder with single-cell fastq.gz files is selected, nothing further is necessary before running quantification.
              Quantification uses",
          tags$a(href = 'https://pachterlab.github.io/kallisto/about', target = '_blank', 'kallisto'),
          "for both bulk and single-cell RNA-Seq data."),
        p("After quantification completes, bulk datasets can be selected in the Datasets tab for differential expression and
              exploratory analyses whereas single-cell datasets will show up in the Single Cell tab for exploratory analyses and sample
              integration."),
        div(class = 'bs-callout bs-callout-info',
            h4('Single cell processing'),
            p('For single cell datasets,',
              tags$a(href = 'https://satijalab.org/seurat/', target='_blank', 'Seurat'),
              'version 3.1 is used for normalization, dimensionality reduction, clustering, and to find marker genes for each cluster.',
              tags$code('Seurat::SCTransform'),
              'is used for normalization. The standard Seurat workflow is used for all other steps.'
            )
        )
      )),
    # Differential expression ----
    docsSubsection(
      id = 'ds-differential-expression',
      name = 'Differential expression',
      content = tagList(
        HTML("<p>To run differential expression analyses between two groups of samples in a dataset, first select a bulk dataset
                  and ensure that the <b>differential expression</b> analysis type toggle is selected. Next, type a name for your new analysis in the <b>Analysis name</b>
                  input and type enter.</p>"),
        img(src = 'docs-ds-diff_expr.png', class = 'bs-docs-img'),
        HTML("This will bring up additional buttons that allow you to label samples in your dataset as either test or control. For example, to
                  label a group of samples as the test group, select the rows that correspond to test samples and click the <b>test</b> button. After
                  indicating control and test samples click <b>Run Analysis</b>."),
        img(src = 'docs-ds-label_rows.png', class = 'bs-docs-img'),
        div(class = 'bs-callout bs-callout-danger',
            h4('Minimum number of samples'),
            p('Differential expression analysis requires at least three samples total between control and test groups.')
        ),
        div(class = 'bs-callout bs-callout-info',
            h4('Differential expression analysis'),
            p('Differential expression analysis is performed with',
              tags$a(href = 'https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29', target='_blank', 'limma-voom.'),
              'Counts are first imported in accordance with the',
              tags$a(href = 'https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#limma-voom', target='_blank', 'tximport vignette.'),
              'Surrogate variables are discovered using',
              tags$a(href = 'https://academic.oup.com/nar/article/42/21/e161/2903156', target='_blank', 'svaseq'),
              'and included in the design matrix.')
        ),
        p("After running differential expression analysis, you can re-select the same analysis to view sample labels, download
               a csv of the limma differential expression results, and inspect multidimensional scaling plots with and without accounting
               for surrogate variables.")
      )
    ),
    # Sample-level gene expression ----
    docsSubsection(
      id = 'ds-gene-plots',
      name = 'Sample-level gene expression',
      content = tagList(p("Sample level gene expression plots allow you to inspect normalized expression values for genes of
                              interest split by groups of interest. These plots allow you to see expression values for all
                              samples in order to avoid being misled by summary statistics."),
                        img(src = 'docs-ds-genes_plot.png', class = 'bs-docs-img'),
                        HTML("<p>To view sample level gene expression plots, select a dataset and make sure the <b>exploratory</b>
                                  analysis type toggle is selected. If you have not already done so, label samples (table rows) for
                                  groups of interest. To do this, select rows corresponding to a group of interest,
                                  type a <b>Group name for selected rows</b> and click the plus button to label the
                                  selected samples. Finally, select genes to <b>Show expression for</b>.</p>"),
                        img(src = 'docs-ds-explore_genes.png', class = 'bs-docs-img'),
                        div(class = 'bs-callout bs-callout-danger',
                            h4('Minimum number of groups'),
                            p("Sample level gene expression plots won't show up unless at least two groups of samples are labeled. Previously
                                  labeled groups will be saved between sessions.")
                        ),
                        div(class = 'bs-callout bs-callout-info',
                            h4('Gene expression normalization'),
                            HTML("<p>Gene expression is normalized with the function <code>DESeq2::vst</code> with <code>blind = FALSE</code>.
                                     For further details see
                                     <a href='http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization' target='_blank'>Data transformations and visualization</a>
                                     in the DESeq2 vignette. The same function is used to normalize gene expression values for MDS plots after <a href='#ds-differential-expression'>differential expression</a> and for
                                     <a href='#ds-cell-plots'>cell-type deconvolution</a>.</p>")
                        )
      )
    ),
    # Cell-type deconvolution ----
    docsSubsection(
      id = 'ds-cell-plots',
      name = 'Cell-type deconvolution',
      content = tagList(p("Cell type deconvolution uses gene expression from single-cell datasets to estimate the proportion
                               of different cell types in bulk RNA-Seq samples. Large differences in cell-type proportions may
                               confound differential expression but also point to important biology."),
                        div(class = 'bs-callout bs-callout-danger',
                            h4('Requires comparable single-cell dataset'),
                            HTML("<p>Cell-type deconvolution requires a single-cell dataset that contains atleast the cell types that you expect
                                     to observe in the bulk dataset that you are deconvoluting.</p>")
                        ),
                        img(src = 'docs-ds-cell_deconvolution.png', class = 'bs-docs-img'),
                        HTML("<p>To run cell-type deconvolution, select a dataset and make sure the <b>exploratory</b> analysis type is
                              selected. If you have not already done so, label samples (table rows) for
                                  groups of interest. To do this, select rows corresponding to a group of interest,
                                  type a <b>Group name for selected rows</b> and click the plus button to label the
                                  selected samples. Also make sure to <b>Toggle cell-type deconvolution</b></p>"),
                        img(src = 'docs-ds-cell_deconvolution_toggle.png', class = 'bs-docs-img'),
                        HTML("Next, select a reference single cell dataset that contains the cell types that you expect to observe in the bulk
                                 dataset that you are deconvoluting. Also select the cell types to deconvolute. You should not select cell types
                                 that you don't expect within your bulk dataset."),
                        img(src = 'docs-ds-cell_deconvolution_ref.png', class = 'bs-docs-img'),
                        div(class = 'bs-callout bs-callout-info',
                            h4('Cell-type deconvolution'),
                            HTML("<p>Cell-type deconvolution is performed using <code>dtangle</code> as described in the
                                     <a href='https://github.com/gjhunt/dtangle/blob/master/vign/sc_vignette.md' target='_blank'>single-cell vignette</a>
                                     except using <code>DESeq2::vst</code> to normalize RNA-Seq counts and the <code>sctransform</code> log normalized counts
                                     for the reference single cell dataset. See the
                                     <a href='https://academic.oup.com/bioinformatics/article/35/12/2093/5165376?guestAccessKey=ac40b15d-bec0-48c1-be94-fbef567f63ec' target='_blank'>dtangle manuscript</a>
                                     for further details.</p>")
                        )

      ))
  ))



# functions ----

#' UI for docs page (add new sections here)
#' @export
#' @keywords internal
#'
#'
#' UI for docs page (add new sections here)
#' @export
#' @keywords internal
docsPageUI <- function(id, tab, active) {
  docsSections <- list(
    # Datasets documentation
    datasetsSection,
    # Single cell documentation
    docsSection(id = 'sc-docs', name = 'Single Cell',
                content = list(
                  docsSubsection(id = 'sc-label-clusters',
                                 name = 'Labeling cell clusters',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))

                  ),
                  docsSubsection(id = 'sc-integrate-samples',
                                 name = 'Integrate datasets',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  ),
                  docsSubsection(id = 'sc-differential-expression',
                                 name = 'Differential expression analysis',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  )
                )),

    # Pathways documentation
    docsSection(id = 'path-docs', name = 'Pathways',
                content = list(
                  docsSubsection(id = 'path-overview',
                                 name = 'Overview',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))

                  ),
                  docsSubsection(id = 'path-query',
                                 name = 'CMAP02/L1000 query genes',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  )
                )),

    # Drugs documentation
    docsSection(id = 'drugs-docs', name = 'Drugs',
                content = list(
                  docsSubsection(id = 'drugs-overview',
                                 name = 'Overview',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))

                  ),
                  docsSubsection(id = 'drugs-custom-queries',
                                 name = 'Custom queries',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  ),
                  docsSubsection(id = 'drugs-advanced',
                                 name = 'Advanced options',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  )
                ))
  )


  section_info <- extract_section_info(docsSections)


  withTags({
    tabPane(tab, active,
            div(class = "bs-docs-container",
                div(class = 'row',
                    div(class = 'col-md-6 col-md-offset-1', role="main", docsSections),
                    div(class = "col-md-4 col-md-offset-1 bd-toc", role="complementary", docsSideNav(section_info))
                )
            )
    )
  })
}



#' A primary section in the Docs page
#' @export
#' @keywords internal
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
#' @export
#' @keywords internal
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
#' @export
#' @keywords internal
docsSideNav <- function(section_info) {

  section_id_names <- section_info$section_id_names
  subsection_id_names <- section_info$subsection_id_names

  nsections <- length(section_id_names)

  withTags({
    nav(class = "bs-docs-sidebar well-form well-bg",
        ul(class="nav bs-docs-sidenav",
           lapply(seq_len(nsections), function(i) {

             # get id and name for current section
             section_id_name <- section_id_names[[i]]

             # get id and names list of subsections for current section
             subsectioni_id_names <- subsection_id_names[[i]]

             li(
               # link for current section
               a(href=paste0("#", section_id_name[1]), section_id_name[2]),
               ul(class="nav",
                  # link for each subsection of current section
                  lapply(subsectioni_id_names, function(subsectioni_id_name) {
                    li(
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
#' @export
#' @keywords internal
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

