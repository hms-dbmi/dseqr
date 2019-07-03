## ui.R ##
shiny::htmlTemplate("template.html",
                    # buttons for dataset integration
                    show_integration   = shiny::actionButton('show_integration', '',
                                                              icon = shiny::icon('object-group', 'far fa-fw'),
                                                              title = 'Toggle dataset integration'),
                    submit_integration = shiny::actionButton('submit_integration', 'Integrate Datasets',
                                                              icon = shiny::icon('object-group', 'fa-fw'),
                                                              title = 'Integrate datasets',  
                                                              class = 'btn-block btn-default'),
                                                            
                    # button for cluster input
                    show_rename    = shiny::actionButton('show_rename', '',
                                                         icon = shiny::icon('tag', 'fa-fw'),
                                                         title = 'Toggle rename cluster'),
                    show_contrasts = shiny::actionButton('show_contrasts', '',
                                                         icon = shiny::icon('chevron-right', 'fa-fw'),
                                                         title = 'Toggle single group comparisons',
                                                         class = 'squashed-btn'),
                    rename_cluster = shiny::actionButton('rename_cluster', '',
                                                         icon = shiny::icon('plus', 'fa-fw'),
                                                         title = 'Rename cluster'),

                    # buttons for GeneCards
                    genecards = shiny::actionButton('genecards', '',
                                                    icon = shiny::icon('external-link-alt', 'fa-fw'),
                                                    title = 'Go to GeneCards'),

                    # plots
                    cluster_plot = shiny::plotOutput('cluster_plot'),
                    marker_plot  = shiny::plotOutput('marker_plot'),
                    biogps_plot  = shiny::plotOutput('biogps')
)

