## ui.R ##
shiny::htmlTemplate("template.html",
                    # styling related
                    show_settings = shiny::actionButton('show_settings', '',
                                                          icon = shiny::icon('cog', 'fa-fw'),
                                                          title = 'Toggle settings'),
                    point_size = shiny::sliderInput('point_size', 'Point size:',
                                                    width = '100%', ticks = FALSE,
                                                    min = 0.5, max = 4, value = 2.5, step = 0.5),
                    point_jitter = shiny::sliderInput('point_jitter', 'Point jitter:',
                                                    width = '100%', ticks = FALSE,
                                                    min = 0, max = 3, value = 0, step = 0.5),

                    # toggle to show dataset integration panel
                    show_integration   = shiny::actionButton('show_integration', '',
                                                             icon = shiny::icon('object-group', 'far fa-fw'),
                                                             title = 'Toggle dataset integration',
                                                             class = 'squashed-btn'),
                    # dataset integration panel
                    test_integration = shiny::selectizeInput('test_integration', 'Test datasets:', multiple = TRUE, choices = '', width = '100%'),
                    ctrl_integration = shiny::selectizeInput('ctrl_integration', 'Control datasets:', multiple = TRUE, choices = '', width = '100%'),
                    integration_name = shiny::textInput('integration_name', 'Name for new integrated analysis:', width = '100%'),
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

