shiny::htmlTemplate("templates/contrasts.html", 
                     # styling related
                    show_settings = shinyWidgets::dropdownButton(
                      shiny::sliderInput('point_size_contrasts', 'Point size:',
                                         width = '100%', ticks = FALSE,
                                         min = 0.5, max = 4, value = 2.5, step = 0.5),
                      shiny::sliderInput('point_jitter', 'Point jitter:',
                                         width = '100%', ticks = FALSE,
                                         min = 0, max = 3, value = 0, step = 0.5),
                      circle = FALSE, right = TRUE, icon = shiny::icon('cog', 'fa-fw')),
                      cluster_plot = shiny::plotOutput('cluster_plot_contrasts')
)
