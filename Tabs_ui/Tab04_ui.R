Tab04 <- tabPanel(title = "Heatmap",
                  sidebarLayout(
                          sidebarPanel(
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Visualize phyloseq object in a heatmap, and visualize sparsity because zero counts are highlighted by a user-defined color")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_4')
                                  ),
                                  
                                  
                                  wellPanel(
                                          radioButtons(inputId = "tcaType4", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered4", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%"),
                                          
                                          
                                          textInput(inputId = "max_abundance_for_color", label = "Abundance of max. color", value = "", width = "200px"),
                                          
                                          tags$h5("If Abundance of max. color is blank, the max abundance value of the data is used."),
                                          
                                          textInput(inputId = "zero_color", label = "Color for zero counts", value = "gray", width = "200px"),
                                          
                                          textInput(inputId = "taxa_index_range", label = "taxa index range, e.g. 5:34, leave blank to see all", value = "", width = "200px"),
                                          
                                          checkboxInput(inputId = "restrictto2compare_4", label = "only show the compare levels in heatmap", value = FALSE),
                                          
                                          checkboxInput(inputId = "ra4", label = "turn to rel. ab.", value = FALSE),
                                          
                                          checkboxInput(inputId = "log4", label = "log abundance", value = FALSE),
                                          
                                          actionButton(inputId = "calculateHM", label = "Plot heatmap"),
                                          
                                          tags$h5("NB: Please be patient, calculation of plot can take 20 seconds.")
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  # tableOutput(outputId = "explorationViews_3"),
                                  tableOutput(outputId = "overviewView4"),
                                  
                                  plotOutput(outputId = "plotHM")
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)