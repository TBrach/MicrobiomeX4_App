Tab15 <- tabPanel(title = "Tax Assignment",
                  sidebarLayout(
                          sidebarPanel(
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Check taxonomic assignment niveau and whether non-zero counts associate with prevalence.")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_15')
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("The objects to look at here:"),
                                          
                                          radioButtons(inputId = "tcaType15", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered15", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%")
                                  ),
                                  
                                  
                                  wellPanel(
                                          tags$h4("Taxonomic Assignment"),
                                          
                                          actionButton(inputId = "calcAssignment", label = "Calc. taxonomic assignment")
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Visualize association of abundance to prevalence"),
                                          
                                          actionButton(inputId = "plot_ab_prev", label = "Visualise prevalence to non-zero counts."),
                                          
                                          tags$h5("NB: Please be patient, calculation of plot can take a bit.")
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  tableOutput(outputId = "overviewView15"),
                                  
                                  tableOutput(outputId = "table_Assignment"),
                                  
                                  plotOutput(outputId = "plot_Assignment", height = "700px"), #, height = "700px"
                                  
                                  plotOutput(outputId = "plotabPrev", height = "700px")
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)