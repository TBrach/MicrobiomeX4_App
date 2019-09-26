Tab07 <- tabPanel(title = "Fisher test",
                  sidebarLayout(
                          sidebarPanel(
                                  
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Find taxa that are more prevalent in one group than in the other. NB: you can only compare two
                                                  groups, so if you have more than two levels in group_var_levels you have to define two levels under 'compare' in the Parameters tab!")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_7')
                                  ),
                                  
                                  
                                  wellPanel(
                                          
                                          tags$h5("NB: the object type only matters for the heatmap here not for the fisher test results."), 
                                          
                                          radioButtons(inputId = "tcaType7", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered7", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%"),
                                          
                                          actionButton(inputId = "fisher", label = "Calculate prevalence differences using fisher.test"),
                                          tags$h5(),
                                          textInput(inputId = "max_abundance_for_colorF", label = "Abundance of max. color", value = "", width = "200px"),
                                          
                                          tags$h5("If Abundance of max. color is blank, the max abundance value of the data is used."),
                                          
                                          textInput(inputId = "zero_color7", label = "Color for zero counts", value = "gray", width = "200px"),
                                          
                                          textInput(inputId = "maxShown7", label = "max number of hits shown in heatmap", value = "40", width = "200px"),
                                          
                                          checkboxInput(inputId = "restrictto2levels", label = "only show the two compare levels in heatmap", value = FALSE),
                                          
                                          checkboxInput(inputId = "ra7", label = "turn to rel. ab.", value = FALSE),
                                          
                                          checkboxInput(inputId = "log7", label = "log abundance in heatmap", value = FALSE)
                                          
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  tableOutput(outputId = "overviewView7"),
                                  
                                  plotOutput(outputId = "plotFisher", height = "750px"), # height = "1000px"
                                  
                                  tableOutput(outputId = "tableFisher")
                                  
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)