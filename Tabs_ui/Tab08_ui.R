Tab08 <- tabPanel(title = "DESeq2 test",
                  sidebarLayout(
                          sidebarPanel(
                                  
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Find taxa that are differentially abundant between the groups using DESeq2. NB: you can only compare two
                                                  groups, so if you have more than two levels in group_var_levels you have to define two levels under 'compare' in the Parameters tab!")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_8')
                                  ),
                                  
                                  
                                  wellPanel(
                                          
                                          radioButtons(inputId = "tcaType8", label = "Select size factor correction",
                                                       choices = c("DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered8", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%"),
                                          
                                          
                                          actionButton(inputId = "DESeq2", label = "Calculate abundance differences using DESeq2"),
                                          
                                          tags$h5(),
                                          
                                          textInput(inputId = "max_abundance_for_colorD", label = "Abundance of max. color", value = "", width = "200px"),
                                          
                                          tags$h5("If Abundance of max. color is blank, the max abundance value of the data is used."),
                                          
                                          textInput(inputId = "zero_color8", label = "Color for zero counts", value = "gray", width = "200px"),
                                          
                                          textInput(inputId = "maxShown8", label = "max number of hits shown in heatmap", value = "40", width = "200px"),
                                          
                                          checkboxInput(inputId = "restrictto2levels_8", label = "only show the two compare levels in heatmap", value = FALSE),
                                          
                                          checkboxInput(inputId = "ra8", label = "turn to rel. ab.", value = FALSE),
                                          
                                          checkboxInput(inputId = "log8", label = "log abundance in heatmap", value = FALSE)
                                          
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  tableOutput(outputId = "overviewView8"),
                                  
                                  plotOutput(outputId = "plotDESeq2", height = "750px"), # height = "1000px"
                                  
                                  tableOutput(outputId = "tableDESeq2")
                                  
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)