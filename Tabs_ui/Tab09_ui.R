Tab09 <- tabPanel(title = "Wilcoxon test",
                  sidebarLayout(
                          sidebarPanel(
                                  
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Find taxa that are differentially abundant between the groups using wilcoxon test. NB: you can only compare two
                                                  groups, so if you have more than two levels in group_var_levels you have to define two levels under 'compare' in the Parameters tab!")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_9')
                                  ),
                                  
                                  
                                  wellPanel(
                                          
                                          radioButtons(inputId = "tcaType9", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered9", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "wilcoxonZeros", label = "Include or exlcude zero counts.",
                                                       choices = c("include zero counts",
                                                                   "exclude zero counts"),
                                                       width = "100%"),
                                          
                                          
                                          actionButton(inputId = "wilcoxon", label = "Calculate abundance differences using wilcoxon"),
                                          tags$h5(),
                                          textInput(inputId = "max_abundance_for_colorW", label = "Abundance of max. color", value = "", width = "200px"),
                                          
                                          tags$h5("If Abundance of max. color is blank, the max abundance value of the data is used."),
                                          
                                          textInput(inputId = "zero_color9", label = "Color for zero counts", value = "gray", width = "200px"),
                                          
                                          textInput(inputId = "maxShown9", label = "max number of hits shown in heatmap", value = "40", width = "200px"),
                                          
                                          checkboxInput(inputId = "restrictto2levels_9", label = "only show the two compare levels in heatmap", value = FALSE),
                                          
                                          checkboxInput(inputId = "ra9", label = "turn to rel. ab.", value = FALSE),
                                          
                                          checkboxInput(inputId = "log9", label = "log abundance in heatmap", value = FALSE)
                                          
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  tableOutput(outputId = "overviewView9"),
                                  
                                  plotOutput(outputId = "plotWilcoxon", height = "750px"), # height = "1000px"
                                  
                                  tableOutput(outputId = "tableWilcoxon")
                                  
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)