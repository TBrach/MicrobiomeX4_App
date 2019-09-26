Tab12 <- tabPanel(title = "Alpha diversity",
                  sidebarLayout(
                          sidebarPanel(
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Calculate and visualize alpha-diversity.")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_12')
                                  ),
                                  
                                  
                                  wellPanel(
                                          radioButtons(inputId = "tcaType12", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered12", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%"),
                                          
                                          selectInput(inputId = "alpha_measure", label = "Choose the distance measure to use:", choices = c("Observed", "Shannon", "Chao1", "ACE", "Simpson", "InvSimpson", "Fisher"), selected = "Observed"),
                                          
                                          checkboxInput(inputId = "rarify", label = "rarify", value = FALSE),
                                          
                                          textInput(inputId = 'ref_12', label = "Ref (if blank comparisons will be used)"),
                                          
                                          actionButton(inputId = "calculatealphDiv", label = "calculate alpha diversity"),
                                          
                                          tags$h5("NB: Please be patient, calculation can take a bit of time.")
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  tableOutput(outputId = "overviewView12"),
                                  
                                  tableOutput(outputId = "tableAlpha"),
                                  
                                  
                                  plotOutput(outputId = "alphaDiv")
                                  
                          )
                          
                  )
)