Tab13 <- tabPanel(title = "Taxa finder",
                  sidebarLayout(
                          sidebarPanel(
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Find specific taxa and compare abundances of those taxa.")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_13')
                                  ),
                                  
                                  wellPanel(
                                          
                                          tags$h4("Select the object to work with on this panel."),
                                          
                                          radioButtons(inputId = "tcaType13", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered13", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%")
                                          
                                  ),
                                  
                                  
                                  wellPanel(
                                          tags$h4("Find taxa"),
                                          
                                          radioButtons(inputId = "taxLevelSearch", label = "Select taxonomic level to search in",
                                                       choices = c("Species",
                                                                   "Genus",
                                                                   "Family",
                                                                   "Order",
                                                                   "Class",
                                                                   "Phylum"),
                                                       width = "100%"),
                                          
                                          textInput(inputId = "searchWord", label = "Search word", value = "Bacteroides"),
                                          
                                          actionButton(inputId = "taxaSearch", label = "Perform taxa search")
                                          
                                          # tags$h5("NB: Please be patient, tax_glom can take several minutes")
                                  ),
                                  
                                  
                                  wellPanel(
                                          tags$h4("Plot taxa abundances"),
                                          
                                          textInput(inputId = "userTaxa", label = "comma-separated Indexes of Taxa to print", value = ""),
                                          
                                          actionButton(inputId = "plotTaxa", label = "Plot abundances of given taxa"),
                                          
                                          textInput(inputId = 'ref_13', label = "Ref (if blank: 'all comparisons' will be used)"),
                                          
                                          checkboxInput(inputId = "log", label = "log abundance", value = FALSE),
                                          
                                          checkboxInput(inputId = "pool", label = "pool abundances", value = FALSE),
                                          
                                          checkboxInput(inputId = "RA", label = "turn to rel. ab.", value = FALSE)
                                          
                                          # tags$h5("NB: Please be patient, tax_glom can take several minutes")
                                  )
                                  
                                  
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  tableOutput(outputId = "overviewView13"),
                                  
                                  tableOutput(outputId = "foundTaxa"),
                                  
                                  tableOutput(outputId = "pValsTaxa"),
                                  
                                  plotOutput(outputId = "plotTaxaAbundances")
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)