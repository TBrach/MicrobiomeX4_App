Tab14 <- tabPanel(title = "Taxa ratio",
                  sidebarLayout(
                          sidebarPanel(
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Compare the ratios of two user defined taxa.")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_14')
                                  ),
                                  
                                  
                                  wellPanel(
                                          tags$h4("plot taxa Abundance ratios"),
                                          
                                          radioButtons(inputId = "tcaType14", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered14", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%"),
                                          
                                          
                                          textInput(inputId = "numerator", label = "Index of numerator", value = ""),
                                          
                                          textInput(inputId = "denominator", label = "Index of denominator", value = ""),
                                          
                                          checkboxInput(inputId = "logRatio", label = "log Ratios", value = FALSE),
                                          
                                          actionButton(inputId = "plotRatioTaxa", label = "plot abundance ratios of the taxa")
                                          
                                          # tags$h5("NB: Please be patient, tax_glom can take several minutes")
                                  )
                                  
                                  
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  tableOutput(outputId = "overviewView14"),
                                  
                                  tableOutput(outputId = "pValsRatioTaxa"),
                                  
                                  plotOutput(outputId = "plotTaxaRatios")
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)