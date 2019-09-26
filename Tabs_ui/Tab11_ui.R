Tab11 <- tabPanel(title = "Tax_glom",
                  sidebarLayout(
                          sidebarPanel(
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Option to tax_glom to the chosen taxonomic level")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_11')
                                  ),
                                  
                                  
                                  wellPanel(
                                          tags$h4("Do tax_glom"),
                                          
                                          tags$h5("NB: this action will also refresh all library size adjusted objects if SFs exist!! So re-do the size-factor calculation if you want to have SFs based on the tax_glomed object.
                                                  It also resets the filtered taxa."),
                                          
                                          radioButtons(inputId = "taxLevel", label = "Select taxonomic level to tax_glom to",
                                                       choices = c("Species",
                                                                   "Genus",
                                                                   "Family",
                                                                   "Order",
                                                                   "Class",
                                                                   "Phylum"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "NArm", label = "Do you want to remove taxa that are NA at chosen taxonomic level?",
                                                       choices = c("FALSE",
                                                                   "TRUE"),
                                                       width = "100%"),
                                          
                                          actionButton(inputId = "taxGlom", label = "Perform tax_glom"),
                                          
                                          tags$h5("NB: Please be patient, tax_glom can take several minutes")
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  tableOutput(outputId = "overviewView11")
                                  
                                  # tableOutput(outputId = "explorationViews_11")
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)