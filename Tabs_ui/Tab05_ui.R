Tab05 <- tabPanel(title = "Filtering",
                  sidebarLayout(
                          sidebarPanel(
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Option to filter the taxa in your phyloseq object for subsequent analyses. Filtering is based on given prevalence and taxa_sums quantile in PC. Taxa remain if their prevalence is higher than the given one
                                                one, or if their taxa_sums is higher than the given quantile. You can choose the object to do the filtering on, only taxa_sums quantiles are affected by previous library size adjustment, 
                                                  prevalence is not. The taxa that remain the filtering are stored and can be used in subsequent analysis.")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_5')
                                  ),
                                  
                                  
                                  wellPanel(
                                          
                                          radioButtons(inputId = "tcaType5", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          textInput(inputId = "filt_prevalence", label = "prevalence", value = "20"),
                                          
                                          
                                          textInput(inputId = "taxa_sums_quantile", label = "taxa_sums quantile in PC", value = "100"),
                                          
                                          tags$h5("taxa whose taxa_sums are above this threshold will be kept even if they do not pass prevalence filter"),
                                          
                                          wellPanel(
                                                  actionButton(inputId = "prevalenceDistribution", label = "Visualize prevalence distribution"),
                                                  
                                                  actionButton(inputId = "filter", label = "Filter")
                                          ),
                                          
                                          tags$h5("NB: Please be patient, may take ~ 10 seconds.")
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  # tableOutput(outputId = "explorationViews_4"),
                                  tableOutput(outputId = "overviewView5"),
                                  
                                  
                                  plotOutput(outputId = "plotFilter") # height = "1000px"
                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)