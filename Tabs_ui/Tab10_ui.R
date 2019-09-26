Tab10 <- tabPanel(title = "Phylum analysis",
                  sidebarLayout(
                          sidebarPanel(
                                  
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Compare phylum to phylum ratios. These are independent of compositionality and especially the Firmicutes to Bacteroides ratio has been reported
                                                  extensively in the gut microbiome field in relation to obesity."),
                                          
                                          
                                          tags$h5("NB: in both cases ratios where either the nominator or denominator phylum is absent are set to NA and ignored!
                                                  Wilcoxon test is used for determining whether the ratios are significantly different between the sample groups.")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_10')
                                  ),
                                  
                                  
                                  wellPanel(
                                          
                                          tags$h5("Since phylum/phylum analyses are independent of size factor corrections, always the raw counts will be used. But you can decide on whether
                                                  all taxa or only the filtered taxa are included."),
                                          
                                          radioButtons(inputId = "tcaType10", label = "Select phyloseq object",
                                                       choices = c("raw counts"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered10", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%"),
                                          
                                          wellPanel(
                                                  
                                                  tags$h5("Step 1: Compare phylum/phylum ratios with boxplots"),
                                                  
                                                  textInput(inputId = "numerator_phylum", label = "Numerator Phylum", value = "Firmicutes", width = "200px"),
                                                  
                                                  textInput(inputId = "denominator_phylum", label = "Denominator Phylum", value = "", width = "200px"),
                                                  
                                                  tags$h5("Leave empty when you want to compare to all other phyla. Otherwise give a comma-separated list without typos."),
                                                  
                                                  textInput(inputId = 'ref_10', label = "Ref (if blank: 'all comparisons' will be used)"),
                                                  
                                                  actionButton(inputId = "firmicutes", label = "Calculate P/P ratios")
                                          ),
                                          
                                          
                                          wellPanel(
                                                  tags$h5("Step 2: compare all phylum/phylum ratios"),
                                                  
                                                  actionButton(inputId = "tile", label = "Calculate phylum to phylum tile plot")
                                          ),
                                          
                                          tags$h5("NB: Please be patient, plots can take 10 - 20 seconds.")
                                          
                                  )
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  tableOutput(outputId = "overviewView10"),
                                  
                                  plotOutput(outputId = "plotPhylum", height = "750px"),
                                  
                                  tableOutput(outputId = "tableFirmicutes")
                          )
                          
                  )
)