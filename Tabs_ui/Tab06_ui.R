Tab06 <- tabPanel(title = "Beta diversity",
                  sidebarLayout(
                          sidebarPanel(
                                  
                                  wellPanel(
                                          tags$h4("Purpose of this tab:"),
                                          
                                          tags$h5("Test for overall microbiome differences between the sample groups. 
                                                  A distance (option to choose different distance measures!) is calculated between each two samples in your data. 
                                                  Are samples within the groups closer to each other than between the groups? Do the groups form clusters in an ordination plot?
                                                  (NB: calculations are done on relative abundances of filtered phyloseq object.)")
                                          
                                  ),
                                  
                                  wellPanel(
                                          tags$h4("Info box"),
                                          
                                          textOutput(outputId = 'infoText_6')
                                  ),
                                  
                                  
                                  wellPanel(
                                          
                                          radioButtons(inputId = "tcaType6", label = "Select phyloseq object",
                                                       choices = c("raw counts",
                                                                   "DESeq_gm_exclZero",
                                                                   "DESeq_poscounts",
                                                                   "library size (RA)"),
                                                       width = "100%"),
                                          
                                          radioButtons(inputId = "filtered6", label = "only filtered taxa?",
                                                       choices = c("all taxa",
                                                                   "filtered taxa"),
                                                       width = "100%")
                                  ),
                                  
                                  wellPanel(
                                          tags$h5("Calculate beta diversity distances"),
                                          
                                          selectInput(inputId = "beta_measure", label = "Choose the distance measure to use:", choices = c("jsd", "bray", "jaccard", "euclidean", "manhattan", "canberra"), selected = "jsd"),
                                          
                                          actionButton(inputId = "calcBetaDiversity", label = "Calculate beta diversity distances"),
                                          
                                          tags$h5("Patience, can take a while if you have many samples. Only includes samples comprised by group_var_levels.")
                                  ),
                                  
                                  wellPanel(
                                          tags$h5("Plot beta diversity distances and calculate significance (adonis)"),
                                          
                                          checkboxInput(inputId = "ellipse", label = "add ellipse shadows", value = FALSE),
                                          
                                          numericInput(inputId = "ellipse_level", label = "ellipse level (0-1)", value = 0.95),
                                          
                                          radioButtons(inputId = "compareUse", label = "Decide on how you want to restrict to the compare levels",
                                                       choices = c("remove other levels",
                                                                   "keep other levels as gray dots"),
                                                       width = "100%"),
                                          
                                          
                                          radioButtons(inputId = "pcoaCorrection", label = "Adjust axes ratio in pcoa?",
                                                       choices = c("yes",
                                                                   "no"),
                                                       width = "100%"),
                                          
                                          
                                          actionButton(inputId = "plotPCoA", label = "Visualize beta diversity")
                                  )
                                  
                                  
                                  
                          ),
                          mainPanel(
                                  # tags$h2("The current phyloseq object"),
                                  
                                  tableOutput(outputId = "overviewView6"),
                                  
                                  tableOutput(outputId = "adonis"),
                                  
                                  plotOutput(outputId = "plotPCoA") # height = "1000px"
                                  
        
                                  

                                  
                          )
                          # textOutput(outputId = 'explorationViews'))
                  )
)