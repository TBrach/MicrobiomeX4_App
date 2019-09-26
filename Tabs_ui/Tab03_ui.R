Tab03 <- tabPanel(title = "Library size adjustment",
                 sidebarLayout(
                         sidebarPanel(
                                 wellPanel(
                                         tags$h4("Purpose of this tab:"),
                                         
                                         tags$h5("Adjust differences in library size of the samples, so samples can be compared to each other. Three variations of calculating size factors are implemented.
                                                 NB: They all start from the raw count object. You could tax_glom first if you think this is better."),
                                         tags$ul(
                                                 tags$li(tags$b("DESeq_gm_exclZero:"), "DESeq2 ratio method using a prevalence filtered object and a geometric mean that excludes zero counts."),
                                                 tags$li(tags$b("DESeq_poscounts:"), "DESeq2 ratio method with a geometric mean that includes zero counts using the same prevalence filtered object."),
                                                 tags$li(tags$b("library size (RA):"), "Putting all samples to same count level (basically relative abundance)")
                                         )
                                         
                                 ),
                                 
                                 wellPanel(
                                         tags$h4("Info box"),
                                         
                                         textOutput(outputId = 'infoText_3')
                                 ),
                                 
                                 
                                 wellPanel(
                                         tags$h4("Calculate size factors"),
                                         
                                         textInput(inputId = "prev_SFs", label = "Prevalence for DESeq_gm_exclZero and DESeq_poscounts", value = "60", width = "250px"),
                                         
                                         actionButton(inputId = "calcSFs", label = "Calc. size factors")
                                 ),
                                 
                                 wellPanel(
                                         tags$h4("Visualize size factor adjustment"),
                                         
                                         radioButtons(inputId = "tcaType3", label = "Select type of library size adjustment",
                                                      choices = c("DESeq_gm_exclZero",
                                                                  "DESeq_poscounts",
                                                                  "library size (RA)"),
                                                      width = "100%"),
                                         
                                         radioButtons(inputId = "filtered3", label = "only filtered taxa?",
                                                      choices = c("all taxa",
                                                                  "filtered taxa"),
                                                      width = "100%"),
                                         
                                         actionButton(inputId = "barplotSFs", label = "Visualise library size adjustment"),
                                         
                                         tags$h5("NB: Please be patient, calculation of plot can take 20 seconds.")
                                 )
                                 
                         ),
                         mainPanel(
                                 # tags$h2("The current phyloseq object"),
                                 tableOutput(outputId = "overviewView3"),
                                 
                                 #tableOutput(outputId = "explorationViews_2"),
                                 
                                 
                                 plotOutput(outputId = "plotSFs")
                                 
                         )
                         # textOutput(outputId = 'explorationViews'))
                 )
)