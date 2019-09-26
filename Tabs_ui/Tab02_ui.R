Tab02 <- tabPanel(title = "Load and explore phyloseq object",
                 sidebarLayout(
                         sidebarPanel(
                                 
                                 wellPanel(
                                         tags$h4("Purpose of this tab:"),
                                         
                                         tags$h5("Get an overview of the data in the phyloseq object:"),
                                         tags$ul(
                                                 tags$li("Step 1: Load phyloseq object"),
                                                 tags$li("Step 2: Get an idea of the phyla and their abundance distribution (in the raw data!). This step sets phylum_colors for later plots."),
                                                 tags$li("Step 3: Visualize sample sizes (library sizes/total counts).")
                                         )
                                         
                                 ),
                                 
                                 wellPanel(
                                         tags$h4("Info box"),
                                         
                                         textOutput(outputId = 'infoText_2')
                                 ),
                                 
                                 wellPanel(
                                         tags$h4("Step1: Load phyloseq object"),
                                         
                                         fileInput(inputId = "loadPhyseq", label = ""), # accept = "text/rds"
                                         
                                         tags$h5(),
                                         wellPanel(
                                                 tags$h4("Look at the components of the phyloseq object"),
                                                 radioButtons(inputId = "component", label = "Select component of phyloseq object",
                                                              choices = c("otu_table",
                                                                          "sample_data",
                                                                          "tax_table"),
                                                              width = "100%"),
                                                 textInput(inputId = "sampleRange", label = "Sample range", value = "1:100"),
                                                 textInput(inputId = "taxaRange", label = "Taxa range", value = "1:100"),
                                                 actionButton(inputId = "checkPhyseq", label = "See component of phyloseq object")
                                         )
                                         
                                 ),
                                 
                                 wellPanel(
                                         tags$h4("Select data to use for steps 2 and 3"),
                                         
                                         radioButtons(inputId = "tcaType2", label = "Select phyloseq object",
                                                      choices = c("raw counts",
                                                                  "DESeq_gm_exclZero",
                                                                  "DESeq_poscounts",
                                                                  "library size (RA)"),
                                                      width = "100%"),
                                         
                                         radioButtons(inputId = "filtered2", label = "only filtered taxa?",
                                                      choices = c("all taxa",
                                                                  "filtered taxa"),
                                                      width = "100%")
                                         
                                 ),
                                 
                                 
                                 wellPanel(
                                         tags$h4("Step 2: Phyla distribution"),
                                         
                                         
                                         actionButton(inputId = "calcPhyla", label = "Calc. phyla distribution and phylum_colors"),
                                         
                                         tags$h5("NB: Step is necessary to calculate phylum_colors for later plots! Can be repeated after size factor adjustment.")
                                 ),
                                 
                                 wellPanel(
                                         tags$h4("Step 3: Overview of sample sizes"),
                                         
                                         actionButton(inputId = "barplotPS", label = "Generate abundance barplot of the samples"),
                                         
                                         checkboxInput(inputId = "RA_2", label = "turn to rel. ab.", value = FALSE),
                                         
                                         tags$h5("NB: Please be patient, plot calculation takes 20 to 30 seconds.")
                                 )
                                 
                         ),
                         mainPanel(
                                 # tags$h2("The current phyloseq object"),
                                 
                                 tableOutput(outputId = "overviewView2"),
                                 
                                 # tableOutput(outputId = "explorationViews"),
                                 
                                 tableOutput(outputId = "phylaViews"),
                                 
                                 plotOutput(outputId = "samplesPS")
                                 
                         )
                         # textOutput(outputId = 'explorationViews'))
                 )
)