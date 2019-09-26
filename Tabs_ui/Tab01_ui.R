Tab01 <- tabPanel(title = "Parameters",
                 sidebarLayout(
                         sidebarPanel(
                                 
                                 wellPanel(
                                         tags$h4("Purpose of this tab:"),
                                         
                                         tags$h5("Set general input parameters")
                                         
                                 ),
                                 
                                 wellPanel(
                                         tags$h4("Info box"),
                                         
                                         textOutput(outputId = 'infoText_1')
                                 ),
                                 
                                 wellPanel(
                                         
                                         tags$h4("Group variables and colors"),
                                         textInput(inputId = "group_var", label = "group_var", value = ""), # , width = "150px"
                                         textInput(inputId = "group_var_levels", label = "group_var_levels (comma-separated list)", value = ""),
                                         textInput(inputId = "color_levels", label = "colors (comma-separated list)", value = ""),
                                         textInput(inputId = "compare", label = "compare (comma-separated list)", value = "")
                                         # splitLayout(
                                         #         textInput(inputId = "grp1", label = "group 1", value = "DK", width = "120px"),
                                         #         textInput(inputId = "grp2", label = "group 2", value = "SL", width = "120px")),
                                         
                                 ),
                                 
                                 wellPanel(
                                         
                                         tags$h4("Plot your colors"),
                                         
                                         actionButton(inputId = "plotColors", label = "Just to see the colors you have chosen")
                                         
                                         
                                 )
                                 
                         ),
                         mainPanel(
                                 
                                 tableOutput(outputId = "overviewView1"),
                                 
                                 plotOutput(outputId = "plotColors")
                                 
                         )
                         # textOutput(outputId = 'explorationViews'))
                 )
)