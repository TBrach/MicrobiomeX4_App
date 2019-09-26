html_styles <- tags$head(
        # set color of the info texts
        tags$style(HTML("
                        .shiny-text-output {
                        color: red;
                        font-size: 15px;
                        }
                        ")),
        # set colors of checkboxes (not in currently)
        tags$style(HTML("
                        .checkbox {
                        color: black;
                        font-size: 22px;
                        line-height: 35px;
                        }
                        ")),
        # next set's text of the radio buttons
        tags$style(HTML("
                        .radio {
                        color: black;
                        font-size: 15px;
                        }
                        ")),
        # sets text of navbarPage title
        # rgb(17, 119, 85), red #cc3f3f
        tags$style(type = 'text/css',
                   '.navbar-default .navbar-brand {
                   color: white;
                   font-weight: 800;
                   font-size: 20px;
                   }',
                   # sets background-color of navbarPage bar and text size of the tabs, not sure how to change the font color
                   '.navbar {
                   background-color: #23373B;
                   font-size: 15px;
                   font-weight: 700;
                   }')
        
        
)
# #814681