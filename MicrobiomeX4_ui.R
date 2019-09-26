# - source all the tabs -
tabpath <- "./Tabs_ui"
source(file.path(tabpath, "html_styles_ui.R"))
source(file.path(tabpath, "Tab01_ui.R"))
source(file.path(tabpath, "Tab02_ui.R"))
source(file.path(tabpath, "Tab03_ui.R"))
source(file.path(tabpath, "Tab04_ui.R"))
source(file.path(tabpath, "Tab05_ui.R"))
source(file.path(tabpath, "Tab06_ui.R"))
source(file.path(tabpath, "Tab07_ui.R"))
source(file.path(tabpath, "Tab08_ui.R"))
source(file.path(tabpath, "Tab09_ui.R"))
source(file.path(tabpath, "Tab10_ui.R"))
source(file.path(tabpath, "Tab11_ui.R"))
source(file.path(tabpath, "Tab12N_ui.R"))
source(file.path(tabpath, "Tab13_ui.R"))
source(file.path(tabpath, "Tab14_ui.R"))
source(file.path(tabpath, "Tab15_ui.R"))
# --



ui <- fluidPage(
        html_styles,
        navbarPage(title = "MicrobiomeX4", Tab01, Tab02, Tab03, Tab15, Tab11,
                   Tab04, Tab12, Tab05, Tab06, Tab07, Tab08, Tab09, Tab10,
                   Tab13, Tab14)
)