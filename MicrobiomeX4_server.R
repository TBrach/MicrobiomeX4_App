server_tab_path <- "./Tabs_server"


server <- function(input, output, session){
        
        # - generate reactive Values for entire app -
        rv <- reactiveValues(viewItem = NULL, overviewView = NULL, group_var = NULL, group_var_levels = NULL,
                             color_levels = NULL, compare = NULL, ps = NULL, infoText = NULL, infoText_2 = NULL, Tr_colors = NULL, phyla = NULL, phylum_colors = NULL, Tr = NULL, Tr_SFs = NULL,
                             ps_tca = NULL, ps_tca_poscounts = NULL, SFs = NULL, SFs_poscounts = NULL, SFs_RA = NULL, ps_RA = NULL, Tr_HM = NULL, 
                             filteredKeepTaxa = NULL, tableAssignment = NULL, Tr_Assignment = NULL, Tr_ab_prev = NULL, TrL_filter = NULL, betaDivDist = NULL, Tr_PCoA = NULL, adonis = NULL,
                             Tr_fisher = NULL, fisherTable = NULL, Tr_DESeq2 = NULL, DESeq2Table = NULL, Tr_wilcoxon = NULL, wilcoxonTable = NULL,
                             firmicutesTable = NULL, Tr_phylum = NULL, Tr_alpha = NULL, alphaTable = NULL, taxaFindTable = NULL, infoText_13 = NULL, Tr_taxaAbundance = NULL,
                             tablepValsTaxa = NULL, Tr_taxaRatios = NULL, tablepValsTaxaRatios = NULL,
                             symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
        # --

        # - Tab 01 -
        source(file.path(server_tab_path, "Tab01_server.R"), local = TRUE)
        # --
        
        # - Tab 02 -
        source(file.path(server_tab_path, "Tab02_server.R"), local = TRUE)
        # --
        
        # - Tab 03 -
        source(file.path(server_tab_path, "Tab03_server.R"), local = TRUE)
        # --
        
        # - Tab 04 -
        source(file.path(server_tab_path, "Tab04_server.R"), local = TRUE)
        # --
        
        # - Tab 05 -
        source(file.path(server_tab_path, "Tab05_server.R"), local = TRUE)
        # --
                        
        # - Tab 06 -
        source(file.path(server_tab_path, "Tab06_server.R"), local = TRUE)
        # --
        
        # - Tab 07 -
        source(file.path(server_tab_path, "Tab07_server.R"), local = TRUE)
        # --
        
        # - Tab 08 -
        source(file.path(server_tab_path, "Tab08_server.R"), local = TRUE)
        # --
        
        # - Tab 09 -
        source(file.path(server_tab_path, "Tab09_server.R"), local = TRUE)
        # --
        
        # - Tab 10 -
        source(file.path(server_tab_path, "Tab10_server.R"), local = TRUE)
        # --
        
        # - Tab 11 -
        source(file.path(server_tab_path, "Tab11_server.R"), local = TRUE)
        # --
        
        # - Tab 12 -
        source(file.path(server_tab_path, "Tab12N_server.R"), local = TRUE)
        # --
        
        # - Tab 13 -
        source(file.path(server_tab_path, "Tab13_server.R"), local = TRUE)
        # --
        
        # - Tab 14 -
        source(file.path(server_tab_path, "Tab14_server.R"), local = TRUE)
        # --
        
        # - Tab 15 -
        source(file.path(server_tab_path, "Tab15_server.R"), local = TRUE)
        # --
        
}