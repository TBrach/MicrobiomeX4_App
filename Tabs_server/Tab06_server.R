
# - Tab6: beta diversity -
# -- output infoText --
output$infoText_6 <- renderText({
        rv$infoText_6
})
# ----

# -- render the overviewViews table --
output$overviewView6 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- calc BetaDiversity --
observeEvent(input$calcBetaDiversity, {
        
        rv$infoText_6 <- NULL
        
        rv$Tr_PCoA <- NULL
        rv$adonis <- NULL
        rv$betaDivDist <- NULL
        
        
        if (input$tcaType6 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType6 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType6 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType6 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_6 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_6 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered6 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_6 <- "Sorry. You have to do the filtering first."
                        return()
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
        }
        
        measure <- input$beta_measure
        
        # see unlist(phyloseq::distanceMethodList)
        
        if (! measure %in% c("jsd", "manhattan", "euclidean", "bray", "jaccard", "canberra")){
                rv$infoText_6 <- "Sorry. Allowed beta diversity distance measures are: jsd, manhattan, euclidean, bray, caberra, and jaccard today."
                return()
                
        }
        
        # --- the input parameter check ---
        group_var <- rv$group_var
        
        group_var_levels <- rv$group_var_levels 
        
        color_levels <- rv$color_levels
        
        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                               color_levels = color_levels, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_6 <- error_message
                return()
        }
        # ------
        
        
        # - for relative abundance I want it to be relative abundances, for Manhatten it has to be relative abundances otherwise there is an error: -
        if (input$tcaType6 == "library size (RA)" || measure == "manhattan"){
                ps <- phyloseq::transform_sample_counts(ps, function(x){x/sum(x)})
        }
        # --
        
        
        dist_list <- calc_beta_div_distances(ps, dist_methods = measure, group_var = group_var, compare = group_var_levels)
        
        rv$betaDivDist <- dist_list
        
        rv$infoText_6 <- "Beta diversity distances have been calculated for all samples defined by group_var_levels."
        
})
# ----





# -- calc BetaDiversity --
observeEvent(input$plotPCoA, {
        
        rv$infoText_6 <- NULL
        rv$Tr_PCoA <- NULL
        rv$adonis <- NULL
        
        
        if (is.null(rv$phylum_colors)) {
                rv$infoText_6 <- "Sorry. You need to calculate phylum colors first."
                return()
        }
        
        if (input$tcaType6 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType6 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType6 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType6 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_6 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_6 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered6 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_6 <- "Sorry. You have to do the filtering first."
                        return()
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
        }
        
        
        measure <- input$beta_measure
        
        # see unlist(phyloseq::distanceMethodList)
        
        if (! measure %in% c("jsd", "manhattan", "euclidean", "bray", "jaccard", "canberra")){
                rv$infoText_6 <- "Sorry. Allowed beta diversity distance measures are: jsd, manhattan, euclidean, bray, caberra, and jaccard today."
                return()
                
        }
        
        
        phylum_colors <- rv$phylum_colors
        
        
        # --- the input parameter check ---
        group_var <- rv$group_var
        
        group_var_levels <- rv$group_var_levels 
        
        color_levels <- rv$color_levels
        
        compare <- rv$compare
        
        error_message <- check_user_parameters_4(group_var = group_var, group_var_levels = group_var_levels, 
                                                 color_levels = color_levels, compare = compare, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_6 <- error_message
                return()
        }
        # ------
        
        
        
        # - for relative abundance I want it to be relative abundances, for Manhatten it has to be relative abundances otherwise there is an error: -
        # NB: probably only relevant for calculating the distances but does not harm
        if (input$tcaType6 == "library size (RA)" || measure == "manhattan"){
                ps <- phyloseq::transform_sample_counts(ps, function(x){x/sum(x)})
        }
        # --

        
        # - calculate adonis Table -
        dist_list <- rv$betaDivDist
        
        if (is.null(dist_list)){
                rv$infoText_6 <- "You have to calculate the distances first"
                return()
        }
        
        # -- check that same group_var_levels were used for dist_list calculation --
        # NB: in  calc_beta_div_distances you use group_var_levels for compare, so you calculate always distances for all levels even if compare is smaller, could be changed!
        group_factor <- sample_data(ps)[[group_var]]
        # make sure group_factor fits to dist_list, i.e. only keep samples covered by compare = group_var_levels!
        group_factor <- factor(group_factor[group_factor %in% group_var_levels], levels = group_var_levels, ordered = T)
        
        dist_objj <- dist_list[[1]]
        
        if (length(dist_objj) != length(group_factor)*(length(group_factor)-1)/2) {
                rv$infoText_6 <- "You must use the same group_var_levels as was used for calculation of the distances."
                return()
        }
        # ----
        
        # -- adjust dist_list and group_factor based on compare for adonis (nb: makes sure adonis is only calculated for compare levels while stil allowing the gray option in the plot) --
        
        dist_list_adonis <- dist_list
        
        if (!all(group_var_levels %in% compare)){ # all(compare %in% group_var_levels) is tested in check_user_parameters_4
                
                group_factor <- sample_data(ps)[[group_var]]
                # make sure group_factor fits to dist_list, i.e. only keep samples covered by compare = group_var_levels!
                group_factor <- factor(group_factor[group_factor %in% compare], levels = compare, ordered = T)
                
                keepSamples <- sample_names(ps)[sample_data(ps)[[group_var]] %in% compare]
                
                dist_list_adonis <- lapply(dist_list_adonis, filter_dist_obj, keepNames = keepSamples)
        }
        # ----
        
        if (length(levels(group_factor)) > 1){ # no adonis test if only one factor level left because that would cause an error
                
                adonis_list <- lapply(dist_list_adonis, function(dist_obj){
                        loop_vegan_adonis(dist_obj = dist_obj, group_fac = group_factor)
                })
                
        } else {
                adonis_list <- NULL
        }
        # --
        
        
        # - generate a pcoa plot -
        gray_levels <- NULL
        
        if (!all(group_var_levels %in% compare)){
                
                if (input$compareUse == "remove other levels"){
                        dist_list <- dist_list_adonis
                        color_levels <- color_levels[names(color_levels) %in% compare]
                } else if (input$compareUse == "keep other levels as gray dots"){
                        gray_levels <- compare
                }
        }
        
        
        
        if (input$pcoaCorrection == "yes"){
                coordCor <- TRUE
        } else {
                coordCor <- FALSE
        }
        
        
        if (input$ellipse) {
                ellipse <- TRUE
                ellipse_level <- input$ellipse_level
                if (ellipse_level < 0 || ellipse_level > 1){
                        rv$infoText_6 <- "ellipse_level must be between 0 and 1"
                        return()
                }
        } else {
                ellipse <- FALSE
                ellipse_level <- 0.95
        }
        
        # - the pcoa needs categorial data as group_var, therefore it throws an error if group_var is numerical -
        if (is.numeric(sample_data(ps)[[group_var]])){
                sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]])
        }
        # --
        
        pcoas <- calc_ordination_from_distances(physeq = ps, group_var = group_var, dist_list = dist_list, color_levels = color_levels, ordination_type = "PCoA",
                                                shape = NULL, coord_cor = coordCor, phylum_colors = phylum_colors, ellipses = ellipse, ellipse_level = ellipse_level,
                                                gray_levels = gray_levels, paired_var = NULL)
        
        rv$Tr_PCoA <- pcoas[["ordination_Tr_samples"]]
        
        rv$adonis <- adonis_list[[1]]
        
        rv$infoText_6 <- "Beta diversity analysis has been performed."
        
})
# ----





# -- render the adonis table --
output$adonis <- renderTable({
        if(is.null(rv$adonis)){
                NULL
        } else {
                adonisShow <- rv$adonis
                adonisShow
        }
}, digits = 4, sanitize.text.function = function(x) x, caption = "adonis test result", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- output plot PCoA --
output$plotPCoA <- renderPlot({
        if (is.null(rv$Tr_PCoA)) {
                NULL
        } else {
                rv$Tr_PCoA
        }
},
height = 500)
# ----

