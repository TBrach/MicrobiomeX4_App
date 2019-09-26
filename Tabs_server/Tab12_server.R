
# - Tab12: alpha diversity -
# -- output infoText --
output$infoText_12 <- renderText({
        rv$infoText_12
})
# ----


# -- render the overviewViews table --
output$overviewView12 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- calc alpha Diversity --
observeEvent(input$calculatealphDiv, {
        
        rv$infoText_12 <- NULL
        
        rv$Tr_alpha <- NULL
        
        rv$alphaTable <- NULL
        
        
        if (input$tcaType12 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType12 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType12 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType12 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_12 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_12 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered12 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_12 <- "Sorry. You have to do the filtering first."
                        return()
                        
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
                
        }
        
        
        # --- the input parameter check ---
        group_var <- rv$group_var
        
        group_var_levels <- rv$group_var_levels 
        
        color_levels <- rv$color_levels
        
        compare <- rv$compare
        
        
        error_message <- check_user_parameters_4(group_var = group_var, group_var_levels = group_var_levels, compare = compare, 
                                                 color_levels = color_levels, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_12 <- error_message
                return()
        }
        # ------
        
        
        # --- adjust color_levels based on compare ---
        color_levels <- color_levels[names(color_levels) %in% compare]
        # ------
        
        # --- check Ref ---
        Ref <- input$ref_12
        
        if (Ref == ""){
                Ref <- NULL
        } else {
                if (! Ref %in% names(color_levels)){
                        rv$infoText_12 <- "Ref is not a level in compare/group_var_levels, needs to be changed."
                        return()
                        
                }
        }
        # ------
        
        
        # --- NB: estimateR.default, in phyloseq::estimate_richness accepts only integers therefore: ---
        OTU <- as(otu_table(ps), "matrix")
        if (mode(OTU) != "integer"){
                OTU <- round(OTU, 1)
                mode(OTU) <- "integer"
                otu_table(ps) <- otu_table(OTU, taxa_are_rows = taxa_are_rows(ps))
        }
        # -------
        
        
        if (input$rarify) {
                rare_level <- min(sample_sums(ps))
                count_table_rare <- vegan::rrarefy(as(otu_table(ps), "matrix"), sample = rare_level)
                
                otu_table(ps) <- otu_table(count_table_rare, taxa_are_rows = taxa_are_rows(ps))
        }
        
        
        alpha_div_measures <- input$alpha_measure
        
        DF_alpha_list <- calc_alphadiv_plusLmResids(physeq = ps, measures = alpha_div_measures, group_var = group_var,
                                                    compare = names(color_levels))
        
        lm_fitlist <- DF_alpha_list[[2]]
        DF_alpha <- DF_alpha_list[[1]]
        
        alpha_div_pVals <- calc_pVals_alphdiv(DF_alpha = DF_alpha, measures = alpha_div_measures, Ref = Ref, group_var = group_var, compare = names(color_levels), test = "t.test")
        
        
        alpha_div_boxplots <- boxplots_alphdiv(DF_alpha = DF_alpha, measures = alpha_div_measures, Ref = Ref, group_var = group_var, shape = NULL, color_levels = color_levels, test = "t.test", hide.ns = FALSE)
        
        alpha_div_lmPlots <- lmPlots_alphdiv(DF_alpha = DF_alpha, lm_fitlist = lm_fitlist, measures = alpha_div_measures, group_var = group_var, shape = NULL, color_levels = color_levels, test = "t.test")
        
        
        TrList <- c(alpha_div_boxplots, alpha_div_lmPlots)
        TrList <- TrList[order(names(TrList))] 
        
        
        rv$Tr_alpha <- TrList
        
        rv$alphaTable <- alpha_div_pVals
        
        
        rv$infoText_12 <- "alpha diversity analysis has been performed."
        
})
# ----



# -- render the alpha table --
output$tableAlpha <- renderTable({
        if(!is.null(rv$alphaTable)){
                alphaShow <- rv$alphaTable
                alphaShow
        } else {
                NULL
        }
}, digits = 4, sanitize.text.function = function(x) x, caption = "alpha diversity test", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- output alpha diversity plot --
output$alphaDiv <- renderPlot({
        if (is.null(rv$Tr_alpha)) {
                NULL
        } else {
                do.call("grid.arrange", c(rv$Tr_alpha, ncol = 3))
        }
},
height = 500)
# ----
