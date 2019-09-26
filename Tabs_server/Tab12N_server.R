
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
        
        # - set all items to NULL -
        rv$infoText_12 <- NULL
        
        rv$Tr_alpha <- NULL
        
        rv$alphaTable <- NULL
        # --
        
        # - get the ps based on which alpha diversity values should be calculated -
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
        # --
        
        
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
                        rv$infoText_12 <- "Ref is not a level in compare/group_var_levels, you need to change it."
                        return()
                        
                }
        }
        # ------
        
        
        # --- NB: estimateR.default, in phyloseq::estimate_richness accepts only integers therefore: ---
        # NB: ps_tca and ps_RA are not integers, but if you do this here, you change the sample_sums a bit, in general
        # for these choices in alpha diversity  you should ignore the Total corrections and only look at the uncorrected
        # the corrected versions are only for raw counts with unadjusted library sizes
        OTU <- as(otu_table(ps), "matrix")
        if (mode(OTU) != "integer"){
                OTU <- round(OTU, 1)
                mode(OTU) <- "integer"
                otu_table(ps) <- otu_table(OTU, taxa_are_rows = taxa_are_rows(ps))
        }
        # -------
        
        
        
        alpha_div_measures <- input$alpha_measure
        
        
        DF_alpha <- calc_alphadiv(physeq = ps, measures = alpha_div_measures, group_var = group_var,
                                  compare = names(color_levels))
        
        
        # - calculate correction factor directly from Total -
        DF_alpha <- dplyr::mutate(DF_alpha, 
                                  alphaCorrFac_Total = Total/gm_own(Total, zeros.count = FALSE))
        # --
        
        # - calculate correction factors from linear fits with Filtered reads and Total -
        # -- Here you have the option to do the linear fit only on Control samples if one group stands out drastically --
        # DF_alpha_ref <- dplyr::filter(DF_alpha, Diet == "Chow")
        DF_alpha_ref <- DF_alpha
        # ----
        
        # -- do the linear fit to Total --
        formulaa <- paste0(alpha_div_measures, " ~ Total")
        fit <- lm(formula = formulaa, data = DF_alpha_ref)
        # ----
        
        # -- predict the alpha diversity measure from the fit --
        prediction <- predict(fit, newdata = DF_alpha)
        # ----
        
        # -- Calculate and add the correcetion factors--
        corFactor <- prediction/gm_own(prediction, zeros.count = FALSE)
 
        corFactorDF <- as.data.frame(corFactor)
        DF_alpha <- cbind(DF_alpha, corFactorDF)
        
        DF_alpha[[paste0(alpha_div_measures[1], "_Totalcorrected")]] <- DF_alpha[[alpha_div_measures[1]]]/DF_alpha$alphaCorrFac_Total
        DF_alpha[[paste0(alpha_div_measures[1], "_lmTotalcorrected")]] <- DF_alpha[[alpha_div_measures[1]]]/DF_alpha$corFactor
        # ----
        # --
        
        # - plot the linear fit used for lmTotal correction -
        fitlists <- list(fit)
        names(fitlists) <- alpha_div_measures
        lmPlots_Total <- lmPlots_alphdiv2(DF_alpha = DF_alpha, lm_fitlist = fitlists, measures = alpha_div_measures, xval = "Total", group_var, shape = NULL, color_levels, test = "t.test", alpha = 1)
        # --
        
        # - plot the boxplots for uncorrected and corrected alpha div measure -
        AddOns <- c("", "_Totalcorrected", "_lmTotalcorrected")
        Combis <- expand.grid(alpha_div_measures, AddOns)
        PlotVariables <- apply(Combis, 1, paste, collapse = "")
        alpha_div_boxplots2 <- boxplots_alphdiv2(DF_alpha = DF_alpha,
                                                 PlotVariables = PlotVariables, 
                                                 group_var = group_var, shape = NULL,
                                                 facet_x = NULL,
                                                 facet_y =  NULL, color_levels = color_levels, Ref = Ref, test = "t.test", hide.ns = TRUE)
        # --
        
        
        TrList <- c(lmPlots_Total, alpha_div_boxplots2)
        
        rv$Tr_alpha <- TrList
        
        rv$alphaTable <- NULL
        
        
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
                do.call("grid.arrange", c(rv$Tr_alpha, ncol = 2))
        }
},
height = 700)
# ----
