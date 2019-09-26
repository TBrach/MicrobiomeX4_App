# - Tab14: Taxa ratio -
# -- output infoText --
output$infoText_14 <- renderText({
        rv$infoText_14
})
# ----


# -- render the overviewViews table --
output$overviewView14 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- plot Ratio Taxa --
observeEvent(input$plotRatioTaxa, {
        
        rv$infoText_14 <- NULL
        
        rv$Tr_taxaRatios <- NULL
        
        rv$tablepValsTaxaRatios <- NULL
        
        if (input$tcaType14 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType14 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType14 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType14 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_14 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_14 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                return()
        }
        
        
        if (input$filtered14 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_14 <- "Sorry. You have to do the filtering first."
                        return()
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
        }
        
        
        # --- the input parameter check ---
        group_var <- rv$group_var
        
        group_var_levels <- rv$group_var_levels 
        
        color_levels <- rv$color_levels
        
        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                               color_levels = color_levels, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_14 <- error_message
                return()
        }
        # ------
        
        
        numerator <- as.numeric(input$numerator)
        
        if (length(numerator) == 0 || is.na(numerator) || !(numerator %in% 1:ntaxa(ps))){
                rv$infoText_14 <- "Couldn't find numerator. Change input."
                return()
        }
        
        
        denominator <- as.numeric(input$denominator)
        
        if (length(denominator) == 0 || is.na(denominator) || !(denominator %in% 1:ntaxa(ps))){
                rv$infoText_14 <- "Couldn't find denominator. Change input."
                return()
        }
        
        
        keepTaxa <- taxa_names(ps)[c(numerator, denominator)]
        
        ps.pruned <- phyloseq::prune_taxa(taxa = keepTaxa, ps)
        
        mdf <- psmelt(ps.pruned)
        
        TT <- as.data.frame(unclass(tax_table(ps.pruned)))
        TT$Annotation <- get_taxon_names_plusTL(TT)
        LookUp <- data.frame(Index = c(numerator, denominator), OTU = keepTaxa, Annotation = TT$Annotation[match(keepTaxa, rownames(TT))])
        
        mdf$Index <- LookUp$Index[match(mdf$OTU, LookUp$OTU)]
        
        mdf <- mdf[mdf[[group_var]] %in% group_var_levels,]
        
        mdf_ratio <- group_by_(mdf, "Sample", group_var)
        
        mdf_ratio <- summarise(mdf_ratio, Ratio = Abundance[Index == numerator]/Abundance[Index == denominator])
        
        mdf_ratio$Labeller <- paste0(numerator, "/", denominator, "_", paste(LookUp$Annotation, collapse = "/"))
        
        mdf_ratio <- mdf_ratio[is.finite(mdf_ratio$Ratio),]
        
        
        mdf_ratio[[group_var]] <- factor(mdf_ratio[[group_var]], levels = names(color_levels), ordered = TRUE)
        
        group_fac <- factor(mdf_ratio[[group_var]])
        fac_levels <- levels(group_fac)
        
        comparisonList <- get_unique_facLevel_combinations(fac_levels)
        
        Tr <- ggplot(mdf_ratio, aes_string(x = group_var, y = "Ratio", col = group_var))
        
        Tr <- Tr +
                geom_boxplot(outlier.colour = NA) +
                geom_jitter(height = 0) +
                facet_wrap(~ Labeller, scales = "free_y") +
                scale_color_manual(values = color_levels) +
                theme_bw() +
                xlab("")
        
        if (input$logRatio){
                Tr <- Tr + scale_y_log10()
                mdf_ratio$Ratio <- log(mdf_ratio$Ratio)
                mdf_ratio <- mdf_ratio[is.finite(mdf_ratio$Ratio),]
        }
        
        
        
        
        if (all(table(mdf_ratio[[group_var]]) > 2)) {
                
                Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = "t.test", hide.ns = FALSE, symnum.args = rv$symnum.args)
                
                formulaa <- as.formula(paste("Ratio ~", group_var, sep = " "))
                
                pVals <- compare_means(formula = formulaa, data = mdf_ratio, group.by = "Labeller", method = "t.test", p.adjust.method = "fdr", symnum.args = rv$symnum.args)
        } else {
                pVals <- NULL
        }
        
        rv$Tr_taxaRatios <- Tr
        
        rv$tablepValsTaxaRatios <- pVals
        
        
        rv$infoText_14 <- "Selected abundance ratios have been plotted and statistical test has been performed if possible."
        
        # ------
        
})
# ----



# -- render the pVals of the ratio comparison --
output$pValsRatioTaxa <- renderTable({
        if(!is.null(rv$tablepValsTaxaRatios)){
                rv$tablepValsTaxaRatios
        } else {
                NULL
        }
}, rownames = TRUE, digits = 0, sanitize.text.function = function(x) x, caption = "p Value table", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- output plot Taxa Abundances --
output$plotTaxaRatios <- renderPlot({
        if (is.null(rv$Tr_taxaRatios)) {
                rv$Tr_taxaRatios
        } else {
                rv$Tr_taxaRatios
        }
}, 
height = 500)
# ----
