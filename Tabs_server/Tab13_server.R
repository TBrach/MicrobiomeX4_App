
# - Tab13: Taxa finder -
# -- output infoText --
output$infoText_13 <- renderText({
        rv$infoText_13
})
# ----

# -- render the overviewViews table --
output$overviewView13 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----


# -- search Taxa --
observeEvent(input$taxaSearch, {
        
        rv$infoText_13 <- NULL
        
        rv$taxaFindTable <- NULL
        
        rv$Tr_taxaAbundance <- NULL
        
        rv$tablepValsTaxa <- NULL
        
        
        if (input$tcaType13 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType13 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType13 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType13 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_13 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_13 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                return()
        }
        
        
        if (input$filtered13 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_13 <- "Sorry. You have to do the filtering first."
                        return()
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
        }
        
        tax_level <- input$taxLevelSearch
        
        TT <- get_tax_table(ps)
        
        searchWord <- input$searchWord
        
        indexes <- grep(pattern = searchWord, TT[[tax_level]], ignore.case = TRUE)
        
        if (length(indexes) == 0){
                rv$infoText_13 <- "No taxon matched your search."
                return()
        }
        
        TT_show <- TT[indexes,]
        
        TT_show <- cbind(data.frame(Index = indexes, Annotation = get_taxon_names(TT_show)), TT_show)
        
        rv$taxaFindTable <- TT_show
        
        rv$infoText_13 <- "Taxa have been found"
        
        # ------
        
})
# ----



# -- search Taxa --
observeEvent(input$plotTaxa, {
        
        rv$infoText_13 <- NULL
        
        rv$Tr_taxaAbundance <- NULL
        
        rv$tablepValsTaxa <- NULL
        
        if (is.null(rv$phylum_colors)) {
                rv$infoText_13 <- "Sorry. You need to calculate phylum colors first."
                return()
        }
        
        
        if (input$tcaType13 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType13 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType13 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType13 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_13 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_13 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                return()
        }
        
        
        if (input$filtered13 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_13 <- "Sorry. You have to do the filtering first."
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
        
        error_message <- check_user_parameters_4(group_var = group_var, group_var_levels = group_var_levels, 
                                                 color_levels = color_levels, compare = compare, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_13 <- error_message
                return()
        }
        # ------
        
        
        userTaxa <- input$userTaxa
        
        userTaxa <- strsplit(x = userTaxa, split = ",")
        
        userTaxa <- as.numeric(unlist(userTaxa))
        
        userTaxa <- userTaxa[!is.na(userTaxa)]
        
        if (length(userTaxa) == 0){
                rv$infoText_13 <- "No correct indexes were specified. Please give a comma-separated list of the indexes of the taxa you want to plot."
                return()
        }
        
        userTaxa <- sort(userTaxa)
        
        if (!all(userTaxa %in% 1:ntaxa(ps))){
                rv$infoText_13 <- "Some indexes were outside of the range of the chosen phyloseq object so please correct your input and try again."
                return()
        }
        
        keepTaxa <- taxa_names(ps)[userTaxa]
        
        
        if (input$RA){ 
                ps <- phyloseq::transform_sample_counts(ps, function(x){x/sum(x)})
        }
        
        ps.pruned <- phyloseq::prune_taxa(taxa = keepTaxa, ps)
        
        
        if (!all(group_var_levels %in% compare)){ # all(compare %in% group_var_levels) is tested in check_user_parameters_4
                
                keepSamples <- sample_names(ps)[sample_data(ps)[[group_var]] %in% compare]
                
                ps.pruned <- phyloseq::prune_samples(samples = keepSamples, ps.pruned)
        }
        
        
        
        mdf <- psmelt(ps.pruned)
        
        TT <- get_tax_table(ps.pruned)
        TT$Annotation <- get_taxon_names_plusTL(TT)
        LookUp <- data.frame(Index = userTaxa, OTU = keepTaxa, Annotation = TT$Annotation[match(keepTaxa, rownames(TT))])
        
        LookUp <- dplyr::mutate(LookUp, Labeller = paste(Index, OTU, Annotation, sep = "_"))
        
        mdf$Index <- LookUp$Index[match(mdf$OTU, LookUp$OTU)]
        mdf$Annotation <- LookUp$Annotation[match(mdf$OTU, LookUp$OTU)]
        mdf$Labeller <- paste(mdf$Index, mdf$OTU, mdf$Annotation, sep = "_")
        
        mdf$Labeller <- factor(mdf$Labeller, levels = LookUp$Labeller, ordered = TRUE) 
        
        mdf <- mdf[mdf[[group_var]] %in% group_var_levels,]
        
        if (input$pool){
                
                mdf$Labeller <- paste0("Indexes_", paste(userTaxa, collapse = "_"))
                
                mdf <- group_by_(mdf, "Sample", "Phylum", group_var, "Labeller")
                
                mdf <- summarise(mdf, Abundance = sum(Abundance))
        }
        
        
        
        mdf[[group_var]] <- factor(mdf[[group_var]], levels = names(color_levels), ordered = TRUE)
        
        group_fac <- factor(mdf[[group_var]])
        fac_levels <- levels(group_fac)
        
        # --- check Ref ---
        Ref <- input$ref_13
        
        if (Ref == ""){
                Ref <- NULL
        } else {
                if (! Ref %in% names(color_levels)){
                        rv$infoText_13 <- "Ref is not a level in compare/group_var_levels, you need to change it."
                        return()
                        
                }
        }
        # ------
        
        if (is.null(Ref)){
                comparisonList <- get_unique_facLevel_combinations(fac_levels)
        }
        
        # phylum_colors <- rv$phylum_colors
        
        # Tr <- ggplot(mdf, aes_string(x = group_var, y = "Abundance", fill = "Phylum", col = group_var))
        Tr <- ggplot(mdf, aes_string(x = group_var, y = "Abundance", col = group_var))
        
        
        Tr <- Tr +
                geom_boxplot(outlier.colour = NA) +
                geom_jitter(height = 0, width = 0.2) +
                facet_wrap(~ Labeller, scales = "free_y") +
                scale_color_manual(values = color_levels) +
                # scale_fill_manual(values = phylum_colors) +
                theme_bw() +
                xlab("")
        
        if (input$log){
                Tr <- Tr + scale_y_log10()
                mdf$Abundance <- log(mdf$Abundance)
                mdf <- mdf[is.finite(mdf$Abundance),]
        }
        
        
        
        if (all(table(mdf[[group_var]]) > 2)) {
                
                if (is.null(Ref)){
                        Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = "t.test", hide.ns = FALSE, symnum.args = rv$symnum.args)
                } else {
                        Tr <- Tr + ggpubr::stat_compare_means(ref.group = Ref, label = "p.signif", method = "t.test", hide.ns = FALSE, symnum.args = rv$symnum.args)
                }
                
                
                formulaa <- as.formula(paste("Abundance ~", group_var, sep = " "))
                
                # - make sure mdf[[group_var]] is unordered factor otherwise this throws an error, not sure why: -
                mdf2 <- mdf
                mdf2[[group_var]] <- factor(mdf2[[group_var]], ordered = F)
                
                pVals <- try(compare_means(formula = formulaa, ref.group = Ref, data = mdf2, group.by = "Labeller", method = "t.test", p.adjust.method = "fdr", symnum.args = rv$symnum.args))
                if (inherits(pVals, "try-error")){
                        pVals <- NULL
                }
                # --
        } else {
                pVals <- NULL
        }
        
        rv$Tr_taxaAbundance <- Tr
        
        rv$tablepValsTaxa <- pVals
        
        rv$infoText_13 <- "selected Taxa have been plotted."
        
        # ------
        
})
# ----




# -- render the foundTaxa table --
output$foundTaxa <- renderTable({
        if(!is.null(rv$taxaFindTable)){
                taxaShow <- rv$taxaFindTable
                taxaShow
        } else {
                NULL
        }
}, rownames = TRUE, digits = 0, sanitize.text.function = function(x) x, caption = "Taxa found in your search", caption.placement = getOption("xtable.caption.placement", "top"))
# ----


# -- render the pVals of the taxa comparison --
output$pValsTaxa <- renderTable({
        if(!is.null(rv$tablepValsTaxa)){
                rv$tablepValsTaxa
        } else {
                NULL
        }
}, rownames = TRUE, digits = 0, sanitize.text.function = function(x) x, caption = "p Value table", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- output plot Taxa Abundances --
output$plotTaxaAbundances <- renderPlot({
        if (is.null(rv$Tr_taxaAbundance)) {
                rv$Tr_taxaAbundance
        } else {
                rv$Tr_taxaAbundance
        }
}, 
height = 700)
# ----

