# - Tab7: fisher Test -
# -- output infoText --
output$infoText_7 <- renderText({
        rv$infoText_7
})
# ----



# -- render the overviewViews table --
output$overviewView7 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- calc fisher test --
observeEvent(input$fisher, {
        
        rv$infoText_7 <- NULL
        
        rv$Tr_fisher <- NULL
        rv$fisherTable <- NULL
        
        if (is.null(rv$phylum_colors)) {
                rv$infoText_7 <- "Sorry. You need to calculate phylum colors first."
                return()
        }
        
        if (input$tcaType7 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType7 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType7 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType7 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_7 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_7 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                return()
        }
        
        
        if (input$filtered7 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_7 <- "Sorry. You have to do the filtering first."
                        return()
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
        }
        
        
        phylum_colors <- rv$phylum_colors
        
        
        # --- the input parameter check ---
        group_var <- rv$group_var
        
        group_var_levels <- rv$group_var_levels 
        
        color_levels <- rv$color_levels
        
        compare <- rv$compare
        
        
        error_message <- check_user_parameters_4(group_var = group_var, group_var_levels = group_var_levels, compare = compare, 
                                                 color_levels = color_levels, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_7 <- error_message
                return()
        }
        # ------
        
        # - check zero color input -
        zero_color <- input$zero_color7
        
        if (!areColors(zero_color)){
                rv$infoText_7 <- "Sorry. The color for the zero counts is not an R color, please change."
                return()
                
        }
        # --
        
        # - Check that compare defines exactly two levels! -
        if (length(compare) != 2){
                rv$infoText_7 <- "Sorry, you have to compare exactly two levels here, please change the compare field in the Parameters tab."
                return()
        }
        # --
        
        
        physeq_to_test <- ps
        
        
        
        # - prune the ps if there are more than 2 group_var_levels and the user wants to see only 2 -
        if (input$restrictto2levels && length(group_var_levels) > 2) {
                SS <- get_sample_data(physeq_to_test)
                SS <- SS[SS[[group_var]] %in% compare,]
                keepSamples <- rownames(SS)
                physeq_to_test <- prune_samples(samples = keepSamples, physeq_to_test)
                # physeq_to_test <- phyloseq::subset_taxa(physeq_to_test, taxa_sums(physeq_to_test) != 0) # Don't understand but throws error therefore:
                keepTaxa <- taxa_names(physeq_to_test)[taxa_sums(physeq_to_test) != 0]
                physeq_to_test <- prune_taxa(keepTaxa, physeq_to_test)
        }
        # --
        
        
        diff_ab_df <- test_diffs_in_prevalence_single(physeq = physeq_to_test, group_var = group_var, compare = compare, p.adj.method = "fdr", minCount = 0L, symnum.args = rv$symnum.args)
        
        hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = "fdr")
        
        taxa_hit_df <- hit_list[["hit_table"]]
        
        # - for the heatmap, see if the user wants to see it in rel. abundance -
        if (input$ra7){
                physeq_to_test <- phyloseq::transform_sample_counts(physeq_to_test, function(x){x/sum(x)})
        }
        # --
        
        
        significance_colors <- brewer.pal(4, "Reds")
        significance_colors <- c(rev(significance_colors), "gray", "violet")
        names(significance_colors) = c("****", "***", "**", "*", "ns", "?")
        taxa_colors <- list("signi_adj" = significance_colors, "Phylum" = phylum_colors)
        sample_colors <- list(color_levels)
        names(sample_colors) <- group_var
        taxa_annotation <- taxa_hit_df$Annotation
        
        
        max_abundance_for_color <- input$max_abundance_for_colorF
        
        if (max_abundance_for_color == ""){
                max_abundance_for_color <- max(otu_table(physeq_to_test))
        } else {
                max_abundance_for_color <- as.numeric(max_abundance_for_color)
                
        }
        
        
        OTU <- as(otu_table(physeq_to_test), "matrix")
        
        if (is.na(max_abundance_for_color) || max_abundance_for_color < min(OTU[OTU > 0]) || max_abundance_for_color > max(OTU, na.rm = TRUE)) {
                rv$infoText_7 <- "Sorry. max_abundance_for_color either not numeric, too small or too high:). Change it."
                return()
                
        }
        
        if (input$log7){
                logger <- TRUE 
        } else {
                logger <- FALSE
        }
        
        max_shown <- as.numeric(input$maxShown7)
        
        if (is.na(max_shown) || max_shown < 1) {
                max_shown <- 40
        }
        
        
        rv$Tr_fisher <- plot_heatmap_physeq(physeq_to_test, sample_colors = sample_colors, taxa_info_df = head(taxa_hit_df, max_shown), taxa_colors = taxa_colors, 
                                            taxa_annotation = head(taxa_annotation, max_shown), max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                                            zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = logger, drop_color_levels = TRUE,
                                            border_color = NA, 
                                            cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                                            annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 22, 
                                            fontsize_row = 12, fontsize_col = 12, fontsize_number = 12)
        
        rv$fisherTable <- head(taxa_hit_df, max_shown)
        
        rv$infoText_7 <- "Fisher test done"
        
})
# ----



# -- render the fisher result table --
output$tableFisher <- renderTable({
        if(!is.null(rv$fisherTable)){
                fisherShow <- rv$fisherTable
                fisherShow
        } else {
                NULL
        }
}, digits = 4, sanitize.text.function = function(x) x, caption = "Results of fisher test result", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- output plot fisher --
output$plotFisher <- renderPlot({
        if (is.null(rv$Tr_fisher)) {
                NULL
        } else {
                grid::grid.newpage()
                grid::grid.draw(rv$Tr_fisher)
        }
},
height = 750)
# ----
# --
