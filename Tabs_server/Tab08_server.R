# - Tab8: DESeq2 Test -
# -- output infoText --
output$infoText_8 <- renderText({
        rv$infoText_8
})
# ----


# -- render the overviewViews table --
output$overviewView8 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- calc DESeq2 --
observeEvent(input$DESeq2, {
        
        rv$infoText_8 <- NULL
        
        rv$Tr_DESeq2 <- NULL
        rv$DESeq2Table <- NULL
        
        
        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                rv$infoText_8 <- "Sorry. You need a loaded raw count phyloseq object, and you need phylum_colors. Make sure these steps are done."
                return()
        }
        
        ps <- rv$ps
        
        if (input$tcaType8 == "DESeq_gm_exclZero") {
                SFs <- rv$SFs
        } else if (input$tcaType8 == "DESeq_poscounts") {
                SFs <- rv$SFs_poscounts
        } else if (input$tcaType8 == "library size (RA)"){
                SFs <- rv$SFs_RA
        }
        
        if (is.null(SFs)){
                rv$infoText_8 <- "Sorry. You need to calculate size factors first!"
                return()
                
        }
        
        
        if (input$filtered8 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_8 <- "Sorry. You have to do the filtering first."
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
                rv$infoText_8 <- error_message
                return()
        }
        # ------
        
        # - test the zero_color input -
        zero_color <- input$zero_color8
        
        if (!areColors(zero_color)){
                rv$infoText_8 <- "Sorry. The color for the zero counts is not an R color, please change."
                return()
                
        }
        # --
        
        # - Check that compare defines exactly two levels! -
        if (length(compare) != 2){
                rv$infoText_8 <- "Sorry, you have to compare exactly two levels here, please change the compare field in the Parameters tab."
                return()
        }
        # --
        
        
        physeq_to_test <- ps
        
        # - prune the ps and SFs if there are more than 2 group_var_levels and the user wants to see only 2 -
        if (input$restrictto2levels_8 && length(group_var_levels) > 2) {
                SS <- get_sample_data(physeq_to_test)
                SS <- SS[SS[[group_var]] %in% compare,]
                keepSamples <- rownames(SS)
                SFs <- SFs[names(SFs) %in% keepSamples]
                physeq_to_test <- prune_samples(samples = keepSamples, physeq_to_test)
                # physeq_to_test <- phyloseq::subset_taxa(physeq_to_test, taxa_sums(physeq_to_test) != 0) # Don't understand but throws error therefore:
                keepTaxa <- taxa_names(physeq_to_test)[taxa_sums(physeq_to_test) != 0]
                physeq_to_test <- prune_taxa(keepTaxa, physeq_to_test)
        }
        # --
        
        
        res_list <- test_differential_abundance_DESeq2single(physeq = physeq_to_test, group_var = group_var, compare = compare, SFs = SFs, p.adjust.method = "fdr", symnum.args = rv$symnum.args, cooksCutoff = TRUE)
        
        diff_ab_df <- res_list[[1]]
        physeq_to_test <- res_list[[2]]
        
        
        hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = "fdr")
        
        taxa_hit_df <- hit_list[["hit_table"]]
        
        # - for the heatmap, see if the user wants to see it in rel. abundance -
        if (input$ra8){
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
        
        
        max_abundance_for_color <- input$max_abundance_for_colorD
        
        if (max_abundance_for_color == ""){
                max_abundance_for_color <- max(otu_table(physeq_to_test))
        } else {
                max_abundance_for_color <- as.numeric(max_abundance_for_color)
                
        }
        
        OTU <- as(otu_table(physeq_to_test), "matrix")
        
        if  (is.na(max_abundance_for_color) || max_abundance_for_color < min(OTU[OTU > 0]) || max_abundance_for_color > max(OTU, na.rm = TRUE)) {
                rv$infoText_8 <- "Sorry. max_abundance_for_color either not numeric, too small or too high:). Change it."
                return()
                
        }
        
        
        
        
        if (input$log8){
                logger <- TRUE 
        } else {
                logger <- FALSE
        }
        
        max_shown <- as.numeric(input$maxShown8)
        
        if (is.na(max_shown) || max_shown < 1) {
                max_shown <- 40
        }
        
        rv$Tr_DESeq2 <- plot_heatmap_physeq(physeq_to_test, sample_colors = sample_colors, taxa_info_df = head(taxa_hit_df, max_shown), taxa_colors = taxa_colors, 
                                            taxa_annotation = head(taxa_annotation, max_shown), max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                                            zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = logger, drop_color_levels = TRUE,
                                            border_color = NA, 
                                            cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                                            annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 22, 
                                            fontsize_row = 12, fontsize_col = 12, fontsize_number = 12)
        
        rv$DESeq2Table <- head(taxa_hit_df, max_shown)
        
        
        rv$infoText_8 <- "DESeq2 test done."
        
})
# ----



# -- render the DESeq2 result table --
output$tableDESeq2 <- renderTable({
        if(!is.null(rv$DESeq2Table)){
                DESeq2Show <- rv$DESeq2Table
                DESeq2Show
        } else {
                NULL
        }
}, digits = 4, sanitize.text.function = function(x) x, caption = "Results of DESeq2 analysis", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- output plot DESeq2 --
output$plotDESeq2 <- renderPlot({
        if (is.null(rv$Tr_DESeq2)) {
                NULL
        } else {
                grid::grid.newpage()
                grid::grid.draw(rv$Tr_DESeq2)
        }
},
height = 750)
# ----
# --    