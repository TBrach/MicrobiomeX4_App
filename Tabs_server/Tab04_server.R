
# -- output infoText --
output$infoText_4 <- renderText({
        rv$infoText_4
})
# ----


# -- render the overviewViews table --
output$overviewView4 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- calc heatmap --
observeEvent(input$calculateHM, {
        
        rv$infoText_4 <- NULL
        rv$Tr_HM <- NULL
        
        
        if (input$tcaType4 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType4 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType4 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType4 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_4 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_4 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered4 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_4 <- "Sorry. You have to do the filtering first."
                        return()
                        
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
                
        }
        
        # - check zero_color input -
        zero_color <- input$zero_color
        
        if (!areColors(zero_color)){
                rv$infoText_4 <- "Sorry. The color for the zero counts is not an R color, please change."
                return()
                
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
                rv$infoText_4 <- error_message
                return()
        }
        # ------
        
        # - prune the ps if the user wants to see only the ones compare says -
        if (input$restrictto2compare_4) {
                SS <- get_sample_data(ps)
                SS <- SS[SS[[group_var]] %in% compare,]
                keepSamples <- rownames(SS)
                ps <- prune_samples(samples = keepSamples, ps)
                # keepTaxa <- taxa_names(ps)[taxa_sums(ps) != 0]
                # ps <- prune_taxa(keepTaxa, ps)
        }
        # --
        
        
        # --- taxa_range option ---
        taxa_index_range <- input$taxa_index_range
        
        if (taxa_index_range != ""){
                
                taxa_index_range <- strsplit(x = taxa_index_range, split = ":")
                taxa_index_range <- round(as.numeric(unlist(taxa_index_range))[1:2])
                if (any(is.na(taxa_index_range))){
                        rv$infoText_4 <- "Sorry couldn't read taxa_index_range, please change input."
                        return()
                }
                start_index <- taxa_index_range[1]
                end_index <- taxa_index_range[2]
                
                if (start_index > end_index || start_index < 1 || end_index > ntaxa(ps)){
                        rv$infoText_4 <- "Sorry taxa_index_range did not fit to ps or start index was bigger than end index, please change input."
                        return()
                }
                
        } else {
                start_index <- 1
                end_index <- ntaxa(ps)
        }
        
        keepTaxa <- taxa_names(ps)[start_index:end_index]
        ps <- prune_taxa(keepTaxa, ps)
        # ------
        
        # - for the heatmap, see if the user wants to see it in rel. abundance -
        if (input$ra4){
                ps <- phyloseq::transform_sample_counts(ps, function(x){x/sum(x)})
        }
        # --
        
        
        # df <- as.data.frame(as(tax_table(ps), "matrix"))
        df <- get_tax_table(ps)
        taxa_annotation <- get_taxon_names_plusTL(df)
        taxa_annotation <- strsplit(taxa_annotation, split = "/")
        taxa_annotation <- sapply(taxa_annotation, `[`, 1)
        taxa_annotation <- make.unique(taxa_annotation)
        
        sample_colors <- list(color_levels)
        names(sample_colors) <- group_var
        
        max_abundance_for_color <- input$max_abundance_for_color
        
        if (max_abundance_for_color == ""){
                max_abundance_for_color <- max(otu_table(ps))
        } else {
                max_abundance_for_color <- as.numeric(max_abundance_for_color)
                
        }
        
        if (is.na(max_abundance_for_color) || max_abundance_for_color < 0 || max_abundance_for_color > max(otu_table(ps), na.rm = TRUE)) {
                rv$infoText_4 <- "Sorry. max_abundance_for_color either not numeric, too small or too high:). Change it."
                return()
                
        }
        
        if (input$log4){
                logger <- TRUE 
        } else {
                logger <- FALSE
        }
        
        rv$Tr_HM <- plot_heatmap_physeq(physeq = ps, sample_colors = sample_colors, taxa_info_df = NULL, taxa_colors = NULL,
                                        taxa_annotation = taxa_annotation, max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1),
                                        zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = logger, drop_color_levels = TRUE,
                                        border_color = NA,
                                        cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE,
                                        annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 16,
                                        fontsize_row = 8, fontsize_col = 8, fontsize_number = 12)
        
        
        sparsity <- round(100*sum(otu_table(ps) == 0)/length(otu_table(ps)), 2)
        
        rv$infoText_4 <- paste("Heatmap has been generated. There are ", ntaxa(ps), " taxa in the plotted data. Overall sparsity of the count table is ", sparsity, "%.", sep = "") 
})
# ----



# -- output plot samplesPS --
output$plotHM <- renderPlot({
        if (is.null(rv$Tr_HM)) {
                NULL
        } else {
                grid::grid.newpage()
                grid::grid.draw(rv$Tr_HM)
        }
},
height = 1200)
# ----
# --
