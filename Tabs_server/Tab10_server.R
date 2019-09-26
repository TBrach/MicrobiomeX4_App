
# - Tab10: Phylum analysis -
# -- output infoText --
output$infoText_10 <- renderText({
        rv$infoText_10
})
# ----

# -- render the overviewViews table --
output$overviewView10 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----


# -- calc Firmicutes to other phyla ratios --
observeEvent(input$firmicutes, {
        
        rv$infoText_10 <- NULL
        
        rv$Tr_phylum <- NULL
        
        rv$firmicutesTable <- NULL
        
        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                rv$infoText_10 <- "Sorry. A phyloseq raw count object must be present and phylum_colors have to be determined. Do these steps first."
                return()
        }
        
        if (input$tcaType10 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else {
                rv$infoText_10 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        
        if (input$filtered10 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_10 <- "Sorry. You have to do the filtering first."
                        return()
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
        }
        
        
        ps_phylum <- phyloseq::tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
        
        taxonomic_level <- "Phylum"
        
        phylum_colors <- rv$phylum_colors
        
        # --- the input parameter check ---
        group_var <- rv$group_var
        
        group_var_levels <- rv$group_var_levels 
        
        color_levels <- rv$color_levels
        
        compare <- rv$compare
        
        error_message <- check_user_parameters_4(group_var = group_var, group_var_levels = group_var_levels, compare = compare, 
                                                 color_levels = color_levels, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_10 <- error_message
                return()
        }
        # ------
        
        # --- adjust color_levels based on compare ---
        color_levels <- color_levels[names(color_levels) %in% compare]
        # ------
        
        # --- check Ref ---
        Ref <- input$ref_10
        
        if (Ref == ""){
                Ref <- NULL
        } else {
                if (! Ref %in% names(color_levels)){
                        rv$infoText_10 <- "Ref is not a level in compare/group_var_levels, you need to change it."
                        return()
                        
                }
        }
        # ------
        
        
        df <- get_tax_table(ps_phylum)
        taxa_annotation <- get_taxon_names(df[, 2:7])
        # taxa_annotation <- strsplit(taxa_annotation, split = "/")
        # taxa_annotation <- sapply(taxa_annotation, `[`, 1)
        taxa_annotation <- make.unique(taxa_annotation)
        
        taxa_order <- c(names(phylum_colors), taxa_annotation[!taxa_annotation %in% names(phylum_colors)])
        
        
        numerator_phylum <- input$numerator_phylum
        
        if (! numerator_phylum %in% df$Phylum){
                rv$infoText_10 <- "Sorry, numerator_phylum was not a phylum in the phyloseq object, please change the input."
                return()
        }
        
        
        denominator_phylum <- input$denominator_phylum
        
        if (denominator_phylum == ""){
                denominator_phylum <- NULL
        } else {
                denominator_phylum <- strsplit(x = denominator_phylum, split = ",")
                denominator_phylum <- unlist(denominator_phylum)
                denominator_phylum <- gsub(" ", "", denominator_phylum)
                
                if (!all(denominator_phylum %in% df$Phylum)){
                        rv$infoText_10 <- "Sorry, not all given denominator_phyla were a phylum in the phyloseq object, please change the input."
                        return()
                }
        }
        
        
        # NB: if you would like the order based on pValues/significance choose tax_order = NULL
        FirmicutesRatioPlots <- plot_taxa_ratios_AllLevels(physeq = ps_phylum, group_var = group_var, color_levels = color_levels, tax_names = taxa_annotation, taxa_nom = numerator_phylum, taxa_den = denominator_phylum, test = "wilcox.test", p_adjust_method = "fdr",
                                                           tax_order = taxa_order, Ref = Ref, symnum.args = rv$symnum.args, hide.ns = FALSE)
        
        rv$Tr_phylum <- FirmicutesRatioPlots[["Tr3"]]
        
        rv$firmicutesTable <- FirmicutesRatioPlots[["pValsLog"]]
        
        
        rv$infoText_10 <- "Ratios of Firmicutes to all other phyla have been calculated."
        
})
# ----


# -- calc phylum to phylum tile plots --
observeEvent(input$tile, {
        
        rv$infoText_10 <- NULL
        
        rv$Tr_phylum <- NULL
        
        rv$firmicutesTable <- NULL
        
        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                rv$infoText_10 <- "Sorry. A phyloseq raw count object must be present and phylum_colors have to be determined. Do these steps first."
                return()
        }
        
        if (input$tcaType10 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else {
                rv$infoText_10 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        
        if (input$filtered10 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_10 <- "Sorry. You have to do the filtering first."
                        return()
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
        }
        
        
        ps_phylum <- phyloseq::tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
        
        
        taxonomic_level <- "Phylum"
        
        
        phylum_colors <- rv$phylum_colors
        
        # --- the input parameter check ---
        group_var <- rv$group_var
        
        group_var_levels <- rv$group_var_levels 
        
        color_levels <- rv$color_levels
        
        compare <- rv$compare
        
        error_message <- check_user_parameters_4(group_var = group_var, group_var_levels = group_var_levels, compare = compare, 
                                                 color_levels = color_levels, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_10 <- error_message
                return()
        }
        # ------
        
        # --- adjust color_levels based on compare, NB: only two levels allowed! ---
        if (length(compare) != 2){
                rv$infoText_10 <- "For tile plots you have to define exactly two levels via compare!"
                return()
        }
        
        color_levels <- color_levels[names(color_levels) %in% compare]
        # ------
        
        
        
        
        raw_TbTmatrixes <- calculate_raw_TbTmatrixes(physeq = ps_phylum)
        
        
        df <- as.data.frame(as(tax_table(ps_phylum), "matrix"))
        taxa_annotation <- get_taxon_names(df[, 2:7])
        # taxa_annotation <- strsplit(taxa_annotation, split = "/")
        # taxa_annotation <- sapply(taxa_annotation, `[`, 1)
        taxa_annotation <- make.unique(taxa_annotation)
        
        taxa_order <- c(names(phylum_colors), taxa_annotation[!taxa_annotation %in% names(phylum_colors)])
        
        TbT_tile <- create_raw_TbT_TilePlot(TbTmatrixes = raw_TbTmatrixes, physeq = ps_phylum, group_var = group_var, color_levels = color_levels, signi_level = 0.05, tax_names = taxa_annotation, tax_order = taxa_order, test = "wilcoxon", p_adjust_method = "none")
        
        rv$Tr_phylum <- TbT_tile
        
        
        rv$infoText_10 <- "Tile plot of phylum/phylum ratios has been generated. If a tile is colored with the color of a sample group, this group has higher (row phylum)/(column phylum) ratios."
        
})
# ----



# -- render the Firmicutes result table --
output$tableFirmicutes <- renderTable({
        if(!is.null(rv$firmicutesTable)){
                firmicutesShow <- rv$firmicutesTable
                firmicutesShow
        } else {
                NULL
        }
}, digits = 5, sanitize.text.function = function(x) x, caption = "Results of Firmicutes ratios analysis", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- output plot phylum analysis --
output$plotPhylum <- renderPlot({
        if (is.null(rv$Tr_phylum)) {
                rv$Tr_phylum
        } else {
                rv$Tr_phylum
        }
},
height = 750)
# ----
# --