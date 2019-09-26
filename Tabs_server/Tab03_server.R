
# -- output infoText --
output$infoText_3 <- renderText({
        rv$infoText_3
})
# ----


# -- render the overviewViews table --
output$overviewView3 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----


# # -- render the explorationViews table --
# output$explorationViews_2 <- renderTable({
#         if(!is.null(rv$ps)){
#                 psShow <- capture.output(rv$ps)
#                 psShow <- as.data.frame(psShow, nrow = 4)
#                 psShow <- psShow[2:4, , drop = FALSE]
#                 colnames(psShow) <- "object component"
#                 psShow <- tidyr::separate(psShow, col = "object component", into = c("object component", "dimension"), sep = ":")
#                 psShow
#         } else {
#                 NULL
#         }
# }, sanitize.text.function = function(x) x, caption = "Loaded phyloseq object", caption.placement = getOption("xtable.caption.placement", "top"))
# # ----



# -- calc size factors --
observeEvent(input$calcSFs, {
        
        rv$infoText_3 <- NULL
        rv$ps_tca <- NULL
        rv$ps_tca_poscounts <- NULL
        rv$ps_RA <- NULL
        rv$SFs <- NULL
        rv$SFs_poscounts <- NULL
        rv$SFs_RA <- NULL
        
        
        if (is.null(rv$ps)) {
                rv$infoText_3 <- "Sorry. You need to load a phyloseq object first (Tab 2)."
                return()
        }
        
        prevalence_for_sf<- as.numeric(input$prev_SFs)
        
        if (is.na(prevalence_for_sf) || prevalence_for_sf > 100 || prevalence_for_sf < 0){
                rv$infoText_3 <- "Sorry. Prevalence for DESeq_gm_exclZero must be between 0 and 100. Please change your input."
                return()
                
        }
        
        min_obs <- 0L
        
        ps_sf_filt <- phyloseq::filter_taxa(rv$ps, function(x){(sum(x > min_obs) > (prevalence_for_sf/100)*length(x))}, prune = TRUE)
        
        taxa_used_for_SFs <- ntaxa(ps_sf_filt)
        
        SFs <- calc_SFs(physeq = ps_sf_filt)
        
        
        group_var <- rv$group_var
        
        if(! group_var %in% colnames(sample_data(rv$ps))) {
                rv$infoText_3 <- "The given group_var is not a variable in the sample data of the loaded phyloseq object. A correct group variable is needed for size factor calculation using poscount method"
                return()
        }
        
        SFs_poscounts <- calc_SFs_DESeq(ps_sf_filt, type = "poscounts", group_var = group_var)
        
        SFs_RA <- sample_sums(rv$ps)/gm_own(sample_sums(rv$ps), zeros.count = FALSE)
        # If you want real relative abundances instead use:
        # SFs_RA <- sample_sums(rv$ps) 
        # But NB: will cause trouble with alpha diversity and also DESeq2!
        
        rv$SFs <- SFs
        rv$SFs_poscounts <- SFs_poscounts
        rv$SFs_RA <- SFs_RA
        
        library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs) 
        rv$ps_tca <- library_size_adjust_list[[1]]
        
        library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs_poscounts) 
        rv$ps_tca_poscounts <- library_size_adjust_list[[1]]
        
        library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs_RA) 
        rv$ps_RA <- library_size_adjust_list[[1]]
        
        rv$infoText_3 <- paste("Library size adjusted phyloseq objects have been calculated. ", taxa_used_for_SFs,  " taxa were used for the DESeq2 based SF calculations.", sep = "") 
})
# ----



# -- generate a barplot before after of selected library size adjustment --
observeEvent(input$barplotSFs, {
        
        rv$infoText_3 <- NULL
        rv$Tr_SFs <- NULL
        
        if (is.null(rv$ps)) {
                rv$infoText_3 <- "Sorry. No phyloseq object is loaded."
                return()
        }
        
        ps <- rv$ps
        
        if (is.null(rv$phylum_colors)) {
                rv$infoText_3 <- "Sorry. phylum_colors have not yet been calculated."
                return()
                
        }
        
        phylum_colors <- rv$phylum_colors
        
        
        if (input$tcaType3 == "raw counts"){
                plot_ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType3 == "DESeq_gm_exclZero"){
                plot_ps <- rv$ps_tca
        } else if (input$tcaType3 == "DESeq_poscounts") {
                plot_ps <- rv$ps_tca_poscounts
        } else if (input$tcaType3 == "library size (RA)") {
                plot_ps <- rv$ps_RA
        } else {
                rv$infoText_3 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(plot_ps)) {
                rv$infoText_3 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered3 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_3 <- "Sorry. You have to do the filtering first."
                        return()
                        
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
                ps_plot <- prune_taxa(keepTaxa, ps)
                
        }
        
        
        # --- the input parameter check ---
        group_var <- rv$group_var
        
        group_var_levels <- rv$group_var_levels 
        
        color_levels <- rv$color_levels
        
        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                               color_levels = color_levels, ps = ps)
        
        if (!is.null(error_message)){
                rv$infoText_3 <- error_message
                return()
        }
        # ------
        
        
        
        psP <- phyloseq::tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
        ps_tcaP <- phyloseq::tax_glom(plot_ps, taxrank = "Phylum", NArm = FALSE)
        rv$Tr_SFs <- plot_sample_bars_compare(physeq = psP, physeq2 = ps_tcaP, x = "Sample", y = "Abundance", group_var = group_var, color_levels = color_levels, color_sample_names = TRUE, fill = "Phylum", col_vec = phylum_colors, order_by_raw_counts = TRUE)
        
        
        rv$infoText_3 <- "Plot illustrating library size adjustment hast been generated." 
})
# ----



# -- output plot samplesPS --
output$plotSFs <- renderPlot({
        if (is.null(rv$Tr_SFs)) {
                rv$Tr_SFs
        } else {
                rv$Tr_SFs
        }
},
height = 700)
# ----
# --
