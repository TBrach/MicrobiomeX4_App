# - Tab5: filtering -
# -- output infoText --
output$infoText_5 <- renderText({
        rv$infoText_5
})
# ----



# -- render the overviewViews table --
output$overviewView5 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- show prevalence distribution --
observeEvent(input$prevalenceDistribution, {
        
        rv$infoText_5 <- NULL
        
        rv$TrL_filter <- NULL
        
        if (input$tcaType5 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType5 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType5 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType5 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_5 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_5 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        prevalence <- as.numeric(input$filt_prevalence)
        
        if (is.na(prevalence) || prevalence > 100 || prevalence < 0){
                rv$infoText_3 <- "Sorry. Prevalence for filtering must be between 0 and 100. Please change your input."
                return()
                
        }
        
        TrList <- plot_ab_pev_distributions(ps, prevalence = prevalence)
        
        rv$TrL_filter <- list(TrList[[2]], TrList[[4]])
        
        rv$infoText_5 <- "Filter test has been performed."
        
})
# ----



# -- do the filtering --
observeEvent(input$filter, {
        
        rv$infoText_5 <- NULL
        rv$TrL_filter <- NULL
        rv$filteredKeepTaxa <- NULL
        
        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                rv$infoText_5 <- "Sorry. You need a loaded phyloseq object and phylum_colors defined."
                return()
                
        }
        
        if (input$tcaType5 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType5 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType5 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType5 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_5 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_5 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        prevalence <- as.numeric(input$filt_prevalence)
        
        if (is.na(prevalence) || prevalence > 100 || prevalence < 0){
                rv$infoText_5 <- "Sorry. Prevalence for filtering must be between 0 and 100. Please change your input."
                return()
                
        }
        
        taxa_sums_quantile <- as.numeric(input$taxa_sums_quantile)
        
        if (is.na(taxa_sums_quantile) || taxa_sums_quantile > 100 || taxa_sums_quantile < 0){
                rv$infoText_5 <- "Sorry. Prevalence for taxa_sums_quantile must be between 0 and 100. Please change your input."
                return()
                
        }
        
        
        min_obs <- 0L
        
        ps_Filt <- phyloseq::filter_taxa(ps, function(x){
                (sum(x > min_obs) > (prevalence/100)*length(x)) || 
                        (sum(x) > quantile(taxa_sums(ps), probs = taxa_sums_quantile/100))
        }, prune = TRUE)
        
        keepTaxa <- taxa_names(ps_Filt)
        
        rv$filteredKeepTaxa <- keepTaxa
        
        
        rv$TrL_filter <- visualize_filtering(physeq = ps, prevalence = prevalence, taxa_sums_quantile = taxa_sums_quantile, phylum_colors = rv$phylum_colors)
        
        
        rv$infoText_5 <- "ps_filt has been generated. Filtering done." 
})
# ----



# -- output plot filter --
output$plotFilter <- renderPlot({
        if (is.null(rv$TrL_filter)) {
                rv$TrL_filter
        } else {
                grid.arrange(rv$TrL_filter[[1]], rv$TrL_filter[[2]], nrow = 2)
        }
},
height = 800)
# --
