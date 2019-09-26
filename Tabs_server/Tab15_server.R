
# -- output infoText --
output$infoText_15 <- renderText({
        rv$infoText_15
})
# ----


# -- render the overviewViews table --
output$overviewView15 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----


# -- calc taxonomic assignment --
observeEvent(input$calcAssignment, {
        
        rv$infoText_15 <- NULL
        rv$tableAssignment <- NULL
        rv$Tr_Assignment <- NULL
        rv$Tr_ab_prev <- NULL
        
        if (input$tcaType15 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType15 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType15 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType15 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_15 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_15 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered15 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_15 <- "Sorry. You have to do the filtering first."
                        return()
                        
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
                
        }
        
        
        assignment_distribution <- get_assignemnt_distribution(ps)
        
        assignment_distribution <- cbind(Level = rownames(assignment_distribution), assignment_distribution)
        
        rv$tableAssignment <- assignment_distribution
        
        assign_vs_ab <- check_assignment_vs_abundance(ps)
        assign_vs_prev <- check_assignment_vs_prevalence(ps)
        
        rv$Tr_Assignment <- list(assign_vs_prev[[2]], assign_vs_ab[[2]])
        
        rv$infoText_15 <- "Taxonomic assignment niveaus have been calculated."
        
})
# ----


# -- render the assignment table --
output$table_Assignment <- renderTable({
        if(!is.null(rv$tableAssignment)){
                rv$tableAssignment
        } else {
                NULL
        }
}, digits = 2, sanitize.text.function = function(x) x, caption = "Taxonomic assignment levels:", caption.placement = getOption("xtable.caption.placement", "top"))
# ----


# -- output plot of Assignment curves --
output$plot_Assignment <- renderPlot({
        if (is.null(rv$Tr_Assignment)) {
                rv$Tr_Assignment
        } else {
                do.call("grid.arrange", c(rv$Tr_Assignment, nrow = 2))
        }
})
# ----


# -- plot prevalence to non-zero abundance --
observeEvent(input$plot_ab_prev, {
        
        rv$infoText_15 <- NULL
        rv$tableAssignment <- NULL
        rv$Tr_Assignment <- NULL
        rv$Tr_ab_prev <- NULL
        
        if (input$tcaType15 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType15 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType15 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType15 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_15 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_15 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered15 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_15 <- "Sorry. You have to do the filtering first."
                        return()
                        
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
                
        }
        
        
        if (is.null(rv$phylum_colors)) {
                rv$infoText_15 <- "Sorry. You need to calculate phylum_colors first."
                return()
                
        }
        
        phylum_colors <- rv$phylum_colors
        
        TrrList <- plot_correlations_abundance_prev_sparsity(physeq = ps, col = "Phylum", col_vec = phylum_colors)
        
        rv$Tr_ab_prev <- list(TrrList[[3]], TrrList[[4]])
        
        rv$infoText_15 <- "Plots of prevalence to non-zero abundances have been calculated."
        
})
# ----



# -- output plot of Assignment curves --
output$plotabPrev <- renderPlot({
        if (is.null(rv$Tr_ab_prev)) {
                rv$Tr_ab_prev
        } else {
                do.call("grid.arrange", c(rv$Tr_ab_prev, ncol = 2))
        }
})
# ----

# --
