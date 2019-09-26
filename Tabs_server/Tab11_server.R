
# -- output infoText --
output$infoText_11 <- renderText({
        rv$infoText_11
})
# ----


# -- render the overviewViews table --
output$overviewView11 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----


# -- tax_glom --
observeEvent(input$taxGlom, {
        
        rv$infoText_11 <- NULL
        rv$filteredKeepTaxa <- NULL
        
        if (is.null(rv$ps)) {
                rv$infoText_11 <- "Sorry. No phyloseq object is loaded."
                return()
        }
        
        ps <- rv$ps
        
        tax_level <- input$taxLevel
        NArmal <- as.logical(input$NArm)
        
        ps <- phyloseq::tax_glom(ps, taxrank = tax_level, NArm = NArmal)
        
        rv$ps <- ps
        
        # --- adjust size factor corrected objects if SFs have already been calculated ---
        
        if (!is.null(rv$SFs) && !is.null(rv$SFs_poscounts) && !is.null(rv$SFs_RA)) {
                
                library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = rv$SFs) 
                rv$ps_tca <- library_size_adjust_list[[1]]
                
                library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = rv$SFs_poscounts) 
                rv$ps_tca_poscounts <- library_size_adjust_list[[1]]
                
                library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = rv$SFs_RA) 
                rv$ps_RA <- library_size_adjust_list[[1]]
                
                rv$infoText_11 <- paste("Tax_glom has been performed on original object and size factor adjusted objects!", sep = "")
                
        } else {
                
                rv$infoText_11 <- paste("Tax_glom has been performed on original object. No size factors had been calculated yet.", sep = "")
                
        }
        
        # ------
        
})
# ----