observe({
        rv$group_var <- input$group_var
        group_var_levels <- input$group_var_levels
        group_var_levels <- unlist(strsplit(group_var_levels, split = ","))
        rv$group_var_levels <- trim_whitespace(group_var_levels)
        color_levels <- input$color_levels
        color_levels <- trim_whitespace(unlist(strsplit(color_levels, split = ",")))
        
        # - just so you can use it as a color viewer:) -
        if (length(color_levels) == 0){
                color_levels <- NULL
        }
        rv$color_levels <- color_levels
        # --
        
        if (length(rv$group_var_levels) != length(color_levels)){
                rv$infoText_1 <- "Length group_var_levels did not fit to length color_levels, you need to change this"
                return()
        }
        
        if (length(color_levels) != 0){
                names(color_levels) <- rv$group_var_levels
        }
        rv$color_levels <- color_levels
        
        compare <- input$compare
        if (compare == ""){
                rv$compare <- rv$group_var_levels
        } else {
                compare <- unlist(strsplit(compare, split = ","))
                rv$compare <- trim_whitespace(compare)
        }
        rv$infoText_1 <- "parameters have been adjusted"
})



observe({
        ps <- rv$ps
        ps_tca <- rv$ps_tca 
        ps_tca_poscounts <- rv$ps_tca_poscounts
        ps_RA <- rv$ps_RA
        filteredTaxa <- rv$filteredKeepTaxa
        if (is.null(filteredTaxa)){
                filteredTaxa <- "NULL"
        } else {
                filteredTaxa <- paste0(length(filteredTaxa), " taxa")
        }
        
        Items <- c("ps: raw_counts",
                   "ps: DESeq_gm_exclZero",
                   "ps: DESeq_poscounts",
                   "ps: relative Abundance",
                   "filtered taxa")
        Status <- c(output_ps(ps), output_ps(ps_tca), output_ps(ps_tca_poscounts),
                    output_ps(ps_RA), filteredTaxa)
        
        rv$overviewView <- data.frame(Item = Items, Status = Status)
        
})



# -- render the overviewViews table --
output$overviewView1 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- output infoText --
output$infoText_1 <- renderText({
        rv$infoText_1
})
# ----



# -- visualise the colors --
observeEvent(input$plotColors, {
        
        rv$Tr_colors <- NULL
        
        color_levels <- rv$color_levels
        
        if (is.null(color_levels)){
                rv$infoText_1 <- "Sorry. No colors are given, nothing to plot."
                return()
        }
        
        if (!all(areColors(color_levels))) {
                rv$infoText_1 <- "Sorry. At least one of the given colors is not an R color. Please change color_levels."
                return()
        }
        
        # pal(color_levels)
        
        Tr <- pal_ggplot(color_levels)
        rv$Tr_colors <- Tr
        
        rv$infoText_1 <- "You might see your color levels now."
        
})
# ----



# -- output color plot --
output$plotColors <- renderPlot({
        if (is.null(rv$Tr_colors)) {
                NULL
        } else {
                rv$Tr_colors
        }
},
height = 500)
# ----
