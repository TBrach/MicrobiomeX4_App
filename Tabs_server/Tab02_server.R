# - output infoText -
output$infoText_2 <- renderText({
        rv$infoText_2
})
# --



# - render the overviewViews table -
output$overviewView2 <- renderTable({
        rv$overviewView
}, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
# --



# -- Load phyloseq object --
observeEvent(input$loadPhyseq, {
        
        rv$ps <- NULL
        rv$ps_tca <- NULL
        rv$ps_tca_poscounts <- NULL
        rv$ps_RA <- NULL
        rv$Tr <- NULL
        rv$infoText_2 <- NULL
        rv$phyla <- NULL
        rv$phylum_colors <- NULL
        rv$filteredKeepTaxa <- NULL
        
        inFile <- input$loadPhyseq
        
        if(is.null(inFile)){
                rv$infoText_2 <- "No correct rds file was selected. Please try again."
                return()
        }
        
        if (!grepl(pattern = ".rds$", inFile$datapath)) {
                rv$infoText_2 <- "Sorry. The chosen file was not a .rds file."
                return()
                
        }
        
        
        ps <- readRDS(file = inFile$datapath)
        
        
        if (class(ps) != "phyloseq"){
                rv$infoText_2 <- "The chosen .rds file did not contain a phyloseq object."
                return()
                
        }
        
        
        
        # - added to remove taxa that are not present in a single sample directly -
        # ps <- phyloseq::subset_taxa(ps, taxa_sums(ps) != 0) # not sure why but this gave me an error
        keepTaxa <- taxa_names(ps)[taxa_sums(ps) > 0]
        ps <- phyloseq::prune_taxa(keepTaxa, ps)
        # --
        
        
        
        rv$ps <- ps
        rv$infoText_2 <- paste("Uploaded phyloseq object with ", nsamples(ps), " samples and ", ntaxa(ps), " taxa.", sep = "")
        
})
# ----



# # -- render the explorationViews table --
# output$explorationViews <- renderTable({
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



# -- visualise physeq components --
observeEvent(input$checkPhyseq, {
        
        rv$phyla <- NULL
        rv$infoText_2 <- NULL
        rv$Tr <- NULL
        
        if (is.null(rv$ps)) {
                rv$infoText_2 <- "Sorry. There is no phyloseq object loaded."
                return()
        }
        
        # -- evaluate sample and taxaRanges --
        sampleRange <- input$sampleRange
        if (is.character(sampleRange)){
                sampleRange <- strsplit(x = input$sampleRange, split = ":")
        } else {
                sampleRange <- list(paste0(1, ":", nsamples(rv$ps)))
        }
        sampleRange <- strsplit(x = input$sampleRange, split = ":")
        firstSample <- as.numeric(sapply(sampleRange, `[`, 1))
        lastSample <- as.numeric(sapply(sampleRange, `[`, 2))
        taxaRange <- strsplit(x = input$taxaRange, split = ":")
        firstTaxon <- as.numeric(sapply(taxaRange, `[`, 1))
        lastTaxon <- as.numeric(sapply(taxaRange, `[`, 2))
        
        if (is.na(firstSample) || is.na(lastSample) || !is.numeric(firstSample) || !is.numeric(lastSample) ||
            firstSample < 1 || lastSample > nsamples(rv$ps) || firstSample > lastSample) {
                
                firstSample <- 1
                lastSample <- nsamples(rv$ps)
        }
        
        if (is.na(firstTaxon) || is.na(lastTaxon) || !is.numeric(firstTaxon) || !is.numeric(lastTaxon) ||
            firstTaxon < 1 || lastTaxon > ntaxa(rv$ps) || firstTaxon > lastTaxon) {
                
                firstTaxon <- 1
                lastTaxon <- ntaxa(rv$ps)
        }
        # ----
        
        
        if (input$component == "sample_data"){
                
                # SD <- as(sample_data(rv$ps), "data.frame")
                SD <- get_sample_data(rv$ps)
                
                SD <- SD[firstSample:lastSample, ]
                
                # - make sure POSIX times are shown as character - 
                Classes <- sapply(SD, class)
                PosixIndexes <- grep(pattern = "POSIX", Classes)
                SD[PosixIndexes] <- lapply(SD[PosixIndexes], as.character)
                # --
                
                rv$phyla <- SD
                rv$infoText_2 <- "Sample data has been extracted and will be shown."
                
        } else if (input$component == "otu_table"){
                
                if (taxa_are_rows(rv$ps)){
                        OTU <- round(as(otu_table(rv$ps), "matrix"))
                } else {
                        OTU <- round(t(as(otu_table(rv$ps), "matrix")))
                }
                
                OTU <- as.data.frame(OTU)
                
                OTU <- OTU[firstTaxon:lastTaxon, firstSample:lastSample]
                
                
                rv$phyla <- OTU
                rv$infoText_2 <- "count table has been extracted and will be shown."
                
        } else if (input$component == "tax_table"){
                TT <- as.data.frame(unclass(tax_table(rv$ps)))
                
                TT <- TT[firstTaxon:lastTaxon,]
                rv$phyla <- TT
                rv$infoText_2 <- "taxa table has been extracted and will be shown."
                
        }
        
        
})
# ----



# -- calc phyla distribution and colors --
observeEvent(input$calcPhyla, {
        
        rv$phylum_colors <- NULL
        rv$phyla <- NULL
        rv$infoText_2 <- NULL
        rv$Tr <- NULL
        
        if (is.null(rv$ps)) {
                rv$infoText_2 <- "Sorry. There is no phyloseq object loaded."
                return()
        }
        
        
        if (input$tcaType2 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType2 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType2 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType2 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_2 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_2 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered2 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_2 <- "Sorry. You have to do the filtering first."
                        return()
                        
                }
                
                keepTaxa <- rv$filteredKeepTaxa
                ps <- prune_taxa(keepTaxa, ps)
                
        }
        
        Phyla <- check_phyla_distribution_NA(ps)
        PhylaForColor <- check_phyla_distribution(ps)
        
        if (nrow(PhylaForColor) < 16) {
                phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), QuantColors15)
        } else {
                phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), c(QuantColors15, viridis(nrow(PhylaForColor)-15)))
        }
        
        rv$phylum_colors <- phylum_colors
        Phyla <- dplyr::select(Phyla, Phylum:PC_of_counts, mean_prevalence_in_PC)
        rv$phyla <- Phyla
        rv$infoText_2 <- "Phylum distribution has been calculated and phylum_colors have been determined." 
})
# ----



# -- render the phylaViews table --
output$phylaViews <- renderTable({
        if(!is.null(rv$phyla)){
                phylaShow <- rv$phyla
                phylaShow
        } else {
                NULL
        }
}, rownames = TRUE, digits = 0, sanitize.text.function = function(x) x, caption = "Phyla distribution or phyloseq component", caption.placement = getOption("xtable.caption.placement", "top"))
# ----



# -- generate a barplot of ps --
observeEvent(input$barplotPS, {
        
        rv$Tr <- NULL
        rv$infoText_2 <- NULL
        rv$phyla <- NULL
        
        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                rv$infoText_2 <- "Sorry. Missing phyloseq object or phylum_colors. Please do these steps first."
                return()
                
        }
        
        phylum_colors <- rv$phylum_colors
        
        if (input$tcaType2 == "raw counts"){
                ps <- rv$ps
                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                # rv$ps <- ps
        } else if (input$tcaType2 == "DESeq_gm_exclZero"){
                ps <- rv$ps_tca
        } else if (input$tcaType2 == "DESeq_poscounts") {
                ps <- rv$ps_tca_poscounts
        } else if (input$tcaType2 == "library size (RA)") {
                ps <- rv$ps_RA
        } else {
                rv$infoText_2 <- "Weird how could you choose non of the given options?"
                return()
        }
        
        if (is.null(ps)) {
                rv$infoText_2 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                return()
        }
        
        
        if (input$filtered2 == "filtered taxa") {
                
                if (is.null(rv$filteredKeepTaxa)){
                        rv$infoText_2 <- "Sorry. You have to do the filtering first."
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
                rv$infoText_2 <- error_message
                return()
        }
        # ------
        
        if (input$RA_2){ 
                ps <- phyloseq::transform_sample_counts(ps, function(x){x/sum(x)})
        }
        
        if (!all(group_var_levels %in% compare)){ # all(compare %in% group_var_levels) is tested in check_user_parameters_4
                
                keepSamples <- sample_names(ps)[sample_data(ps)[[group_var]] %in% compare]
                
                ps <- phyloseq::prune_samples(samples = keepSamples, ps)
                
                color_levels <- color_levels[names(color_levels) %in% compare]
        }
        
        
        psP <- phyloseq::tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
        
        
        Tr <- plot_sample_bars(physeq = psP, x = "Sample", y = "Abundance", group_var, color_levels, fill = "Phylum",
                               color_sample_names = TRUE, col_vec = phylum_colors, facet_grid = NULL, order_by_firmicutes = FALSE)
        
        
        rv$Tr <- Tr
        rv$infoText_2 <- "Barplot of samples has been generated." 
})
# ----


# -- output plot samplesPS --
output$samplesPS <- renderPlot({
        if (is.null(rv$Tr)) {
                rv$Tr
        } else {
                rv$Tr
        }
}, 
height = 610)
# ----
# --
