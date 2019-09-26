# --
#######################################
### FUNCTION: test_diffs_in_prevalence_single
#######################################
# Function performs fisher exact test on prevalence (absence presence) between the levels
# in a grouping factor
# INPUT:
# physeq: phyloseq
# group_var: name of the column in sample_data(physeq) that defines the groups
# p.adj.method, used in p.adjust
# minCount: present are taxa in species with more counts than minCount
# OUTPUT:
# list of data.frames, one data frame for each combi of levels in your grouping factor
# The data frames are ordered by p_value, and the tax_table has been cbound:)

test_diffs_in_prevalence_single <- function(physeq, group_var, compare = NULL, p.adj.method = "fdr", minCount = 0L, symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        
        CT <- as(otu_table(physeq), "matrix")
        CT <- CT > minCount
        
        
        prev_list <- lapply(group_var_levels, function(level){
                data.frame(Present = colSums(CT[group_fac == level, ]),
                           Absent = colSums(!CT[group_fac == level, ]))
        })
        
        pr_ab_gr1 <- prev_list[[1]] #pr_ab = presence absence
        pr_ab_gr2 <- prev_list[[2]]
        
        rowwise_compare_matrix <- cbind(pr_ab_gr1, pr_ab_gr2)
        
        FisherTests <- lapply(1:nrow(rowwise_compare_matrix), function(e){
                mat_fisher <- matrix(c(rowwise_compare_matrix[e, 1],
                                       rowwise_compare_matrix[e, 3],
                                       rowwise_compare_matrix[e, 2],
                                       rowwise_compare_matrix[e, 4]), ncol = 2)
                fisher.test(mat_fisher, conf.int = TRUE, conf.level = 0.95)
                # fisher.test(mat_fisher, conf.int = FALSE)
        })
        
        # - take out oddsRatios -
        oddRatios <- sapply(FisherTests, function(TestResult){TestResult$estimate})
        oddRatios_lb <- sapply(FisherTests, function(TestResult){TestResult$conf.int[1]})
        oddRatios_ub <- sapply(FisherTests, function(TestResult){TestResult$conf.int[2]})
        
        direction <- rep(group_var_levels[1], length(oddRatios))
        direction[oddRatios < 1] <- group_var_levels[2]
        
        oddRatios_lb[oddRatios < 1] <- 1/oddRatios_lb[oddRatios < 1]
        oddRatios_ub[oddRatios < 1] <- 1/oddRatios_ub[oddRatios < 1]
        oddRatios_lb_final <- pmin(oddRatios_lb, oddRatios_ub)
        oddRatios_ub_final <- pmax(oddRatios_lb, oddRatios_ub)
        oddRatios[oddRatios < 1] <- 1/oddRatios[oddRatios < 1]
        # --
        
        
        # - take out p_vals and assign significance levels -
        p_vals <- sapply(FisherTests, function(TestResult){TestResult$p.value})
        p_vals_adj <- p.adjust(p_vals, p.adj.method)
        
        symnum.args$x <- p_vals
        significance <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- p_vals_adj
        significance_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
        # --
        
        # - take out prevalence percentages -
        prev_PC_gr1 <- round(100*(pr_ab_gr1[, "Present"]/sum(group_fac == group_var_levels[1])), 1)
        prev_PC_gr2 <- round(100*(pr_ab_gr2[, "Present"]/sum(group_fac == group_var_levels[2])), 1)
        # --
        
        df <- data.frame(p_val = p_vals, p_val_adj = p_vals_adj, signi = significance, signi_adj = significance_adj, oddsRatio = round(oddRatios, 2), oddsRatio_lb = round(oddRatios_lb_final, 2),
                         oddsRatio_ub = round(oddRatios_ub_final, 2), direction = direction, comparison = paste(group_var_levels, collapse = " vs "),  prev_PC_gr1 = prev_PC_gr1,  prev_PC_gr2 = prev_PC_gr2)
        
        df <- cbind(as.data.frame(df), as.data.frame(unclass(tax_table(physeq))))
        df$Taxon <- colnames(CT)
        df <- arrange(df, p_val) %>% select(Taxon, 1:(ncol(df)-1))
        df
        
}
# --



# --
#######################################
### FUNCTION: get_taxon_names
#######################################


get_taxon_names <- function(df) {
        
        df1 <- df[, colnames(df) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", 
                                        "Genus", "Species")]
        
        if (ncol(df1) == 0) {stop("the provided data frame did not contain the expected taxonomy columns such as Phylum, Class etc.")}
        
        df1[] <- lapply(df1, as.character)
        
        Last_NotNA_Position <- apply(df1, 1, function(x){length(which(!is.na(x)))})
        Last_NotNA_Position[Last_NotNA_Position == 0] <- 1
        
        Names <- vector(mode = "character", length = nrow(df1))
        
        for (i in 1:nrow(df1)) {
                if (Last_NotNA_Position[i] == 7){
                        Names[i] <- paste(df1[i, "Genus"], df1[i, "Species"], sep = " ")
                } else {
                        Names[i] <- df1[i, Last_NotNA_Position[i]] 
                }
        }
        
        Names[is.na(Names)] <- "NA"
        Names
        
}
# --




# --
#######################################
### FUNCTION: get_taxon_names_plusTL
#######################################


get_taxon_names_plusTL <- function(df) {
        
        df1 <- df[, colnames(df) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", 
                                        "Genus", "Species")]
        
        if (ncol(df1) == 0) {stop("the provided data frame did not contain the expected taxonomy columns such as Phylum, Class etc.")}
        
        df1[] <- lapply(df1, as.character)
        
        Last_NotNA_Position <- apply(df1, 1, function(x){length(which(!is.na(x)))})
        Last_NotNA_Position[Last_NotNA_Position == 0] <- 1
        
        Names <- vector(mode = "character", length = nrow(df1))
        
        Prefixes <- c("k", "p", "c", "o", "f", "g", "s")
        
        for (i in 1:nrow(df1)) {
                if (Last_NotNA_Position[i] == 7){
                        Names[i] <- paste0("gs_", df1[i, "Genus"], " ", df1[i, "Species"])
                } else {
                        Names[i] <- paste0(Prefixes[Last_NotNA_Position[i]], "_", df1[i, Last_NotNA_Position[i]])
                }
        }
        
        Names[is.na(Names)] <- "NA"
        Names
        
}
# --






# --
#######################################
### FUNCTION: format_hit_table
#######################################

format_hit_table <- function (result_df, p.adjust.threshold = 0.1, p.adjust.method = NULL) {
        
        if (!all(c("p_val", "p_val_adj", "Taxon", "direction", "signi", "signi_adj") %in% colnames(result_df))) {
                stop("result_df should be a data.frame generated by one of the differential abundance tests and contain all corresponding columns.")
        }

        if (!is.null(p.adjust.method)) {
                result_df$p_val_adj = p.adjust(result_df$p_val, method = p.adjust.method)
        }
        
        
        result_df <- arrange(result_df, p_val_adj, p_val)
        
        no_hits <- sum(result_df$p_val_adj <= p.adjust.threshold, na.rm = TRUE)
        
        keepTaxa <- no_hits
        
        if (keepTaxa < 10 && nrow(result_df) >= 10) {
                keepTaxa <- 10
        } else if (keepTaxa < 10 && nrow(result_df < 10)){
                keepTaxa <- nrow(result_df)
        }
        
        
        df <- result_df[1:keepTaxa,]
        
        taxa_annotation <- get_taxon_names_plusTL(df)
        taxa_annotation <- strsplit(taxa_annotation, split = "/")
        taxa_annotation <- sapply(taxa_annotation, `[`, 1)
        df$Annotation <- taxa_annotation
        
        df <- select(df, Taxon, Annotation, p_val, p_val_adj, signi, signi_adj, direction, colnames(df)[!(colnames(df) %in% c("Taxon", "Annotation", "p_val", "p_val_adj", "signi", "signi_adj", "direction"))])
        
        rownames(df) <- df$Taxon
        
        list(hit_table = df, no_hits = no_hits)
}
# --





# --
#######################################
### FUNCTION: format_summary_hit_table
#######################################

format_summary_hit_table <- function (taxa_all_DESeq, taxa_all_Wilcox, taxa_all_fisher, p.adjust.threshold = 0.1,
                                      symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        # - combine the result tables -
        da_summary <- dplyr::select(taxa_all_DESeq, Taxon, p_val_adj_DESeq = p_val_adj, direction_DESeq = direction, log2FoldChange, Mean_grp1:Species)
        da_summary$p_val_adj_Wilcox = taxa_all_Wilcox$p_val_adj[match(da_summary$Taxon, taxa_all_Wilcox$Taxon)]
        da_summary$p_val_adj_Fisher = taxa_all_fisher$p_val_adj[match(da_summary$Taxon, taxa_all_fisher$Taxon)]
        da_summary <- rowwise(da_summary) %>% dplyr::mutate(p_val_adj_median = median(c(p_val_adj_DESeq, p_val_adj_Wilcox, p_val_adj_Fisher), na.rm =TRUE),
                                                            p_val_adj_minimum = min(c(p_val_adj_DESeq, p_val_adj_Wilcox, p_val_adj_Fisher), na.rm = TRUE))
        da_summary <- dplyr::arrange(da_summary, p_val_adj_median, p_val_adj_minimum)
        # -- 
        
        # - catch the hits -
        result_df <- dplyr::filter(da_summary, p_val_adj_DESeq <= p.adjust.threshold |  p_val_adj_Wilcox <= p.adjust.threshold | p_val_adj_Fisher <= p.adjust.threshold)
        
        if (nrow(result_df) > 0){
                taxa_annotation <- get_taxon_names_plusTL(result_df)
                taxa_annotation <- strsplit(taxa_annotation, split = "/")
                taxa_annotation <- sapply(taxa_annotation, `[`, 1)       
        } else {
                taxa_annotation <- character(0)
        }
        result_df$Annotation <- taxa_annotation
        
        symnum.args$x <- result_df$p_val_adj_minimum
        result_df$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        
        result_df <- dplyr::select(result_df, Taxon, Annotation, p_val_adj_median, p_val_adj_minimum, p_val_adj_DESeq, p_val_adj_Fisher, p_val_adj_Wilcox,  signi, direction_DESeq:Species)
        result_df <- data.frame(result_df)
        rownames(result_df) <- result_df$Taxon # used for plots later
        result_df
}
# --




# --
#########################
## plot_heatmap_physeq ##
#########################
# function that uses pheatmap to draw a heatmap of the physeq object, or a pruned version of it determined by taxa_info_df
# INPUTs
# physeq = physeq object
# sample_colors = named list of named color character vectors. The sample_data(physeq) data.frame is used as sample_info_df to color the samples (columns) in the heatmap.
    # Only those variables/columns in sample_info_df that are names in sample_colors are considered for coloring the samples. I.e. if sample_colors = NULL no coloring occurs.
    # NB: if a variable is a name in sample_colors but the color_character is NA, default colors will be used.
    # NB!!: The named character vectors in sample_colors are also used to prune the samples!! This works as follows:
    # The character vectors in sample_colors are checked in order. If the first character vector (sample_colors[[1]]) is for example "Group" and it contains only colors for "A" and "B", 
    # but sample_info_df[["Group"]] contains also "C", the C samples will be removed! The procedure continues for all entries in sample_colors
    # finally the order of names(sample_colors) defines the order of the color bars, the first entry being drawn closest to the heatmap and so on. 
# taxa_info_df = a data frame similar to sample_info_df (see sample_colors) that together with taxa_colors will be used to label the taxa.
    # NB: MORE importantly, it is also used to restrict the included taxa and to order the taxa!! The rownames(taxa_info_df) must be the taxa_names in physeq 
    # of the taxa you want to include in the heatmap in the order you want the taxa to be shown! IF NULL all taxa will be plotted in the order given by physeq
# taxa_colors: see sample_colors
# taxa_annotation: if not NULL will be used as labels_row in pheatmap, i.e. a different labelling of the taxa instead of taxa_names(physeq) or rownames(taxa_info_df)
# max_abundance_for_color: the maximum abundance value to which the maximum color is assigned, all abundances above get also this color, if NULL max(otu_table(physeq))
    # is used
# gradient_steps:  gradient_steps: vector of numbers between 0 and 1, default c(0.15, 0.3, 0.45, 1). Defines how the main colors in the gradient will be positioned between min_in_data and
    # max_abundance_for_color. min_in_data (= lowest non-zero count/abundance) will be added to first position after gradient_steps have been normalized with max_abundance_for_color
# zero_color: defines the color of the 0 values (minimum Value in otu_table(physeq).
# color_function: a function such as the default viridis function that generates color strings when given a number.
# color_steps_bw_markers: numeric value defining how many equally distributed breaks will be introduced between the markers defined by gradient_steps.
# log_transform: option to log10 transform the counts
# drop_color_levels: if TRUE, color levels in sample_colors and taxa_colors are removed (not shown in legend so) when not present in the data
# rest pheatmap arguments

plot_heatmap_physeq <- function (physeq, sample_colors = NULL, taxa_info_df = NULL, taxa_colors = NULL, taxa_annotation = NULL, 
                                 max_abundance_for_color = NULL, gradient_steps = c(0.15, 0.3, 0.45, 1), zero_color = "white",
                                 color_function = viridis,
                                 color_steps_bw_markers = 10, log_transform = FALSE, drop_color_levels = FALSE,
                                 border_color = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = TRUE, annotation_names_col = TRUE, annotation_legend = TRUE,
                                 legend = TRUE, font_size = 10, fontsize_row = 8, fontsize_col = 8, fontsize_number = 6, filename = NA, ...) {
        
        # - keep only taxa in taxa_info_df in physeq -
        if (!is.null(taxa_info_df)) {
                
                if (!all(rownames(taxa_info_df) %in% taxa_names(physeq))) {
                        stop("not all taxa in taxa_info_df are taxa in physeq!")
                }
                pruned_ps <- prune_taxa(rownames(taxa_info_df), physeq)
        } else {
                pruned_ps <- physeq
                taxa_colors <- NULL
        }
        # --
        
        # - test sample_colors input and prepare sample_info_df -
        if (!is.null(sample_colors)) {
                
                if (!is.list(sample_colors)) {
                        stop("sample_colors must be a named list.")
                }
                
                
                if (!all(names(sample_colors) %in% colnames(sample_data(pruned_ps)))) {
                        stop("all names of the sample_colors list must be variables/columns in sample_data(physeq).")
                }
                
                sample_info_df <- as(sample_data(pruned_ps), "data.frame")
                sample_info_df <- sample_info_df[, match(names(sample_colors), colnames(sample_info_df)), drop = FALSE] # the match construct makes sure that sample_info_df columns are ordered after names(sample_colors)
                # which makes sure that the column coloring occurs in the whished for order
                
                # -- test all entries in sample_colors, prune pruned_ps and order sample_info_df, and set default color where necessary --
                for (i in 1:length(sample_colors)){
                        variable_name <- names(sample_colors)[i]
                        color_character <- sample_colors[[i]]
                        
                        if (!all(is.na(color_character))){
                                
                                # to not remove all samples from pruned_ps:
                                if (!(sum(names(color_character) %in% unique(sample_info_df[[variable_name]])) > 0)) {
                                        stop(paste("Not a single name of a color for ", variable_name, " matched to an entry in the corresponding column in sample_data(physeq)", sep = ""))
                                }
                                
                                # remove all samples that no color was assigned to in sample_colors, i.e.option to restrict the comparison to defined samples
                                if (!all(unique(sample_info_df[[variable_name]]) %in% names(color_character))) {
                                        
                                        keepSamples <- sample_names(pruned_ps)[sample_info_df[[variable_name]] %in% names(color_character)]
                                        pruned_ps <- prune_samples(keepSamples, pruned_ps)
                                        sample_info_df <- as(sample_data(pruned_ps), "data.frame")
                                        sample_info_df <- sample_info_df[, colnames(sample_info_df) %in% names(sample_colors), drop = FALSE]
                                }
                                
                                # remove unused color levels if asked for it
                                if (drop_color_levels && !all(names(color_character) %in% unique(sample_info_df[[variable_name]]))) {
                                        color_character <- color_character[names(color_character) %in% unique(sample_info_df[[variable_name]])]
                                        sample_colors[[i]] <- color_character
                                }
                                
                                if (!all(areColors(color_character))) {
                                        warning(paste("Not all entries in sample_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(sample_info_df, variable_name)
                                        sample_colors[[i]] <- color_character
                                }
                                
                                # define the order of the samples based on the sample_color vectors (actual ordering below with dplyr::arrange)
                                sample_info_df[[variable_name]] <- factor(sample_info_df[[variable_name]], levels = names(color_character), ordered = TRUE)
                                
                        } else {
                                color_character <- assign_default_colors(sample_info_df, variable_name)
                                sample_colors[[i]] <- color_character
                        }
                }
                # --- order the samples based on the sample_colors ---
                # NB: when you give several sample_color variables the first one given will be the lowest bar in the heatmap, i.e. the bar directly on top of the columns
                # You want to order the samples thus first by the last variable given, then break ties by the next and so on
                # with help from the following two stackoverflow I finally came with a solution that still appears a bit clumsy to me
                # <https://stackoverflow.com/questions/31681894/ordering-data-frame-using-variable-as-column-name>
                # https://stackoverflow.com/questions/1743698/evaluate-expression-given-as-a-string
                
                orderCommand <- paste0("sample_info_df[,", ncol(sample_info_df):1, "]")
                orderCommand <- paste0(orderCommand, collapse = ", ")
                orderCommand <- paste0("order(", orderCommand, ")")
                orderCommand <- parse(text = orderCommand)
                
                # --- you had to add here the rownames (or whatever) again to avoid having a 1 column data frame which would be changed into a factor
                sample_info_df$ID_Saver <- rownames(sample_info_df)
                # ------
                
                sample_info_df <- sample_info_df[eval(orderCommand), ]
                
                sample_info_df <- sample_info_df[-ncol(sample_info_df)]
                # ----
        } else {
                sample_info_df <- NULL
        }
        # --
        
        # - set the taxa colors (and pool sample_colors and taxa_colors into annotation_colors) -
        if (!is.null(taxa_colors)) { # see above, guarantees also that taxa_info_id is not NULL
                
                if (!is.list(taxa_colors)) {
                        stop("taxa_colors must be a named list.")
                }
                
                if (!all(names(taxa_colors) %in% colnames(taxa_info_df))) {
                        stop("all names of the taxa_colors list must be variables/columns in taxa_info_df.")
                }
                
                taxa_info_df <- taxa_info_df[, match(names(taxa_colors), colnames(taxa_info_df)), drop = FALSE] # again match to get coloring ordered after taxa_colors
                
                # -- test all entries in taxa_colors --
                for (i in 1:length(taxa_colors)){
                        variable_name <- names(taxa_colors)[i]
                        color_character <- taxa_colors[[i]]
                        
                        if (!all(is.na(color_character))){
                                
                                
                                if (!all(unique(taxa_info_df[[variable_name]]) %in% names(color_character))) {
                                        warning(paste("There were levels in taxa_info_df at ", variable_name, " for which no color was assigned in taxa_colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(taxa_info_df, variable_name)
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                # remove unused color levels if asked for it
                                if (drop_color_levels && !all(names(color_character) %in% unique(taxa_info_df[[variable_name]]))) {
                                        color_character <- color_character[names(color_character) %in% unique(taxa_info_df[[variable_name]])]
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                
                                if (!all(areColors(color_character))) {
                                        warning(paste("Not all entries in taxa_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(taxa_info_df, variable_name)
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                # define the order of the taxa coloring based on the taxa_color vectors
                                # taxa_info_df[[variable_name]] <- factor(taxa_info_df[[variable_name]], levels = names(color_character), ordered = TRUE) # was unnecessary, pheatmap orders based on annotation_colors
                                
                        } else {
                                color_character <- assign_default_colors(taxa_info_df, variable_name)
                                taxa_colors[[i]] <- color_character
                        }
                }
                # ----
        } 
        
        annotation_colors <- c(sample_colors, taxa_colors)
        # -- 
        
        # - check or set taxa_annotation -
        if (is.null(taxa_annotation)){
                taxa_annotation <- taxa_names(pruned_ps) # or you could use 
                # taxa_annotation <- get_taxon_names(as.data.frame(tax_table(pruned_ps)))
        }
        
        if (length(taxa_annotation) != ntaxa(pruned_ps)) {
                warning("taxa_annotation did not fit in length to nrow(taxa_info_df) or ntaxa(physeq), used taxa_names")
                taxa_annotation <- taxa_names(pruned_ps)
        }
        
        taxa_annotation <- as.character(taxa_annotation)
        
        taxa_annotation[is.na(taxa_annotation)] <- "NA"
        
        taxa_annotation <- make.unique(taxa_annotation)
        # --
        
        # - generate count data frame in which taxa are rows in the order determined by taxa_info_df -
        if (taxa_are_rows(pruned_ps)) {
                pruned_ps <- t(pruned_ps)
        }
        
        DF_CT <- as(otu_table(pruned_ps), "matrix")
        DF_CT <- as.data.frame(t(DF_CT))
        
        if (!is.null(taxa_info_df)){
                DF_CT <- DF_CT[rownames(taxa_info_df), ]
        }
        # --
        
        # - order the samples in DF_CT based on sample_colors list (see above for sample_info_df) -
        if (!is.null(sample_info_df)){
                DF_CT <- DF_CT[, rownames(sample_info_df)]
        }
        # --
        
        # - test and adjust max_abundance_for_color -
        if (is.null(max_abundance_for_color)) {
                max_abundance_for_color <- max(DF_CT)
        }
        if (max_abundance_for_color < min(DF_CT) && max_abundance_for_color > max(DF_CT)) {
                max_abundance_for_color <- max(DF_CT)
        }
        # --
        
        # - set breaks and colors for pheatmap and do a possible log transform -
        ZeroValue <- min(DF_CT) # should be 0 in most cases of microbiome data, unless pseudocount had been added. In rare cases of high taxonomy it might be a non-zero count so actually min_in_data 
        min_in_data <- min(DF_CT[DF_CT > ZeroValue]) # the lowest non-zero value
        max_in_data <- max(DF_CT)
        
        # -- make sure the last entry in gradient steps = 1 --
        if (!all(gradient_steps >= 0 & gradient_steps <= 1)) {
                gradient_steps <- c(0.15, 0.3, 0.45, 1)
        }
        
        if (gradient_steps[length(gradient_steps)] != 1) {
                gradient_steps <- c(gradient_steps, 1)
        }
        # ----
        
        # -- added this even though it might be confusing to the user, but otherwise you get an error --
        if (max_abundance_for_color < min_in_data) {
                max_abundance_for_color <- min_in_data # NB: in this case all non-zero values should become maximal color
        }
        # -----
        
        
        # you want that the final color gradient covers the values min to max_abundance_for_color
        # 0 values will be set to a different color (usually red or white), values above max_abundance_for_color should be all max_color
        # normalise gradient steps with max_abundance_for_color
        myBreaks <- gradient_steps * max_abundance_for_color
        myBreaks <- myBreaks[myBreaks >= min_in_data]
        # add break at min_in_data
        if (length(myBreaks) != 1 || myBreaks != max_abundance_for_color){ # in case max_abundance_for_color == min_in_data
                myBreaks <- c(min_in_data, myBreaks)
        }
        
        if (length(myBreaks) > 1){
                
                # now myBreaks goes from min_in_data up to max_abundance_for_color (provided the last gradient_steps was 1)
                myColors = color_function(length(myBreaks)) # see help pheatmap, breaks should be 1 element longer than color, now it is same legnth
                
                # myColors contains now the colors that represent the gradient_steps values (= markers).
                # now we want to introduce breaks between these markers and make linear color gradients between the markers. The number of breaks between the markers is defined
                # by color_steps_bw_markers
                
                myColors <- unlist(lapply(1:(length(myColors)-1), function(i) {
                        colorRampPalette(colors = c(myColors[i], myColors[i+1]))(color_steps_bw_markers)
                }))
                
                
                # -- do a possible log transform using min_in_data/5 as pseudocounts --
                if (log_transform) {
                        if (ZeroValue > 0){ # > 0 if pseudocount had been added or simply no zero counts in data, < 0 if data h
                                pseudocount <- ZeroValue
                        } else if (ZeroValue < 0){
                                stop("you asked for log_transform but the lowest count in your data is already below 0!")
                        } else {
                                pseudocount <- min_in_data/2
                        }
                        DF_CT[DF_CT == ZeroValue] <- pseudocount
                        DF_CT <- log10(DF_CT)
                        
                        myBreaks <- log10(myBreaks)
                        myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
                                seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
                        })
                        myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)]) #length(myBreaks) is now length(myBreaks) * color_steps_bw_markers + 1, the markers are at positions 1, 1+color_steps_bw_markers, 1+2*color_steps_bw_markers 
                        
                        if (log10(max_in_data) > myBreaks[length(myBreaks)]) {
                                myBreaks <- c(log10(pseudocount), myBreaks, log10(max_in_data))
                                myColors <- c(zero_color, myColors, myColors[length(myColors)])
                        } else {
                                myBreaks <- c(log10(pseudocount), myBreaks)
                                myColors <- c(zero_color, myColors)
                        }
                        
                } else {
                        myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
                                seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
                        })
                        myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)])
                        
                        # add max_in_data values to myBreaks, otherwise all values above max_abundance_for_color will be white
                        if (max_in_data > myBreaks[length(myBreaks)]) {
                                myBreaks <- c(ZeroValue, myBreaks, max_in_data)
                                myColors <- c(zero_color, myColors, myColors[length(myColors)])
                        } else {
                                myBreaks <- c(ZeroValue, myBreaks)
                                myColors <- c(zero_color, myColors)
                        }
                }
                
        } else {
                
                myColors = color_function(2)
                myColors <- myColors[2]
                
                if (log_transform) {
                        if (ZeroValue > 0){ # > 0 if pseudocount had been added or simply no zero counts in data, < 0 if data h
                                pseudocount <- ZeroValue
                        } else if (ZeroValue < 0){
                                stop("you asked for log_transform but the lowest count in your data is already below 0!")
                        } else {
                                pseudocount <- min_in_data/2
                        }
                        DF_CT[DF_CT == ZeroValue] <- pseudocount
                        DF_CT <- log10(DF_CT)
                        
                        myBreaks <- log10(myBreaks)
                        
                        
                        if (log10(max_in_data) > myBreaks[length(myBreaks)]) {
                                myBreaks <- c(log10(pseudocount), myBreaks, log10(max_in_data))
                                myColors <- c(zero_color, myColors, myColors[length(myColors)])
                        } else {
                                myBreaks <- c(log10(pseudocount), myBreaks)
                                myColors <- c(zero_color, myColors)
                        }
                } else {
                        
                        if (max_in_data > myBreaks[length(myBreaks)]) {
                                myBreaks <- c(ZeroValue, myBreaks, max_in_data)
                                myColors <- c(zero_color, myColors, myColors[length(myColors)])
                        } else {
                                myBreaks <- c(ZeroValue, myBreaks)
                                myColors <- c(zero_color, myColors)
                        }
                        
                }
                
        }
        
        
        # ----
        # --
        
        # - not necessary, just for clarity -
        if (is.null(annotation_colors)) {annotation_colors <- NA}
        if (is.null(sample_info_df)) {sample_info_df <- NA}
        if (is.null(taxa_info_df)) {taxa_info_df <- NA}
        # --
        
        # --
        hm.parameters <- list(DF_CT, color = myColors, breaks = myBreaks, border_color = border_color, cluster_cols = cluster_cols, cluster_rows = cluster_rows,
                              show_rownames = show_rownames, show_colnames = show_colnames, annotation_col = sample_info_df,
                              annotation_row = taxa_info_df, annotation_colors = annotation_colors, labels_row = taxa_annotation, annotation_names_row = annotation_names_row,
                              annotation_names_col = annotation_names_col, annotation_legend = annotation_legend, legend = legend, font_size = font_size, 
                              fontsize_row = fontsize_row, fontsize_col = fontsize_col, fontsize_number = fontsize_number)
        do.call("pheatmap", hm.parameters)
        
        
}
# --





# --
#######################################
### test_differential_abundance_DESeq2single
#######################################
## Inputs
# physeq: phyloseq object
# group_var: name of column that defines group fac in sample_data
# SFs: often you might want to give the SizeFactors already because you wanted to calculate them on non-filtered data,
# when SFs are not NULL, type is ignored
# type: type in estimateSizeFactors, ignored when Size factors given
## OUTPUT:
# list, first: result df of DESEQ2 analysis, second: the adjusted phyloseq object after size factor correction


# ATTENTION: you could add here block as Mani had in test_differential_abundance_DESeq2

test_differential_abundance_DESeq2single <- function(physeq, group_var, compare = NULL, cooksCutoff = TRUE, SFs = NULL, type = "ratio", p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))){
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        
        # NB: DESeq2 can not deal with ordered factors, sees them somehow as just one level, therefore
        sample_data(physeq)[[group_var]] <- factor(group_fac, levels = c(group_var_levels, setdiff(levels(group_fac), group_var_levels)), ordered = FALSE)
        
        DES = phyloseq::phyloseq_to_deseq2(physeq, formula(paste("~", group_var)))
        
        
        if (is.null(SFs)){
                if(type == "ratio"){
                        GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = FALSE)
                }
                
                dds <- estimateSizeFactors(DES, type = type, geoMeans = GM)
                # NB: geoMeans is ignored when type = "iterate"
                # NB2: "iterate" takes much longer than "ratio", and when using GM via gm_own with "ratio" the size factors 
                # correlate with above 99% with size factors from "iterate"
                # NB3: "poscounts" is the same as "ratio" with GM calculated with zeros.count = TRUE!!
                # SFs2 <- sizeFactors(dds)
                
        } else {
                dds <- DES
                sizeFactors(dds) = SFs
                # identical(sizeFactors(dds), SFs) 
        }
        
        
        dds <-  DESeq(dds, fitType = "parametric", test = "Wald", 
                      quiet = TRUE, minReplicatesForReplace = Inf) 
        
        # to get the size factor adjusted physeq object
        physeq_out <- physeq
        otu_table(physeq_out) <- otu_table(t(counts(dds, normalized = TRUE)), taxa_are_rows = FALSE)
        
        
        
        # - analyse the res -
        res <- as.data.frame(results(dds, contrast = c(group_var, group_var_levels), cooksCutoff = cooksCutoff))
        
        res$p_val_adj <- p.adjust(res$pvalue, method = p.adjust.method) # NB: in case of "fdr" same as default DESeq2
        
        CT <- counts(dds, normalized = TRUE)
        i = 1
        j = 2
        n1 <- sum(group_fac == group_var_levels[i])
        n2 <- sum(group_fac == group_var_levels[j])
        # res$Median_grp1 <- apply(CT[, group_fac == group_var_levels[i]], 1, median)
        # res$Median_grp2 <- apply(CT[, group_fac == group_var_levels[j]], 1, median)
        res$Mean_grp1 <- apply(CT[, group_fac == group_var_levels[i]], 1, mean)
        res$Mean_grp2 <- apply(CT[, group_fac == group_var_levels[j]], 1, mean)
        # res$baseMeanSelf <- apply(CT, 1, mean) # exactly the same as baseMean!
        # res$Zeros_grp1 <- apply(CT[, group_fac == group_var_levels[i]], 1, function(cnts){sum(cnts == 0)})
        # res$Zeros_grp2 <- apply(CT[, group_fac == group_var_levels[j]], 1, function(cnts){sum(cnts == 0)})
        res$Present_grp1 <- apply(CT[, group_fac == group_var_levels[i]], 1, function(cnts){sum(cnts != 0)})
        res$Present_grp2 <- apply(CT[, group_fac == group_var_levels[j]], 1, function(cnts){sum(cnts != 0)})
        res$prev_PC_grp1 <- round(100*(res$Present_grp1/n1),1)
        res$prev_PC_grp2 <- round(100*(res$Present_grp2/n2), 1)
        #res$Sparsity_grp1 <- 100*(res$Zeros_grp1/n1)
        #res$Sparsity_grp2 <- 100*(res$Zeros_grp2/n2)
        
        # - add sginificance and direction -
        symnum.args$x <- res$pvalue
        res$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- res$p_val_adj
        res$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character() # note ? in case of p_val = NA
        
        res$direction <- rep(group_var_levels[2], nrow(res))
        res$direction[res$log2FoldChange > 0] <- group_var_levels[1]
        
        res$comparison <- paste(group_var_levels, collapse = " vs ")
        # --
        
        res$Taxon <- rownames(res)
        
        res <- dplyr::select(res, Taxon, teststat = stat, p_val = pvalue, p_val_adj,
                             signi, signi_adj, direction, comparison,
                             baseMean, log2FoldChange, Mean_grp1,
                             Mean_grp2, prev_PC_grp1, prev_PC_grp2)
        
        # NB: I dropped here padj from DESeq since same as p_val_adj in case of p.adjust.method = "fdr"
        res <- cbind(res, as.data.frame(unclass(tax_table(physeq))))
        res <- dplyr::arrange(res, desc(abs(teststat)))
        list(result_df = res, physeq_out = physeq_out) 
        
}
# --





# --
#######################################
### plot_hittaxa_boxAndviolin ##
#################


plot_hittaxa_boxAndviolin <- function(physeq, group_var, color_levels, taxa_info_df = NULL, taxa_annotation = NULL, 
                                      facet_cols = 5, shape = NULL, excludeZero = FALSE, logTransform = FALSE){
        
        # - make sure group_var and shape are in the sample_data -
        if (!group_var %in% colnames(sample_data(physeq))) {
                stop("group_var must be a variable in sample_data(physeq)")
        }
        
        if (!is.null(shape) && !shape %in% colnames(sample_data(physeq))) {
                shape = NULL
        }
        # --
        
        
        # - show only samples "represented" in color levels -
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
        }
        
        
        # make sure the order of the group levels is as defined by color_levels
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
        # --
        
        # - keep only taxa in taxa_info_df in physeq -
        if (!is.null(taxa_info_df)) {
                
                if (!all(rownames(taxa_info_df) %in% taxa_names(physeq))) {
                        stop("not all taxa in taxa_info_df are taxa in physeq!")
                }
                pruned_ps <- prune_taxa(rownames(taxa_info_df), physeq)
        } else {
                pruned_ps <- physeq
        }
        # --
        
        # - check or set taxa_annotation -
        if (is.null(taxa_annotation)){
                if (!is.null(taxa_info_df)) {
                        taxa_annotation <- paste(taxa_info_df$Taxon, "_", taxa_info_df$p_val_adj, sep = "")
                } else {
                        taxa_annotation <- taxa_names(pruned_ps) # or you could use 
                        # taxa_annotation <- get_taxon_names(as.data.frame(tax_table(pruned_ps)))
                }
        }
        
        if (length(taxa_annotation) != ntaxa(pruned_ps)) {
                warning("taxa_annotation did not fit in length to nrow(taxa_info_df) or ntaxa(physeq), used taxa_names")
                taxa_annotation <- taxa_names(pruned_ps)
        }
        
        taxa_annotation <- as.character(taxa_annotation)
        
        taxa_annotation[is.na(taxa_annotation)] <- "NA"
        
        taxa_annotation <- make.unique(taxa_annotation)
        
        if (!is.null(taxa_info_df)) {
                taxa_lookup <- data.frame(Taxon = taxa_info_df$Taxon, Annotation = taxa_annotation)
        } else {
                taxa_lookup <- data.frame(Taxon = taxa_names(pruned_ps), Annotation = taxa_annotation)
        }
        # --

        # - generate data frame for plotting, add annotation, and make sure order of annotation fits -
        DF_CT <- psmelt(pruned_ps)
        
        DF_CT$Annotation <- taxa_lookup$Annotation[match(DF_CT$OTU, taxa_lookup$Taxon)]
        DF_CT$Annotation <- factor(DF_CT$Annotation, levels = taxa_annotation, ordered = TRUE) #
        # --
        
        # - exclude zeros and log transform if asked for -
        if (excludeZero) {
                DF_CT <- filter(DF_CT, Abundance != 0)
        }
        
        if (logTransform) {
                
                minValue <- min(DF_CT$Abundance)
                
                if (minValue < 0){
                        stop("you asked for log_transform but the lowest count in your data is already below 0!")
                }
                if (minValue == 0) {
                        pseudocount <- min(DF_CT$Abundance[DF_CT$Abundance > 0])/2
                        DF_CT$Abundance[DF_CT$Abundance == 0] <- pseudocount
                }
                DF_CT$Abundance <- log10(DF_CT$Abundance)
                
        }
        # --
        
        
        
                
        Tr <- ggplot(DF_CT, aes_string(x = group_var, y = "Abundance", col = group_var))
        Tr <- Tr +
                geom_boxplot(outlier.color = NA) +
                geom_point(aes_string(shape = shape), size = 1, alpha = 0.6, position = position_jitterdodge(jitter.width = .8)) +
                facet_wrap(~ Annotation, ncol = facet_cols, scales = "free_y") +
                scale_color_manual("", values = color_levels) +
                xlab("") +
                ylab("abundance") +
                theme_bw()
        
        if (is.null(shape)) {
                Tr <- Tr + theme(legend.position = "none")
        }
        
        # if (ttestp == "yes"){
        #         Tr <- Tr + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
        # }
                
        Tr1 <- ggplot(DF_CT, aes_string(x = group_var, y = "Abundance", col = group_var))
        Tr1 <- Tr1 +
                geom_violin(fill = NA) +
                geom_point(aes_string(shape = shape), size = 1, alpha = 0.6, position = position_jitterdodge(jitter.width = 0.8)) +
                facet_wrap(~ Annotation, ncol = facet_cols, scales = "free_y") +
                scale_color_manual("", values = color_levels) +
                xlab("") +
                ylab("abundance") +
                theme_bw()
        
        if (is.null(shape)) {
                Tr1 <- Tr1 + theme(legend.position = "none")
        }    
             
        list(Boxplot = Tr, ViolinPlot = Tr1)   
}             
                
# --






# --
#######################################
### test_differential_abundance_Wilcoxonsingle ##
#################

test_differential_abundance_Wilcoxonsingle <- function(physeq, group_var, compare = NULL, excludeZeros = FALSE, p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))){
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        CT <- as(otu_table(physeq), "matrix")
        
        i <- group_var_levels[1]
        j <- group_var_levels[2]
        
        res_mat <- apply(CT, 2, function(taxon_counts){
                x <- taxon_counts[group_fac == i]
                Zeros_grp1 <- sum(x == 0)
                # # Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                Present_grp1 <- length(x)-Zeros_grp1
                prev_PC_grp1 <- 100*(Present_grp1/length(x))
                if(excludeZeros){
                        x <- x[x != 0]
                }
                Median_grp1 <- median(x, na.rm = T) # NA in case all 0 and excludeZeros was TRUE
                Mean_grp1 <- mean(x, na.rm = T) # NaN in case all 0 and excludeZeros was TRUE
                if (is.na(Mean_grp1)){ Mean_grp1 = NA }
                
                y <- taxon_counts[group_fac == j]
                Zeros_grp2 <- sum(y == 0)
                Present_grp2 <- length(y)-Zeros_grp2
                # # Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                prev_PC_grp2 <- 100*(Present_grp2/length(y))
                if(excludeZeros){
                        y <- y[y != 0]
                }
                Median_grp2 <- median(y, na.rm = T)
                Mean_grp2 <- mean(y, na.rm = T)
                if (is.na(Mean_grp2)){ Mean_grp2 = NA }
                
                if (length(x) != 0 && length(y) != 0){
                        wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                        pValue <- wilcTest$p.value
                        W <- wilcTest$statistic
                        # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                        Ranks <- rank(c(x, y))
                        n1 <- length(x)
                        n2 <- length(y)
                        # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                        # how about the other W?
                        # Wy <- sum(Ranks[(n1+1):(n1+ n2)]) - (n2*(n2+1)/2)
                        standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                        
                        # # if you want to check that multtest::mt.minP would give the same statistic
                        # mati <- matrix(c(x,y), nrow = 1)
                        # grFac <- c(rep(fac_levels[i], n1), rep(fac_levels[j], n2))
                        # grFac <- factor(grFac, levels = c(fac_levels[i], fac_levels[j]))
                        # standStat2 <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                        # # identical(standStat, standStat2) # TRUE
                        # uncomment all with standStat2 to test all the way
                        
                } else {
                        pValue = NA
                        W <- NA
                        standStat = NA
                        n1 <- length(x)
                        n2 <- length(y)
                        # standStat2 = NA
                }
                
                
                c(standStat, pValue, Median_grp1, Median_grp2, 
                  Mean_grp1, Mean_grp2, prev_PC_grp1, prev_PC_grp2, n1, n2, W) 
        })
        
        res_mat <- t(res_mat)
        
        DF <- data.frame(Taxon = rownames(res_mat), res_mat)
        colnames(DF) <- c("Taxon", "teststat", "p_val", "Median_grp1", "Median_grp2", "Mean_grp1", "Mean_grp2", "prev_PC_grp1", "prev_PC_grp2", "n1", "n2", "W")
        DF$p_val_adj <- p.adjust(DF$p_val, method = p.adjust.method)
        
        
        symnum.args$x <- DF$p_val
        DF$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- DF$p_val_adj
        DF$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
        
        DF$direction <- i
        DF$direction[DF$teststat > 0] <- j
        DF$comparison <- paste(group_var_levels, collapse = " vs ")
        
        
        DF <- dplyr::select(DF,  Taxon:p_val, p_val_adj:comparison, Median_grp1:W)
        
        DF <- cbind(DF, as.data.frame(unclass(tax_table(physeq))))
        DF <- dplyr::arrange(DF, desc(abs(teststat)))
        DF
}
# --






# --
#######################################
### test_differential_abundance_WilcoxonsingleManiCoin ##
#################

test_differential_abundance_WilcoxonsingleManiCoin <- function(physeq, group_var, compare = NULL, block = NULL, excludeZeros = FALSE, p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))){
        library(coin)
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        CT <- as(otu_table(physeq), "matrix")
        
        i <- group_var_levels[1]
        j <- group_var_levels[2]
        
        res_mat <- apply(CT, 2, function(taxon_counts){
                x <- taxon_counts[group_fac == i]
                Zeros_grp1 <- sum(x == 0)
                # # Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                Present_grp1 <- length(x)-Zeros_grp1
                prev_PC_grp1 <- 100*(Present_grp1/length(x))
                
                if (excludeZeros){
                        x <- x[x != 0]
                } 
                Median_grp1 <- median(x, na.rm = T) # NA in case all 0 and excludeZeros was TRUE
                Mean_grp1 <- mean(x, na.rm = T) # NaN in case all 0 and excludeZeros was TRUE
                if (is.na(Mean_grp1)){ Mean_grp1 = NA }
                
                y <- taxon_counts[group_fac == j]
                Zeros_grp2 <- sum(y == 0)
                Present_grp2 <- length(y)-Zeros_grp2
                # # Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                prev_PC_grp2 <- 100*(Present_grp2/length(y))
                if(excludeZeros){
                        y <- y[y != 0]
                }
                Median_grp2 <- median(y, na.rm = T)
                Mean_grp2 <- mean(y, na.rm = T)
                if (is.na(Mean_grp2)){ Mean_grp2 = NA }
                
                data <- cbind(RA = taxon_counts, sample_data(physeq))
                # make sure level order is as given by compare (group_var_levels)
                data[[group_var]] <- factor(data[[group_var]], levels = c(group_var_levels, setdiff(levels(group_fac), group_var_levels)), ordered = TRUE)
                
                if (excludeZeros){
                        subset <- group_fac %in% c(i, j) & taxon_counts != 0
                } else {
                        subset <- group_fac %in% c(i, j)
                }
                
                
                formula <- paste0("RA ~ ", group_var)
                if (!is.null(block)) {
                        formula <- paste0(formula, " | ", block)
                }
                formula <- as.formula(formula)
                
                
                if (length(x) != 0 && length(y) != 0){
                        wt <- coin::wilcox_test(formula, data = data, subset = subset, 
                                                conf.int = FALSE, distribution = "asymptotic", 
                                                alternative = "two.sided")
                        pValue <- pvalue(wt)
                        Zz <- statistic(wt)
                        # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                        Ranks <- rank(c(x, y))
                        n1 <- length(x)
                        n2 <- length(y)
                        # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                        # how about the other W?
                        # Wy <- sum(Ranks[(n1+1):(n1+ n2)]) - (n2*(n2+1)/2)
                        standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                        
                        # # if you want to check that multtest::mt.minP would give the same statistic
                        # mati <- matrix(c(x,y), nrow = 1)
                        # grFac <- c(rep(fac_levels[i], n1), rep(fac_levels[j], n2))
                        # grFac <- factor(grFac, levels = c(fac_levels[i], fac_levels[j]))
                        # standStat2 <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                        # # identical(standStat, standStat2) # TRUE
                        # uncomment all with standStat2 to test all the way
                        
                } else {
                        pValue = NA
                        Zz <- NA
                        standStat = NA
                        n1 <- length(x)
                        n2 <- length(y)
                        # standStat2 = NA
                }
                
                
                c(standStat, pValue, Median_grp1, Median_grp2, 
                  Mean_grp1, Mean_grp2, prev_PC_grp1, prev_PC_grp2, n1, n2, Zz) 
        })
        
        res_mat <- t(res_mat)
        
        DF <- data.frame(Taxon = rownames(res_mat), res_mat)
        colnames(DF) <- c("Taxon", "teststat", "p_val", "Median_grp1", "Median_grp2", "Mean_grp1", "Mean_grp2", "prev_PC_grp1", "prev_PC_grp2", "n1", "n2", "Z")
        DF$p_val_adj <- p.adjust(DF$p_val, method = p.adjust.method)
        
        
        symnum.args$x <- DF$p_val
        DF$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- DF$p_val_adj
        DF$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
        
        DF$direction <- i
        DF$direction[DF$teststat > 0] <- j
        DF$comparison <- paste(group_var_levels, collapse = " vs ")
        
        DF <- dplyr::select(DF,  Taxon:p_val, p_val_adj:comparison, Median_grp1:Z)
        
        DF <- cbind(DF, as.data.frame(unclass(tax_table(physeq))))
        DF <- dplyr::arrange(DF, desc(abs(teststat)))
        DF
}


# --



# --
#######################################
### evaluate_fisher_DeSeq_wilcoxon ##
#################
# just a pool function in case you want to evaluate DA on different backgrounds for the same group_var several times used in 00119_DA_Genus_Rifaximin.Rmd of the GalaxyMouse project
# currently very specific for this particular purpose but can of cause be adjusted in the future (more arguments)

evaluate_fisher_DeSeq_wilcoxon <- function(physeq_to_test_fisher, physeq_to_test_DESeq, physeq_to_test_wilcoxon, SFs_all, Trigger, group_var, group_var_levels,
                                           SampleVariable){
        
        
        # - Find the samples to test -
        SS <- as(sample_data(physeq_to_test_fisher), "data.frame")
        SS <- SS[SS$Rel_Time_toTreatment == "After" & SS[[SampleVariable]] == Trigger,] 
        keepSamples <- as.character(SS$Sample_ID)
        # --
        
        
        # - fisher test -
        physeq_to_test <- physeq_to_test_fisher
        physeq_to_test <- prune_samples(samples = keepSamples, physeq_to_test)
        # This gave annoying problems I don't know why:
        # physeq_to_test <- phyloseq::subset_taxa(physeq_to_test, phyloseq::taxa_sums(physeq_to_test) != 0)
        # therefore
        keepTaxa <- names(which(phyloseq::taxa_sums(physeq_to_test) != 0))
        physeq_to_test <- prune_taxa(taxa = keepTaxa, physeq_to_test)
        
        diff_ab_df <- test_diffs_in_prevalence_single(physeq = physeq_to_test, group_var = group_var, compare = group_var_levels, p.adj.method = "fdr", minCount = 0L, symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
        
        hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = NULL)
        taxa_hit_df <- hit_list[["hit_table"]]
        
        taxa_hit_df_fisher <- taxa_hit_df
        taxa_all_fisher <- diff_ab_df
        # --
        
        
        # - DESeq test -
        physeq_to_test <- physeq_to_test_DESeq
        physeq_to_test <- prune_samples(samples = keepSamples, physeq_to_test)
        keepTaxa <- names(which(phyloseq::taxa_sums(physeq_to_test) != 0))
        physeq_to_test <- prune_taxa(taxa = keepTaxa, physeq_to_test)
        SFs_after <- SFs_all[names(SFs_all) %in% keepSamples]
        
        res_list <- test_differential_abundance_DESeq2single(physeq = physeq_to_test, group_var = group_var, compare = group_var_levels, SFs = SFs_after, p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), cooksCutoff = TRUE)
        diff_ab_df <- res_list[[1]]
        physeq_to_test <- res_list[[2]]
        
        hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = NULL)
        taxa_hit_df <- hit_list[["hit_table"]]
        
        taxa_hit_df_DESeq <- taxa_hit_df
        taxa_all_DESeq <- diff_ab_df
        physeq_DESeq <- physeq_to_test
        # --
        
        
        
        
        # - Wilcoxon test -
        physeq_to_test <- physeq_to_test_wilcoxon
        physeq_to_test <- prune_samples(samples = keepSamples, physeq_to_test)
        keepTaxa <- names(which(phyloseq::taxa_sums(physeq_to_test) != 0))
        physeq_to_test <- prune_taxa(taxa = keepTaxa, physeq_to_test)
        
        diff_ab_df <- test_differential_abundance_Wilcoxonsingle(physeq = physeq_to_test, group_var, compare = group_var_levels, excludeZeros = TRUE, p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
        
        hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = NULL)
        taxa_hit_df <- hit_list[["hit_table"]]
        
        taxa_hit_df_Wilcox <- taxa_hit_df
        taxa_all_Wilcox <- diff_ab_df
        # --
        
        # - output the results -
        out <- list(taxa_hit_df_fisher = taxa_hit_df_fisher, taxa_all_fisher = taxa_all_fisher,
                    taxa_hit_df_DESeq = taxa_hit_df_DESeq, taxa_all_DESeq = taxa_all_DESeq, physeq_DESeq = physeq_DESeq,
                    taxa_hit_df_Wilcox = taxa_hit_df_Wilcox, taxa_all_Wilcox = taxa_all_Wilcox)
        out
        # --
        
}
# --





# --
#######################################
### calc_IgAp.Prop ##
#################
# NB: This function is currently only for the DanFunD_IgA study, so it has to be adjusted if you want to use it for more generally.
# But I thought to save it as a function is the best way to save it for later

calc_IgAp.Prop <- function(physeq_to_test, Fraction, Subject){

        # - generate a Data.frame df_calc_IgAp.Prop with Subject, Taxon, IgA-, IgA+, and IgAp.Prop columns -
        OTU <- get_otu_table(physeq_to_test)
        # -- change into a data.frame Taxon Sample rel.Ab --
        OTU <- as.data.frame(OTU)
        OTU$Taxon <- rownames(OTU)
        OTU_w <- tidyr::gather(OTU, key = "Sample", value = "rel.Ab", -Taxon)
        # ----
        # -- add Fraction and Subject information --
        SS <- get_sample_data(physeq_to_test)
        OTU_w$Fraction <- SS[[Fraction]][match(OTU_w$Sample, rownames(SS))]
        OTU_w$Subject <- SS[[Subject]][match(OTU_w$Sample, rownames(SS))]
        OTU_w2 <- dplyr::select(OTU_w, Subject, Taxon, Fraction, rel.Ab)
        # ----
        # -- get separate columns for IgA+ and IgA- and calculate IgAp.Prop --
        df_calc_IgAp.Prop <- tidyr::spread(OTU_w2, key = Fraction, value = rel.Ab)
        df_calc_IgAp.Prop$IgAp.Prop <- df_calc_IgAp.Prop$`IgA+/D2-`/(df_calc_IgAp.Prop$`IgA-` + df_calc_IgAp.Prop$`IgA+/D2-`)
        # ----
        # --
        # - change to wide format also, i.e. ntaxa x subject matrix/data.frame with IgAp.Prop values -
        df_IgAp.Prop <- dplyr::select(df_calc_IgAp.Prop, Taxon, Subject, IgAp.Prop)
        df_IgAp.Prop <- tidyr::spread(df_IgAp.Prop, key = Subject, value = IgAp.Prop)
        rownames(df_IgAp.Prop) <- df_IgAp.Prop$Taxon
        df_IgAp.Prop <- df_IgAp.Prop[, -1]
        # df_IgAp.Prop <- as.matrix(df_IgAp.Prop)
        # --
        list(df_calc_IgAp.Prop = df_calc_IgAp.Prop, df_IgAp.Prop = df_IgAp.Prop, OTU_w = OTU_w)
}
# --





# --
#######################################
### plot_heatmap_IgAp.Prop ##
#################
# NB: THIS function is similar to plot_heatmap_physeq, it is only for the purpose of heatmapping IgAp.Props
# It could surely be simplified when more IgA Seq experiments will be conducted

plot_heatmap_IgAp.Prop <- function (physeq_to_test, DF_IgAp.Prop, Subject, sample_colors = NULL, taxa_info_df = NULL, taxa_colors = NULL, taxa_annotation = NULL, 
                                 max_abundance_for_color = NULL, gradient_steps = c(0.15, 0.3, 0.45, 1), zero_color = "white",
                                 color_function = viridis,
                                 color_steps_bw_markers = 10, drop_color_levels = FALSE,
                                 border_color = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = TRUE, annotation_names_col = TRUE, annotation_legend = TRUE,
                                 legend = TRUE, font_size = 10, fontsize_row = 8, fontsize_col = 8, fontsize_number = 6, filename = NA, ...) {
        
        pruned_ps <- physeq_to_test
        
        
        # - test sample_colors input and prepare sample_info_df -
        if (!is.null(sample_colors)) {
                
                if (!is.list(sample_colors)) {
                        stop("sample_colors must be a named list.")
                }
                
                
                if (!all(names(sample_colors) %in% colnames(sample_data(pruned_ps)))) {
                        stop("all names of the sample_colors list must be variables/columns in sample_data(physeq).")
                }
                
                
                # - here you have to restrict to only the Subject samples! THis is different from plot_heatmap_physeq -
                sample_info_df <- get_sample_data(pruned_ps)
                sample_info_df <- sample_info_df[, match(c(Subject, names(sample_colors)), colnames(sample_info_df)), drop = FALSE] # the match construct makes sure that sample_info_df columns are ordered after names(sample_colors)
                # which makes sure that the column coloring occurs in the whished for order
                colnames(sample_info_df)[1] <- "Subject"
                sample_info_df <- sample_info_df[!duplicated(sample_info_df$Subject),]
                rownames(sample_info_df) <- sample_info_df$Subject
                sample_info_df <- sample_info_df[-1]
                # --
                
                # -- test all entries in sample_colors, prune pruned_ps and order sample_info_df, and set default color where necessary --
                for (i in 1:length(sample_colors)){
                        variable_name <- names(sample_colors)[i]
                        color_character <- sample_colors[[i]]
                        
                        if (!all(is.na(color_character))){
                                
                                # to not remove all samples from pruned_ps:
                                if (!(sum(names(color_character) %in% unique(sample_info_df[[variable_name]])) > 0)) {
                                        stop(paste("Not a single name of a color for ", variable_name, " matched to an entry in the corresponding column in sample_data(physeq)", sep = ""))
                                }
                                
                                # HAS TO BE CHANGED: remove all samples that no color was assigned to in sample_colors, i.e.option to restrict the comparison to defined samples
                                if (!all(unique(sample_info_df[[variable_name]]) %in% names(color_character))) {
                                        
                                        keepSamples <- sample_names(pruned_ps)[sample_info_df[[variable_name]] %in% names(color_character)]
                                        pruned_ps <- prune_samples(keepSamples, pruned_ps)
                                        sample_info_df <- as(sample_data(pruned_ps), "data.frame")
                                        sample_info_df <- sample_info_df[, colnames(sample_info_df) %in% names(sample_colors), drop = FALSE]
                                }
                                
                                # remove unused color levels if asked for it
                                if (drop_color_levels && !all(names(color_character) %in% unique(sample_info_df[[variable_name]]))) {
                                        color_character <- color_character[names(color_character) %in% unique(sample_info_df[[variable_name]])]
                                        sample_colors[[i]] <- color_character
                                }
                                
                                if (!all(areColors(color_character))) {
                                        warning(paste("Not all entries in sample_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(sample_info_df, variable_name)
                                        sample_colors[[i]] <- color_character
                                }
                                
                                # define the order of the samples based on the sample_color vectors (actual ordering below with dplyr::arrange)
                                sample_info_df[[variable_name]] <- factor(sample_info_df[[variable_name]], levels = names(color_character), ordered = TRUE)
                                
                        } else {
                                color_character <- assign_default_colors(sample_info_df, variable_name)
                                sample_colors[[i]] <- color_character
                        }
                }
                # --- order the samples based on the sample_colors ---
                # NB: when you give several sample_color variables the first one given will be the lowest bar in the heatmap, i.e. the bar directly on top of the columns
                # You want to order the samples thus first by the last variable given, then break ties by the next and so on
                # with help from the following two stackoverflow I finally came with a solution that still appears a bit clumsy to me
                # <https://stackoverflow.com/questions/31681894/ordering-data-frame-using-variable-as-column-name>
                # https://stackoverflow.com/questions/1743698/evaluate-expression-given-as-a-string
                
                orderCommand <- paste0("sample_info_df[,", ncol(sample_info_df):1, "]")
                orderCommand <- paste0(orderCommand, collapse = ", ")
                orderCommand <- paste0("order(", orderCommand, ")")
                orderCommand <- parse(text = orderCommand)
                
                # --- you had to add here the rownames (or whatever) again to avoid having a 1 column data frame which would be changed into a factor
                sample_info_df$ID_Saver <- rownames(sample_info_df)
                # ------
                
                sample_info_df <- sample_info_df[eval(orderCommand), ]
                
                sample_info_df <- sample_info_df[-ncol(sample_info_df)]
                # ----
        } else {
                sample_info_df <- NULL
        }
        # --
        
        
        # - set the taxa colors (and pool sample_colors and taxa_colors into annotation_colors) -
        if (!is.null(taxa_colors)) { # see above, guarantees also that taxa_info_id is not NULL
                
                if (!is.list(taxa_colors)) {
                        stop("taxa_colors must be a named list.")
                }
                
                if (!all(names(taxa_colors) %in% colnames(taxa_info_df))) {
                        stop("all names of the taxa_colors list must be variables/columns in taxa_info_df.")
                }
                
                taxa_info_df <- taxa_info_df[, match(names(taxa_colors), colnames(taxa_info_df)), drop = FALSE] # again match to get coloring ordered after taxa_colors
                
                # -- test all entries in taxa_colors --
                for (i in 1:length(taxa_colors)){
                        variable_name <- names(taxa_colors)[i]
                        color_character <- taxa_colors[[i]]
                        
                        if (!all(is.na(color_character))){
                                
                                
                                if (!all(unique(taxa_info_df[[variable_name]]) %in% names(color_character))) {
                                        warning(paste("There were levels in taxa_info_df at ", variable_name, " for which no color was assigned in taxa_colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(taxa_info_df, variable_name)
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                # remove unused color levels if asked for it
                                if (drop_color_levels && !all(names(color_character) %in% unique(taxa_info_df[[variable_name]]))) {
                                        color_character <- color_character[names(color_character) %in% unique(taxa_info_df[[variable_name]])]
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                
                                if (!all(areColors(color_character))) {
                                        warning(paste("Not all entries in taxa_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(taxa_info_df, variable_name)
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                # define the order of the taxa coloring based on the taxa_color vectors
                                # taxa_info_df[[variable_name]] <- factor(taxa_info_df[[variable_name]], levels = names(color_character), ordered = TRUE) # was unnecessary, pheatmap orders based on annotation_colors
                                
                        } else {
                                color_character <- assign_default_colors(taxa_info_df, variable_name)
                                taxa_colors[[i]] <- color_character
                        }
                }
                # ----
        } 
        
        annotation_colors <- c(sample_colors, taxa_colors)
        # -- 
        
        
        # - check or set taxa_annotation -
        if (is.null(taxa_annotation)){
                taxa_annotation <- taxa_names(pruned_ps) # or you could use 
                # taxa_annotation <- get_taxon_names(as.data.frame(tax_table(pruned_ps)))
        }
        
        if (length(taxa_annotation) != ntaxa(pruned_ps)) {
                warning("taxa_annotation did not fit in length to nrow(taxa_info_df) or ntaxa(physeq), used taxa_names")
                taxa_annotation <- taxa_names(pruned_ps)
        }
        
        taxa_annotation <- as.character(taxa_annotation)
        
        taxa_annotation[is.na(taxa_annotation)] <- "NA"
        
        taxa_annotation <- make.unique(taxa_annotation)
        # --
        
        
        # - prepare DF_IgAp.Prop so that taxa are rows in the order determined by taxa_info_df -
        if (!is.null(taxa_info_df)){
                DF_IgAp.Prop <- DF_IgAp.Prop[rownames(taxa_info_df), ]
        }
        # --
        
        # - order the samples in DF_IgAp.Prop based on sample_colors list (see above for sample_info_df) -
        if (!is.null(sample_info_df)){
                DF_IgAp.Prop <- DF_IgAp.Prop[, rownames(sample_info_df)]
        }
        # --
        
        
        # - test and adjust max_abundance_for_color -
        if (is.null(max_abundance_for_color)) {
                max_abundance_for_color <- max(DF_IgAp.Prop, na.rm = TRUE)
        }
        if (max_abundance_for_color < 0 && max_abundance_for_color > 1) {
                max_abundance_for_color <- max(DF_IgAp.Prop, na.rm = TRUE)
        }
        # --
        
        
        # - set breaks and colors for pheatmap  -
        # -- set all NaN, i.e. where taxon was not present to -.1 --
        DF_IgAp.Prop[is.na(DF_IgAp.Prop)] <- -.1
        ZeroValue <- min(DF_IgAp.Prop) # here than obviously -.1
        min_in_data <- min(DF_IgAp.Prop[DF_IgAp.Prop > ZeroValue]) # the non NaN:
        # NB NB NB: I had trouble that the way the function was, 0 values were plotted in Gray, this is an easy ugly fix
        # because I was too lazy to go through the myBreaks code
        min_in_data <- min_in_data - 0.001
        max_in_data <- max(DF_IgAp.Prop)
        
        # -- make sure the last entry in gradient steps = 1 --
        if (!all(gradient_steps >= 0 & gradient_steps <= 1)) {
                gradient_steps <- c(0.15, 0.3, 0.45, 1)
        }
        
        if (gradient_steps[length(gradient_steps)] != 1) {
                gradient_steps <- c(gradient_steps, 1)
        }
        # ----
        
        # -- added this even though it might be confusing to the user, but otherwise you get an error --
        if (max_abundance_for_color < min_in_data) {
                max_abundance_for_color <- min_in_data # NB: in this case all non-zero values should become maximal color
        }
        # -----
        
        
        # you want that the final color gradient covers the values min to max_abundance_for_color
        # 0 values will be set to a different color (usually red or white), values above max_abundance_for_color should be all max_color
        # normalise gradient steps with max_abundance_for_color
        myBreaks <- gradient_steps * max_abundance_for_color
        myBreaks <- myBreaks[myBreaks >= min_in_data]
        # add break at min_in_data
        if (length(myBreaks) != 1 || myBreaks != max_abundance_for_color){ # in case max_abundance_for_color == min_in_data
                myBreaks <- c(min_in_data, myBreaks)
        }
        
        if (length(myBreaks) > 1){
                
                # now myBreaks goes from min_in_data up to max_abundance_for_color (provided the last gradient_steps was 1)
                myColors = color_function(length(myBreaks)) # see help pheatmap, breaks should be 1 element longer than color, now it is same legnth
                
                # myColors contains now the colors that represent the gradient_steps values (= markers).
                # now we want to introduce breaks between these markers and make linear color gradients between the markers. The number of breaks between the markers is defined
                # by color_steps_bw_markers
                
                myColors <- unlist(lapply(1:(length(myColors)-1), function(i) {
                        colorRampPalette(colors = c(myColors[i], myColors[i+1]))(color_steps_bw_markers)
                }))
                
                
                
                myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
                        seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
                })
                myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)])
                
                # add max_in_data values to myBreaks, otherwise all values above max_abundance_for_color will be white
                if (max_in_data > myBreaks[length(myBreaks)]) {
                        myBreaks <- c(ZeroValue, myBreaks, max_in_data)
                        myColors <- c(zero_color, myColors, myColors[length(myColors)])
                } else {
                        myBreaks <- c(ZeroValue, myBreaks)
                        myColors <- c(zero_color, myColors)
                }
                
        } else {
                
                myColors = color_function(2)
                myColors <- myColors[2]
                
                
                
                if (max_in_data > myBreaks[length(myBreaks)]) {
                        myBreaks <- c(ZeroValue, myBreaks, max_in_data)
                        myColors <- c(zero_color, myColors, myColors[length(myColors)])
                } else {
                        myBreaks <- c(ZeroValue, myBreaks)
                        myColors <- c(zero_color, myColors)
                }
                
                
        }
        
        
        # ----
        # --
        
        # - not necessary, just for clarity -
        if (is.null(annotation_colors)) {annotation_colors <- NA}
        if (is.null(sample_info_df)) {sample_info_df <- NA}
        if (is.null(taxa_info_df)) {taxa_info_df <- NA}
        # --
        
        # --
        hm.parameters <- list(DF_IgAp.Prop, color = myColors, breaks = myBreaks, border_color = border_color, cluster_cols = cluster_cols, cluster_rows = cluster_rows,
                              show_rownames = show_rownames, show_colnames = show_colnames, annotation_col = sample_info_df,
                              annotation_row = taxa_info_df, annotation_colors = annotation_colors, labels_row = taxa_annotation, annotation_names_row = annotation_names_row,
                              annotation_names_col = annotation_names_col, annotation_legend = annotation_legend, legend = legend, font_size = font_size, 
                              fontsize_row = fontsize_row, fontsize_col = fontsize_col, fontsize_number = fontsize_number)
        do.call("pheatmap", hm.parameters)
}
# --




# --
#######################################
### test_differential_IgAp.Prop_Wilcoxonsingle ##
#################

test_differential_IgAp.Prop_Wilcoxonsingle <- function(DF_IgAp.Prop, SS, group_var, compare = NULL, 
                                                       excludeZeros = FALSE, p.adjust.method = "fdr", 
                                                       symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), 
                                                                          symbols = c("****", "***", "**", "*", "ns"))){
        
        
        if(! group_var %in% colnames(SS)) {
                stop("The given group_var is not a variable in the sample data SS.")
        }
        
        group_fac <- factor(SS[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        i <- group_var_levels[1]
        j <- group_var_levels[2]
        
        mat_IgAp.Prop <- as.matrix(DF_IgAp.Prop)
        
        res_mat <- apply(mat_IgAp.Prop, 1, function(taxon_IgAp.Props){
                x <- taxon_IgAp.Props[group_fac == i]
                Zeros_grp1 <- sum(is.na(x))
                # # Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                Present_grp1 <- length(x)-Zeros_grp1
                prev_PC_grp1 <- 100*(Present_grp1/length(x))
                if(excludeZeros){
                        x <- x[!is.na(x)]
                }
                Median_grp1 <- median(x, na.rm = T) # NA in case all 0 and excludeZeros was TRUE
                Mean_grp1 <- mean(x, na.rm = T) # NaN in case all 0 and excludeZeros was TRUE
                if (is.na(Mean_grp1)){ Mean_grp1 = NA }
                
                y <- taxon_IgAp.Props[group_fac == j]
                Zeros_grp2 <- sum(is.na(y))
                Present_grp2 <- length(y)-Zeros_grp2
                # # Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                prev_PC_grp2 <- 100*(Present_grp2/length(y))
                if(excludeZeros){
                        y <- y[!is.na(y)]
                }
                Median_grp2 <- median(y, na.rm = T)
                Mean_grp2 <- mean(y, na.rm = T)
                if (is.na(Mean_grp2)){ Mean_grp2 = NA }
                
                if (length(x) != 0 && length(y) != 0){
                        wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                        pValue <- wilcTest$p.value
                        W <- wilcTest$statistic
                        # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                        Ranks <- rank(c(x, y))
                        n1 <- length(x)
                        n2 <- length(y)
                        # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                        # how about the other W?
                        # Wy <- sum(Ranks[(n1+1):(n1+ n2)]) - (n2*(n2+1)/2)
                        standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                        
                        # # if you want to check that multtest::mt.minP would give the same statistic
                        # mati <- matrix(c(x,y), nrow = 1)
                        # grFac <- c(rep(fac_levels[i], n1), rep(fac_levels[j], n2))
                        # grFac <- factor(grFac, levels = c(fac_levels[i], fac_levels[j]))
                        # standStat2 <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                        # # identical(standStat, standStat2) # TRUE
                        # uncomment all with standStat2 to test all the way
                        
                } else {
                        pValue = NA
                        W <- NA
                        standStat = NA
                        n1 <- length(x)
                        n2 <- length(y)
                        # standStat2 = NA
                }
                
                
                c(standStat, pValue, Median_grp1, Median_grp2, 
                  Mean_grp1, Mean_grp2, prev_PC_grp1, prev_PC_grp2, n1, n2, W) 
        })
        
        res_mat <- t(res_mat)
        
        DF <- data.frame(Taxon = rownames(res_mat), res_mat)
        colnames(DF) <- c("Taxon", "teststat", "p_val", "Median_grp1", "Median_grp2", "Mean_grp1", "Mean_grp2", "prev_PC_grp1", "prev_PC_grp2", "n1", "n2", "W")
        DF$p_val_adj <- p.adjust(DF$p_val, method = p.adjust.method)
        
        
        symnum.args$x <- DF$p_val
        DF$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- DF$p_val_adj
        DF$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
        
        DF$direction <- i
        DF$direction[DF$teststat > 0] <- j
        DF$comparison <- paste(group_var_levels, collapse = " vs ")
        
        
        DF <- dplyr::select(DF,  Taxon:p_val, p_val_adj:comparison, Median_grp1:W)
        
        #DF <- cbind(DF, as.data.frame(unclass(tax_table(physeq))))
        DF <- dplyr::arrange(DF, desc(abs(teststat)))
        DF
}
# --

############# WORK here #########








