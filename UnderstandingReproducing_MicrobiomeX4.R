# - Index -
functionpath <- "/Users/jvb740/Coursera_MOOC/20161202_LearningShiny_FantasySports/shinyy/Apps/Teaching_Apps/MicrobiomeX4_App/Functions"
source(file.path(functionpath, "_n_000_helper_functions.R"))
source(file.path(functionpath, "_n_010_explore_ps_functions.R"))
source(file.path(functionpath, "_n_020_alpha_diversity_functions.R"))
source(file.path(functionpath, "_n_030_preprocess_filtering_functions.R"))
source(file.path(functionpath, "_n_040_beta_diversity_functions.R"))
source(file.path(functionpath, "_n_050_diff_abundance_functions.R"))
source(file.path(functionpath, "_n_060_phylum_analysis_functions.R"))
source(file.path(functionpath, "shiny_app_functions.R"))
# --

# - TMicrobiomeX4_server -
rv <- list()
# --

# - Tab01 -
rv$group_var <- "Patient_ID."
group_var_levels <- "301, 302, 303, 304, 305, 306, 307, 308, 309, 310"
group_var_levels <- unlist(strsplit(group_var_levels, split = ","))
rv$group_var_levels <- trim_whitespace(group_var_levels)

color_levels <- "green, red, blue, yellow, black, orange, gray, brown, cyan, darkgreen"
color_levels <- trim_whitespace(unlist(strsplit(color_levels, split = ",")))
if (length(rv$group_var_levels) != length(color_levels)){
        rv$infoText_1 <- "Length group_var_levels did not fit to length color_levels, you need to change this"
        return()
}
names(color_levels) <- rv$group_var_levels
rv$color_levels <- color_levels

compare <- "301, 302"
if (compare == ""){
        rv$compare <- rv$group_var_levels
} else {
        compare <- unlist(strsplit(compare, split = ","))
        rv$compare <- trim_whitespace(compare)
}
# --


# - Tab02 -
datapath <- "/Users/jvb740/Coursera_MOOC/20161202_LearningShiny_FantasySports/shinyy/Apps/Teaching_Apps/MicrobiomeX4_App/Data"
ps <- readRDS(file = file.path(datapath, "ps_healthy_2.rds"))
keepTaxa <- taxa_names(ps)[taxa_sums(ps) > 0]
ps <- phyloseq::prune_taxa(keepTaxa, ps)
rv$ps <- ps

# -- calculate phyla --
Phyla <- check_phyla_distribution_NA(ps)
PhylaForColor <- check_phyla_distribution(ps)

if (nrow(PhylaForColor) < 16) {
        phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), QuantColors15)
} else {
        phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), c(QuantColors15, viridis(nrow(PhylaForColor)-15)))
}

rv$phylum_colors <- phylum_colors
# ----
# --


# - Tab 03 Size factors -
prevalence_for_sf<- 40
min_obs <- 0L
ps_sf_filt <- phyloseq::filter_taxa(rv$ps, function(x){(sum(x > min_obs) > (prevalence_for_sf/100)*length(x))}, prune = TRUE)
taxa_used_for_SFs <- ntaxa(ps_sf_filt)
SFs <- calc_SFs(physeq = ps_sf_filt)
group_var <- rv$group_var
SFs_poscounts <- calc_SFs_DESeq(ps_sf_filt, type = "poscounts", group_var = group_var)
SFs_RA <- sample_sums(rv$ps)/gm_own(sample_sums(rv$ps), zeros.count = FALSE)
rv$SFs <- SFs
rv$SFs_poscounts <- SFs_poscounts
rv$SFs_RA <- SFs_RA
library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs) 
rv$ps_tca <- library_size_adjust_list[[1]]
library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs_poscounts) 
rv$ps_tca_poscounts <- library_size_adjust_list[[1]]
library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs_RA) 
rv$ps_RA <- library_size_adjust_list[[1]]
# --


# - Tab 06 Beta Diversity -
ps <- rv$ps_RA
measure <- "bray"
group_var <- rv$group_var
group_var_levels <- rv$group_var_levels 
color_levels <- rv$color_levels
error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                       color_levels = color_levels, ps = ps)
ps <- phyloseq::transform_sample_counts(ps, function(x){x/sum(x)})
dist_list <- calc_beta_div_distances(ps, dist_methods = measure, group_var = group_var, compare = group_var_levels)
rv$betaDivDist <- dist_list
phylum_colors <- rv$phylum_colors
compare <- rv$compare
error_message <- check_user_parameters_4(group_var = group_var, group_var_levels = group_var_levels, 
                                         color_levels = color_levels, compare = compare, ps = ps)
group_factor <- sample_data(ps)[[group_var]]
# make sure group_factor fits to dist_list, i.e. only keep samples covered by compare = group_var_levels!
group_factor <- factor(group_factor[group_factor %in% group_var_levels], levels = group_var_levels, ordered = T)

dist_objj <- dist_list[[1]]
dist_list_adonis <- dist_list
if (!all(group_var_levels %in% compare)){ # all(compare %in% group_var_levels) is tested in check_user_parameters_4
        
        group_factor <- sample_data(ps)[[group_var]]
        # make sure group_factor fits to dist_list, i.e. only keep samples covered by compare = group_var_levels!
        group_factor <- factor(group_factor[group_factor %in% compare], levels = compare, ordered = T)
        
        keepSamples <- sample_names(ps)[sample_data(ps)[[group_var]] %in% compare]
        
        dist_list_adonis <- lapply(dist_list_adonis, filter_dist_obj, keepNames = keepSamples)
}
if (length(levels(group_factor)) > 1){ # no adonis test if only one factor level left because that would cause an error
        
        adonis_list <- lapply(dist_list_adonis, function(dist_obj){
                loop_vegan_adonis(dist_obj = dist_obj, group_fac = group_factor)
        })
        
} else {
        adonis_list <- NULL
}

gray_levels <- NULL

if (!all(group_var_levels %in% compare)){
        
        if (input$compareUse == "remove other levels"){
                dist_list <- dist_list_adonis
                color_levels <- color_levels[names(color_levels) %in% compare]
        } else if (input$compareUse == "keep other levels as gray dots"){
                gray_levels <- compare
        }
}
coordCor <- TRUE
ellipse <- TRUE
ellipse_level <- 0.95
# --



