# --
###########################
### FUNCTION: check_user_parameters ###
############################
check_user_parameters <- function(group_var, group_var_levels, color_levels, ps){
        

if(!group_var %in% colnames(sample_data(ps))) {
        error_message <- "The given group_var is not a variable in the sample data of the phyloseq object. Please change input parameters."
        return(error_message)
}


if (!all(group_var_levels %in% unique(sample_data(ps)[[group_var]]))) {
        error_message <- "Sorry. Group 1 or Group 2 is not a level in the given group_var column. Please change input parameters."
        return(error_message)
}


if (!all(areColors(color_levels))) {
        error_message <- "Sorry. At least one of the given colors is not an R color. Please change input parameters."
        return(error_message)
}
        NULL
}
# --



# --
###########################
### FUNCTION: check_user_parameters_4 ###
############################
check_user_parameters_4 <- function(group_var, group_var_levels, color_levels, compare, ps){
        
        
        if(!group_var %in% colnames(sample_data(ps))) {
                error_message <- "The given group_var is not a variable in the sample data of the phyloseq object. Please change input parameters."
                return(error_message)
        }
        
        
        if (!all(group_var_levels %in% unique(sample_data(ps)[[group_var]]))) {
                error_message <- "Sorry. One of your group_var_levels is not a level in the given group_var column. Please change input parameters."
                return(error_message)
        }
        
        if (!all(compare %in% unique(sample_data(ps)[[group_var]]))) {
                error_message <- "Sorry. One of your compare levels is not a level in the given group_var column. Please change input parameters."
                return(error_message)
        }
        
        if (!all(compare %in% group_var_levels)) {
                error_message <- "Sorry. Your compare levels must be included in group_var_levels."
                return(error_message)
        }
        

        
        if (!all(areColors(color_levels))) {
                error_message <- "Sorry. At least one of the given colors is not an R color. Please change input parameters."
                return(error_message)
        }
        NULL
}
# --







# --
###########################
### FUNCTION: output_ps ###
############################
output_ps <- function(ps){
        
        
        if(is.null(ps)) {
                return("NULL")
        } else {
                
                TT <- as.data.frame(as(tax_table(ps), "matrix"))
                Assigned <- sapply(TT, function(col){!all(is.na(col))})
                currentLevel <- names(which(Assigned))[length(which(Assigned))]
                if(is.null(ps@phy_tree)){
                        Nnode <- "NULL"
                } else {
                        Nnode <- ps@phy_tree$Nnode
                }
                
                paste0("Object with ", ntaxa(ps), " taxa and ", nsamples(ps), " samples. Assigned to ", currentLevel, ". No of sample_data variables: ", ncol(sample_data(ps)), ". No of nodes in tree: ",
                       Nnode, ".")
        }
        
}
# --


