# MixrobiomeX4_App

The app allows microbiome analyses starting from a phyloseq object.



## To Improve



- in beta_diversity tab: # NB: in  calc_beta_div_distances you use group_var_levels for compare, so you calculate always distances for all levels even if compare is smaller, could be changed but not strictly necessary
- in general there are more plots where you could get more consistent about only showing the compare levels if asked to do so.
- allow do download plots
- add a Phewas tab, should have been made by the students last year





## alpha diversity

- NB: the function **estimateR.default** in phyloseq::estimate_richness accepts only integers, but ps_tca and ps_tca_RA are not integers after size factor adjustment.
- **I therefore put them to integers before calculating alpha diversity, this changes the sample_sums a bit, so that the sample_sums of ps_tca_RA are not all equal, OF COURSE: for these choices you should ignore the corrected versions of the alpha diversity plots, here sample_sums have been corrected beforehand!** In other words:
- **the corrected versions are only for raw counts with unadjusted library sizes**



### Beta Diversity Tab

- I separated the distance calculation from the PCOA plot, because the distance calculation can take a lot of time, and I want to allow the user to make a lot of plots with just calculating the distances once

- Distance calculation: distances are calculated for all samples that are part of group_var_levels, so don't change group_var_levels between distances calculation and plotting. The app will tell you that dist_list does not fit to the group_var_levels when you plot.

- plotting and adonis: 
    - adonis is always only calculated for the levels defined by compare
    -


