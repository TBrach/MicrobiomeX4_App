# MixrobiomeX4_App

The app allows microbiome analyses starting from a phyloseq object.


## To Improve

- allow to download plots
- I guess there are more plots in the app where the plot should only show the compare levels if asked to do so by the user.
- make a tab to analyse metadata variables (simple things like table(SD$variable))
- maybe add a Phewas tab
- a tab to analyse microbiome parameters against a continuous variable
- implement reference-frame compositional analysis?
- in beta_diversity tab: # NB: in  calc_beta_div_distances you use group_var_levels for compare, so you calculate always distances for all levels even if compare is smaller. This can make sense if you want to plot the other levels as gray dots. So I think it makes sense even thought it costs time.


## alpha diversity

- I recommend to calculate alpha-diversity measures on raw counts with unadjusted library sizes, and then consider the "library size adjusted corrected plots."
    - If you start with already library-size adjusted input, ignore the corrected plots.
- NB: the function estimateR.default in phyloseq::estimate_richness accepts only integers, but ps_tca and ps_tca_RA are not integers after size factor adjustment. I therefore put them to integers before calculating alpha diversity, this changes the sample_sums a bit, so that the sample_sums of ps_tca_RA are not all equal.


### Beta Diversity Tab

- I separated the distance calculation from the PCOA plot, because the distance calculation can take a lot of time, and I want to allow the user to play with the plot choices after the distanes have been calculated. 
    - NB: Distances are calculated for all samples that are part of group_var_levels, so don't change group_var_levels between distances calculation and plotting. If you do, the app will tell you that dist_list does not fit to the group_var_levels when you plot.
- adonis significance values are always only calculated for the levels defined by compare in the first tab.



