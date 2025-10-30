# priorRelevance
This public repository contains the analysis code and behavioral data for the following manuscript:
"How relevant is the prior? Bayesian causal inference for dynamic perception in volatile environments"
by David Meijer, Roberto Barumerli and Robert Baumgartner. It was accepted to eLife and is currently available as a reviewed preprint (revision is on the way):
https://elifesciences.org/reviewed-preprints/105385 (https://doi.org/10.7554/eLife.105385.1).

The behavioral data (.mat files in the "final_behav_eye_data_Krishnamurthy" folder) are part of the dataset that was collected by 
Kamesh Krishnamurthy, Matthew R. Nassar, Shilpa Sarode & Joshua I. Gold, and it was published in Nature Human Behaviour (2017):
"Arousal-related adjustments of perceptual biases optimize perception in dynamic environments", https://doi.org/10.1038/s41562-017-0107. 
The authors gave us permission to share their behavioral data. 
Despite the file names, we did not include the eyetracker data here to reduce file sizes, and because it is not relevant to our study.
Note that the estimation responses and single-sound trials (condition number 3) were included, although they are not used in our analysis. 
The data from the first session of each participant is not included here, because those were considered as training sessions. 

The fitting function ("krishnamurthy_reanalysis_fitfun_2024.m") in the fitting_functions folder loads the relevant behavioral data and calls the model code to optimize the parameters. 
All model code is found in the "BCPmodel2024" folder. The model code can be called with different settings, such that a different model is fit for different settings. 
For example usage, have a look at "modelPlay.m" in the various folder, or have a look at one of the 21 scripts in the fitting_functions folder, which we used for the paper. 
Note that model fitting (i.e., optimizing parameters) can take a long time (minutes to hours), which is why we ran the fitting functions on a computing server.

To optimize the model parameters we used the Bayesian Adaptive Direct Search (BADS) algorithm, by Luigi Acerbi and Wei Ji Ma. 
Please download the toolbox from their repository directly, https://github.com/acerbilab/bads, and place it in the "bads-master" folder.

The output files of the model fitting procedure (i.e., "BCPfitResults.mat" files) are not included in this repository because of their file sizes. 
Instead, the fitted parameter values were extracted with "extract_fits.m" and saved as "fitted_params_2024.mat" in the "fitted_data_Krishna2017_2024" folder. 

The raw manuscript figures can be created by calling each of the "Create_fig" scripts. These use the helper functions in the "utilities" folder. 
"computeSurprisal.m" was used to create "surprisal_surface_plot_data.mat", which is loaded in "Create_fig_5.m". 

The "Create_fig" scripts also run the basic statistical tests that are reported in the manuscript. 
We used JASP (https://jasp-stats.org/) for the repeated-measures ANOVA and accompanying post-hoc tests (see "stats" folder).

For any questions, please feel free to send an email to: MeijerDavid1@gmail.com
