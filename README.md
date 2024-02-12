# crass_bgc

<img align="left" width="15%" src=figures/CRASS_logo.png> 

Scripts, data, and figures for Collaborative for Research in Arid Stream Systems Biogeochemistry working group. Note, these resources were created primarily to support the following manuscript, but they also contain data and analysis scripts for other projects:

H.E. Lowman, J.R. Blaszczak, A. Cale, X. Dong, S. Earl, J. Grabow, N.B. Grimm, T.K. Harms, J.M. Melack, A.M. Reinhold, B. Summers, A. Webster. *In Review at Biogeochemistry Letters*. Persistent and lagged effects of fire on stream solutes linked to intermittent precipitation in arid lands.

**Please visit the Zenodo webpage for this project to access the most relevant and annotated versions of data and code - link forthcoming!**

- The `data_raw` folder contains raw chemistry, discharge, precipitation, and watershed characteristic data from the Mediterranean (Santa Barbara, CA) and monsoonal (Valles Caldera, NM) stream locations.
- The `data_working` folder contains the aggregated dataset used to fit the MARSS models as well as the MARSS model outputs.
- The `figures` folder contains all versions of figures created in this project from initial data visualization and examination to final manuscript figures generated.
- The `R` folder contains additional folders with scripts for initial data visualization (`data_vis`, `large_rivers`), sample collaborator scripts (`others_scripts`), and the MARSS model prep and fit scripts (`marss`) which are presented in more detail here:
  - Specific conductance data was prepped for modeling using `cond_marss_model_data_prep.R` and fit using `cond_marss_models.R`.
  - Nutrient (i.e., nitrate, ammonium, and phosphate) data was prepped for modeling using `nutrient_marss_model_data_prep.R` and fit using `no3_marss_models.R`, `nh4_marss_models.R`, and `po4_marss_models.R`. The final manuscript figure was assembled using `nutrient_marss_figure.R`.
  - Supplementary manuscript figures were created using `chem_figures.R`.
 
Data for this project were downloaded from multiple published data sources including:

LTER, SBC, Melack, J.M. (2020) SBC LTER: Land: Stream chemistry in the Santa Barbara Coastal drainage area, 2000 -2018. [https://doi.org/10.6073/pasta/ae8784dcba16a65aa3c8e44321b7c1bb](https://doi.org/10.6073/pasta/ae8784dcba16a65aa3c8e44321b7c1bb)

At sites monitored by the USGS, discharge data was downloaded using the `dataRetrieval` package in R.

For additional information regarding this project, please contact Heili at heili.lowman _at_ duke.edu or Tamara at tharms _at_ ucr.edu
