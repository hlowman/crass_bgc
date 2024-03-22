# crass_bgc

<img align="left" width="15%" src=figures/CRASS_logo.png> 

Scripts, data, and figures for Collaborative for Research in Arid Stream Systems Biogeochemistry working group. Note, these resources were created primarily to support the following manuscript, but they also contain data and analysis scripts for other projects:

H.E. Lowman, J.R. Blaszczak, A. Cale, X. Dong, S. Earl, J. Grabow, N.B. Grimm, T.K. Harms, J.M. Melack, A.M. Reinhold, B. Summers, A. Webster. *In Review at Biogeochemistry Letters*. Persistent and lagged effects of fire on stream solutes linked to intermittent precipitation in arid lands.

**Please visit the Zenodo webpage for this project to access the most relevant and annotated versions of data and code - link forthcoming!**

- The `data_raw` folder contains raw chemistry, discharge, precipitation, and watershed characteristic data from the Mediterranean (Santa Barbara, CA) and monsoonal (Valles Caldera, NM) stream locations.
- The `data_working` folder contains the aggregated datasets used to fit the MARSS models as well as the MARSS model outputs.
- The `figures` folder contains all versions of figures created in this project from initial data visualization and examination to final manuscript figures generated.
- The `R` folder contains additional folders with scripts for initial data visualization (`data_vis`, `large_rivers`), sample collaborator scripts (`others_scripts`), and the MARSS model prep and fit scripts (`marss`) which are presented in more detail below:
  - Specific conductance data was prepped for modeling using `cond_marss_model_data_prep.R` and fit using `cond_marss_models.R`.
  - Nutrient (i.e., nitrate, ammonium, and phosphate) data was prepped for modeling using `nutrient_marss_model_data_prep.R` and fit using `no3_marss_models.R`, `nh4_marss_models.R`, and `po4_marss_models.R`. The final manuscript figure was assembled using `nutrient_marss_figure.R`.
  - Supplementary manuscript figures were created using `chem_figures.R`.
  - Additional MARSS models fit to satisfy reviewer comments may be found in `supp_marss_models.R`.
 
Data for this project were downloaded from multiple published data sources including:

(1) LTER SBC, Melack J.M. (2019) Precipitation at Arroyo Hondo at Upper Outlaw Trail (HO202) ver 10. Environmental Data Initiative. [https://doi.org/10.6073/pasta/128b6fc45f7d840a673b3c3b2f2eaa32](https://doi.org/10.6073/pasta/128b6fc45f7d840a673b3c3b2f2eaa32)
(2) LTER SBC, Melack J.M. (2019) SBC LTER: Land: Hydrology: Santa Barbara County Flood Control District - Precipitation at Cater Water Treatment Plant (CaterWTP229). [https://doi.org/10.6073/PASTA/8E12EFBC31D84DBED0C105F8C8A2486F](https://doi.org/10.6073/PASTA/8E12EFBC31D84DBED0C105F8C8A2486F)
(3) LTER SBC, Melack J.M. (2019) SBC LTER: Land: Hydrology: Santa Barbara County Flood Control District - Precipitation at El Deseo (ElDeseo255). [https://doi.org/10.6073/PASTA/C5819A9D15A98696E65F0EFD8B45C14A](https://doi.org/10.6073/PASTA/C5819A9D15A98696E65F0EFD8B45C14A)
(4) LTER SBC, Melack J.M. (2019) SBC LTER: Land: Hydrology: Stream discharge and associated parameters at Arroyo Burro Creek, Cliff Drive (AB00). [https://doi.org/10.6073/PASTA/66F46E14C73B5CE07EF104B91C09360B](https://doi.org/10.6073/PASTA/66F46E14C73B5CE07EF104B91C09360B)
(5) LTER SBC, Melack J.M. (2019) SBC LTER: Land: Hydrology: Stream discharge and associated parameters at Gaviota Creek, Hwy 101 South Rest Stop (GV01). [https://doi.org/10.6073/PASTA/3CDF7FD60F88F0822ECFAD5A92EDF9D6](https://doi.org/10.6073/PASTA/3CDF7FD60F88F0822ECFAD5A92EDF9D6)
(6) LTER SBC, Melack J.M. (2019) SBC LTER: Land: Hydrology: Stream discharge and associated parameters at Arroyo Hondo Creek , Upstream Side of 101 Bridge (HO00). [https://doi.org/10.6073/PASTA/85F4639AC0AC76C1EB4B3A82842CE171](https://doi.org/10.6073/PASTA/85F4639AC0AC76C1EB4B3A82842CE171)
(7) LTER SBC, Melack J.M. (2019) SBC LTER: Land: Hydrology: Stream discharge and associated parameters at Rattlesnake Creek, Las Canoas Rd (RS02). [https://doi.org/10.6073/PASTA/00AD8507A9EED8E923867C67D907AF12](https://doi.org/10.6073/PASTA/00AD8507A9EED8E923867C67D907AF12)
(8) LTER, SBC, Melack, J.M. (2020) SBC LTER: Land: Stream chemistry in the Santa Barbara Coastal drainage area, 2000-2018. Environmental Data Initiative. [https://doi.org/10.6073/pasta/ae8784dcba16a65aa3c8e44321b7c1bb](https://doi.org/10.6073/pasta/ae8784dcba16a65aa3c8e44321b7c1bb)
(9) LTER SBC, Melack J.M. (2022) SBC LTER: Land: Hydrology: Precipitation at Gaviota at Las Cruces School (GV202). [https://doi.org/10.6073/PASTA/9E655E16C8386E8A0CD475F9EB4A8FCD](https://doi.org/10.6073/PASTA/9E655E16C8386E8A0CD475F9EB4A8FCD)
(10) Western Regional Climate Center (2021) Valles Caldera National Preserve Climate Stations. [https://wrcc.dri.edu/vallescaldera](https://wrcc.dri.edu/vallescaldera)

Stream chemistry data for the monsoonal sites located in the Valles Caldera National Preserve were provided by L. Crossey, C. McGibbon, and M. Albonico; it is in the process of being made publicly available by the National Parks Service at [https://irma.nps.gov/aqwebportal](https://irma.nps.gov/aqwebportal).

At sites monitored by the USGS, discharge data was downloaded using the `dataRetrieval` package in R.

For additional information regarding this project, please contact Heili at heili.lowman _at_ duke.edu or Tamara at tharms _at_ ucr.edu
