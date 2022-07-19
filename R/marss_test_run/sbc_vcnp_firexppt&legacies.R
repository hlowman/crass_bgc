# SB & VC MARSS models with fire x ppt interactions and legacy effects
# Script started July 19, 2022
# Heili Lowman, Alex Webster

# The following script will run a MARSS analysis at SBC LTER & NM Valles Caldera sites for the CRASS project.

# A few notes about the process below:
# Data was tidied in sbc_marss_model.R script and exported to file:
# marss_data_sb_vc_060622.rds (as of July 19, 2022)
# It is imported from that file in this script.

# The first goal of this script is to demo a MARSS model with 12 states (each watershed as a unique state) and the following vars:
# - cum. monthly ppt
# - fire p/a  (1=fire ignited in this month, 0=no fire ignitions, one column for all fires)
# - fire p/a 1 yr legacy
# - fire p/a 2 yr legacy
# - fire p/a 3 yr legacy
# - fire p/a 4 yr legacy
# - fire p/a 5 yr legacy
# - fire p/a x ppt 1 yr legacy
# - fire p/a x ppt 2 yr legacy
# - fire p/a x ppt 3 yr legacy
# - fire p/a x ppt 4 yr legacy
# - fire p/a x ppt 5 yr legacy

# Other data notes (copied from the sbc_marss_model.R script):
# Following a discussion with John Melack, the following sites have both enough chemistry
# and precip data, as well as fires that occurred within that timeframe to be included in the
# MARSS analysis performed below.

# Arroyo Burro - AB00
# Atascadero - AT07
# Gaviota - GV01
# Arroyo Hondo - HO00
# Mission Creek (at Rocky Nook) - MC06
# Refugio - RG01
# Rattlesnake - RS02
# San Pedro - SP02

# The following site may be added later, but a longer precipitation record needs to be
# identified for it - Bell Canyon (BC02).

# See "data_raw/VCNP_sonde_site_codes_names.csv" for other possible site names used across other files (E.g., sonde data, GIS, etc.)
# "Redondo Creek" ~ "RED",
# "East Fork Jemez River" ~ "EFJ",
# "San Antonio - West" ~ "RSAW",
# "San Antonio Creek - Toledo" ~ "RSA",
# "Indios Creek" ~ "IND",
# "Indios Creek - Post Fire (Below Burn)" ~ "IND_BB",
# "Indios Creek - above burn" ~ "IND_AB",
# "Sulfur Creek" ~ "SULF",

#### Import data ####

dat1 = readRDS("data_working/marss_data_sb_vc_060622.rds")
