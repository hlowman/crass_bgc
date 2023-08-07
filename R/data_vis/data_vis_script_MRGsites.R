## Plot time series of water quality variables at sites along the Middle Rio Grande, New Mexico
## Created by: Betsy Summers
## Date: 9/13/21
## Note 1: Time stamp of imported data imported include daily averages but 15 min time interval are available.
## Note 2: Each site has 1 data file

### Load libraries ----
library(tidyverse)
library(timetk)
library(lubridate)
library(dplyr)
library(DataExplorer)
library(naniar)
library(viridis)

### Load site data ----
# Data files listed in upstream to downstream order
Cochiti_tbl <- read_csv("data_raw/Cochiti_WQ_2012-2017_DAvg.csv") # Bernalillo, NM
Bern_tbl <- read_csv("data_raw/Bernalillo_DailyAvg_2007-2019_WQ.csv") # Bernalillo, NM
Alameda_tbl <- read_csv("data_raw/Alameda_DailyAvg_2006-2019_WQ.csv") # Alameda site, NM
Bravo_tbl <- read_csv("data_raw/RioBravo_DailyAvg_2006-2019_WQ.csv") # Rio Bravo on Rio Grande site, NM
I25_tbl <- read_csv("data_raw/I25_DailyAvg_2006-2019_WQ.csv") # I-25 on Rio Grande site, NM


## Create figures of daily average water quality for each site
### Cochiti ----
# Located upstream of Albuquerque and upstream Cochiti reservoir
Cochiti_tbl <- Cochiti_tbl %>%
  select(day, DO.obs, SpCond, pH, tempwater, turbidity)

Cochiti_tbl_long <- Cochiti_tbl %>%
  pivot_longer(
    cols = DO.obs:turbidity, 
    names_to = "Variable",
    values_to = "value")

fig_Coch <- Cochiti_tbl_long %>% # use pivoted dataset
  mutate(day = as.Date(day, format="%m/%d/%Y")) %>%
  # mutate(Variable = factor(as.character(Variable))) %>% # recreate column
  ggplot(aes(x = day, y = value)) +
  geom_point()+ aes(color = Variable) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Date",
       y = "Units", 
       title = "Cochiti") +
  facet_grid(facets = vars(Variable), 
             scales = "free") +
  theme(legend.position = "none")+ 
  theme_bw()
fig_Coch

ggsave("figures/Fig_MRG_Cochiti_WQ_Davg.png", fig_panel,
       height = 10, width = 10, units = "in")

### Bernalillo ----  
# Located upstream of Albuquerque
Bern_tbl <- Bern_tbl %>%
  select(date, DO.obs, SpCond, pH, tempwater, turbidity.capped)

Bern_tbl_long <- Bern_tbl %>%
  pivot_longer(
    cols = DO.obs:turbidity.capped, 
    names_to = "Variable",
    values_to = "value")

fig_Bern <- Bern_tbl_long %>% # use pivoted dataset
  mutate(date = as.Date(date, format="%m/%d/%Y")) %>%
 # mutate(Variable = factor(as.character(Variable))) %>% # recreate column
  ggplot(aes(x = date, y = value)) +
  geom_point()+ aes(color = Variable) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Date",
       y = "Units", 
       title = "Bernalillo") +
  facet_grid(facets = vars(Variable), 
             scales = "free") +
  theme(legend.position = "none")+ 
  theme_bw()
fig_Bern

ggsave("figures/Fig_MRG_Bernalillo_WQ_Davg.png", fig_panel,
       height = 10, width = 10, units = "in")

### Alameda ----  
# Located in Albuquerque
Alameda_tbl <- read_csv("data_raw/Alameda_DailyAvg_2006-2019_WQ.csv") # Alameda site, NM

Alameda_tbl <- Alameda_tbl %>%
  select(day, DO.obs, SpCond, pH, tempwater, turbidity.capped)

Alameda_tbl_long <- Alameda_tbl %>%
  pivot_longer(
    cols = DO.obs:turbidity.capped, 
    names_to = "Variable",
    values_to = "value")

fig_Alameda <- Alameda_tbl_long %>% # use pivoted dataset
  mutate(day = as.Date(day, format="%m/%d/%Y")) %>%
  # mutate(Variable = factor(as.character(Variable))) %>% # recreate column
  ggplot(aes(x = day, y = value)) +
  geom_point()+ aes(color = Variable) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Date",
       y = "Units", 
       title = "Rio Grande - Alameda") +
  facet_grid(facets = vars(Variable), 
             scales = "free") +
  theme(legend.position = "none")+ 
  theme_bw()
fig_Alameda

ggsave("figures/Fig_MRG_Alameda_WQ_Davg.png", fig_panel,
       height = 10, width = 10, units = "in")

### Rio Bravo ---- 
# Located in Albuquerque
Bravo_tbl <- read_csv("data_raw/RioBravo_DailyAvg_2006-2019_WQ.csv") # Rio Bravo on Rio Grande site, NM

Bravo_tbl <- Bravo_tbl %>%
  select(day, DO.obs, SpCond, pH, tempwater, turbidity.capped)

Bravo_tbl_long <- Bravo_tbl %>%
  pivot_longer(
    cols = DO.obs:turbidity.capped, 
    names_to = "Variable",
    values_to = "value")

fig_Bravo <- Bravo_tbl_long %>% # use pivoted dataset
  mutate(day = as.Date(day, format="%m/%d/%Y")) %>%
  # mutate(Variable = factor(as.character(Variable))) %>% # recreate column
  ggplot(aes(x = day, y = value)) +
  geom_point()+ aes(color = Variable) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Date",
       y = "Units", 
       title = "Rio Grande - Rio Bravo") +
  facet_grid(facets = vars(Variable), 
             scales = "free") +
  theme(legend.position = "none")+ 
  theme_bw()
fig_Bravo

ggsave("figures/Fig_MRG_RioBravo_WQ_Davg.png", fig_panel,
       height = 8, width = 10, units = "in")

## I-25 site 
# Located at end of Albuquerque
I25_tbl <- read_csv("data_raw/I25_DailyAvg_2006-2019_WQ.csv") # I-25 on Rio Grande site, NM

I25_tbl <- I25_tbl %>%
  select(day, DO.obs, SpCond, pH, tempwater, turbidity)

I25_tbl_long <- I25_tbl %>%
  pivot_longer(
    cols = DO.obs:turbidity, 
    names_to = "Variable",
    values_to = "value")

fig_I25 <- I25_tbl_long %>% # use pivoted dataset
  mutate(day = as.Date(day, format="%m/%d/%Y")) %>%
  # mutate(Variable = factor(as.character(Variable))) %>% # recreate column
  ggplot(aes(x = day, y = value)) +
  geom_point()+ aes(color = Variable) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Date",
       y = "Units", 
       title = "Rio Grande - Rio Bravo") +
  facet_grid(facets = vars(Variable), 
             scales = "free") +
  theme(legend.position = "none")+ 
  theme_bw()
fig_I25

ggsave("figures/Fig_MRG_I25_WQ_Davg.png", fig_panel,
       height = 9, width = 10, units = "in")
