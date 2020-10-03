
library(readxl)
library(tidyverse)
library(sf)
library(ggplot2)
#Load the data.
tow = read_excel(path = "data/DFW/wdfw_vTow.xlsx")
effort = read_excel("data/DFW/wdfw_2_ Calc_Area_Swept.xlsx")
biomass = read_excel("data/DFW/wdfw_3_ Calc_Weights_Nums.xlsx")
#region = read_csv(paste0(localfile, "data/wdfw_spatial_regions.txt"))
#```


### Preliminary tidying 
###Convert Lat/Lon into a single column each by using decimal degrees. Just use the first set of locations. 
###```{r}
tow$LAT = tow$LatDeg1 + tow$LatMin1/60
tow$LON = -(tow$LonDeg1 + tow$LonMin1/60)
###```

###Convert Lat/Lon to UTM (might come in handy for later?)
###```{r}
tmp = tow %>% st_as_sf(coords = c("LON", "LAT"), crs = 4326) %>%
  st_transform(crs = 32610)

tow$UTM_E = st_coordinates(tmp)[,1]
tow$UTM_N = st_coordinates(tmp)[,2]
#```

#Remove unsatisfactory tows (0 and 1 are good) and ROV trawls.  
#```{r}
# Delete "unsatisfactory" tows
tow = filter(tow, PerformanceID <= 1) %>%
      filter(SurveyType != "ROVTRAWL")
#```

#Remove Canadian sampling.
#```{r}
tow = filter(tow, RegionID != "GC")
#```

### Add new location designation
#Based on the residuals from intial fits, split up some of the areas.  Replace #the RegionID with Zone. 
#```{r}
subbasins_new = st_read("data/DFW/shapefiles/subbasins_custom.shp") %>%
  st_transform(4326)

# Check map
test = st_as_sf(tow, coords = c("LON", "LAT"), crs = 4326) %>%
  st_join(subbasins_new)

#ggplot() + 
#  geom_sf(data = subbasins_new, aes(color = ZONE, alpha = 0.5, fill = NULL)) + 
#  geom_sf(data = test, aes(color = ZONE, fill = ZONE), alpha = 0.5)

# Do merge
tow$RegionID = st_as_sf(tow, coords = c("LON", "LAT"), crs = 4326) %>%
   st_join(subbasins_new) %>%
   pull(Name)
 
# Remove NAs
tow = tow %>% filter(!is.na(RegionID))
#```




### Create dataframe with shared tow-level information
#We want a data frame `dat` where each row is a unique haul, with the following columns:

#* Tow$ID: tow$TowID (unique ID for each tow)
#* LOC: tow$RegionID (WDFW sub basin, categorical)
# * YEAR: tow$SurveyYear
# * DATE: tow$HaulDate
# * LAT: latitude, calculate decimal degrees from `tow$LatDeg1` and `tow$LatMin2`
# * LON: longitude, calculate decimal degrees from `tow$LonDeg1` and `tow$LonMin2`
# * DEPTH: average depth (in m), convert from fathoms `tow$`AvgDepth(Fa)``
# * logDEPTH: log depth centered
# * EFFORT: WDFW-calculated area swept ($m^2$) `effort$`AreaSwept(m2)``
# * BIOMASS: WDFW-calculated total catch in kg `biomass$`TotalWt(kg)``
# * NUMBERS: WDFW-calculated total catch in kg `biomass$`TotalInd``

# Create empty data frame where each row will be a unique haul/tow. 
# ```{r}
cols = c("TowID","LOC", "YEAR","DATE","LAT","LON", "UTM_N", "UTM_E", "DEPTH", "logDEPTH", "EFFORT", "BIOMASS","NUMBERS")

tows = unique(tow$TowID)
ntows = length(tows)
dat = matrix(NA, nrow = ntows, ncol = length(cols))
dat = as.data.frame(dat)
names(dat) = cols
#```

#Fill in columns from `tow`.  
#```{r}
dat$TowID = tow$TowID
dat$LOC = tow$RegionID
dat$YEAR = tow$SurveyYear
dat$DATE = as.Date(tow$HaulDate,format = "%Y%m%d")
dat$LAT = tow$LAT
dat$LON = tow$LON
dat$UTM_N = tow$UTM_N
dat$UTM_E = tow$UTM_E
dat$DEPTH = tow$`AvgDepth(Fa)`*1.8288 # convert fathoms to meters
dat$logDEPTH = log(dat$DEPTH) - mean(log(dat$DEPTH))
#```

#Add the calculated effort data. 
#```{r}
dat$EFFORT = left_join(dat, effort, by = "TowID") %>% pull(`AreaSwept(m2)`)
dat$logEFFORT = log(dat$EFFORT) - mean(log(dat$EFFORT))
#```


### Add the species-specific biomass data
#The list of species we're interested in:
#```{r}
focal_species = c("English sole", "Spiny dogfish", "Spotted ratfish", "Speckled sanddab", "Pacific cod", "Walleye pollock", "Pacific whiting (hake)", "Pacific tomcod", "Shiner perch", "Starry flounder", "Pacific sanddab","Plainfin midshipman","Blackbelly eelpout", "Lingcod","Longnose skate", "Big skate", "Rock sole")
#```

#Make a list of dataframes. 
#```{r}
species_data = rep(list(dat), length(focal_species))
names(species_data) = focal_species
#```

# Add the catch of each species. 
# ```{r}
for(i in seq(species_data)){
  species_data[[i]]$BIOMASS = filter(biomass, SpeciesID == focal_species[i]) %>%
           select(c("TowID", "TotalWt(kg)")) %>%
           right_join(species_data[[i]], by = "TowID") %>%
           pull(`TotalWt(kg)`)
  
  species_data[[i]]$BIOMASS = ifelse(is.na(species_data[[i]]$BIOMASS), 0, species_data[[i]]$BIOMASS)
  species_data[[i]]$NUMBERS = filter(biomass, SpeciesID == focal_species[i]) %>%
    select(c("TowID", "TotalInd")) %>%
    right_join(species_data[[i]], by = "TowID") %>%
    pull(`TotalInd`)
  
  species_data[[i]]$NUMBERS = ifelse(is.na(species_data[[i]]$NUMBERS), 0, species_data[[i]]$NUMBERS)
  species_data[[i]]$PRESENT = ifelse(species_data[[i]]$BIOMASS > 0, 1, 0)
  species_data[[i]]$CPUE = species_data[[i]]$BIOMASS/species_data[[i]]$EFFORT
}


### Add Rock Sole, sum all species
i <- length(focal_species)
species_data[[i]]$BIOMASS = filter(biomass, SpeciesID %in% c("Rock sole uniden.","Northern rock sole", "Southern rock sole"))  %>%
  group_by(TowID) %>%
  summarize(`TotalWt(kg)` = sum(`TotalWt(kg)`)) %>%
  right_join(species_data[[length(focal_species)]], by = "TowID") %>%
  pull(`TotalWt(kg)`)
species_data[[i]]$BIOMASS = ifelse(is.na(species_data[[i]]$BIOMASS), 0, species_data[[i]]$BIOMASS)

species_data[[i]]$NUMBERS = filter(biomass, SpeciesID %in% c("Rock sole uniden.","Northern rock sole", "Southern rock sole"))  %>%
  group_by(TowID) %>%
  summarize(`TotalInd` = sum(`TotalInd`)) %>%
  right_join(species_data[[i]], by = "TowID") %>%
  pull(`TotalInd`)

species_data[[i]]$NUMBERS = ifelse(is.na(species_data[[i]]$NUMBERS), 0, species_data[[i]]$NUMBERS)



species_data[[i]]$PRESENT = ifelse(species_data[[i]]$BIOMASS > 0, 1, 0)
species_data[[i]]$CPUE = species_data[[i]]$BIOMASS/species_data[[i]]$EFFORT

### Save the data
#```{r}
saveRDS(species_data, file = "outputs/DFW/data_formatted/wdfw_processed.RData")
readRDS("outputs/DFW/data_formatted/wdfw_processed.RData")
#```
# 
# Comment: I can see how having all the steps in separate scripts (or functions), but calling all from one could be usefull.  Would allow me to change options easily, like which subbasin set I was using. 