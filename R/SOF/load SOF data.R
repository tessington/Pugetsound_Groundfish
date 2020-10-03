rm(list = ls())

wd <- "~/Dropbox/Desktop/Active Projects/Puget Sound Baseline/Database Cleanup/Project Data/Groundfish/"

require(dplyr)
require(reshape2)
require(readxl)
setwd(wd)

# get tow data
tow_data <- read_xlsx("SOF/data/1 Coordinate Cleaning/Revised_location_table.xlsx")

# find the station IDs that are in WDFW data
wdfw.stations<-unique(tow_data$`WDFW StationID`)

# need to create new station IDs for all other stations.  Essentially assume each is unknown
na.index <- which(is.na(tow_data$`WDFW StationID`))

for (i in 1:length(na.index)){
  tow_data$`WDFW StationID`[na.index[i]] <- "Unid"
}

# get the catch data
catch_data <- read.csv("SOF/data/4 Catch Cleaning/catch_data.csv", header = T, stringsAsFactors = F)

# get the effort data
effort_data <- read.csv("SOF/data/5 Effort Cleaning/effort_data.csv", header = T, stringsAsFactors = F)

# I have all the data now.  Turn this into a wide table, with columns for HaulID, depth, duration (when available), and catch for each species.
# 
species.2.use <- c("English sole","Spiny dogfish", "Pacific whiting", "Big skate", "Pacific cod", "C-O sole", "Rock sole", "Flathead sole", "Sablefish", "Cabezon", "Spotted ratfish", "Pacific tomcod", "Blackbelly eelpout", "Pacific sanddab","Slender sole","Starry flounder", "Petrale sole","Dover sole", "Plainfin midshipman", "Rex sole","Longnose skate","Lingcod", "Quillback rockfish","Copper rockfish", "Greenstriped rockfish","Yelloweye rockfish","Canary rockfish","Redstripe rockfish","Bocaccio rockfish","Tiger rockfish","Vermillion rockfish","Rougheye rockfish","Brown rockfish","Walleye pollock", "Shiner perch")

red_catch_data <- catch_data %>%
  filter(species %in% species.2.use)

red_catch_data$species <- as.factor(red_catch_data$species)
red_catch_data$count <- as.numeric(red_catch_data$count)

wide.catch <- dcast(red_catch_data, haul_id~species, na.rm = T, fun.aggregate = sum, value.var = "count")
wide.pa <- dcast(red_catch_data, haul_id~species, fun.aggregate = length, value.var = "species")


# create a cool list, kind of like what Elizabeh created, to make it easier to combine stuff
# 
thedata.catch <- merge(tow_data, wide.catch, by = "haul_id" )
thedata.pa <- merge(tow_data, wide.pa, by = "haul_id" )
# now link to effort data tow ids.  I want to combine rows that are the same tow but have multiple catch entries


thedata.catch <- merge(thedata.catch, effort_data, by = "haul_id")
thedata.pa <- merge(thedata.pa, effort_data, by = "haul_id")

# need to do a bit of shifting samplying type from other category
gears.2.move <-c("Semi-Baloon Shrimp Trawl", "Gulf Shrimp Trawl", "Bottom Sampler", "Dredge")
move.index <- which(thedata.catch$Other %in% gears.2.move)
thedata.catch$`Sampling type`[move.index]<- thedata.catch$Other[move.index]

# now move over WDFW station ID codes into subregion codes, where applicable
move.index <- which(!thedata.catch$`WDFW StationID`=="Unid")
thedata.catch$`Sub region`[move.index] <- thedata.catch $`WDFW StationID`[move.index]
# create list 
species.data <- list()

for (i in 1:length(species.2.use)){
  tmp.species.data <- thedata.catch %>%
    select(haul_id, Year, Month, Day, Basin, Region, `Sub region`, `WDFW StationID`, `Sampling type`, duration_calc_hr, min_depth_fms, max_depth_fms)
  tmp.species.data$present <- pmin(unlist(thedata.pa %>% select(species.2.use[i])), 1)
  
  tmp.species.data$count <- unlist(thedata.catch %>% select(species.2.use[i]))
  
  # final thing, combine catch into summed catch, useful for when there are multiple 
  # entries per species and tow (e.g. when there are different counts by sex)
  species.data[[i]] <-tmp.species.data
}
names(species.data)<-species.2.use

save(file = "SOF/data/SOFdata.Rdata", species.data)

# some meta data analysis of species data.
# How many samples have depth information
# 
with.depth <- which(!is.na(species.data$`English sole`$min_depth_fms))
with.site <-  which(!is.na(species.data$`English sole`$`Sub region`))


length(with.depth)
length(with.site)
length(intersect(with.depth, with.site))
