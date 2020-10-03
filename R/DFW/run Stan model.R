### Need to evaluate model selection for Stan model and possibly diagnose high level of divergences that seem to happen with the presence / absenc
### Part of the model
rm(list = ls())
require(extraDistr)
require(rstan)
require(loo)
require(dplyr)
source("R/project_functions.R")
species.names <- load_species_names()

for (spc in 1:length(species.names)){
# Step one, load data for one species,  just to get the data frame set right
species <- species.names[spc]

data <- format_data(species, min.sample.size = 3)
tmpdata <- data$tmpdata
tmpdata$y <- with(tmpdata, C / EFFORT)
data$tmpdata <- tmpdata

# Run Stan for Model Selection --------------------------------------------
model.tests <- list()
species <- species.names[spc]
### Extract posterior medians ###

model.tests <- get_models()

model.name <- "Stan/delta_ar1_model_v2.stan"


for (j in 1:5){
  model.stan <-
    run_stan_model(
      data,
      formula.text = model.tests[[j]],
      model_name = model.name,
      iters = 2000
    )
  filename <- paste("outputs/DFW/modelselection/fit_",species, "_", j, "mcmc.Rdata", sep = "")
  save(file= filename, model.stan)
}
}