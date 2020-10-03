### Need to evaluate model selection for Stan model and possibly diagnose high level of divergences that seem to happen with the presence / absenc
### Part of the model
rm(list = ls())
require(rstan)
require(loo)
require(dplyr)
source("R/project_functions.R")
species.names <- load_species_names()
model.tests <- get_models()
best.models <- get.best.models()


for (spc in 1:length(species.names)){
# Step one, load data for one species,  just to get the data frame set right
species <- species.names[spc]

data <- format_data(species, min.sample.size = 3)
tmpdata <- data$tmpdata
tmpdata$y <- with(tmpdata, C / EFFORT)
data$tmpdata <- tmpdata

# Run Stan for Model Selection --------------------------------------------
species <- species.names[spc]
### Extract posterior medians ###


model.name <- "Stan/delta_ar1_model_v2.stan"


  model.stan <-
    run_stan_model(
      data,
      formula.text = model.tests[[best.models$pmodel[spc]]],
      model_name = model.name,
      iters = 2000
    )
  filename <- paste("outputs/DFW/modelselection/fit_",species, "_", best.models$pmodel[spc], "mcmc.Rdata", sep = "")
  save(file= filename, model.stan)
  # run the following if the best model for y is different from the best model for p
  if (!best.models$pmodel[spc]==best.models$ymodel[spc]) {
  model.stan <-
    run_stan_model(
      data,
      formula.text = model.tests[[best.models$ymodel[spc]]],
      model_name = model.name,
      iters = 2000
    )
  filename <- paste("outputs/DFW/modelselection/fit_",species, "_", best.models$ymodel[spc], "mcmc.Rdata", sep = "")
  save(file= filename, model.stan)
  }
}