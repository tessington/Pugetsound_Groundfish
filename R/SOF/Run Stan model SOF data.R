# code to run all sorts of tweaks to the data to get it in the same format / structure as the WDFW data
# core data file is the .Rdata file called species.data that is a list by species that contains
# trawl ID information along with effort, catch, presence/ absence.  Note that catch may be 0 but presences =1
#  This just means count (catch) was not recorded
# catch by weight is also included.  Ultimately I may have to etstimate catch by number for those
# becuase I'm worried that the counts were not taken for large catches
rm(list = ls())

model.name <- "Stan/SOF_delta_ar1_model_v2.stan"
n_chains <- 3
iters = 3000

require(dplyr)
require(reshape2)
require(extraDistr)
require(readxl)
require(rstan)
require(shinystan)

# load project functions
source("R/project_functions.R")

load("data/SOF/SOFdata.Rdata")
# load data for any species.  Calculate mean depth and calculate sin and cosine of month
species.list <- load_species_names()
model.tests <- get_models()
best.models <- get.best.models()
pmodel <- best.models$pmodel
ymodel <- best.models$ymodel

# get basin, regions, and mean/sd depths from DFW data

dfw.data <- format_data("English sole", min.sample.size =3)
dfw.basin.list <- levels(dfw.data$tmpdata$Basin)
dfw.region.list <- levels(dfw.data$tmpdata$Region)
mean.depth <- mean(dfw.data$tmpdata$DEPTH)
sd.depth  <- sd(dfw.data$tmpdata$DEPTH)
mean.depth2 <- mean(dfw.data$tmpdata$DEPTH2)
sd.depth2  <- sd(dfw.data$tmpdata$DEPTH2)

depth.std <- rbind(c(mean.depth, sd.depth), c(mean.depth2, sd.depth2))

dfw.basin.names <-  levels(dfw.data$tmpdata$Basin)
dfw.region.names <-  levels(dfw.data$tmpdata$Region)



# cycle through species
for (spc in 1:15){
  species.2.use <- species.list[spc]
  
  pmodel.index <- pmodel[spc]
  ymodel.index <- ymodel[spc]
  
  # get information on model covariates (beta) to use and set up names 
  p.model.text <- paste(model.tests[[pmodel[spc]]], "+SamplingType",sep = "")
  y.model.text <- paste(model.tests[[ymodel[spc]]], "+SamplingType",sep = "")
  dfw_beta_p_names <- colnames(model.matrix(object =eval(parse(text = model.tests[pmodel[spc]])), data = dfw.data$tmpdata))
  dfw_beta_y_names <- colnames(model.matrix(object =eval(parse(text = model.tests[ymodel[spc]])), data = dfw.data$tmpdata))
  all.dfw.names = list(beta_p =dfw_beta_p_names,
                       beta_y=dfw_beta_y_names,
                       basin=dfw.basin.names,
                       region=dfw.region.names)
  
 
  
  # load species data
  extract.expression <- paste('thedata <- species.data$`', species.2.use, '`', sep = "")
  eval(parse(text = extract.expression))
  ### Modify Data ####
  
  thedata <- format.sof.data(thedata, depth.std)
  thedata$SamplingType <- relevel(thedata$SamplingType, "Otter trawl") # make Otter trawl the reference level (b/c it is most common)
  
  
  # now get priors
  priors <- assign.priors(thedata, species.2.use,pmodel.index, ymodel.index, all.dfw.names)
  #### Set up STAN run ####
  xp <- model.matrix(object = eval(parse(
    text = p.model.text)), data = thedata)
  xy <- model.matrix(object = eval(parse(
    text = y.model.text)), data = thedata)
  t <- model.matrix(~-1+Year, data = thedata)
  b <- model.matrix(~-1+Basin, data = thedata)
  r<- model.matrix(~-1 + Region, data = thedata) 
  
  ndata <- nrow(xp)
  nbetap <- ncol(xp)
  nbetay <- ncol(xy)
  ngamma <- ncol(b)
  ntheta <- ncol(r)
  nyear <- ncol(t)
  
  region.list <- levels(thedata$Region)
  # for each location, return basin index
  basin.list <- levels(thedata$Basin)
  basin.index <- rep(NA, length(region.list))
  
  for (i in 1:length(region.list)) {
    tmp.region.index <- which(thedata$Region == region.list[i])[1]
    basin.index[i] <- which(basin.list == thedata$Basin[tmp.region.index])
  }
  
  y <- thedata$count 
  pa <- thedata$present
  present_index <- which(y>0)
  npresent <- length(present_index) #note, not that same as which  pa ==1, because some catch amounts not recorded
  stan.data <- list(
    ndata = ndata,
    nbetap = nbetap,
    nbetay = nbetay,
    nyear = nyear,
    ngamma = ngamma,
    ntheta = ntheta,
    npresent = npresent,
    xp = xp,
    xy = xy,
    b = b,
    r = r,
    t = t,
    y = as.vector(y),
    pa = pa,
    present_index = present_index,
    basin_index = basin.index,
    beta_p_prior = priors$beta_p,
    beta_y_prior = priors$beta_y,
    gamma_p_raw_prior = priors$gamma_p_raw,
    gamma_y_raw_prior = priors$gamma_y_raw,
    theta_p_raw_prior = priors$theta_p_raw,
    theta_y_raw_prior = priors$theta_y_raw,
    shape_prior = priors$shape,
    logitrho_p_prior = priors$logitrho_p,
    logitrho_y_prior = priors$logitrho_y,
    sigma_psi_p_prior = priors$sigma_psi_p,
    sigma_psi_y_prior = priors$sigma_psi_y,
    sigma_basin_p_prior = priors$sigma_basin_p,
    sigma_basin_y_prior = priors$sigma_basin_y,
    sigma_region_p_prior = priors$sigma_region_p,
    sigma_region_y_prior = priors$sigma_region_y
  )
  
  params = c(
    "beta_p",
    "beta_y",
    "psi_p",
    "psi_y",
    "shape",
    "sigma_psi_p",
    "sigma_psi_y",
    "rho_p",
    "rho_y",
    "logitrho_p",
    "logitrho_y",
    "gamma_pa",
    "gamma_y",
    "gamma_p_raw",
    "gamma_y_raw",
    "theta_p",
    "theta_y",
    "theta_p_raw",
    "theta_y_raw",
    "sigma_basin_p",
    "sigma_basin_y",
    "sigma_region_p",
    "sigma_region_y",
    'log_lik_pa',
    'log_lik_y',
    'logit_phat',
    'log_yhat'
  )
  
  
  init_ll <- lapply(1:n_chains, function(id)
    initf3(chain_id = id, nbetap = nbetap, nbetay = nbetay, ngamma = ngamma, ntheta = ntheta, nyear = nyear))
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  
  #### run STAN model #####
  # write stan model
  m <- stan_model(model.name)
  model.stan <- sampling(m,
                         data = stan.data,
                         iter = iters,
                         pars = params,
                         warmup = floor(iters / 2),
                         chains = n_chains,
                         thin = 1,
                         algorithm = 'NUTS',
                         init = init_ll,
                         verbose = FALSE,
                         control = list(adapt_engaged = TRUE, adapt_delta = 0.99)
  )
  
  
  ## Save STAN model output
  save.text <- paste("outputs/SOF/MCMC/",species.2.use, "_fitmcmc.Rdata", sep = "")
  save(file = save.text, model.stan)
}

# run simple diagnostics on each fitted model, testing for divergences

count_diverge <-rep(NA, 15)

# Run diagnostics - number of divergences
for (spc in 1:15) {
  
  species <- species.list[spc]
  # diagnostics for each model

    filename <- paste("outputs/SOF/MCMC/",species,"_fitmcmc.Rdata", sep = "")
    load(file= filename)
    sampler_params <- get_sampler_params(model.stan, inc_warmup = FALSE)
    count_diverge[spc]<-sum( sapply(sampler_params, function(x) sum(x[, "divergent__"])))
}
