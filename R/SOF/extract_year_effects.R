require(rstan)
source("R/project_functions.R")

species.list <-  load_species_names()
# list best fitting modelfor each species / model part
best.models <- get.best.models()
pmodel <- best.models$pmodel
ymodel <- best.models$ymodel

# create a base x matrix to get intercept right: assume may 15 sample date, average depth
model.tests <- get_models()

# create list to save
p.list <-y.list <-  rep(list(NA), 15)
names(p.list) <- names(y.list) <- species.list
load("data/SOF/SOFdata.Rdata")

for(spc in 1:15){
  
  species <- species.list[spc]
  # load species data
  extract.expression <- paste('thedata <- species.data$`', species, '`', sep = "")
  eval(parse(text = extract.expression))
  ### Modify Data ####
  
  thedata <- format.sof.data(thedata, rbind(c(85,54), c(10281,12323)))
  thedata$SamplingType <- relevel(thedata$SamplingType, "Otter trawl") # make Otter trawl the reference level (b/c it is most common)
  
  # extract best presence-absence parameters
  load(
    file = paste(
      "outputs/SOF/MCMC/",
      species,
      "_fitmcmc.Rdata",
      sep = ""
    )
  )
  p.fit <- y.fit <- model.stan
  beta_p <- rstan::extract(object = p.fit, pars = "beta_p")[[1]]
  nbeta_p <- ncol(beta_p)
  psi_p <-  rstan::extract(object = p.fit, pars = "psi_p")[[1]] # these are the year effects
  n.years <- ncol(psi_p)
  
  formula.txt <- model.tests[[pmodel[spc]]]
  tmp.x <- model.matrix(object = eval(parse(text =formula.txt)), data = thedata)
  beta_names <- colnames(tmp.x)
  # id columns for date
  date.cols <- grep("Date", beta_names)
  xp <- rep(0, nbeta_p);xp[1] <- 1; xp[date.cols] <- c(cos(2 * (5.5 - 1) / 12 * pi), sin(2 * (5.5 -
                                                                                                1) / 12 * pi))
  xp <- matrix(xp, nrow = nbeta_p, ncol = 1)
  
  
  y.fit <- model.stan
  psi_y <-  rstan::extract(object = y.fit, pars = "psi_y")[[1]]
  beta_y <- rstan::extract(object = y.fit, pars = "beta_y")[[1]]
  nbeta_y <- ncol(beta_y)
  
  
  # same as above to get the y model
  formula.txt <- model.tests[[ymodel[spc]]]
  tmp.x <- model.matrix(object = eval(parse(text =formula.txt)), data = thedata)
  beta_names <- colnames(tmp.x)
  # id columns for date
  date.cols <- grep("Date", beta_names)
  xy <- rep(0, nbeta_y);xy[1] <- 1; xy[date.cols] <- c(cos(2 * (5.5 - 1) / 12 * pi), sin(2 * (5.5 -
                                                                                                1) / 12 * pi))
  xy <- matrix(xy, nrow = nbeta_y, ncol = 1)
  #calculate linear terms of model
  logit.phat <-
    psi_p + (beta_p %*% xp) %*% matrix(1, nrow = 1, ncol = n.years)
  phat <- inv_logit(logit.phat)
  
  #log y = log(phat)* linear terms
  logyhat <-
    log(phat) + psi_y+ (beta_y %*% xy) %*% matrix(1, nrow = 1, ncol = n.years)
  # calculate ci on linear scales
  logit.output <-
    calc_ci(psi_p + (beta_p %*% xp) %*% matrix(1, nrow = 1, ncol = n.years), 0.5)
  logy.output <- calc_ci(logyhat, 0.5)
  
  # apply inverse link functions and save to lists
  p.list[[spc]] <- inv_logit(logit.output)
  y.list[[spc]] <- exp(logy.output)
}

save(file= "outputs/SOF/year_effects.Rdata", p.list, y.list)
