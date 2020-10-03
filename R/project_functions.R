### Project source files

load_species_names<- function() c("English sole",
                       "Spiny dogfish",
                       "Spotted ratfish",
                       "Pacific cod",
                       "Pacific whiting",
                       "Pacific sanddab",
                       "Pacific tomcod",
                       "Walleye pollock",
                       "Plainfin midshipman",
                       "Blackbelly eelpout",
                       "Lingcod",
                       "Shiner perch",
                       "Longnose skate",
                       "Big skate",
                       "Rock sole")

format_data <- function(species.2.use, reduce.data = F, min.sample.size =1) {
  require(dplyr)
  require(readxl)
  # load data
  # 
  
  species_data = readRDS("outputs/DFW/data_formatted/wdfw_processed.RData")
  
  # load tow data to get stationIDs
  tow_data <- read.csv("outputs/DFW/wdfw_vTow_final.csv", header = T)
  # extract out species that you want
  extract.expression <- paste('thedata <- species_data$`', species.2.use, '`', sep = "")
  eval(parse(text = extract.expression))
  
  
  # remove first few years of the data
  thedata <- thedata %>%
    filter(YEAR >=1990)
  
  # code to add tow ID to thedata
  tow.indices <- match(thedata$TowID,  tow_data$TowID, nomatch = NA)
  thedata <- thedata[which(thedata$TowID %in% tow_data$TowID),] # remove the tows in thedata that do not appear in tow_data.  These are largely georgia basin tows that were removed
  tow.indices <- tow.indices[which(!is.na(tow.indices))] # get rid of the NAs (no matches) in the tow index
  thedata$StationID <- as.character(tow_data$StationID[tow.indices])
  
  
  thedata$Basin <- as.character(tow_data$Basin[tow.indices])
  thedata$Region <- as.character(tow_data$Region[tow.indices])
  thedata$Month <- as.numeric(tow_data$Month[tow.indices])
  thedata$Day <- as.numeric(tow_data$Day[tow.indices])
  
  
  # Try this: remove regions if there are too few data (fewer than designated # of samples)
  region.summary <- thedata %>%
    group_by(Region) %>%
    summarise(N = n())
  regions.2.rem <- which(region.summary$N <= min.sample.size)
  # Remove these from the data
  thedata <- thedata %>%
    filter(!Region%in% region.summary$Region[regions.2.rem])
  
  
  # 
  ######## Little bit here to reduce the data set by some amount for diagnostics
  if (reduce.data) {
    ndata <- nrow(thedata)
    sample.id <- sample(1:ndata, floor(ndata/2))
    
    thedata <- thedata[sample.id,]
  }

  thedata$C <- thedata$NUMBERS
  
  tmpdata <- thedata %>%
    mutate(EFFORT = EFFORT / 10000, DEPTH2 = DEPTH * DEPTH) %>%
    select(YEAR, DEPTH,DEPTH2,PRESENT,CPUE, StationID, Basin, Region, EFFORT, C, Month, Day)
  
  tmpdata$Basin <- factor(tmpdata$Basin)
  tmpdata$Region <- factor(tmpdata$Region)
  tmpdata$StationID <- factor(tmpdata$StationID)
  tmpdata$YEAR <- factor(tmpdata$YEAR, levels = 1990:2016)
  #scale year and depth
  tmpdata <- tmpdata %>% 
    mutate(sDEPTH = scale(DEPTH), sDEPTH2 = scale(DEPTH2))
  wdfw.depth <- c(mean(tmpdata$DEPTH), sd(tmpdata$DEPTH))
  wdfw.depth2 <- c(mean(tmpdata$DEPTH2), sd(tmpdata$DEPTH2))
  # convert scaled variables to numeric
  
  tmpdata$sDEPTH <- as.numeric(tmpdata$sDEPTH)
  tmpdata$sDEPTH2 <- as.numeric(tmpdata$sDEPTH2)
  
  tmpdata$cosDate <- cos(2*(tmpdata$Month -1 + (tmpdata$Day-1)/30) / 12 * pi)
  tmpdata$sinDate <- sin(2*(tmpdata$Month -1 + (tmpdata$Day-1)/30) / 12 * pi)
  # For each site, return location index.
  station.list <- levels(tmpdata$StationID)
  region.list <- levels(tmpdata$Region)
  region.index <- rep(NA, length(station.list))
  
  
  for (i in 1:length(station.list)) {
    station.index <- which(tmpdata$StationID == station.list[i])[1]
    region.index[i] <- which(region.list == tmpdata$Region[station.index])
  }
  
  # for each location, return basin index
  basin.list <- levels(tmpdata$Basin)
  basin.index <- rep(NA, length(region.list))
  
  for (i in 1:length(region.list)) {
    tmp.region.index <- which(tmpdata$Region == region.list[i])[1]
    basin.index[i] <- which(basin.list == tmpdata$Basin[tmp.region.index])
  }
  
  #this is for diagnostic run, redo region index to be a site-level basin index
  region.index <- rep(NA, length(station.list))
  
  for (i in 1:length(station.list)) {
    station.index <- which(tmpdata$StationID == station.list[i])[1]
    region.index[i] <- which(basin.list == tmpdata$Basin[station.index])
  }
  return(list(tmpdata=tmpdata, basin.index = basin.index, region.index = region.index))
}

inv_logit <- function(x) 1 / (1 + exp(-x))
model.select <- function(species, models.2.use = 1:4) {
  loo.outputs_pa <- list()
  loo.outputs_y <- list()
  lpd_point_pa <- c()
  for (i in 1:length(models.2.use)) {
    model <- models.2.use[i]
    ## load fitted model
    filename <- paste("outputs/DFW/modelselection/fit_",species, "_", model, "mcmc.Rdata", sep = "")
    
    load(filename)
    # get log likelihoods and effective sample size
    log_lik1 <- extract_log_lik(model.stan, parameter_name = "log_lik_pa", merge_chains = FALSE)
    
    rel_n_eff <- relative_eff(exp(log_lik1))
    # LOO cross validation
    x<-loo(log_lik1, r_eff = rel_n_eff, cores = 2, K = 10)
    # save to list and then extract elpd, save in lpd_point object
    loo.outputs_pa[[i]]<-x
    if (i==1) { 
    lpd_point_tmp <-   loo.outputs_pa[[i]]$pointwise[,"elpd_loo"]
    ndata <- length(lpd_point_tmp)
    lpd_point_pa<- array(NA, c(ndata, length(models.2.use)))
    lpd_point_pa[,i] <-lpd_point_tmp
    }
    if (i>1)  lpd_point_pa[,i] <- loo.outputs_pa[[i]]$pointwise[,"elpd_loo"]
  
    # repeat for density:
    log_lik1 <- extract_log_lik(model.stan, parameter_name = "log_lik_y", merge_chains = FALSE)
    rel_n_eff <- relative_eff(exp(log_lik1))
    x<-loo(log_lik1, r_eff = rel_n_eff, cores = 2, K = 10)
    loo.outputs_y[[i]]<-x
    if (i==1) { 
      lpd_point_tmp <-   loo.outputs_y[[i]]$pointwise[,"elpd_loo"]
      ndata <- length(lpd_point_tmp)
      lpd_point_y<- array(NA, c(ndata, length(models.2.use)))
      lpd_point_y[,i] <-lpd_point_tmp
    }
    if (i>1)  lpd_point_y[,i] <- loo.outputs_y[[i]]$pointwise[,"elpd_loo"]
    
    # flag signifies whether function has gone through loop yet
  }
  # get bayesian model weights
  stacking_wts_pa <- stacking_weights(lpd_point_pa)
  stacking_wts_y <- stacking_weights(lpd_point_y)
  return(list(pa = stacking_wts_pa, y = stacking_wts_y))
}  

initf2 <- function(chain_id = 1, nyear, nbeta, ngamma, ntheta) {
  list(
    beta_p = as.array(runif(nbeta, -1.0, 1.0)),
    beta_y = as.array(runif(nbeta, -1.0, 1.0)),
    shape = runif(n = 1, 1, 10),
    gamma_p_raw = rnorm(ngamma, 0, 1),
    gamma_y_raw = rnorm(ngamma, 0, 1),
    theta_p_raw = rnorm(ntheta, 0,1),
    theta_y_raw = rnorm(ntheta,0,1),
    logitrho_p = runif(1, -1, 1),
    logitrho_y = runif(1, -1, 1),
    sigma_psi_p = runif(1, 0, 1),
    sigma_psi_y = runif(1, 0, 1),
    eps_p_raw = rnorm(nyear,0, 1 ),
    eps_y_raw = rnorm(nyear,0, 1 ),
    sigma_basin_p = runif(n = 1, 0, 2),
    sigma_basin_y = runif(n = 1, 0.0, 2),
    sigma_region_p = runif(n = ngamma, 0, 1.5),
    sigma_region_y = runif(n = ngamma, 0, 1.5)
  )
}

## Divergent Graphical Analysis: Parallel chains
parallel_chain_plot <- function(model.stan, pars.2.use) {
  #pars.2.use <- -c(1:15,31:55,81, 82, 84,85:90,97:131, 167, 169:173)
  
  sampler_params <- get_sampler_params(model.stan, inc_warmup = FALSE)
  mcmc <- as.array(model.stan)
  
  parnames <- dimnames(mcmc)$parameters[pars.2.use]
  
  # standardize using normal cdf for easier visualization
  
  par(mar = c(3, 5,1,5), mfrow = c(3,1))
  
  for (chain in 1:3) {
    divergent.index <- which(sampler_params[[chain]][,5]==1)
    non.divergent.index <- which(sampler_params[[chain]][,5]==0)
    
    ylims <- c(-10, 10)
    plot(1:length(parnames), c(),
         type = "n", 
         col = "gray70", 
         axes = F, 
         xlab = "", 
         ylab = "", 
         ylim = ylims)
    
    
    axis(1, at = 1:length(parnames), labels = parnames, las = 2)
    axis(2, las =1)
    for (i in 1:length(non.divergent.index)) {
      lines(1:length(parnames), mcmc[non.divergent.index[i], chain, pars.2.use],
            col = "gray 70"
      )
    }
    for (i in 1:length(divergent.index)) {
      lines(1:length(parnames), mcmc[divergent.index[i], chain, pars.2.use],
            col = "red"
      )
    }
  }
  # standardize using normal cdf for easier visualization
  par.means <- apply(X=mcmc[,,],MARGIN = 3, FUN= mean)
  par.sigmas <- apply(X=mcmc[,,],MARGIN = 3, FUN= sd)
  
  std.mcmc <- pnorm(mcmc[,,], mean = par.means, sd = par.sigmas)
  par(mar = c(7, 5,1,5), mfrow = c(3,1))
  for (chain in 1:3) {
    std.chain <- matrix( NA, nrow = length(mcmc[,chain,1]), ncol = length(pars.2.use))
    for (i in 1:length(pars.2.use)) {
      # by parameter, standardize based on cumulative normal 
      std.chain[,i] <- pnorm(mcmc[,chain,pars.2.use[i]], 
                             mean = par.means[pars.2.use[i]],
                             sd = par.sigmas[pars.2.use[i]])
    }
    
    #  std.chain <- pnorm(mcmc[,chain,], mean = par.means, sd = par.sigmas)
    
    divergent.index <- which(sampler_params[[chain]][,5]==1)
    non.divergent.index <- which(sampler_params[[chain]][,5]==0)
    
    ylims <- c(-10, 10)
    plot(1:length(parnames), c(),
         type = "n", 
         col = "gray70", 
         axes = F, 
         xlab = "", 
         ylab = "", 
         ylim = c(0,1))
    
    
    axis(1, at = 1:length(parnames), labels = parnames, las = 2)
    axis(2, las =1)
    for (i in 1:length(non.divergent.index)) {
      lines(1:length(parnames), std.chain[non.divergent.index[i], ],
            col = "gray 70"
      )
    }
    for (i in 1:length(divergent.index)) {
      lines(1:length(parnames), std.chain[divergent.index[i], ],
            col = "red"
      )
    }
  }
}


run_stan_model <- function(data, formula.text, model_name, n_chains = 3, iters = 2000) {
  tmpdata <- data$tmpdata
  ndata <- nrow(tmpdata)
  nyear <-ncol(model.matrix(~-1 + YEAR, data = tmpdata)) 
  x <- model.matrix(object = eval(parse(text =formula.text)), data = tmpdata)
  t <- model.matrix(object = ~ -1 + YEAR, data = tmpdata)
 
  nbeta <- ncol(x)
  b <- model.matrix(~-1+Basin, data = tmpdata)
  r<- model.matrix(~-1 + Region, data = tmpdata) 
  ndata <- nrow(x)
  ngamma <- ncol(b)
  ntheta <- ncol(r)
  
  pa <- tmpdata$PRESENT
  y <- tmpdata$y
  
  npresent <- sum(pa)
  present_index <- which(pa ==1)
  basin_index <- data$basin.index
  
  
  stan.data <- list(
    ndata = ndata,
    nbeta = nbeta,
    nyear = nyear,
    ngamma = ngamma,
    ntheta = ntheta,
    x = x,
    t = t,
    b = b,
    r = r,
    y = as.vector(y),
    pa = pa,
    npresent = npresent,
    present_index = present_index,
    basin_index = basin_index
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
  
  # get initial values for mcmc
  init_ll <- lapply(1:n_chains, function(id)
    initf2(chain_id = id, nyear = nyear, nbeta = nbeta, ngamma = ngamma, ntheta = ntheta))
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # run STAN model
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
    control = list(adapt_engaged = TRUE, adapt_delta = 0.9)
  )
  return(model.stan)
}


# Model Diagnostics -------------------------------------------------------
model.diagnostic <- function(model.stan) {
sampler_params <- get_sampler_params(model.stan, inc_warmup = FALSE)
n_divergent <- sum(sampler_params[[1]][,5])+ sum(sampler_params[[2]][,5]) + sum(sampler_params[[3]][,5])
return(n_divergent)
}

# Calculate Credible Interval ---------------------------------------------
calc_ci <- function(fit, per.2.use = 0.8) {
  require(KernSmooth)
n.years <- ncol(fit)
per <- c( (1- per.2.use) / 2, 0.5, 0.5 + 0.5 * per.2.use)
# a place to save median and upper and lower bound
output <- matrix(NA, nrow = n.years, ncol = 3)
colnames(output) <- c("lower","median","upper")

# loop through columns, fit smoother, and report percentiles
for (i in 1:n.years) {
 
  t.smooth <- bkde(fit[,i],
                   kernel = "normal",
                   canonical = FALSE)
  delta.t <- t.smooth$x[2] - t.smooth$x[1]
  cumprob <- cumsum(t.smooth$y)* delta.t
  
  # find lower percentile

  output[i,] <- approx(cumprob, t.smooth$x, 
                          xout = per,
                          method = "linear")$y
}
return(output)
}


# Generate list of alternative models -------------------------------------

get_models <- function() {
  model.tests <-  list()
  model.tests[[1]] <- "~ sDEPTH*Basin - Basin + sDEPTH2*Basin - Basin + cosDate + sinDate "
  model.tests[[2]] <- "~ sDEPTH*Basin - Basin +  cosDate + sinDate "
  model.tests[[3]] <- "~ sDEPTH + sDEPTH2 + cosDate + sinDate "
  model.tests[[4]] <- "~ sDEPTH + cosDate + sinDate"
  model.tests[[5]] <- "~cosDate + sinDate"
  return(model.tests)
}


# Get best fitting model indices ------------------------------------------
get.best.models <- function() {
  pmodel <-  c(1,
               3,
               1,
               3,
               1,
               1,
               1,
               1,
               1,
               1,
               5,
               1,
               3,
               1,
               1)
  
  ymodel <- c(1,
              1,
              1,
              1,
              1,
              2,
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              2,
              1)
  return(list(pmodel = pmodel, ymodel = ymodel))
}

# Format SOF data ---------------------------------------------------------

format.sof.data <- function(thedata, depth.std) {
 
  thedata$Month <- as.numeric(thedata$Month)
  thedata$Day <- as.numeric(thedata$Day)
  thedata$cosDate <- cos(2*(thedata$Month -1 + (thedata$Day-1)/30) / 12 * pi)
  thedata$sinDate <- sin(2*(thedata$Month -1 + (thedata$Day-1)/30) / 12 * pi)
  thedata$Depth <- (thedata$min_depth_fms + thedata$max_depth_fms) * 1.88 / 2
  # keep only those tows with bottom contact gear and with depth and Month data.
  thedata <- thedata %>% 
    filter(`Sampling type` %in% c("Otter trawl", "Bottom trawl","Gulf Shrimp Trawl", "Semi-Balloon Shrimp Trawl"), 
           Depth >0, Month >0)

  
  # turn basins, regions into factors
  thedata$Basin <- as.factor(thedata$Basin)
  thedata$Region <- as.factor(thedata$Region)
  thedata$SamplingType <- as.factor(thedata$`Sampling type`)
  thedata$Year <- factor(thedata$Year, levels = 1948:1977)
  
  # scale depth the same way it was scaled in WDFW run
  thedata$Depth2 <- thedata$Depth^2
  thedata$sDEPTH <- (thedata$Depth - depth.std[1,1])/depth.std[1,2]
  thedata$sDEPTH2 <- (thedata$Depth2 - depth.std[2,1])/ depth.std[2,2]
  return(thedata)
}

assign.priors <- function(thedata, species.2.use,pmodel.index, ymodel.index, all.dfw.names) {
 # subfunctions
  # fuction to return matches is two  vectors in both directions
  match.both <- function(x,y) {
    x.in.y <- which(x %in% y)
    y.in.x <- which(y %in% x)
    output <- cbind(x.in.y, y.in.x)
    colnames(output) <- c("xiny","yinx")
    return(output)
  }

# Assign priors -----------------------------------------------------------

  
  # fit priors
  fit.priors <- function(fit, sof.names= NULL, dfw.names=NULL, parname, fit.type = "cauchy", default.prior) {
    # function to fit caucy parameters
    fit.cauchy <- function(pars, dat) sum(-dcauchy(dat, pars[1], pars[2], log = T))
    fitpars <- rstan::extract(fit, pars = parname)[[1]]
    if(!is.null(sof.names)) {
    indices <- match.both(sof.names, dfw.names)
    prior <- matrix(default.prior,nrow = length(sof.names), ncol = 2, byrow = TRUE)
    for (i in 1:nrow(indices)) {
      if (fit.type=="cauchy") prior[indices[i,1],1:2] <- optim(par = c(0, 2.5), fn = fit.cauchy, dat = fitpars[,indices[i,2]])$par
      if (fit.type=="normal") prior[indices[i,1],1:2] <- c(mean(fitpars[,indices[i,2]]), sd(fitpars[,indices[i,2]]))
    }
    }
    if(is.null(sof.names)) {
        if (fit.type=="cauchy") prior <- optim(par = c(0, 2.5), fn = fit.cauchy, dat = fitpars)$par
        if (fit.type=="normal") prior <- c(mean(fitpars), sd(fitpars))
    } 
    return(prior)
  }
   model.tests <- get_models() 
 # load wdfw fits
  mcmc.filename <- paste("outputs/DFW/modelselection/fit_",species.2.use,"_", pmodel.index,"mcmc.Rdata", sep="")
  load(file = mcmc.filename)
  pfit <- model.stan
  
  mcmc.filename <- paste("outputs/DFW/modelselection/fit_",species.2.use,"_", ymodel.index,"mcmc.Rdata", sep="")
  load(file = mcmc.filename)
  yfit <- model.stan
  ###Get names of parameters / categories in sof data
  sof.basin.names <- levels(thedata$Basin)
  sof.region.names <- levels(thedata$Region)
  p.model.text <-paste(model.tests[[pmodel.index]], "+SamplingType", sep = "")
  y.model.text <-paste(model.tests[[ymodel.index]], "+SamplingType", sep = "")
  
  x.p <- model.matrix(object = eval(parse(text = p.model.text)), data = thedata)
  x.y <- model.matrix(object = eval(parse(text = y.model.text)), data = thedata)
  
  # indices of betas that have priors from dfw data
  beta_p_names <- colnames(x.p)
  beta_y_names <- colnames(x.y)
  
 
  priors <- list()
    
    #beta parameters
  priors$beta_p <-
    fit.priors(
      fit = pfit,
      sof.names = beta_p_names,
      dfw.names = all.dfw.names$beta_p,
      parname = "beta_p",
      fit.type = "cauchy",
      default.prior = c(0, 2.5)
    )
  priors$beta_y <-
    fit.priors(
      fit = yfit,
      sof.names = beta_y_names,
      dfw.names = all.dfw.names$beta_y,
      parname = "beta_y",
      fit.type = "cauchy",
      default.prior = c(0, 2.5)
    )
  # extract variances
  priors$shape <- 
    fit.priors(
      fit = yfit,
      parname = "shape",
      fit.type = "normal",
      default.prior = c(0, 2.5)
    )
  
  priors$sigma_basin_p <-
    fit.priors(
      fit = pfit,
      parname = "sigma_basin_p",
      fit.type = "cauchy",
      default.prior = c(0, 2.5)
    )
  priors$sigma_basin_y <-
    fit.priors(
      fit = yfit,
      parname = "sigma_basin_y",
      fit.type = "cauchy",
      default.prior = c(0, 2.5)
    )
  
  priors$sigma_psi_p <-
    fit.priors(
      fit = pfit,
      parname = "sigma_psi_p",
      fit.type = "cauchy",
      default.prior = c(0, 2.5)
    )
  priors$sigma_psi_y <-
    fit.priors(
      fit = yfit,
      parname = "sigma_psi_y",
      fit.type = "cauchy"
    )
  priors$logitrho_p <-
    fit.priors(
      fit = pfit,
      parname = "logitrho_p",
      fit.type = "normal"
    )
  priors$logitrho_y <-
    fit.priors(
      fit = yfit,
      parname = "logitrho_y",
      fit.type = "normal"
    )
  
  # basins
  
  priors$gamma_p_raw <-fit.priors(
    fit = pfit,
    sof.names = sof.basin.names,
    dfw.names = all.dfw.names$basin,
    parname = "gamma_p_raw",
    fit.type = "normal",
    default.prior = c(0, 1)
  )
  
  priors$gamma_y_raw <-fit.priors(
    fit = yfit,
    sof.names = sof.basin.names,
    dfw.names = all.dfw.names$basin,
    parname = "gamma_y_raw",
    fit.type = "normal",
    default.prior = c(0, 1)
  ) 
  
  priors$sigma_region_p <-fit.priors(
    fit = pfit,
    sof.names = sof.basin.names,
    dfw.names = all.dfw.names$basin,
    parname = "sigma_region_p",
    fit.type = "normal",
    default.prior = c(0, 1)
  )
  
  priors$sigma_region_y <-fit.priors(
    fit = yfit,
    sof.names = sof.basin.names,
    dfw.names = all.dfw.names$basin,
    parname = "sigma_region_y",
    fit.type = "normal",
    default.prior = c(0, 1)
  ) 
  # region
  priors$theta_p_raw <-fit.priors(
    fit = pfit,
    sof.names = sof.region.names,
    dfw.names = all.dfw.names$region,
    parname = "theta_p_raw",
    fit.type = "normal",
    default.prior = c(0, 1)
  )
  priors$theta_y_raw <-fit.priors(
    fit = yfit,
    sof.names = sof.region.names,
    dfw.names = all.dfw.names$region,
    parname = "theta_y_raw",
    fit.type = "normal",
    default.prior = c(0, 1)
  )
  

  
  return(priors)
}

initf3 <- function(chain_id = 1, nyear, nbetap, nbetay, ngamma, ntheta) {
  list(
    beta_p = as.array(runif(nbetap, -1.0, 1.0)),
    beta_y = as.array(runif(nbetay, -1.0, 5.0)),
    shape = runif(n = 1, 1, 10),
    gamma_p_raw = rnorm(ngamma, 0, 1),
    gamma_y_raw = rnorm(ngamma, 0, 1),
    theta_p_raw = rnorm(ntheta, 0,1),
    theta_y_raw = rnorm(ntheta,0,1),
    logitrho_p = runif(1, -1, 1),
    logitrho_y = runif(1, -1, 1),
    sigma_psi_p = runif(1, 0, 1),
    sigma_psi_y = runif(1, 0, 1),
    eps_p_raw = rnorm(nyear,0, 1 ),
    eps_y_raw = rnorm(nyear,0, 1 ),
    sigma_basin_p = runif(n = 1, 0, 2),
    sigma_basin_y = runif(n = 1, 0.0, 2),
    sigma_region_p = runif(n = ngamma, 0, 1.5),
    sigma_region_y = runif(n = ngamma, 0, 1.5)
  )
 
}

# Plot DFA Output ---------------------------------------------------------


plot_loadings2 <- function (rotated_modelfit, names = NULL, facet = TRUE, violin = TRUE, 
                            conf_level = 0.95, threshold = NULL) 
{
  v <- reshape2::melt(rotated_modelfit$Z_rot, varnames = c("iter", 
                                                           "name", "trend"))
  v$trend <- paste0("Trend ", v$trend)
  v$trend <- as.factor(v$trend)
  if (!is.null(names)) 
    v$name <- names[v$name]
  v$name <- as.factor(v$name)
  v <- dplyr::group_by(v, .data$name, .data$trend)
  v <- dplyr::mutate(v, q_lower = sum(.data$value < 0)/length(.data$value), 
                     q_upper = 1 - .data$q_lower, prob_diff0 = max(.data$q_lower, 
                                                                   .data$q_upper))
  v <- dplyr::ungroup(v)
  vsum <- dplyr::group_by(v, .data$name, .data$trend)
  vsum <- dplyr::summarize(vsum, lower = quantile(.data$value, 
                                                  probs = (1 - conf_level)/2), upper = quantile(.data$value, 
                                                                                                probs = 1 - (1 - conf_level)/2), median = median(.data$value), 
                           q_lower = sum(.data$value < 0)/length(.data$value), q_upper = 1 - 
                             .data$q_lower, prob_diff0 = max(.data$q_lower, .data$q_upper))
  df <- dplyr::ungroup(vsum)
  if (!is.null(threshold)) {
    df <- df[df$prob_diff0 >= threshold, ]
    v <- v[v$prob_diff0 >= threshold, ]
  }
  if (!violin) {
    p1 <- ggplot(df, aes_string(x = "name", y = "median", 
                                col = 2, alpha = "prob_diff0")) + geom_point(size = 3, 
                                                                             position = position_dodge(0.3)) + geom_errorbar(aes_string(ymin = "lower", 
                                                                                                                                        ymax = "upper"), position = position_dodge(0.3), 
                                                                                                                             width = 0) + geom_hline(yintercept = 0, lty = 2) + 
      coord_flip() + xlab("Time Series") + ylab("Loading") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                panel.background = element_blank())
  }
  if (violin) {
    p1 <- ggplot(v, aes_string(x = "name", y = "value", fill = 2, 
                               alpha = "prob_diff0")) + geom_violin(color = "black") + 
      geom_hline(yintercept = 0, lty = 2) + coord_flip() + 
      xlab("Time Series") + ylab("Loading")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                panel.background = element_blank())
  }
  if (facet) {
    p1 <- p1 + facet_wrap(~trend, scales = "free_x")
  }
  p1
}

plot_community <- function(filename, p = F, yearlist, sort.order = NULL) {
require(reshape2)
  load(file = filename)
nyear <- nrow(p.list[[1]])

# Create a matrix that has Year, Median Abundance, and Species. Might be easiest to make a large df and then use reshape to turn it into short form
output <- data.frame(matrix(NA, nrow = nyear, ncol = 16))
species.list <- load_species_names()
names(output) <- c("Year", species.list)
output$Year <- yearlist

for (i in 1:15) {
if(p)  output[,i+1] <- p.list[[i]][,2]
if(!p) output[,i+1] <- y.list[[i]][,2]
}


# now get mean by species
species.means <- colMeans(output[,2:16])
if (is.null(sort.order)) sort.order <- names(sort(species.means, decreasing = TRUE))
long.output <- melt(output, id.vars = "Year", variable.name = "Species", value.name= "Catch")
long.output$Species <- factor(long.output$Species, levels= sort.order)

p1 <- ggplot(long.output, aes(x=Year, y=Catch, fill=Species)) + 
  geom_area(alpha = 1, size = 2) +
  scale_fill_viridis(discrete = T, option = "plasma", direction = -1) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=14),
        axis.title= element_text(size = 16)) +
  ylab("Catch Rate") +
  scale_y_continuous(expand = c(0,0))
return(list(p1, sort.order, output)) 
}