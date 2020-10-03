# Script to plot fit to SOF data
require(rstan)

source("R/project_functions.R")
spc <- 9
load("data/SOF/SOFdata.Rdata")
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
species <- species.list[spc]
extract.expression <- paste('thedata <- species.data$`', species, '`', sep = "")
eval(parse(text = extract.expression))
### Modify Data ####

thedata <- format.sof.data(thedata, depth.std)
thedata$SamplingType <- relevel(thedata$SamplingType, "Otter trawl") # make Otter trawl the reference level (b/c it is most common)

# for presence - absences
formula.txt <- paste(model.tests[[pmodel[spc]]],"+SamplingType", sep = "")
tmp.x <- model.matrix(object = eval(parse(text =formula.txt)), data = thedata)
nbetap <- ncol(tmp.x)
tmp.t <- model.matrix(object = ~ -1 + Year, data = thedata)
tmp.b<- model.matrix(~-1 + Basin, data = thedata) 
formula.txt <- paste0(model.tests[[ymodel[spc]]], "+SamplingType")
tmp.xy <- model.matrix(object = eval(parse(text =formula.txt)), data = thedata)
nbetay <- ncol(tmp.xy)

ngamma <- ncol(tmp.b)
nyear <- ncol(tmp.t)

filename <- paste("outputs/SOF/MCMC/",species, "_fitmcmc.Rdata", sep = "")
load(file= filename)
rslt <- summary(model.stan)$summary[, "50%"]
parnames <- names(rslt)
psi_p <- rstan::extract(model.stan, pars = "psi_p")[[1]]
nyear <- ncol(psi_p)
est.psip <- rslt[grep("psi_p", parnames)][1:nyear]
est.psiy <- rslt[grep("psi_y", parnames)][1:nyear]

beta_p <- rstan::extract(model.stan, pars = "beta_p")[[1]]
nbetap <- ncol(beta_p)
beta_y <- rstan::extract(model.stan, pars = "beta_y")[[1]]
nbetay <- ncol(beta_y)
gammay <- rstan::extract(model.stan,pars = "gamma_y")[[1]]
est.betap <- rslt[grep("beta_p", parnames)][1:nbetap]
est.betay <- rslt[grep("beta_y", parnames)][1:nbetay]
est.gammap <- rslt[grep("gamma_pa", parnames)]
est.gammay <- apply(X = gammay, MARGIN = 2, median)

p.predict <- inv_logit(tmp.x %*% est.betap + tmp.b %*% est.gammap + tmp.t %*% est.psip)
y.predict <- exp(tmp.xy %*% est.betay + tmp.b %*% est.gammay )
thedata$ppredict <- p.predict
thedata$ypredict <- y.predict

basins <- levels(thedata$Basin)

for ( i in 1:length(basins)){
  basin.data <- thedata %>%
    filter(Basin == basins[i])
  par(mfrow = c(2,2))
  with(basin.data, plot(sDEPTH, present,
                        type = "p",
                        pch = 21,
                        bg = "black",
                        ylab = "present",
                        xlab = "smoothed depth"
  )
  )
  points(basin.data$sDEPTH, basin.data$ppredict,
         pch = 21,
         bg = "red"
  )
  
  with(basin.data, plot(sDEPTH2, present,
                        type = "p",
                        pch = 21,
                        bg = "black",
                        ylab = "present",
                        xlab = "smoothed depth^2"
  )
  )
  with(basin.data, points(sDEPTH2, ppredict,
                          type = "p",
                          pch = 21,
                          bg = "red"
  )
  )
  
  basin.data <- basin.data %>%
    filter(present==1)
  if(nrow(basin.data)>0) {
    
    with(basin.data, plot(sDEPTH, count,
                          type = "p",
                          pch = 21,
                          bg = "black",
                          ylab = "cpue",
                          xlab = "smoothed depth"
    )
    )
    with(basin.data, points(sDEPTH, ypredict,
                            pch = 21,
                            bg = "red"
    )
    )
    
    with(basin.data, plot(sDEPTH2, count,
                          type = "p",
                          pch = 21,
                          bg = "black",
                          ylab = "cpue",
                          xlab = "smoothed depth^2"
    )
    )
    with(basin.data, points(sDEPTH2, ypredict,
                            type = "p",
                            pch = 21,
                            bg = "red"
    )
    )
  } else {
    plot(c(), c(), type = "n", xlim = c(0,1), ylim = c(0,1))
    plot(c(), c(), type = "n", xlim = c(0,1), ylim = c(0,1))
  }
  mtext(text = basins[i], side = 3, outer= T, line = -1)
}



