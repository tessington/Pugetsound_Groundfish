rm(list = ls())
require(extraDistr)
require(rstan)
require(loo)
require(dplyr)
source("R/project_functions.R")
species.names <- load_species_names()

# Get Sampler Diagnostics -------------------------------------------------
select.output <- matrix(NA, nrow = 15, ncol = 10)
rownames(select.output) <- species.names
colnames(select.output) <- c("p1","p2","p3","p4","p5","y1","y2","y3","y4","y5")

for (spc in 1:15) {
  
species <- species.names[spc]
# diagnostics for each model
count_diverge<- rep(NA, 4)
for (j in 1:5) {
  filename <- paste("outputs/DFW/modelselection/fit_",species, "_", j, "mcmc.Rdata", sep = "")
  load(file= filename)
  sampler_params <- get_sampler_params(model.stan, inc_warmup = FALSE)
  count_diverge[j]<-sum( sapply(sampler_params, function(x) sum(x[, "divergent__"])))
}

# Perform Model Selection -------------------------------------------------

models.2.use <- 1:5
model.compare<- model.select(species, models.2.use)
print(count_diverge)
print(model.compare)
select.output[spc,1:5] <- model.compare[[1]]
select.output[spc,6:10] <- model.compare[[2]]
}

save(file = "outputs/DFW/modelselection/modelselect.Rdata",select.output)
write.csv( x= select.output, file = "outputs/DFW/modelselection/modelselect.csv")

# View Model Fits ---------------------------------------------------------
species <- species.names[spc]
### Extract posterior medians ###
model <- 1
model.tests <- get_models()
formula.txt <- model.tests[[model]]

###Load Data###
data <- format_data(species)
thedata <-data$tmpdata
tmp.x <- model.matrix(object = eval(parse(text =formula.txt)), data = thedata)
nbeta <- ncol(tmp.x)
tmp.t <- model.matrix(object = ~ -1 + YEAR, data = thedata)
tmp.b<- model.matrix(~-1 + Basin, data = thedata) 
ndata <- nrow(tmp.x)
ngamma <- ncol(tmp.b)
nyear <- ncol(tmp.t)


filename <- paste("outputs/DFW/modelselection/fit_",species, "_", model, "mcmc.Rdata", sep = "")
load(file= filename)
rslt <- summary(model.stan)$summary[, "50%"]
parnames <- names(rslt)

est.psip <- rslt[grep("psi_p", parnames)][1:nyear]
est.psiy <- rslt[grep("psi_y", parnames)][1:nyear]

est.betap <- rslt[grep("beta_p", parnames)][1:nbeta]
est.betay <- rslt[grep("beta_y", parnames)][1:nbeta]
est.gammap <- rslt[grep("gamma_pa", parnames)]
est.gammay <- rslt[grep("gamma_y", parnames)]

p.predict <- inv_logit(tmp.x %*% est.betap + tmp.b %*% est.gammap + tmp.t %*% est.psip)
y.predict <- exp(tmp.x %*% est.betay + tmp.b %*% est.gammay )
thedata$ppredict <- p.predict
thedata$ypredict <- y.predict

basins <- levels(thedata$Basin)

for ( i in 1:length(basins)){
  basin.data <- thedata %>%
    filter(Basin == basins[i])
  par(mfrow = c(2,2))
  with(basin.data, plot(sDEPTH, PRESENT,
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

  with(basin.data, plot(sDEPTH2, PRESENT,
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
    filter(PRESENT==1)
  if(nrow(basin.data)>0) {
  
  with(basin.data, plot(sDEPTH, C/EFFORT,
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

  with(basin.data, plot(sDEPTH2, C / EFFORT,
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
  
  
  
