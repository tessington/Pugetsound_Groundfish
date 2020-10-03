##rm(list = ls())
## Need to decide how to standardize.  Best to standardize to maximum

require(rstan)
require(RColorBrewer)
require(viridis)
source("R/project_functions.R")

species.list <- load_species_names()

# get best models
group.names <- c("flatfish", 
                "gadids",
                "chondrichthys",
                "other")

# species group indices
species.groups <- list(flatfish = c(1,6,15),
                       gadids = c(4,5,7,8),
                       chondrichthys = c(2,3,13,14),
                       other = c(9,10,11,12)
)

group.2.run <- 4

spc.2.run <-species.groups[[group.2.run]]
plotfilename <- paste("Graphics/sparkline_",group.names[group.2.run],".pdf", sep = "")
pdf(file = plotfilename, width = 8, height = 6)

xlim1 <- c(1945,1976)
xlim2 <- c(1990, 2016)
par(xpd =NA, mfrow = c(4+1,2), omi = c(0.5, 0.5, 0.5, 0.5), mar = c(2,3,1,1))

color.list <- plasma(n=16)[c(1,12)]
span=0.15
degree=2
n.points=41
family="gaussian"


# add transparency to a list of colors
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


add.2.plot.new <- function(p.output, y.output, first.year, last.year, col.2.use, span = 0.2, degree = 2, n.points = 41, family = "gaussian", pa = T){
  # use this fuction for pre-extracted year effects
 
  
  if (pa) abundance <- p.output
  if(!pa) abundance <- y.output
 
  
  if(pa) { 
    median.abundance <- abundance[,2]
    lower.abundance <- abundance[,1]
    upper.abundance <- abundance[,3]
  }
  
  if(!pa) {
    median.abundance <- abundance[,2] / max(abundance[,2]) # scale by max to preserve sense of varation
    lower.abundance <-abundance[,1] / max(abundance[,2])
    upper.abundance <- abundance[,3] / max(abundance[,2])
  }
  
  col.2.shade<-add.alpha(col.2.use,0.5)
  year.list <- first.year: last.year
  # make polygon
  bottom.poly.smooth<-loess.smooth(first.year:last.year,lower.abundance, span = span, degree = degree,
                                   family = family, evaluation = n.points*2)
  
  top.poly.smooth<-loess.smooth(first.year:last.year,upper.abundance, span = span, degree = degree,
                                family = family, evaluation = n.points*2)
  
  x.poly<-c(bottom.poly.smooth$x,rev(top.poly.smooth$x))
  y.poly<-c(bottom.poly.smooth$y,rev(top.poly.smooth$y))
  # x.poly <- c(year.list, rev(year.list))
  #  y.poly <- c(lower.abundance, rev(upper.abundance))
  polygon(x.poly,y.poly,col=col.2.shade,border=col.2.shade)
  
  lines(first.year : last.year, median.abundance, 
        lwd = 2, 
        col = col.2.use)
  
  ## remove WDFW stan data
} 
## Establish Plot

for (spc in spc.2.run) {
  species.2.use <- species.list[spc]
  plot(0,0,
       type = "n",
       xlim = xlim1,
       ylim = c(0,1),
       ylab = "",
       xlab = "", 
       axes = F)
  
# add species name on two lines
  pos <- regexpr(" ", species.2.use)
  
  if(pos[1]>0) {
    first.name <- substr(species.2.use, 1, pos[1]-1)
    last.name <- substr(species.2.use, pos[1]+1, nchar(species.2.use))
    expression.string <- bquote()
    text(x = c(1940, 1940), y = c(0.6, 0.3), labels =  c(first.name,last.name), pos = 3, cex = 1.5) 
  } else {
    text(x = c(1940), y = c(0.45), labels =  species.2.use, pos = 3, cex = 1.5)
  }
  
  #mtext(side = 2, text = species.2.use, las = 1, line = -3, adj = 1)
  axis(side = 2, at = c(0,1), las =1, pos = c(1948,0) )


## Get year effects = beta[1] + 
# first year = 1978, equals
## load in WDFW stan data for species
  ## load SOF stan data
  
  load(file = "outputs/SOF/year_effects.Rdata")
  col.2.use <- color.list[100]
  first.year <- 1948
  last.year <- 1977
  
  add.2.plot.new(p.list[[spc]], y.list[[spc]],first.year, last.year, color.list[1], span = span, degree = 2, pa = T)
  add.2.plot.new(p.list[[spc]], y.list[[spc]], first.year, last.year, color.list[2], span = span, degree = 2, pa = F)
  
  
  # make second plot for WDFW data
  plot(0,0,
       type = "n",
       xlim = xlim2,
       ylim = c(0,1),
       ylab = "",
       xlab = "", 
       axes = F)
  axis(side = 2, at = c(0,1), las =1, pos = c(1990,0) )
  
  
  # load fitted year effects for DFW data
  load(file = "outputs/DFW/year_effects.Rdata")

first.year <- 1990
last.year <- 2016
add.2.plot.new(p.list[[spc]], y.list[[spc]],first.year, last.year, color.list[1], span = span, degree = 2, pa = T)
add.2.plot.new(p.list[[spc]], y.list[[spc]], first.year, last.year, color.list[2], span = span, degree = 2, pa = F)



#abline(h = 0, col = "gray", xpd  = FALSE)
}

# add final plots just to have the x axis ticks
plot(0,0,
          type = "n",
          xlim = xlim1,
          ylab = "",
          xlab = "", 
          axes = F)
     
     axis(3, at = seq(1950, 1975, by = 5))
     
     plot(0,0,
          type = "n",
          xlim = xlim2,
          ylab = "",
          xlab = "", 
          axes = F)
     
     axis(3, at = seq(1990, 2015, by = 5))
dev.off()
system2("open", args = c("-a Skim.app", plotfilename))