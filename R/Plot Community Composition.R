##rm(list = ls())
## Need to decide how to standardize.  Best to standardize to maximum

require(rstan)
require(RColorBrewer)
require(viridis)
require(reshape2)
source("R/project_functions.R")

species.list <- load_species_names()

# get best models

# species group indices

  ## load SOF stan data
plotfilename <- "graphics/community.pdf"
#pdf(file = plotfilename, height = 5, width = 11, useDingbats= F)
par(mfrow = c(1,2))
  p1 <- plot_community(filename = "outputs/SOF/year_effects.Rdata", yearlist = 1948:1977)
  p2 <- plot_community(filename = "outputs/DFW/year_effects.Rdata", yearlist = 1990:2016, sort.order = p1[[2]])
  
  gridExtra::grid.arrange(p1[[1]]+ theme(legend.title = element_blank(), legend.text=element_text(size = 12)),
                          p2[[1]]+ theme(legend.position = "none"), nrow = 1, widths = c(3,2))
#dev.off()
#system2("open", args = c("-a Skim.app", plotfilename))
                          