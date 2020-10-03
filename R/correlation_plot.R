library(ggplot2)
library(ggcorrplot)
require(viridis)
source("R/project_functions.R")
dataset = "SOF" # DFW or SOF
colors <- plasma(4)
colors <- colors[1:3]
colors[2] <- "white"
colors <- rev(colors)
# SOF Data
sof_data <- plot_community(filename = "outputs/SOF/year_effects.Rdata", yearlist = 1948:1977)[[3]][,-1]
sof_data <- scale(sof_data)

for(i in 1:ncol(sof_data)) {
  sof_data[-1,i] = diff(sof_data[,i])
}
sof_data = sof_data[-1,]
sof_cor = cor(sof_data)
n <- nrow(sof_data)
# get pvalues
sof.t.cor <- sof_cor / sqrt((1-sof_cor^2)/(n-2))
sof.p.val <- 1-pt(abs(sof.t.cor), lower.tail = TRUE, df = n-2) 
sof.results <- matrix(0, nrow = nrow(sof.t.cor), ncol = nrow(sof.t.cor))
colnames(sof.results) <- colnames(sof_cor); rownames(sof.results) <- rownames(sof.results)
sof.results[which(sof.p.val<=0.05)] <- sof_cor[which(sof.p.val<=0.05)]
diag(sof.results) <- 0; sof.results[lower.tri(sof.results)]<-0
c(length(which(sof.results<0)), length(which(sof.results>0))) # counts the number of negative and positive significant correlations
sof.neg <- which(sof.results <0)
sof.pos <- which(sof.results >0)
# DFW Data
dfw_data <- plot_community(filename = "outputs/DFW/year_effects.Rdata", yearlist = 1990:2016)[[3]][,-1]
dfw_data <- scale(dfw_data)
for(i in 1:ncol(dfw_data)) {
  dfw_data[-1,i] = diff(dfw_data[,i])
}
dfw_data = dfw_data[-1,]
dfw_cor = cor(dfw_data)

n <- nrow(dfw_data)
# get pvalues
dfw.t.cor <- dfw_cor / sqrt((1-dfw_cor^2)/(n-2))
dfw.p.val <- 1-pt(abs(dfw.t.cor), lower.tail = TRUE, df = n-2) 

dfw.results <- matrix(0, nrow = nrow(dfw.t.cor), ncol = nrow(dfw.t.cor))
colnames(dfw.results) <- colnames(dfw_cor); rownames(dfw.results) <- rownames(dfw.results)
dfw.results[which(dfw.p.val<=0.05)] <- dfw_cor[which(dfw.p.val<=0.05)]
diag(dfw.results) <- 0; dfw.results[lower.tri(dfw.results)] <- 0
c(length(which((dfw.results)<0)), length(which((dfw.results)>0))) # counts the number of negative and positive significant correlations
dfw.neg <- which(dfw.results <0)
dfw.pos <- which(dfw.results>0)

# how many negative correlations remained
length(intersect(sof.neg, dfw.neg))
#how many positive correlations remained
length(intersect(sof.pos, dfw.pos))

# how many switched!
length(intersect(sof.neg, dfw.pos)) + length(intersect(sof.pos, dfw.neg))

### Plot Stuff
g1 = ggcorrplot(sof_cor,outline.col = "white",
  type = "lower", colors = colors) + ggtitle("Historical") + 
  theme(axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10))
g2 = ggcorrplot(dfw_cor,outline.col = "white",
  type = "lower",colors = colors) + ggtitle("Contemporary") + 
  theme(axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10))
pdf("graphics/Correlation_plot_1.pdf",height=5,width=10)
gridExtra::grid.arrange(g1,g2,nrow=1)
dev.off()


plotfilename <- "graphics/correllation_plot_2.pdf"

g2 = ggcorrplot(dfw_cor,outline.col = "white",
  type = "lower",hc.order = TRUE, colors = colors) + ggtitle("Contemporary") + 
  theme(axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10))
species.order <- levels(g2$data$Var2)
species.indices <- match(species.order, colnames(sof_cor))
# rearrange SOF rows
sof_cor_sorted <- sof_cor[species.indices,]
sof_cor_sorted[,1:15] <- sof_cor_sorted[, species.indices]
colnames(sof_cor_sorted) <- species.order; rownames(sof_cor_sorted) <- species.order
# rearrange SOF data
g1 = ggcorrplot(sof_cor_sorted,outline.col = "white",
                type = "lower", colors = colors) + ggtitle("Historical") + 
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

# get order of species in contemporary plot
pdf(plotfilename,height=5,width=10)
gridExtra::grid.arrange(g1,g2,nrow=1)
dev.off()
system2("open", args = c("-a Skim.app", plotfilename))