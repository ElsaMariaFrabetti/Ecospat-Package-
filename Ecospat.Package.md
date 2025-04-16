# Ecospat Package for species distribution modelling

1 LOAD DATA

```install.packages("ecospat", dependencies = TRUE)
library(ecospat)
citation("ecospat")
```
Test data for ecospat library
```data("ecospat.testData")
names(ecospat.testData)
```
Test data for Niche Overlap Analysis
```data("ecospat.testNiche.inv")
names(ecospat.testNiche.inv)
data("ecospat.testNiche.nat")
names(ecospat.testNiche.nat)
```
Test tree for Phylogenetic Diversity Analysis
```install.packages("ape")
if(requireNamespace("ape"))
fpath <- system.file("extdata", "ecospat.testTree.tre", package = "ecospat")
tree <- ape::read.tree(fpath)
tree$tip.label
plot(tree, cex = 0.6)
```

2 PRE-MODELLING ANALYSIS

2.1 Spatial Auto-correlation 

Mantel Correlogram 
```ecospat.mantel.correlogram(dfvar = ecospat.testData[c(2:16)], colxy = 1:2, n = 100, colvar = 3:7, max = 1000, nclass = 10, nperm = 100)
```
2.2 Predictor Variable Selection

Number of predictors with Pearson Correlation 
```colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method = "pearson")
ecospat.npred(x, th=0.75)
```
Number of predictors with Spearman Correlation
```x <- cor(colvar, method = "spearman")
ecospat.npred(x, th=0.75)
```
2.3 Climate Analogy Tools
```x <- ecospat.testData[c(4:8)]
p <- x[1:90,] #A projection dataset
ref <- x[91:300,] #A reference dataset
ecospat.climan(ref,p)
```
Extrapolation detection, creating a MESS object 
```x <- ecospat.testData[c(2,3,4:8)]
proj <- x[1:90,] #A projection dataset
cal <- x[91:300,] #A calibration dataset
mess.object <- ecospat.mess(proj, cal, w = "default")
ecospat.plot.mess(mess.object, cex = 1, pch = 15)
```
2.4 Phylogenetic Diversity Measures
```if(requireNamespace("ape"))
fpath <- system.file("extdata", "ecospat.testTree.tre", package = "ecospat")
tree <- ape::read.tree(fpath)
data <- ecospat.testData[9:52]
pd <- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = TRUE)
plot(pd)
```
2.5 Niche Quantification and Comparison with Ordination techniques
```install.packages("sf") #necessary to carry the grid
library(ade4)
inv <- ecospat.testNiche.inv
nat <- ecospat.testNiche.nat
pca.env <- ade4::dudi.pca(rbind(nat,inv)[,3:10], scannf = F, nf=2)
ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig) #The correlation circle indicates the contribution of original predictors to the PCA axes
```
Now we can predict the scores on the axes
```scores.globclim <- pca.env$li #PCA scores for the whole study area
scores.sp.nat <- ade4::suprow(pca.env, nat[which(nat[,11]==1),3:10])$li #PCA scores for the sp native distr
scores.sp.inv <- ade4::suprow(pca.env, inv[which(inv[,11]==1),3:10])$li #PCA scores for the sp invasive distr
scores.clim.nat <- ade4::suprow(pca.env, nat[,3:10])$li #PCA scores for the whole native study area
scores.clim.inv <- ade4::suprow(pca.env, inv[,3:10])$li #PCA scores for the whole invaded study area
```
Calculate the occurrence densities grid
```grid.clim.nat <- ecospat.grid.clim.dyn(glob = scores.globclim, glob1 = scores.clim.nat, sp = scores.sp.nat, R = 100, th.sp = 0) #gridding the native niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob = scores.globclim, glob1 = scores.clim.inv, sp = scores.sp.inv, R = 100, th.sp = 0) #gridding the invasive niche
```
Calculate niche overlap 
```D.overlap <- ecospat.niche.overlap(grid.clim.nat, grid.clim.inv, cor = TRUE)$D
D.overlap #Compute Schoener's D, index of niche overlap = 22%
```
Perform the niche equivalency test according to Warren et al.
```eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv, rep = 10, intersection = 0.1, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower")
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
```
Perform the niche similarity test
```sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv, rep = 10, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower", intersection = 0.1, rand.type = 1)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
```
Delimiting niche categories and quantifying niche dynamics in analogue climates
```niche.dyn <- ecospat.niche.dyn.index(grid.clim.nat, grid.clim.inv)
```
Visualizing niche categories, niche dynamics and climate analogy between ranges
```ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.25, interest=2, title = "Niche Overlap", name.axis1 = "PC1", name.axis2 = "PC2")
ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv) #Plot niche overlap
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity") #Plot similarity test for niche expansion 
ecospat.plot.overlap.test(sim.test, "stability", "Similarity") 
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob = as.data.frame(rbind(nat,inv)[,10]), glob1 = as.data.frame(nat[,10]), sp=as.data.frame(nat[which(nat[,11]==1),10]), R=1000, th.sp = 0) #Gridding the native niche
grid.clim.t.inv <-ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat,inv)[,10]),glob1=as.data.frame(inv[,10]),sp=as.data.frame(inv[which(inv[,11]==1),10]), R=1000, th.sp=0) #Gridding the invaded niche
t.dyn <- ecospat.niche.dyn.index(grid.clim.t.nat, grid.clim.t.inv)                                       
ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0, interest=2, title = "Niche Overlap", name.axis1 = "Average temperature") #Plot the niche dynamics along one gradient (T)
```
2.6 Biotic Interactions

Species co-occurrence analysis with a presence-absence matrix
```data <- ecospat.testData[c(9:16, 54:57)]
ecospat.co_occurrences(data)
```
Pairwise co-occurrence analysis with calculation of the C-score index
```data <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
nperm <- 100
outpath <- getwd()
ecospat.Cscore(data, nperm, outpath)
```
2.7 Data Preparation 

Correlation plot of variables 
```data <- ecospat.testData[,4:8]
ecospat.cor.plot(data)
```
Calibration and evaluation dataset
```data <- ecospat.testData
caleval <- ecospat.caleval(data = ecospat.testData[53], xy=data[2:3], row.num = 1:nrow(data), nrep = 2, ratio = 0.7, disaggregate = 0.2, pseudoabs = 100, npres = 10, replace = FALSE)
head(caleval) #Obtaining an evaluation and calibration dataset with a desired ratio of disaggregation
```









