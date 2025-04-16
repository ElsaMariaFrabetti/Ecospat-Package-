# Ecospat Package for species distribution modelling

1 LOAD DATA

```
install.packages("ecospat", dependencies = TRUE)
library(ecospat)
citation("ecospat")
```
Test data for ecospat library
```
data("ecospat.testData")
names(ecospat.testData)
```
Test data for Niche Overlap Analysis
```
data("ecospat.testNiche.inv")
names(ecospat.testNiche.inv)
data("ecospat.testNiche.nat")
names(ecospat.testNiche.nat)
```
Test tree for Phylogenetic Diversity Analysis
```
install.packages("ape")
if(requireNamespace("ape"))
fpath <- system.file("extdata", "ecospat.testTree.tre", package = "ecospat")
tree <- ape::read.tree(fpath)
tree$tip.label
plot(tree, cex = 0.6)
```
![Rplot1](https://github.com/user-attachments/assets/a004c030-7b66-4be0-b16c-85f3077d857f)

2 PRE-MODELLING ANALYSIS

2.1 Spatial Auto-correlation 

Mantel Correlogram 
```
ecospat.mantel.correlogram(dfvar = ecospat.testData[c(2:16)], colxy = 1:2, n = 100, colvar = 3:7, max = 1000, nclass = 10, nperm = 100)
```
![Rplot2 1](https://github.com/user-attachments/assets/c830eae8-ccec-4c0b-a5f6-1f56e7d8e18b)

2.2 Predictor Variable Selection

Number of predictors with Pearson Correlation 
```
colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method = "pearson")
ecospat.npred(x, th=0.75)
```
Number of predictors with Spearman Correlation
```
x <- cor(colvar, method = "spearman")
ecospat.npred(x, th=0.75)
```
2.3 Climate Analogy Tools
```
x <- ecospat.testData[c(4:8)]
p <- x[1:90,] #A projection dataset
ref <- x[91:300,] #A reference dataset
ecospat.climan(ref,p)
```
Extrapolation detection, creating a MESS object 
```
x <- ecospat.testData[c(2,3,4:8)]
proj <- x[1:90,] #A projection dataset
cal <- x[91:300,] #A calibration dataset
mess.object <- ecospat.mess(proj, cal, w = "default")
ecospat.plot.mess(mess.object, cex = 1, pch = 15)
```
![Rplot2 3](https://github.com/user-attachments/assets/d81c612b-6531-4977-9ab6-887109158a99)

2.4 Phylogenetic Diversity Measures
```
if(requireNamespace("ape"))
fpath <- system.file("extdata", "ecospat.testTree.tre", package = "ecospat")
tree <- ape::read.tree(fpath)
data <- ecospat.testData[9:52]
pd <- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = TRUE)
plot(pd)
```
![Rplot2 4](https://github.com/user-attachments/assets/9228acb6-e959-49de-9636-4d33abf5d454)

2.5 Niche Quantification and Comparison with Ordination techniques
```
install.packages("sf") #necessary to carry the grid
library(ade4)
inv <- ecospat.testNiche.inv
nat <- ecospat.testNiche.nat
pca.env <- ade4::dudi.pca(rbind(nat,inv)[,3:10], scannf = F, nf=2)
ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig) #The correlation circle indicates the contribution of original predictors to the PCA axes
```
![Rplot2 5](https://github.com/user-attachments/assets/213c5938-4b05-4689-b9c5-4c3f8d70c5ad)

Now we can predict the scores on the axes
```
scores.globclim <- pca.env$li #PCA scores for the whole study area
scores.sp.nat <- ade4::suprow(pca.env, nat[which(nat[,11]==1),3:10])$li #PCA scores for the sp native distr
scores.sp.inv <- ade4::suprow(pca.env, inv[which(inv[,11]==1),3:10])$li #PCA scores for the sp invasive distr
scores.clim.nat <- ade4::suprow(pca.env, nat[,3:10])$li #PCA scores for the whole native study area
scores.clim.inv <- ade4::suprow(pca.env, inv[,3:10])$li #PCA scores for the whole invaded study area
```
Calculate the occurrence densities grid
```
grid.clim.nat <- ecospat.grid.clim.dyn(glob = scores.globclim, glob1 = scores.clim.nat, sp = scores.sp.nat, R = 100, th.sp = 0) #gridding the native niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob = scores.globclim, glob1 = scores.clim.inv, sp = scores.sp.inv, R = 100, th.sp = 0) #gridding the invasive niche
```
Calculate niche overlap 
```
D.overlap <- ecospat.niche.overlap(grid.clim.nat, grid.clim.inv, cor = TRUE)$D
D.overlap #Compute Schoener's D, index of niche overlap = 22%
```
Perform the niche equivalency test according to Warren et al.
```
eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv, rep = 10, intersection = 0.1, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower")
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
```
![Rplot2 5'](https://github.com/user-attachments/assets/97e0e206-05e0-4cb9-812a-fbb88d7eeb4d)

Perform the niche similarity test
```
sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv, rep = 10, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower", intersection = 0.1, rand.type = 1)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
```
![Rplot2 5''](https://github.com/user-attachments/assets/04fcc4ec-a81b-4330-8299-cf55ad6c1a89)

Delimiting niche categories and quantifying niche dynamics in analogue climates
```
niche.dyn <- ecospat.niche.dyn.index(grid.clim.nat, grid.clim.inv)
```
Visualizing niche categories, niche dynamics and climate analogy between ranges
```
ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.25, interest=2, title = "Niche Overlap", name.axis1 = "PC1", name.axis2 = "PC2")
ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv) #Plot niche overlap
```
![Rplot2 5'''](https://github.com/user-attachments/assets/b8e3a11c-d0a8-4632-b3d5-88180ef62538)

```
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity") #Plot similarity test for niche expansion 
```
![Rplot2 5''''](https://github.com/user-attachments/assets/2a8ca50d-0263-441e-b351-08307f2c1407)

```
ecospat.plot.overlap.test(sim.test, "stability", "Similarity") 
```
![Rplot2 5'''''](https://github.com/user-attachments/assets/f765fbee-a492-4539-97f4-9468d3a7fa49)

```
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
```
![Rplot2 5''''''](https://github.com/user-attachments/assets/50fe3b8c-2377-4154-a1b8-f4ba0fe4356a)

```
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob = as.data.frame(rbind(nat,inv)[,10]), glob1 = as.data.frame(nat[,10]), sp=as.data.frame(nat[which(nat[,11]==1),10]), R=1000, th.sp = 0) #Gridding the native niche
grid.clim.t.inv <-ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat,inv)[,10]),glob1=as.data.frame(inv[,10]),sp=as.data.frame(inv[which(inv[,11]==1),10]), R=1000, th.sp=0) #Gridding the invaded niche
t.dyn <- ecospat.niche.dyn.index(grid.clim.t.nat, grid.clim.t.inv)                                       
ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0, interest=2, title = "Niche Overlap", name.axis1 = "Average temperature") #Plot the niche dynamics along one gradient (T)
```
![Rplot2 5'''''''''''''''](https://github.com/user-attachments/assets/10dea4b4-75bb-48e0-923b-4ec13291bacb)

2.6 Biotic Interactions

Species co-occurrence analysis with a presence-absence matrix
```
data <- ecospat.testData[c(9:16, 54:57)]
ecospat.co_occurrences(data)
```
![Rplot2 6](https://github.com/user-attachments/assets/e64536f4-ce21-42b0-858b-131b4603d5b0)

Pairwise co-occurrence analysis with calculation of the C-score index
```
data <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
nperm <- 100
outpath <- getwd()
ecospat.Cscore(data, nperm, outpath)
```
![Rplot2 6'](https://github.com/user-attachments/assets/dbf2a77b-30a2-44f0-97dd-13fe915b9c98)

2.7 Data Preparation 

Correlation plot of variables 
```
data <- ecospat.testData[,4:8]
ecospat.cor.plot(data)
```
![Rplot2 7](https://github.com/user-attachments/assets/d4dd598b-a5c7-46d0-a011-5b4d93d10bd1)

Calibration and evaluation dataset
```
data <- ecospat.testData
caleval <- ecospat.caleval(data = ecospat.testData[53], xy=data[2:3], row.num = 1:nrow(data), nrep = 2, ratio = 0.7, disaggregate = 0.2, pseudoabs = 100, npres = 10, replace = FALSE)
head(caleval) #Obtaining an evaluation and calibration dataset with a desired ratio of disaggregation
```
3 CORE NICHE MODELLING









