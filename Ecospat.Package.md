# Ecospat Package for species distribution modelling

# 1 LOAD DATA

``` r
install.packages("ecospat", dependencies = TRUE)
library(ecospat)
citation("ecospat")
```
Test data for ecospat library
``` r
data("ecospat.testData")
names(ecospat.testData)
```
Test data for Niche Overlap Analysis
``` r 
data("ecospat.testNiche.inv")
names(ecospat.testNiche.inv)
data("ecospat.testNiche.nat")
names(ecospat.testNiche.nat)
```
Test tree for Phylogenetic Diversity Analysis
``` r
install.packages("ape")
if(requireNamespace("ape"))
fpath <- system.file("extdata", "ecospat.testTree.tre", package = "ecospat")
tree <- ape::read.tree(fpath)
tree$tip.label
plot(tree, cex = 0.6)
```
![Rplot1](https://github.com/user-attachments/assets/1c19b08c-d21a-457f-b00b-d13e767742dc)


## 2 PRE-MODELLING ANALYSIS

2.1 Spatial Auto-correlation 

Mantel Correlogram 
``` r
ecospat.mantel.correlogram(dfvar = ecospat.testData[c(2:16)], colxy = 1:2, n = 100, colvar = 3:7, max = 1000, nclass = 10, nperm = 100)
```
![Rplot2 1](https://github.com/user-attachments/assets/040ab1b5-8647-49ea-b6f4-a7f1b7663b10)


2.2 Predictor Variable Selection

Number of predictors with Pearson Correlation 
``` r
colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method = "pearson")
ecospat.npred(x, th=0.75)
```
Number of predictors with Spearman Correlation
``` r
x <- cor(colvar, method = "spearman")
ecospat.npred(x, th=0.75)
```
## 2.3 Climate Analogy Tools
``` r
x <- ecospat.testData[c(4:8)]
p <- x[1:90,] #A projection dataset
ref <- x[91:300,] #A reference dataset
ecospat.climan(ref,p)
```
Extrapolation detection, creating a MESS object 
``` r
x <- ecospat.testData[c(2,3,4:8)]
proj <- x[1:90,] #A projection dataset
cal <- x[91:300,] #A calibration dataset
mess.object <- ecospat.mess(proj, cal, w = "default")
ecospat.plot.mess(mess.object, cex = 1, pch = 15)
```
![Rplot2 3](https://github.com/user-attachments/assets/7ac803c9-f869-46c2-8016-2613321380c1)

## 2.4 Phylogenetic Diversity Measures
``` r
if(requireNamespace("ape"))
fpath <- system.file("extdata", "ecospat.testTree.tre", package = "ecospat")
tree <- ape::read.tree(fpath)
data <- ecospat.testData[9:52]
pd <- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = TRUE)
plot(pd)
```
![Rplot2 4](https://github.com/user-attachments/assets/7bf48819-a96c-4f78-92e1-be0cdec01f68)

## 2.5 Niche Quantification and Comparison with Ordination techniques
``` r
install.packages("sf") #necessary to carry the grid
library(ade4)
inv <- ecospat.testNiche.inv
nat <- ecospat.testNiche.nat
pca.env <- ade4::dudi.pca(rbind(nat,inv)[,3:10], scannf = F, nf=2)
ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig) #The correlation circle indicates the contribution of original predictors to the PCA axes
```
![Rplot2 5](https://github.com/user-attachments/assets/24af55d5-af40-40a8-96f4-73ce1b51be85)

Now we can predict the scores on the axes
``` r
scores.globclim <- pca.env$li #PCA scores for the whole study area
scores.sp.nat <- ade4::suprow(pca.env, nat[which(nat[,11]==1),3:10])$li #PCA scores for the sp native distr
scores.sp.inv <- ade4::suprow(pca.env, inv[which(inv[,11]==1),3:10])$li #PCA scores for the sp invasive distr
scores.clim.nat <- ade4::suprow(pca.env, nat[,3:10])$li #PCA scores for the whole native study area
scores.clim.inv <- ade4::suprow(pca.env, inv[,3:10])$li #PCA scores for the whole invaded study area
```
Calculate the occurrence densities grid
``` r
grid.clim.nat <- ecospat.grid.clim.dyn(glob = scores.globclim, glob1 = scores.clim.nat, sp = scores.sp.nat, R = 100, th.sp = 0) #gridding the native niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob = scores.globclim, glob1 = scores.clim.inv, sp = scores.sp.inv, R = 100, th.sp = 0) #gridding the invasive niche
```
Calculate niche overlap 
``` r
D.overlap <- ecospat.niche.overlap(grid.clim.nat, grid.clim.inv, cor = TRUE)$D
D.overlap #Compute Schoener's D, index of niche overlap = 22%
```
Perform the niche equivalency test according to Warren et al.
``` r
eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv, rep = 10, intersection = 0.1, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower")
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
```
![Rplot2 5 3](https://github.com/user-attachments/assets/d0f083a4-0970-4ab1-9778-535426c6ed14)

Perform the niche similarity test
``` r
sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv, rep = 10, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower", intersection = 0.1, rand.type = 1)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
```
![Rplot2 5 4](https://github.com/user-attachments/assets/b7761c98-76be-4f51-940a-6bc8358046f1)

Delimiting niche categories and quantifying niche dynamics in analogue climates
``` r
niche.dyn <- ecospat.niche.dyn.index(grid.clim.nat, grid.clim.inv)
```
Visualizing niche categories, niche dynamics and climate analogy between ranges
```
ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.25, interest=2, title = "Niche Overlap", name.axis1 = "PC1", name.axis2 = "PC2")
ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv) #Plot niche overlap
```
![Rplot2 5 5](https://github.com/user-attachments/assets/aa998473-3884-4b16-9f81-142da66d98cd)

``` r
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity") #Plot similarity test for niche expansion 
```
![Rplot2 5 6](https://github.com/user-attachments/assets/ce5f3e00-81c1-4924-aaed-bbf9ec2246d5)

``` r
ecospat.plot.overlap.test(sim.test, "stability", "Similarity") 
```
![Rplot2 5 7](https://github.com/user-attachments/assets/a0e2b93d-7ca2-47b0-a7a9-1e480c3032f9)

``` r
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
```
![Rplot2 5 8](https://github.com/user-attachments/assets/c89c9b95-0f7e-4684-87be-1c1a0737690c)

``` r
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob = as.data.frame(rbind(nat,inv)[,10]), glob1 = as.data.frame(nat[,10]), sp=as.data.frame(nat[which(nat[,11]==1),10]), R=1000, th.sp = 0) #Gridding the native niche
grid.clim.t.inv <-ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat,inv)[,10]),glob1=as.data.frame(inv[,10]),sp=as.data.frame(inv[which(inv[,11]==1),10]), R=1000, th.sp=0) #Gridding the invaded niche
t.dyn <- ecospat.niche.dyn.index(grid.clim.t.nat, grid.clim.t.inv)                                       
ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0, interest=2, title = "Niche Overlap", name.axis1 = "Average temperature") #Plot the niche dynamics along one gradient (T)
```
![Rplot2 5 9](https://github.com/user-attachments/assets/52eec957-f7a4-4374-ad8b-52c5a3013cb1)

## 2.6 Biotic Interactions

Species co-occurrence analysis with a presence-absence matrix
``` r
data <- ecospat.testData[c(9:16, 54:57)]
ecospat.co_occurrences(data)
```
![Rplot2 6](https://github.com/user-attachments/assets/1d8ac9a2-6069-4429-960f-226b0a1715f4)

Pairwise co-occurrence analysis with calculation of the C-score index
``` r
data <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
nperm <- 100
outpath <- getwd()
ecospat.Cscore(data, nperm, outpath)
```
![Rplot2 6 2](https://github.com/user-attachments/assets/fde6024e-9b29-4957-ae6e-403a8991f9f1)

## 2.7 Data Preparation 

Correlation plot of variables 
``` r
data <- ecospat.testData[,4:8]
ecospat.cor.plot(data)
```
![Rplot2 7](https://github.com/user-attachments/assets/b903a5f6-b089-42cf-8aff-9d52e777ca61)

Calibration and evaluation dataset
``` r
data <- ecospat.testData
caleval <- ecospat.caleval(data = ecospat.testData[53], xy=data[2:3], row.num = 1:nrow(data), nrep = 2, ratio = 0.7, disaggregate = 0.2, pseudoabs = 100, npres = 10, replace = FALSE)
head(caleval) #Obtaining an evaluation and calibration dataset with a desired ratio of disaggregation
```
# 3 CORE NICHE MODELLING

## 3.1 Model Evaluation

Presence-only evaluation indices - Boyce Index
``` r
fit <- ecospat.testData$glm_Saxifraga_oppositifolia
obs <- ecospat.testData$glm_Saxifraga_oppositifolia[which(ecospat.testData$Saxifraga_oppositifolia==1)]
ecospat.boyce(fit, obs, nclass = 0, window.w = "default", res = 100, PEplot = TRUE)$cor
```
![Rplot3 1](https://github.com/user-attachments/assets/4003da90-41da-4b2a-960b-d3f84c504407)

Accuracy of community prediction
``` r
eval <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
pred <- ecospat.testData[c(73:92)]
CommunityEval <- ecospat.CommunityEval(eval, pred, proba = TRUE, ntir = 5, verbose = T)
```
## 3.2 Spatial Predictions and Projections

ESM - Ensemble of Small Models
``` r
install.packages("biomod2", dependencies = TRUE) #dependencies added because not automatically installed
library(biomod2) 
xy <- inv[,1:2] #species occurrences
head(xy)
sp_occ <- inv[11]
current <- inv[3:7] #environment
head(current)
```
BIOMOD
``` r
t1 <- Sys.time()
sp <- 1
myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var=as.numeric(sp_occ[,sp]), expl.var=current, resp.xy=xy, resp.name=colnames(sp_occ)[sp]) 
``` 
Calibration of simple bivariate models
``` r
my.ESM <- ecospat.ESM.Modeling(data = myBiomodData, models = c("GLM"), NbRunEval = 2, DataSplit = 70, weighting.score = c("AUC"), parallel=F)
```
Evaluation and average of simple bivariate models to ESMs
``` r
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM, weighting.score = c("SomersD"),threshold = 0)
```
Projection of simple bivariate models into new space
``` r
my.ESM_proj_current <- ecospat.ESM.Projection(ESM.modeling.output = my.ESM, new.env = current)
```
Projection of calibrated ESMs into new space
``` r
my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my.ESM_proj_current, ESM.EnsembleModeling.output = my.ESM_EF)
```
## 3.3 Spatial Prediction of Communities 
``` r
proba <- ecospat.testData[,73:92]
sr <- as.data.frame(rowSums(proba))
```
## 3.4 SESAM Framework 
``` r
prr <- ecospat.SESAM.prr(proba, sr)
head(prr)[,1:4]
```
# 4 POST-MODELLING

## 4.1 Spatial Predictions of Species Assemblages
``` r
presence <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
pred <- ecospat.testData[c(73:92)]
nbpermut <- 100
outpath <- getwd()
ecospat.cons_Cscore(presence, pred, nbpermut, outpath)
```
![Rplot4 1](https://github.com/user-attachments/assets/9f3bca06-428c-4bfc-8381-e2e528bdb6fa)











