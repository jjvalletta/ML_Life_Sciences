#-----------------------------------------------------------------------------#
# Title:    Scripts for lab sessions
# Author:   John Joseph Valletta
# Date:     20/02/2015     
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Preamble
#-----------------------------------------------------------------------------#
rm(list = setdiff(ls(), lsf.str())) # Clean-up workspace

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# 1) Getting started using an artificial dataset
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Data: generate artificial data from a 3D multivariate gaussian distribution
# Task: cluster datapoints
# Method: $k$-means and gaussian mixture model
#-----------------------------------------------------------------------------#
# k-means
library(MASS) # mvrnorm (multivariate normal)
library(rgl) # plot3d
library(RColorBrewer)
N <- 100 # number of data points
sigma2 <- c(1, 2, 3, 4, 5) # variance of data points
kRange <-  seq(from=2, to=20, by=1)
intraClustSS <- matrix(data=NA, nrow=length(kRange), ncol=length(sigma2))  
interClustSS <- intraClustSS
iCnt <- 1
for (iSigma2 in sigma2) {
    covMatrix <- iSigma2*diag(3) # assume the same variance in all three directions
    groupA <- mvrnorm(n=N, mu=c(3, 3, 9), Sigma=covMatrix)
    groupB <- mvrnorm(n=N, mu=c(3, 9, 3), Sigma=covMatrix)
    groupC <- mvrnorm(n=N, mu=c(9, 9, 9), Sigma=covMatrix)
    groupD <- mvrnorm(n=N, mu=c(9, 3, 3), Sigma=covMatrix)
    groupE <- mvrnorm(n=N, mu=c(6, 6, 6), Sigma=covMatrix) 
    xTrain <- rbind(groupA, groupB, groupC, groupD, groupE) # training dataset is 5*N rows by 3 columns
    # k-means, let k vary from 2 to 20
    for (k in kRange){
        fit <- kmeans(x=xTrain, centers=k)
        intraClustSS[(k-1), iCnt] <- fit$tot.withinss # it's (k-1) because k starts from 2
        interClustSS[(k-1), iCnt] <- fit$betweenss
    }
    iCnt <- iCnt + 1
}
minY <- min(c(intraClustSS, interClustSS))
maxY <- max(c(intraClustSS, interClustSS))
pdf("Problem1Plot1.pdf", paper="a4r", width=11.69, height=8.27)
plotColour <- brewer.pal(n=length(sigma2)+1, name="Blues")
plot(kRange, interClustSS[, 1], type="o", pch=1, col=plotColour[2], lwd=2, lty=1, xlim=c(0, 30), 
     ylim=c(minY, maxY), xlab="no. of clusters (k)", ylab="sum-of-squares", main="Sum-of-squares")
text(x=23, y=interClustSS[length(kRange),1], paste("sigma2 =", toString(sigma2[1])), col=plotColour[2], font=2)
for (i in seq(length(sigma2)-1)) {
    lines(kRange, interClustSS[, (i+1)], type="o", pch=1, col=plotColour[i+2], lwd=2, lty=1)  
    text(x=23, y=interClustSS[length(kRange), (i+1)], paste("sigma2 =", toString(sigma2[i+1])), col=plotColour[i+2], font=2)
}
plotColour <- brewer.pal(n=length(sigma2)+1, name="Reds")
for (i in seq(length(sigma2))) {
    lines(kRange, intraClustSS[, i], type="o", pch=1, col=plotColour[i+1], lwd=2, lty=1) 
    text(x=23, y=intraClustSS[length(kRange), i], paste("sigma2 =", toString(sigma2[i])), col=plotColour[i+1], font=2)
}
legend("center", c("inter-cluster (between)", "intra-cluster (within)"), bty="n",
       col=c("blue", "red"), pch=1, lwd=2, lty=1)
dev.off()
# Inter- and intra-cluster SS doesn't change after k=5
fit <- kmeans(x=xTrain, centers=5)
plot3d(xTrain, col=fit$cluster, xlab="x", ylab="y", zlab="z")
# rgl.postscript(filename="test.eps") # to print figure
#-----------------------------------------------------------------------------#
# Gaussian mixture models
library(mclust) # gaussian mixture models
covMatrix <- 5*diag(3) # assume the same variance in all three directions
groupA <- mvrnorm(n=N, mu=c(3, 3, 9), Sigma=covMatrix)
groupB <- mvrnorm(n=N, mu=c(3, 9, 3), Sigma=covMatrix)
groupC <- mvrnorm(n=N, mu=c(9, 9, 9), Sigma=covMatrix)
groupD <- mvrnorm(n=N, mu=c(9, 3, 3), Sigma=covMatrix)
groupE <- mvrnorm(n=N, mu=c(6, 6, 6), Sigma=covMatrix) 
xTrain <- rbind(groupA, groupB, groupC, groupD, groupE) # training dataset is 5*N rows by 3 columns
AIC <- rep(NA, length(kRange))  
BIC <- rep(NA, length(kRange))
for (k in kRange) {
    fit <- Mclust(data=xTrain, G=k, modelNames="EII")    
    BIC[k-1] <- fit$bic
    AIC[k-1] <- 2*fit$df - 2*fit$loglik # Have to compute this not implicitly returned by mclust
}
pdf("infoCriterion.pdf", paper="a4r", width=11.69, height=8.27)
par(mai=c(0.8, 2.2, 0.8, 0.8))
lwd <- 8
kSize <- 3
plot(kRange, -AIC, type="o", pch=1, col="blue", lwd=lwd, lty=1, 
     xlab="no. of clusters (k)", ylab="information criterion", main="Information criterion",
     cex.lab=kSize, cex.axis=kSize, cex.main=kSize, cex.sub=kSize, cex=kSize) # Plot -ve AIC so to keep same scale as BIC
lines(kRange, BIC, type="o", pch=1, col="red", lwd=lwd, cex=kSize)
legend("bottomright", c("Akaike", "Bayesian"), bty="n",
       col=c("blue", "red"), pch=1, lwd=lwd, lty=1, cex=kSize) 
dev.off()
fit <- Mclust(data=xTrain, G=seq(20), modelNames="EII")
summary(fit)
plot(fit, what="density", lwd=2)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# 2) Clustering gene expression data
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Data: gene expression (data of T- and B-cell Acute Lymphocytic Leukemia)    
# Task: deduce the number of distinct phenotypes (unsupervised)
# Method: agglomerative hierarchical clustering
#-----------------------------------------------------------------------------#
# Retrieve the dataset
source("http://bioconductor.org/biocLite.R") # Install Bioconductor
biocLite() # Install core packages (will take a few mins)
biocLite("ALL") # Install the ALL package
library(ALL) # Load the ALL package
data(ALL) # Loads ALL dataset to workspace
xTrain <- exprs(ALL) # Extract the 12625 x 128 dataset
colnames(xTrain) <- ALL$BT # Replace patient ID by B or T to assess clustering
# Hierarchical clustering
pdf("Problem2Plot1.pdf", paper="a4r", width=11.69, height=8.27)
distance <- dist(as.matrix(t(xTrain), method="euclidean")) # Compute distance between patients
fit <- hclust(distance, method="complete") # Perform agglomerative hierarchical clustering using the "complete" linkage function
plot(fit, cex=0.5, main="Distance=Euclidean, Linkage=Complete") # Plot the resultant Dendrogram 
dev.off()
# Hierarchical clustering
pdf("Problem2Plot2.pdf", paper="a4r", width=11.69, height=8.27)
distance <- dist(as.matrix(t(xTrain), method="euclidean")) # Compute distance between patients
fit <- hclust(distance, method="single") # Perform agglomerative hierarchical clustering using the "complete" linkage function
plot(fit, cex=0.5, main="Distance=Euclidean, Linkage=Single") # Plot the resultant Dendrogram 
dev.off()
# Hierarchical clustering
pdf("Problem2Plot3.pdf", paper="a4r", width=11.69, height=8.27)
distance <- dist(as.matrix(t(xTrain), method="euclidean")) # Compute distance between patients
fit <- hclust(distance, method="average") # Perform agglomerative hierarchical clustering using the "complete" linkage function
plot(fit, cex=0.5, main="Distance=Euclidean, Linkage=Average") # Plot the resultant Dendrogram 
dev.off()
# Use correlation distance
distmethod <- function(x) as.dist(1-cor(x))
# Hierarchical clustering
pdf("Problem2Plot4.pdf", paper="a4r", width=11.69, height=8.27)
distance <- distmethod(xTrain)
fit <- hclust(distance, method="complete") # Perform agglomerative hierarchical clustering using the "complete" linkage function
plot(fit, cex=0.5, main="Distance=Correlation, Linkage=Complete") # Plot the resultant Dendrogram 
dev.off()
# Hierarchical clustering
pdf("Problem2Plot5.pdf", paper="a4r", width=11.69, height=8.27)
distance <- distmethod(xTrain)
fit <- hclust(distance, method="single") # Perform agglomerative hierarchical clustering using the "complete" linkage function
plot(fit, cex=0.5, main="Distance=Correlation, Linkage=Single") # Plot the resultant Dendrogram 
dev.off()
# Hierarchical clustering
pdf("Problem2Plot6.pdf", paper="a4r", width=11.69, height=8.27)
distance <- distmethod(xTrain)
fit <- hclust(distance, method="average") # Perform agglomerative hierarchical clustering using the "complete" linkage function
plot(fit, cex=0.5, main="Distance=Correlation, Linkage=Average") # Plot the resultant Dendrogram 
dev.off()

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# 3) Species distribtuion modelling
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Data: Bradypus variegatus (brown-throated slot) geographic distribution data    
# Task: produce a map showing the likely geographic distribution of Bradypus variegatus (unsupervised)
# Method: Gaussian mixture models
#-----------------------------------------------------------------------------#
# Retrieve and plot data
library(raster) # functions for gridded spatial data
library(dismo) # package containing the Bradypus variegatus dataset
library(rworldmap) # access to map of the world
dataPath <- file.path(system.file(package="dismo"), "ex") # path to where data is located
bradypusFilePath <- file.path(dataPath, "bradypus.csv") # absolute path to Bradypus dataset
data <-  read.table(file=bradypusFilePath, header=T, sep=",", skip=0) # read data
worldMap <- getMap(resolution = "low") # access world map
plot(worldMap, xlim = c(-85, -40), ylim = c(-25, 20), asp = 1, axes=T) # xlim = longtitude, ylim = latitude
points(data[, 2], data[, 3], pch=20, col="blue", cex=0.5) # plot a point for every observation
#-----------------------------------------------------------------------------#
# Density estimation problem
#-----------------------------------------------------------------------------#
# Step 1 - Read in and plot environmental variables
grdFiles <- list.files(path=dataPath, pattern='grd', full.names=T) # read path of all .grd files
envData <- stack(grdFiles) # stacks all environmental variables
plot(envData) # plot to confirm data is OK
#-----------------------------------------------------------------------------#
# Step 2 - Extract env data i.e [temperature, precipitation,...] = f(latitude, longitude)
xTrain <- extract(envData, data[, 2:3]) 
#-----------------------------------------------------------------------------#
# Step 3 - Standardise data
meanXTrain <- apply(xTrain, 2, mean)
sdXTrain <- apply(xTrain, 2, sd)
xTrain <- sweep(xTrain, 2, meanXTrain, FUN="-")
xTrain <- sweep(xTrain, 2, sdXTrain, FUN="/")
#-----------------------------------------------------------------------------#
# Step 4 - Fit a distribution
library(mclust)
envVariables <- c("bio1", "bio12", "biome") # chosen environmental variables
fit <- densityMclust(data=xTrain[, envVariables]) # fit distribution
#-----------------------------------------------------------------------------#
# Step 5 - Compute fitted distribution over geographical grid of interest (xTest)
# Standardise test data (xTest is the whole grid not just the 116 obs.)
xTest <- as.data.frame(subset(envData, subset=envVariables)) 
xTest <- sweep(xTest, 2, meanXTrain[envVariables], FUN="-")
xTest <- sweep(xTest, 2, sdXTrain[envVariables], FUN="/")
# Read the distribution at the test points
notNA <- complete.cases(xTest) # ignore sea points (NA = sea) 
pred <- rep(NA, dim(xTest)[1])
pred[notNA] <- predict(fit, newdata=xTest[notNA, ])
# Create a data frame which consists values for the distribution for every grid cell
df <- data.frame(pred=as.matrix(pred, nrow=nrow(envData), ncol=ncol(envData)))
coordinates(df) <- coordinates(envData) # set spatial coordinates (same as envData)
gridded(df) <- TRUE # coerce to SpatialPixelsDataFrame
rasterDF <- raster(df) # coerce to raster
library(RColorBrewer) # package for pretty colour palettes
pdf("Problem3Plot1.pdf", paper="a4r", width=11.69, height=8.27)
plot(rasterDF, col=brewer.pal(n=9, name="Reds"), xlab="longtitude", ylab="latitude",
     main="EnvVariables=Mean Temperature, Precipitation and Habitat Type")
points(data[, 2], data[, 3], pch=20, col="blue", cex=0.5) # plot a point for every observation
dev.off()

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# 4) Predicting forest cover type from cartographic attributes
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Data: forest cover type   
# Task: predict forest cover type for a set of cartographic measurements (supervised)
# Method: decision trees and random forests
#-----------------------------------------------------------------------------#
# Retrieve dataset
library(RCurl) # To compose general HTTP requests 
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE)) # To avoid SSL certificate problem 
dataPath <- getURL("https://dl.dropboxusercontent.com/u/57002389/ML_Life_Sciences/Data/ForestCoverData.csv")
data <- read.table(text=dataPath, header=T, sep=",", skip=0) # Read using a text connection
data <- data[complete.cases(data), ] # Only 2 cases of missing data so ignore
data$Cover_Type <- factor(data$Cover_Type)
NTOTAL <- dim(data)[1] # Total number of observations
NTRAIN <- 12000 # Let us use this many training datapoints
set.seed(101) # Just so we can reproduce results
trainIndices <- sample(x=NTOTAL, size=NTRAIN, replace=FALSE) 
# Fit a simple decision tree
library(rpart)
library(rattle) # For fancyRpartPlot to plot rpart tree nicely
fit <- rpart(Cover_Type ~ ., data=data, subset=trainIndices, method="class")
fancyRpartPlot(fit) # Plot decision tree
predClass <- predict(fit, type="class") # Pred class on training dataset
confusionMatrix <- table(data$Cover_Type[trainIndices], predClass)
classError <- 1 - diag(prop.table(confusionMatrix, 1)) # Compute misclassification rates
# Fit a Random Forest
library(randomForest)
fit <- randomForest(Cover_Type ~ ., data=data, subset=trainIndices, ntree=200, importance=TRUE) # Takes a few secs
fit$confusion # Confusion matrix
varImpPlot(fit) # Variable importance plot
# Compute predictive performance on test set
predClass <- predict(fit, newdata=data[-trainIndices, ])
confusionMatrix <- table(data$Cover_Type[-trainIndices], predClass)
classError <- 1 - diag(prop.table(confusionMatrix, 1)) # Compute misclassification rates

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# 5) Classifying patients by their gene expression signature
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Data: gene expression (data of T- and B-cell Acute Lymphocytic Leukemia)     
# Task: predict a patient's phenotype using their gene expression signature (supervised)
# Method: decision trees
#-----------------------------------------------------------------------------#
# Retrieve the dataset
library(ALL) # Load the ALL package
data(ALL) # Loads ALL dataset to workspace
xTrain <- t(exprs(ALL)) # Extract the 128 x 12625 dataset
yTrain <- factor(substr(ALL$BT,1,1)) # Class is either B or T
df <- data.frame(x=xTrain, y=yTrain) # Create a data frame, x - all inputs, y - class
# Consider first a logistic regression model
# fit <- glm(y ~ ., data=df, family=binomial(link="logit"))
# Fit a decision tree
library(rpart) #
library(rpart.plot) # for prp; plots an rpart tree
library(rattle) # for fancyRpartPlot function plots an rpart tree nicely
xTrain <- xTrain[, seq(50)]
df <- data.frame(x=xTrain, y=yTrain)
fit <- rpart(y ~ ., data=df, method="class")
prp(fit, varlen=30) # Nice looking tree
fancyRpartPlot(fit) # Coloured tree
predClass = predict(fit, type="class") # Pred class on training dataset
confusionMatrix <- table(yTrain, predClass)
classError <- 1 - diag(prop.table(confusionMatrix, 1)) # Compute misclassification rates

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# 6) Identifying patients with Parkinsons's disease from speech recordings
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Data: audio recordings   
# Task: predict presence of Parkinson's disease (supervised)
# Method: support vector machine (SVM)
#-----------------------------------------------------------------------------#
# Retrieve dataset
library(RCurl) # To compose general HTTP requests 
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE)) # To avoid SSL certificate problem 
dataPath <- getURL("https://archive.ics.uci.edu/ml/machine-learning-databases/parkinsons/parkinsons.data")
data <- read.table(text=dataPath, header=T, sep=",", skip=0, row.names=1) # Read using a text connection
data$status <- factor(data$status) # Convert to factor 0/1 absence/presence of Parkinson's
iiPredictors <- which(names(data)!="status") # Column number of all predictors i.e everything apart from "status"
# Scale all features to -1 to +1 (all data columns except "status")
# x' = 2*(x - xmin)/(xmax - xmin) - 1
xMin <- apply(data[, iiPredictors], 2, min)
xMax <- apply(data[, iiPredictors], 2, max)
data[, iiPredictors] <- sweep(data[, iiPredictors], 2, xMin, FUN="-")
data[, iiPredictors] <- 2*data[, iiPredictors]
data[, iiPredictors] <- sweep(data[, iiPredictors], 2, (xMax-xMin), FUN="/")
data[, iiPredictors] <- data[, iiPredictors] - 1
#-----------------------------------------------------------------------------#
# Split into training/testing 60% - 40%
set.seed(101) # Just so we can reproduce results
trainIndices <- sample(x=dim(data)[1], size=round(dim(data)[1]*0.6), replace=FALSE)
#-----------------------------------------------------------------------------#
# Train SVM
library(e1071) # An SVM library
fit <- svm(status ~., data=data, subset=trainIndices, type="C-classification", kernel="linear", probability=TRUE)
#-----------------------------------------------------------------------------#
# Compute misclassification rates on testing dataset
predClass <- predict(fit, newdata=data[-trainIndices, ], type="class")
confusionMatrix <- table(data$status[-trainIndices], predClass)
classError <- 1 - diag(prop.table(confusionMatrix, 1))
#-----------------------------------------------------------------------------#
# Plot Receiver operating characteristic (ROC)
library(ROCR) # ROC curves library
predProb <- predict(fit, newdata=data[-trainIndices, ], probability=TRUE) # Compute prediction probabilities rather than just class
# Lower level - Negative class (0), Upper level - Positive class (1)
predObj <- prediction(attr(predProb, "probabilities")[, 1], data$status[-trainIndices]) # prediction(probability of class=positive class, actual class label or factor)
ROCObj <- performance(predObj, "tpr", "fpr")
plot(ROCObj)
abline(a=0, b=1, col="red", lty=2, lwd=3) # Randomly guessing line
#-----------------------------------------------------------------------------#
# Compute area under curve (AUC)
# AUC = 0.5 (randomly guessing class), AUC = 1 (perfect predictor)
areaUnderCurve <- performance(predObj, "auc") # Compute AUC
areaUnderCurve <- unlist(slot(areaUnderCurve, "y.values")) # Extract AUC (converting S4 class to vector)
#-----------------------------------------------------------------------------#
library("corrgram")
corrgram(data[, iiPredictors], lower.panel=panel.shade, upper.panel=panel.conf)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# The End
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#