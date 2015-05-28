#-----------------------------------------------------------------------------#
# Title:    Scripts for lab sessions
# Author:   John Joseph Valletta
# Date:     20/02/2015     
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Preamble
rm(list = setdiff(ls(), lsf.str()))

#-----------------------------------------------------------------------------#
# 2) Getting started using an artificial dataset
#-----------------------------------------------------------------------------#
#Data: generate artificial data from a 3D multivariate gaussian distribution
#Task: cluster datapoints
#Method: $k$-means and gaussian mixture model
#-----------------------------------------------------------------------------#
# k-means
library(MASS) # mvrnorm (multivariate normal)
library(rgl) # plot3d
N <- 100 # number of data points
sigma2 <- 1 # variance of data points
covMatrix <- sigma2*diag(3) # assume the same variance in all three directions
groupA <- mvrnorm(n=N, mu=c(3, 3, 9), Sigma=covMatrix)
groupB <- mvrnorm(n=N, mu=c(3, 9, 3), Sigma=covMatrix)
groupC <- mvrnorm(n=N, mu=c(9, 9, 9), Sigma=covMatrix)
groupD <- mvrnorm(n=N, mu=c(9, 3, 3), Sigma=covMatrix)
groupE <- mvrnorm(n=N, mu=c(6, 6, 6), Sigma=covMatrix)
xTrain <- rbind(groupA, groupB, groupC, groupD, groupE) # training dataset is 5*N rows by 3 columns 
plot3d(xTrain, col="grey", xlab="x", ylab="y", zlab="z")

# k-means, let k vary from 2 to 10
kRange <-  seq(from=2, to=10, by=1)
intraClustSS <- rep(NA, length(kRange))  
interClustSS <- rep(NA, length(kRange))
for (k in kRange){
    fit <- kmeans(x=xTrain, centers=k)
    intraClustSS[k-1] <- fit$tot.withinss # it's (k-1) because k starts from 2
    interClustSS[k-1] <- fit$betweenss
}
minY <- min(c(intraClustSS, interClustSS))
maxY <- max(c(intraClustSS, interClustSS))
pdf("sumOfSquares.pdf", paper="a4r")
plot(kRange, interClustSS, type="o", pch=1, col="blue", lwd=3, lty=1, 
     ylim=c(minY, maxY), xlab="no. of clusters (k)", ylab="sum-of-squares")
points(kRange, intraClustSS, type="o", pch=1, col="red", lwd=3)
legend("right", c("inter-cluster (between)", "intra-cluster (within)"), bty="n",
       col=c("blue", "red"), pch=1, lwd=3, lty=1)
dev.off()
# Inter- and intra-cluster SS doesn't change after k=5
fit <- kmeans(x=xTrain, centers=5)
plot3d(xTrain, col=fit$cluster, xlab="x", ylab="y", zlab="z")
# rgl.postscript(filename="test.eps") # to print figure

#-----------------------------------------------------------------------------#
# gaussian mixture models
library(mclust) # gaussian mixture models
AIC <- rep(NA, length(kRange))  
BIC <- rep(NA, length(kRange))
for (k in kRange) {
    fit <- Mclust(data=xTrain, G=k, modelNames="EII")    
    BIC[k-1] <- fit$bic
    AIC[k-1] <- 2*fit$df - 2*fit$loglik # Have to compute this not implicitly returned by mclust
}
pdf("infoCriterion.pdf", paper="a4r")
plot(kRange, -AIC, type="o", pch=1, col="blue", lwd=3, lty=1, 
     xlab="no. of clusters (k)", ylab="information criterion") # Plot -ve AIC so to keep same scale as BIC
lines(kRange, BIC, type="o", pch=1, col="red", lwd=3)
legend("topleft", c("Akaike", "Bayesian"), bty="n",
       col=c("blue", "red"), pch=1, lwd=3, lty=1) 
dev.off()
# Could have done this instead of looping...
fit <- Mclust(data=xTrain, G=seq(10), modelNames="EII") 
summary(fit)
plot(fit, what="density", lwd=2)

#-----------------------------------------------------------------------------#
# 3) Clustering gene expression data
#-----------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R") # Install Bioconductor
biocLite() # Install core packages (will take a few mins)
biocLite("ALL") # Install the ALL package
library(ALL) # Load the ALL package
data(ALL) # Loads ALL dataset to workspace
xTrain <- exprs(ALL) # Extract the 12625 x 128 dataset
colnames(xTrain) <- ALL$BT # Replace patient ID by B or T to assess clustering

distance <- dist(as.matrix(t(xTrain), method="euclidean")) # Compute distance between patients
fit <- hclust(distance, method="complete") # Perform agglomerative hierarchical clustering using the "complete" linkage function
plot(fit, cex=0.5) # Plot the resultant Dendrogram 

distmethod <- function(x) as.dist(1-cor(x))
distance <- distmethod(xTrain)

#-----------------------------------------------------------------------------#
# 4) Species distribtuion modelling
#-----------------------------------------------------------------------------#
# Let us start by plotting where we have observed 
library(raster) # functions for gridded spatial data
library(dismo) # package containing the Bradypus variegatus dataset
library(rworldmap) # access to map of the world
dataPath <- file.path(system.file(package="dismo"), "ex") # path to where data is located
bradypusFilePath <- file.path(dataPath, "bradypus.csv") # absolute path to Bradypus dataset
data <-  read.table(file=bradypusFilePath, header=T, sep=",", skip=0) # read data
worldMap <- getMap(resolution = "low") # access world map
plot(worldMap, xlim = c(-85, -40), ylim = c(-25, 20), asp = 1, axes=T) # xlim = longtitude, ylim = latitude
points(data[, 2], data[, 3], pch=20, col="blue", cex=0.5) # plot a point for every observation

# Density estimation problem
# Step 1
grdFiles <- list.files(path=dataPath, pattern='grd', full.names=T) # read path of all .grd files
envData <- stack(grdFiles) # stacks all environmental variables
plot(envData) # plot to confirm data is OK
# Step 2
# Extract env data i.e [temperature, precipitation,...] = f(latitude, longitude)
xTrain <- extract(envData, data[, 2:3]) 
# Step 3
# Standardise data
meanXTrain <- apply(xTrain, 2, mean)
sdXTrain <- apply(xTrain, 2, sd)
xTrain <- sweep(xTrain, 2, meanXTrain, FUN="-")
xTrain <- sweep(xTrain, 2, sdXTrain, FUN="/")
# Step 4
# Fit data
library(mclust)
envVariables <- c("bio1","bio12", "biome") # chosen environmental variables
fit <- densityMclust(data=xTrain[, envVariables]) # fit distribution
# Step 5
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
plot(rasterDF, col=brewer.pal(n=9, name="Reds"), xlab="longtitude", ylab="latitude")
points(data[, 2], data[, 3], pch=20, col="blue", cex=0.5) # plot a point for every observation

#-----------------------------------------------------------------------------#
# 5) Predicting forest cover type from cartographic attributes
#-----------------------------------------------------------------------------#
library(RCurl) # To compose general HTTP requests 
dataPath <- getURL("https://dl.dropboxusercontent.com/u/57002389/ML_Life_Sciences/Data/ForestCoverData.csv")
data <- read.table(text=dataPath, header=T, sep=",", skip=0) # Read using a text connection
data <- data[complete.cases(data), ] # Only 2 cases of missing data so ignore
data$Cover_Type <- factor(data$Cover_Type)
NTOTAL <- dim(data)[1] # Total number of observations
NTRAIN <- 12000 # Let us use this many training datapoints
set.seed(101) # Just so we can reproduce results
trainIndices <- sample(x=NTOTAL, size=NTRAIN, replace=FALSE) 
# Let us try and tree first
library(rpart)
library(rattle) # For fancyRpartPlot to plot rpart tree nicely
fit <- rpart(Cover_Type ~ ., data=data, subset=trainIndices, method="class")
fancyRpartPlot(fit) # Plot decision tree
predClass = predict(fit, type="class") # Pred class on training dataset
confusionMatrix <- table(data$Cover_Type[trainIndices], predClass)
classError <- 1 - diag(prop.table(confusionMatrix, 1)) # Compute misclassification rates

library(randomForest)
fit <- randomForest(Cover_Type ~ ., data=data, subset=trainIndices, ntree=200, importance=TRUE) # Takes a few secs
fit$confusion # Confusion matrix
varImpPlot(fit) # Variable importance plot
# Testing set
predClass <- predict(fit, newdata=data[-trainIndices, ])
confusionMatrix <- table(data$Cover_Type[-trainIndices], predClass)
classError <- 1 - diag(prop.table(confusionMatrix, 1)) # Compute misclassification rates
#-----------------------------------------------------------------------------#
# 6) Classifying gene expression data
#-----------------------------------------------------------------------------#
library(ALL) # Load the ALL package
data(ALL) # Loads ALL dataset to workspace
xTrain <- t(exprs(ALL)) # Extract the 128 x 12625 dataset
yTrain <- factor(substr(ALL$BT,1,1)) # Class is either B or T
df <- data.frame(x=xTrain, y=yTrain) # Create a data frame, x - all inputs, y - class

# Consider first a logistic regression model
fit <- glm(y ~ ., data=df, family=binomial(link="logit"))

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
# 7) Identifying patients with Parkinsons's disease using speech recordings
#-----------------------------------------------------------------------------#
library(RCurl) # To compose general HTTP requests 
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
# Split into training/testing 60% - 40%
set.seed(101) # Just so we can reproduce results
trainIndices <- sample(x=dim(data)[1], size=round(dim(data)[1]*0.6), replace=FALSE)
# Train SVM
library(e1071) # An SVM library
fit <- svm(status ~., data=data, subset=trainIndices, type="C-classification", kernel="linear", probability=TRUE)
# Compute misclassification rates on testing dataset
predClass <- predict(fit, newdata=data[-trainIndices, ], type="class")
confusionMatrix <- table(data$status[-trainIndices], predClass)
classError <- 1 - diag(prop.table(confusionMatrix, 1))
# Plot Receiver operating characteristic (ROC)
library(ROCR) # ROC curves library
predProb <- predict(fit, newdata=data[-trainIndices, ], probability=TRUE) # Compute prediction probabilities rather than just class
# Lower level - Negative class (0), Upper level - Positive class (1)
predObj <- prediction(attr(predProb, "probabilities")[, 1], data$status[-trainIndices]) # prediction(probability of class=positive class, actual class label or factor)
ROCObj <- performance(predObj, "tpr", "fpr")
plot(ROCObj)
abline(a=0, b=1, col="red", lty=2, lwd=3) # Randomly guessing line
# Compute area under curve (AUC)
# AUC = 0.5 (randomly guessing class), AUC = 1 (perfect predictor)
areaUnderCurve <- performance(predObj, "auc") # Compute AUC
areaUnderCurve <- unlist(slot(areaUnderCurve, "y.values")) # Extract AUC (converting S4 class to vector)

library("corrgram")
corrgram(data[, iiPredictors], lower.panel=panel.shade, upper.panel=panel.conf)

#------------------------------------------------------------------------------#
# Plot to show true positive, false positive etc...
library(scales) # alpha
xTest <- seq(from=-12, to=12, by=0.01)
yA <- dnorm(xTest, mean=-2, sd=2)
yB <- dnorm(xTest, mean=+2, sd=4)
yMin <- 0
yMax <- 0.22
lwd <- 5
kSize <- 1.5 # factor to increase size of labels etc...
fontSize <- 2 # 2 is bold + italic
#------------------------------------------------------------------------------#
# The plot
pdf("DefinitionsPlot.pdf", paper="a4r", width=11.69, height=8.27)
plot(xTest, yA, col="skyblue", type="l", lwd=lwd, 
     xlab="predictor", ylab="frequency",
     cex.lab=kSize, cex.axis=kSize, cex.main=kSize, cex.sub=kSize, ylim=c(yMin, yMax))
lines(xTest, yB, col="tomato2", lwd=lwd)
text(x=-2, y=0.21, "no disease", col="skyblue", font=4, cex=kSize) # Bold + italic
text(x=+2, y=0.11, "disease", col="tomato2", font=4, cex=kSize) # Bold + italic
#legend("topleft", c("not protected", "protected"), col=c("tomato2", "skyblue"), lty=1, cex=kSize, lwd=lwd, bty="n")
segments(x0=-2.5, y0=0, x1=-2.5, y1=0.2, col="black", lwd=lwd, lty=2) # threshold line
# True Negative
polygon(c(xTest[xTest < -2.5], rev(xTest[xTest < -2.5])), c(yA[xTest < -2.5], rev(yB[xTest < -2.5])), col=alpha("skyblue", 0.2))
text(x=-4, y=0.06, "True \n Negative", col="skyblue", font=fontSize, cex=kSize)
# False Negative
polygon(c(xTest[xTest < -2.5], rev(xTest[xTest < -2.5])), c(yB[xTest < -2.5], 0*rev(yB[xTest < -2.5])), col=alpha("grey", 0.2))
text(x=-4, y=0.014, "False \n Negative", col="grey50", cex=kSize, font=fontSize)
# False Positive
polygon(c(xTest[xTest > -2.5], rev(xTest[xTest > -2.5])), c(yA[xTest > -2.5], 0*rev(yB[xTest > -2.5])), col=alpha("black", 0.2))
text(x=-0.5, y=0.06, "False \n Positive", col="black", cex=kSize, font=fontSize)
# True Positive
polygon(c(xTest[xTest > 0.5], rev(xTest[xTest > 0.5])), c(yB[xTest > 0.5], rev(yA[xTest > 0.5])), col=alpha("tomato2", 0.2))
text(x=3.5, y=0.06, "True \n Positive", col="tomato2", font=fontSize, cex=kSize)
# Classification cut-off
#text(x=-2.5, y=0.14, "threshold", col="black")
dev.off()