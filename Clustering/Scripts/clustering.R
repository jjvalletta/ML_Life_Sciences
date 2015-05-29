#-----------------------------------------------------------------------------#
# Title:    Scripts to produce plots for Clustering presentation
# Author:   John Joseph Valletta
# Date:     20/02/2015     
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Preamble
rm(list = setdiff(ls(), lsf.str())) # Remove all objects except functions (lsf.str() would remove only functions)
# rm(list=ls()) # Remove everything
library(scatterplot3d) # For 3D scatter plot
library(fields) # rdist
library(mclust) # gaussian mixture models
library(MASS) # mvrnorm
library(RColorBrewer)
library(cluster) # silhouette diagram
library(pracma) # meshgrid

#-----------------------------------------------------------------------------#
# Plot iris data set without class knowledge
plot(iris$Petal.Width, iris$Sepal.Length, pch=16, col="grey", xlab="petal width (cm)", ylab="sepal length (cm)")
scatterplot3d(x=iris$Sepal.Length, 
              y=iris$Sepal.Width, 
              z=iris$Petal.Length,
              xlab='x1',
              ylab="x2",
              zlab="x3",
              pch=16, color="blue")

#-----------------------------------------------------------------------------#
# Let's cluster data step by step
# 1) Select K centroids (K rows chosen at random)
# 2) Assign data point to closest centroid
# 3) Recalculates the centroids as the average of all data points in a cluster (i.e., the centroids are p-length mean vectors, where p is the number of variables)
# 4) Assigns data points to their closest centroids
# 5) Continues steps 3 and 4 until the observations are not reassigned or the maximum number of iterations (R uses 10 as a default) is reached.

set.seed(110)
x <- iris[, c(4, 1)]
N <- dim(x)[1]
D <- dim(x)[2]
k <- 3
maxIter <- 10
iiCentres <- sample(N, k)
# Normalise features for distance calculation
meanX <- apply(x, 2, mean)
stdX <- apply(x, 2, sd)
normX <- sweep(x, 2, meanX, FUN="-")
normX <- sweep(normX, 2, stdX, FUN="/")
normCentres <- normX[iiCentres, ]
clustAssign <- rep(NA, N)
for (i in seq(maxIter)) {
    xDist <- rdist(normCentres, normX) # Compute distance between centres and every pt.   
    clustAssign <- apply(xDist, 2, which.min) # Assign cluster 
    # Plot results
    pdf(paste("Iteration", i, ".pdf", sep=""), paper="a4r")
    plot(x[, 1], x[, 2], col=clustAssign, main=paste("Iteration", i), xlim=c(0, 2.6), 
         ylim=c(4, 8), xlab="petal width (cm)", ylab="sepal length (cm)")
    # Display unnormalised centres
    centres <- sweep(normCentres, 2, stdX, FUN="*")
    centres <- sweep(centres, 2, meanX, FUN="+")
    points(centres, pch=15, cex=2, col=seq(k))
    dev.off()
    # Recompute centres
    for (j in seq(k)) {
        normCentres[j, ] <- colMeans(normX[clustAssign==j, ])
    }
}

#-----------------------------------------------------------------------------#
# Let's plot cluster regions (not pretty at all but it will do for now)
xx <- seq(from=0, to=2.6, by=0.03) # row
yy <- seq(from=4, to=8, by=0.03) # column
zz <- meshgrid(xx, yy)

xDist <- rdist(centres, cbind(as.vector(zz$X), as.vector(zz$Y))) # Compute distance between centres and every pt.   
allClust <- apply(xDist, 2, which.min) # Assign cluster 
allClust <- matrix(allClust, nrow=length(xx), byrow=T)
filled.contour(xx, yy, allClust, plot.axes=points(x[, 1], x[, 2]), 
               levels=seq(k+1), col=c("grey", "red", "green"), 
               xlab="petal width (cm)", ylab="sepal length (cm)")

#-----------------------------------------------------------------------------#
# Let's try mixture of models
# Generate some ficticious data
set.seed(303)
classA <- mvrnorm(n=100, mu=c(5, 10), Sigma=diag(2))
classB <- mvrnorm(n=100, mu=c(10, 10), Sigma=4*diag(2))
# Plot to check
plot(classA, col="blue", xlim=c(1,12), ylim=c(5,15))
points(classB, col="red")
# 1 D case
k <- 2
x <-c(classA, classB)
fit <- Mclust(data=x, G=k, modelNames="V")
xGrid <- seq(0, 20, length=100)
pdfX <- matrix(data=NA, nrow=length(xGrid), ncol=k)
h<-hist(x, col="lightgrey", xlab="x", ylab="probability distribution", xlim=c(0, 20), ylim=c(0, 0.2), freq=FALSE, main="Density Plot")
for (i in seq(k)) {
    p <- fit$parameters$pro[i]
    mu <- fit$parameters$mean[i]
    sigma <- sqrt(fit$parameters$variance$sigmasq[i])
    pdfX[, i] <- p*dnorm(x=xGrid, mean=mu, sd=sigma)
    lines(xGrid, pdfX[, i], type="o", pch="", lty=2, lwd=2, col="blue")
}
lines(xGrid, rowSums(pdfX), type="o", pch="", lty=1, lwd=4)
# 2 D case
x <- rbind(classA, classB)
fit <- Mclust(data=x, G=k)
plot(fit, what="density", col = brewer.pal(n=7, name="RdPu"), lwd=3, xlab="x1", ylab="x2")
points(x[, 1], x[, 2], pch=20, col="grey")

#-----------------------------------------------------------------------------#
# Plot silhouette diagram
xData <- iris[, -5]
# K = 3
k <- 3
fit <- kmeans(xData, k)
silh <- silhouette(x=fit$cluster, dist=daisy(xData)^2)
pdf(file="silhK3.pdf", paper="a4r")
plot(silh, main="Silhouette plot k=3")
dev.off()
# K = 4
k <- 4
fit <- kmeans(xData, k)
silh <- silhouette(x=fit$cluster, dist=daisy(xData)^2)
pdf(file="silhK4.pdf", paper="a4r")
plot(silh, main="Silhouette plot k=4")
dev.off()

#-----------------------------------------------------------------------------#
# NbClust
library(NbClust)
xData <- iris[, -5]
NbClust(data=xData, distance="euclidean", min.nc=2, max.nc=10, method="kmeans", index="all")

#-----------------------------------------------------------------------------#
# Classification figure for Introduction slides
x <- seq(from=-5, to=5, by=0.01)
kSize <- 2.5
lwd <- 10
pdf("classification.pdf", paper="a4r", width=11.69, height=8.27)
par(mai=c(1, 1, 1, 1)) #c(bottom, left, top, right)
plot(x, dnorm(x, mean=-1.5, sd=1), col="tomato2", ylim=c(0, 0.45), lwd=lwd, type="l", 
     xlab="antibody level", ylab="probability density", cex.lab=kSize, cex.axis=kSize, cex.main=kSize, cex.sub=kSize,)
text(-1.5, 0.42, "susceptible", col="tomato2", cex=kSize, font=2)
lines(x, dnorm(x, mean=1.5, sd=1), col="skyblue", lwd=lwd, cex=kSize, font=2)
text(1.5, 0.42, "protected", col="skyblue", cex=kSize, font=2)
dev.off()