#-----------------------------------------------------------------------------#
# Title:    Scripts to produce plots for Classification presentation
# Author:   John Joseph Valletta
# Date:     20/02/2015     
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Preamble
#-----------------------------------------------------------------------------#
rm(list = setdiff(ls(), lsf.str())) # Remove all objects except functions (lsf.str() would remove only functions)
# rm(list=ls()) # Remove everything
library(tree)
library(rpart) #
library(rpart.plot) # for prp; plots an rpart tree
library(rattle) # for fancyRpartPlot function plots an rpart tree nicely
library(class) # for k-nearest neighbour
library(randomForest) # random Forest
library(MASS) # mvrnorm
library(ellipse) # ellipse
library(e1071) # svm

#-----------------------------------------------------------------------------#
# Plot iris data set without class knowledge
#-----------------------------------------------------------------------------#
speciesCol <- c("red", "black", "green")
speciesSym <- c(21, 22, 24) # col = just border of symbol, bg = set filling
pdf(file="irisPlot.pdf", paper="a4r")
plot(iris$Petal.Width, iris$Sepal.Length, pch=speciesSym[iris$Species], col=speciesCol[iris$Species], xlab="petal width (cm)", ylab="sepal length (cm)", bg=speciesCol[iris$Species], type="p")
legend("topleft", c("setosa", "versicolor", "virginica"), bty="n",
       col=speciesCol, pch=speciesSym, pt.bg=speciesCol)
dev.off()

#-----------------------------------------------------------------------------#
# Overfit/underfit
#-----------------------------------------------------------------------------#
set.seed(10)
xSeq <- seq(from=-pi, to=pi, by=0.01) # test set
xTest <- data.frame(x=xSeq)
x <- seq(from=-pi, to=pi, by=0.5)
y <- sin(x)
y <- y + rnorm(length(y), mean=0, sd=0.3)
df <- data.frame(x, y)
comment <- c("(underfitting)", "(parsimonious)", "(overfitting)")
k <- 1
for (i in c(1, 3, 12)) {
    fit <- lm(y ~ poly(x, i), data=df)   
    mainTitle <- paste(i, "degree polynomial", comment[k]); k<-k+1;
    pdf(file=paste("polyFit", i, ".pdf", sep=""), paper="a4r")
    plot(xSeq, sin(xSeq), col="grey", type="l", lwd=4, main=mainTitle, 
         xlab="x", ylab="y", xlim=c(-3.2, 3.2), ylim=c(-1.3, 1.3))
    points(y ~ x, data=df, pch=16, col="black", cex=1.5)
    lines(xSeq, predict(fit, newdata=xTest), col="blue", lwd=4)
    legend("topleft", c("underlying function (true)", "observed data points", "model fit"), bty="n",
           col=c("grey", "black", "blue"), pch=c(NA, 16, NA), lwd=4, lty=c(1, NA, 1))
    dev.off()
}

#-----------------------------------------------------------------------------#
# Decision trees
#-----------------------------------------------------------------------------#
fit <- tree(Species ~ Petal.Width + Sepal.Length, data=iris)
plot(fit)
text(fit, pretty=0) # Note how the last level is completely redundant
# Fancier trees
fit2 <- rpart(Species ~ Petal.Width + Sepal.Length, data=iris, method="class")
pdf(file="crudeTree.pdf", paper="a4r")
prp(fit2, varlen=30) # Nice looking tree
dev.off()
pdf(file="niceTree.pdf", paper="a4r")
fancyRpartPlot(fit2) # Coloured tree
dev.off()
# Plot the decision boundaries
pdf(file="irisPlotBound.pdf", paper="a4r")
plot(iris$Petal.Width, iris$Sepal.Length, pch=speciesSym[iris$Species], col=speciesCol[iris$Species], xlab="petal width (cm)", ylab="sepal length (cm)", bg=speciesCol[iris$Species], type="p")
legend("topleft", c("setosa", "versicolor", "virginica"), bty="n",
       col=speciesCol, pch=speciesSym, pt.bg=speciesCol)
abline(v=c(0.8, 1.75), col="grey", lwd=4)
dev.off()

#-----------------------------------------------------------------------------#
# Let's do kNN for varying levels of k
#-----------------------------------------------------------------------------#
xTrain <-  iris[, c(4, 1)]
yTrain <- iris$Species
# Create grid of values
dx <- 0.01
dy <- 0.01
eps <- 0.1
xMin <- min(xTrain[ , 1]) - eps
xMax <- max(xTrain[ , 1]) + eps
yMin <- min(xTrain[ , 2]) - eps
yMax <- max(xTrain[ , 2]) + eps
xSeq <- seq(from=xMin, to=xMax, by=dx)
ySeq <- seq(from=yMin, to=yMax, by=dy)
xTest <- expand.grid(x=xSeq, y=ySeq)

# Fit for various k
for (k in c(1, 5, 15, 30)) {
    fit <- knn(train=xTrain, test=xTest, cl=yTrain, k=k, use.all=FALSE)
    classAssign <- matrix(data=as.numeric(fit), nrow=length(xSeq), ncol=length(ySeq), byrow=FALSE)
    pdf(file=paste(k, "nearestNeighbour.pdf", sep=""), paper="a4r")
    image(xSeq, ySeq, classAssign, xlab="petal width (cm)", 
          ylab="sepal length (cm)", main=paste(k, "-nearest neighbour", sep=""), 
          col=c("tomato", "grey70", "lightgreen"))
    points(iris$Petal.Width, iris$Sepal.Length, pch=speciesSym[iris$Species], lwd=2, bg=speciesCol[iris$Species], col="black")
    legend("topleft", c("setosa", "versicolor", "virginica"), bty="n",
           col="black", pch=speciesSym, pt.bg=speciesCol)
    dev.off()
    
}
# Alternate Plot
# I don't really like the faint white grid lines that filled.contour produces, 
# it drives me nuts, and it doesn't seem like there's an easy solution 
# so use image instead
# filled.contour(xSeq, ySeq, classAssign, levels=c(-0.5, 1.5, 2.5, 3.5), xlab="petal width (cm)", 
#         ylab="sepal length (cm)", main=paste(k, "-nearest neighbour", sep=""), 
#         col=c("red", "grey35", "green"),
#         plot.axes=points(iris$Petal.Width, iris$Sepal.Length, pch=21, lwd=3, bg=c("red", "black", "green")[iris$Species]))

#-----------------------------------------------------------------------------#
# Random forests
#-----------------------------------------------------------------------------#
# Let us construct 4 trees, just for demo purposes
# I'm aware that randomised node optimisation occurs at every level of the tree
# but for this e.g with only two levels, it suffice to assume the same node 
# optimisation for every level, again just for demoing the technique
fit <- rpart(Species ~ Petal.Length + Sepal.Width, data=iris, method="class")
pdf(file="treeFit1.pdf", paper="a4r")
fancyRpartPlot(fit)
dev.off()
fit <- rpart(Species ~ Sepal.Length + Sepal.Width, data=iris, method="class")
pdf(file="treeFit2.pdf", paper="a4r")
fancyRpartPlot(fit)
dev.off()

fit <- randomForest(Species ~., data=iris, ntree=100, importance=TRUE)
pdf(file="varImportance.pdf", paper="a4r")
par(pty="s")
varImpPlot(fit, main="")
dev.off()
fit$confusion
par(pty="m")

#-----------------------------------------------------------------------------#
# Support Vector Machines (SVMs)
#-----------------------------------------------------------------------------#
# Which is the best separating line?
set.seed(65)
classA <- mvrnorm(n=10, mu=c(4, 6), Sigma=matrix(data=c(3,1,0.5,1), 2, 2))
classB <- mvrnorm(n=10, mu=c(10, 4), Sigma=matrix(data=c(3,1,0.5,1), 2, 2))
pdf(file="svmSepLine1.pdf", paper="a4r")
plot(classB, col="red", xlim=c(2, 14), ylim=c(2, 9), xlab="x1", ylab="x2", pch=16, cex=2)
points(classA, col="blue", pch=16, cex=2)
abline(a=-15, b=3, col="grey", lwd=4)
dev.off()
pdf(file="svmSepLine2.pdf", paper="a4r")
plot(classB, col="red", xlim=c(2, 14), ylim=c(2, 9), xlab="x1", ylab="x2", pch=16, cex=2)
points(classA, col="blue", pch=16, cex=2)
abline(a=0, b=0.9, col="grey", lwd=4)
dev.off()
pdf(file="svmSepLine3.pdf", paper="a4r")
plot(classB, col="red", xlim=c(2, 14), ylim=c(2, 9), xlab="x1", ylab="x2", pch=16, cex=2)
points(classA, col="blue", pch=16, cex=2)
abline(a=-0.5, b=0.6, col="grey", lwd=4)
dev.off()
pdf(file="svmSepLine4.pdf", paper="a4r")
plot(classB, col="red", xlim=c(2, 14), ylim=c(2, 9), xlab="x1", ylab="x2", pch=16, cex=2)
points(classA, col="blue", pch=16, cex=2)
abline(v=6.4, col="grey", lwd=4)
dev.off()

# Support vectors
pdf(file="svmLinear.pdf", paper="a4r")
plot(classB, col="red", xlim=c(2, 14), ylim=c(2, 9), xlab="x1", ylab="x2", pch=16, cex=2)
points(classA, col="blue", pch=16, cex=2)
points(classB[c(7,10), ], col="black", pch=21, cex=3, bg="red", lwd=5)
points(classA[c(2,7), ], col="black", pch=21, cex=3, bg="blue", lwd=5)
abline(a=1.3, b=0.75, col="grey", lwd=4, lty=2)
abline(a=-2.5, b=0.75, col="grey", lwd=4, lty=2)
abline(a=-0.6, b=0.75, col="black", lwd=4, lty=1)
dev.off()

# Non linearly separable
set.seed(351)
classA <- mvrnorm(n=20, mu=c(8, 8), Sigma=matrix(data=c(0.8,0,0,0.35), 2, 2))
classB <- ellipse(matrix(data=c(1.5,0,0,1), 2, 2), centre=c(8,8), npoints=40)
classB <- classB + rnorm(n=20, mean=0, sd=0.3)
sepLine <- ellipse(matrix(data=c(0.8,0,0,0.35), 2, 2), centre=c(8,8), npoints=200)
pdf(file="svmNonLinear1.pdf", paper="a4r")
plot(classB, col="red", xlim=c(5, 11), ylim=c(5, 11), xlab="x1", ylab="x2", pch=16, cex=2)
points(classA, col="blue", pch=16, cex=2)
dev.off()
pdf(file="svmNonLinear2.pdf", paper="a4r")
plot(classB, col="red", xlim=c(5, 11), ylim=c(5, 11), xlab="x1", ylab="x2", pch=16, cex=2)
lines(sepLine, col="black", lwd=4)
points(classA, col="blue", pch=16, cex=2)
dev.off()

# Non linearly separable - simple example
set.seed(50)
classA1 <- rnorm(10, mean=-5, sd=1)
classA2 <- rnorm(10, mean=5, sd=1)
classA <- c(classA1, classA2)
classB <- rnorm(20, mean=0, sd=1)
onesA <- rep(1, length(classA))
onesB <- rep(1, length(classB))
pdf(file="svmSimple1.pdf", paper="a4r")
plot(x=classB, y=onesB, col="red", xlab="x", ylab="", xlim=c(-7, 7), ylim=c(0, 2), pch=16, cex=2, 
     main="1D (original) data is not linearly separable")
points(x=classA, y=onesA, col="blue", pch=16, cex=2)
dev.off()
pdf(file="svmSimple2.pdf", paper="a4r")
plot(x=classB, y=classB^2, col="red", xlab="x", ylab="x^2", xlim=c(-7, 7), ylim=c(0, 50), pch=16, cex=2, 
     main="2D (transformed) data is now linearly separable")
points(x=classA, y=classA^2, col="blue", pch=16, cex=2)
abline(h=11, col="black", lwd=4)
dev.off()

# Iris example
# Fit using linear kernel
fitLinear <- svm(Species ~ Petal.Width + Sepal.Length, data=iris, type="C-classification", kernel="linear")
fitRBF <- svm(Species ~ Petal.Width + Sepal.Length, data=iris, type="C-classification", kernel="radial")
# Create grid of values to plot
dx <- 0.01
dy <- 0.01
eps <- 0.1
xMin <- min(xTrain[ , 1]) - eps
xMax <- max(xTrain[ , 1]) + eps
yMin <- min(xTrain[ , 2]) - eps
yMax <- max(xTrain[ , 2]) + eps
xSeq <- seq(from=xMin, to=xMax, by=dx)
ySeq <- seq(from=yMin, to=yMax, by=dy)
xTest <- expand.grid(x=xSeq, y=ySeq)
colnames(xTest) <- c("Petal.Width", "Sepal.Length")
# Linear kernel
classAssign <- predict(fitLinear, newdata=xTest)
classAssign <- matrix(data=as.numeric(classAssign), nrow=length(xSeq), ncol=length(ySeq), byrow=FALSE)
pdf(file="svmLinear.pdf", paper="a4r")
image(xSeq, ySeq, classAssign, xlab="petal width (cm)", 
      ylab="sepal length (cm)", main="Linear kernel", 
      col=c("tomato", "grey70", "lightgreen"))
points(iris$Petal.Width, iris$Sepal.Length, pch=speciesSym[iris$Species], lwd=2, bg=speciesCol[iris$Species], col="black")
legend("topleft", c("setosa", "versicolor", "virginica"), bty="n",
       col="black", pch=speciesSym, pt.bg=speciesCol)
dev.off()
# RBF kernel
classAssign <- predict(fitRBF, newdata=xTest)
classAssign <- matrix(data=as.numeric(classAssign), nrow=length(xSeq), ncol=length(ySeq), byrow=FALSE)
pdf(file="svmRBF.pdf", paper="a4r")
image(xSeq, ySeq, classAssign, xlab="petal width (cm)", 
      ylab="sepal length (cm)", main="Radial Basis Function (Gaussian) kernel", 
      col=c("tomato", "grey70", "lightgreen"))
points(iris$Petal.Width, iris$Sepal.Length, pch=speciesSym[iris$Species], lwd=2, bg=speciesCol[iris$Species], col="black")
legend("topleft", c("setosa", "versicolor", "virginica"), bty="n",
       col="black", pch=speciesSym, pt.bg=speciesCol)
dev.off()
