#-----------------------------------------------------------------------------#
# Title:    Scripts to produce plots for Practical Issues presentation
# Author:   John Joseph Valletta
# Date:     05/03/2015     
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Preamble
#-----------------------------------------------------------------------------#
rm(list = setdiff(ls(), lsf.str())) # Remove all objects except functions (lsf.str() would remove only functions)

#-----------------------------------------------------------------------------#
# Overfit/underfit
#-----------------------------------------------------------------------------#
set.seed(10)
# Training set
xTrain <- seq(from=-pi, to=pi, by=0.5)
yTrain <- sin(xTrain)
yTrain <- yTrain + rnorm(length(yTrain), mean=0, sd=0.3)
# Validation set
xVal <- runif(n=length(xTrain), min=-pi, max=pi)
xVal <- data.frame(x=xVal)
yVal <- sin(xVal)
yVal <- yVal + rnorm(length(yVal), mean=0, sd=0.3)
# Testing set
xSeq <- seq(from=-pi, to=pi, by=0.01) # test set
xTest <- data.frame(x=xSeq)
df <- data.frame(x=xTrain, y=yTrain)
# Compute RMSE for train and val dataset
rmseTrain <- rep(NA, 10)
rmseVal <- rep(NA, 10)
for (i in seq(10)) {
    fit <- lm(y ~ poly(x, i), data=df)   
    rmseTrain[i] <- sqrt(mean((df$y - fitted(fit))^2))
    rmseVal[i] <- sqrt(mean((yVal - predict(fit, newdata=xVal))^2))
}
pdf(file="trainvsval.pdf", paper="a4r")
lwd <- 8
kSize <- 2
plot(rmseTrain, type="o", pch=1, col="blue", lwd=lwd, xlim=c(0.8,10.2), ylim=c(-0.01,0.55),
     xlab="polynomial degree (model complexity)", ylab="root mean squared errror (rmse)",
     cex.lab=kSize, cex.axis=kSize, cex.main=kSize, cex.sub=kSize, cex=kSize)
points(rmseVal, type="o", pch=1, col="red", lwd=lwd, cex=kSize)
abline(h=0, col="grey", lwd=lwd, lty=2, cex=kSize)
text(9, 0.13, "training", col="blue", cex=kSize, font=2)
text(9, 0.295, "validation", col="red", cex=kSize, font=2)
legend("topright", c("training dataset", "validation dataset"), bty="n",
       col=c("blue", "red"), pch=1, lwd=lwd, lty=1, cex=kSize)
dev.off()
