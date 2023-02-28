# A regression model with p >> n
# The data set consists of 240 samples from patients with diffuse large B-cell lymphoma (DLBCL),
# With gene expression measurements for 7399 genes. The outcome is survived time, either observed or right censored.

# We propose a multiple linear regression model for the n=138 patients died before the end of the study, to explain the logarithm of the survival time as a linear function of the expression of the p=7399 genes.
#n.total = 240
#n.train = 138
#n.val = n.total - n.train
require(glmnet)
express <- read.csv("journal.pbio.0020108.sd012.CSV",header=FALSE)
surv <- read.csv("journal.pbio.0020108.sd013.CSV",header=FALSE)
death <-(surv[,2]==1)
log.surv.death <- log(surv[death,1]+.05)
expr.death <- as.matrix(t(express[,death]))

live <- (surv[,2]==0)
log.surv.live <- log(surv[live,1]+.05)
expr.live <- as.matrix(t(express[,live]))
X.val.noscale <- expr.live
Y.val.noscale <- log.surv.live
X.val.scale <- scale(X.val.noscale,center=TRUE,scale=TRUE)
Y.val.scale <- scale(Y.val.noscale,center=TRUE,scale=TRUE)
n.val <- dim(X.val.scale)
# Use glmnet and cv.glmnet to obtain the Lasso estimation for regressing log.surv against expr. How many coefficient different from zero are in the Lasso estimator?
# Illustrate the result with two graphics

X.noscale = expr.death
Y.noscale = log.surv.death
X.scale = scale(expr.death,center=TRUE,scale=TRUE)
Y.scale = scale(log.surv.death,center=TRUE,scale=TRUE)
n <- dim(X.scale)[1]
p <- dim(X.scale)[2]
lasso.1 <- glmnet(X.scale,Y.scale,standandize = FALSE, intercept = FALSE)# Lasso.1 corresponds to the scaled data
cv.lasso.1 <- cv.glmnet(X.scale,Y.scale,standandize= FALSE, intercept = FALSE) # by default nfolds = 10
print(cv.lasso.1$lambda.min)
print(cv.lasso.1$lambda.1se)

lasso.2 <- glmnet(X.noscale,Y.noscale, standardize=TRUE, intercept=TRUE)
cv.lasso.2 <- cv.glmnet(X.noscale,Y.noscale, standardize=TRUE, intercept=TRUE)
#print(coef(lasso.2,s=cv.lasso.2$lambda.1se))
print(cv.lasso.2$lambda.min)
print(cv.lasso.2$lambda.1se)
mean.Y <- mean(Y.noscale)
mean.X <- apply(X.noscale,2,mean)
sd.X <- apply(X.noscale,2,sd)

# intercept:
# a) No centering, no scaling
coef(lasso.2,s=cv.lasso.2$lambda.1se)[1]
# b) from the fitted model centering and scaling
mean.Y - sum((mean.X/sd.X) * coef(lasso.1,s=cv.lasso.2$lambda.1se)[-1])

# coefficients:
# a) No centering, no scaling
coef(lasso.2,s=cv.lasso.2$lambda.1se)[-1]
# b) from the fitted model centering and scaling
coef(lasso.1,s=cv.lasso.2$lambda.1se)[-1] / sd.X

# cv.glment returns an object 
# lambda: the values of lambda used in the fits
# cvm: the mean cross-validated erro- a vector of length(lambda)
# lambda.min: value of lambda that gives minimum cvm
# lambda.1se: largest value of lambda such that the error is within 1 standard error of minimum
op <- par(mfrow=c(2,1))
plot(cv.lasso.1)
plot(lasso.1,xvar="lambda")
abline(v=log(cv.lasso.1$lambda.min),col=2,lty=2) # lambda min = 0.3650
abline(v=log(cv.lasso.1$lambda.1se),col=2,lty=2) # lambda.1se = 0.6088
print(coef(lasso.1,s=cv.lasso.1$lambda.min))
print(coef(lasso.1,s=cv.lasso.1$lambda.1se))
par(op)

mean.Y <- mean(Y.noscale)
mean.X <- apply(X.noscale,2,mean)
sd.X <- apply(X.noscale,2,sd)
# Compute the fitted values with the Lasso estimated model(you can use predict). Plot the observed values for the response variable against the Lasso fitted values

# intercept:
# a) No centering, no scaling
coef(lasso.2,s=cv.lasso.2$lambda.1se)[1]
# b) from the fitted model centering and scaling
mean.Y - sum((mean.X/sd.X) * coef(lasso.1,s=cv.lasso.1$lambda.1se)[-1])

# coefficients:
# a) No centering, no scaling
coef(lasso.2,s=cv.lasso.2$lambda.1se)[-1]
coef(lasso.2,s=cv.lasso.2$lambda.min)[-1]
# b) from the fitted model centering and scaling
coef(lasso.1,s=cv.lasso.1$lambda.1se)[-1] / sd.X
coef(lasso.1,s=cv.lasso.1$lambda.min)[-1] / sd.X

# From the graph, there are 4 coefficient different from 0 when $\log(\lambda) = -1$
coefficients.1.lambdamin = coef(lasso.1,s=cv.lasso.1$lambda.min)
coefficients.1.lambda1se = coef(lasso.1,s=cv.lasso.1$lambda.1se)
coefficients.1.lambdamin[coefficients.1.lambdamin[,1]!=0,]
numberOfCoefficientsNonzero.1.lambdamin = length(coefficients.1.lambdamin[coefficients.1.lambdamin[,1]!=0,]) # 3
numberOfCoefficientsNonzero.1.lambda1se = length(coefficients.1.lambdamin[coefficients.1.lambda1se[,1]!=0,]) # 0
# From the graph, there are 24 coefficient different from 0 when $\log(\lambda) = -1$
coefficients.2.lambdamin = coef(lasso.2,s=cv.lasso.2$lambda.min)
coefficients.2.lambda1se = coef(lasso.2,s=cv.lasso.2$lambda.1se)
coefficients.2.lambdamin[coefficients.2.lambdamin[,1]!=0,]
numberOfCoefficientsNonzero.2.lambdamin = length(coefficients.2.lambdamin[coefficients.2.lambdamin[,1]!=0,])-1 # 3
numberOfCoefficientsNonzero.2.lambda1se = length(coefficients.2.lambdamin[coefficients.2.lambda1se[,1]!=0,])-1 # 0

# prediction in the prediction set to predict the survival time of people still alive
# a) No centering, no scaling
Y.val.noscale <- log.surv.live
Y.val.noscale.hat <- predict(lasso.2,newx=X.val.noscale,type = "response",s=cv.lasso.2$lambda.1se)
plot(Y.val.noscale.hat,Y.val.noscale,asp=1)
abline(a=0,b=1,col=2)
# b) from the fitted model centering and scaling
Y.val.hat <- predict(lasso.1,newx=X.val.scale,type = "response",s=cv.lasso.1$lambda.1se)
Y.val.noscale.hat.1 <- mean.Y + Y.val.hat 
plot(Y.val.noscale.hat.1,Y.val.noscale,asp=1)
abline(a=0,b=1,col=2)
plot(Y.val.noscale.hat,Y.val.noscale.hat.1,asp=1)
abline(a=0,b=1,col=2)

# prediction in the prediction set to predict the survival time of people still alive
# a) No centering, no scaling
Y.val.noscale <- log.surv.death
Y.val.noscale.hat <- predict(lasso.2,newx=X.scale,type = "response",s=cv.lasso.2$lambda.1se)
plot(Y.val.noscale.hat,Y.val.noscale,asp=1)
abline(a=0,b=1,col=2)
# b) from the fitted model centering and scaling
Y.val.hat <- predict(lasso.1,newx=X.val.scale,type = "response",s=cv.lasso.1$lambda.1se)
Y.val.noscale.hat.1 <- mean.Y + Y.val.hat 
plot(Y.val.noscale.hat.1,Y.val.noscale,asp=1)
abline(a=0,b=1,col=2)
plot(Y.val.noscale.hat,Y.val.noscale.hat.1,asp=1)
abline(a=0,b=1,col=2)


