# Lasso for Boston Housing Data
# The Boston House price dataset contians housing values in 56 surburbs of Boston corresponding to year 1978
# This is the list of the available variables:
#1. CRIM: per capita crime rate by town
#2  ZN: proportion of the residential land zone for lots over 25,000 sq.ft.
#3. INDUS: proportion of non-retail business ycres per town
#4 CHAS: Charles River dummy variable(=1 if tract bounds river; 0 otherwise)
#5 NOX: Nitric Oxides Concentration(parts per 10 million)
#6 RM: average number of rooms per dwelling
#7 AGE: proportion of owner-occupied units built prior to 1940
#8 DIS: weighted distances to five Boston employment centres
#9 RAD index of accessibilitz to radial highways
#10: TAXfull-value property-tax rate per $10,000
#11 PTRATIO. Pupil-teacher ratio by town
#12 B: 1000(Bk-0.63)^2where Bk is the proportion of blacks by town
#13 LSTAT. % LOWER STATUS OF THE POPULATION
#14 MEDV: Median value of owner-occupied dhomes in $1000´s
#15 UTM coordinats of the geographical centers of each neighbourhood

# For the Boston House-price corrreced dataset using Lasso estimation in glmnet to fit the regression model where the response is MEDV and the explanatory variables are the remaining 14 variables in the previous list
# Try to provide an inerpretation to the estimated model
library(Hmisc)
library(glmnet)
load("boston.Rdata")
#select only the variable of interest to the model, in the order desired:
boston = boston.c[c("CRIM","ZN","INDUS","CHAS","NOX","RM","AGE","DIS","RAD","TAX","PTRATIO","B","LSTAT","MEDV")]
labs = list(
  'per capita crime rate by town',
  'proportion of the residential land zone for lots over 25,000 sq.ft.',
  'proportion of non-retail business ycres per town',
  'Charles River dummy variable(=1 if tract bounds river; 0 otherwise)',
  'Nitric Oxides Concentration(parts per 10 million)',
  'average number of rooms per dwelling',
  'proportion of owner-occupied units built prior to 1940',
  'weighted distances to five Boston employment centres',
  'index of accessibilitz to radial highways',
  'full-value property-tax rate per $10,000',
  'Pupil-teacher ratio by town',
  '1000(Bk-0.63)^2where Bk is the proportion of blacks by town',
  '% LOWER STATUS OF THE POPULATION',
  'Median value of owner-occupied dhomes in $1000´s'
)
label(boston,which=NULL) = labs
names(boston.c) # names of the variables
describe(boston,1)# name,label,summary statistics
nobs = dim(boston)[1] # 506 observations
# Scaled data
Y = scale(boston$MEDV,center=TRUE,scale=TRUE)
# since CHAS is a binary variable, we need to use as.vector to transform the numeric to category
boston$CHAS = as.vector(boston$CHAS)
X1 = scale(boston[,c(-4,-14)],center=TRUE,scale=TRUE)
X = cbind(X1,as.numeric(boston$CHAS)-1)
p = dim(X)[2]
n = dim(X)[1]
# Nonscaled data
Y.noscale = scale(boston$MEDV,center = FALSE,scale=FALSE)
X1.noscale = scale(boston[,c(-4,-14)],center=FALSE,scale=FALSE)
X.noscale = cbind(X1.noscale,as.numeric(boston$CHAS)-1)
#X.noscale = scale(boston[,c(-14)],center=FALSE,scale=FALSE)#as.matrix(boston[,c(-14)])

lasso.1 <- glmnet(X,Y,standandize = FALSE, intercept = FALSE)# Lasso.1 corresponds to the scaled data
cv.lasso.1 <- cv.glmnet(X,Y,standandize= FALSE, intercept = FALSE) # by default nfolds = 10

lasso.2 <- glmnet(X.noscale,Y.noscale, standardize=TRUE, intercept=TRUE)
cv.lasso.2 <- cv.glmnet(X.noscale,Y.noscale, standardize=TRUE, intercept=TRUE)
print(coef(lasso.2,s=cv.lasso.2$lambda.1se))
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


# label gives which column (variable) corresponds to each line

# Interpretation:
# The lambda controls the shrinkage content, increases from 0 to +Infty
# When lambda is almost 0, log(lambda) is -Infty, all coefficients are nonzero.
# When lambda is lambda.min, where the MSPE is minimized, there are 7 coefficients being nonzero, the 7th one being zero
# intercept:
# coefficients:
# a) No centering, no scaling
coef(lasso.2,s=cv.lasso.2$lambda.min)[-1]# 11 coefficients nonzero
# b) from the fitted model centering and scaling
coef(lasso.1,s=cv.lasso.2$lambda.min)[-1] / sd.X
# When lambda is lambda.1se, where the MSPE is minimized, there are 7 coefficients being nonzero, the 7th one being zero
# intercept:
# coefficients:
# a) No centering, no scaling
coef(lasso.2,s=cv.lasso.2$lambda.1se)[-1] # 8 coefficients nonzero
# b) from the fitted model centering and scaling
coef(lasso.1,s=cv.lasso.2$lambda.1se)[-1] / sd.X

# Exercise glmnet to fit the model using the ridge regression. Compare the 10-fold cross validation results from function cv.glmnet with your own functions
ridge.1 <- glmnet(X,Y,alpha = 0,standandize = FALSE, intercept = FALSE)# Lasso.1 corresponds to the scaled data
cv.ridge.1 <- cv.glmnet(X,Y,alpha = 0,standandize= FALSE, intercept = FALSE) # by default nfolds = 10

ridge.2 <- glmnet(X.noscale,Y.noscale,alpha = 0, standardize=TRUE, intercept=TRUE)
cv.ridge.2 <- cv.glmnet(X.noscale,Y.noscale,alpha = 0, standardize=TRUE, intercept=TRUE)

# intercept:
# a) No centering, no scaling
coef(ridge.2,s=cv.ridge.2$lambda.1se)[1]
# b) from the fitted model centering and scaling
mean.Y - sum((mean.X/sd.X) * coef(ridge.1,s=cv.ridge.2$lambda.1se)[-1])

# coefficients:
# a) No centering, no scaling
coef(ridge.2,s=cv.ridge.2$lambda.1se)[-1]
# b) from the fitted model centering and scaling
coef(ridge.1,s=cv.ridge.2$lambda.1se)[-1] / sd.X


# Our own ridge regression function
pmse_kcv= function(X,Y,lambda.v,K){
  n.lambdas = length(lambda.v)
  p = dim(X)[2]
  # Split the data into k-folds
  set.seed(123)
  folds = sample(rep(1:k,length.out = dim(X)[1]))
  pmse.lambda = matrix(0,nrow=n.lambdas,ncol=1)
  for (k in 1:K){
    # Get the training and validation sets for this fold
    x.train = X[folds != k,]
    y.train = Y[folds != k]
    x.val = X[folds == k,]
    y.val = Y[folds == k]
    nv = dim(x.val)[1]
    # now fi the model for each lambda
    for (l in  1:n.lambdas){
      beta.lambda = t(solve(t(x.train)%*%x.train + lambda.v[l]*diag(1,p)))%*%t(x.train)%*%y.train
      # Now we can calculate the pmse for lambda
      error.lambda = y.val -x.val %*% beta.lambda
      pmse.lambda[l] = pmse.lambda[l] + t(error.lambda)%*%error.lambda
    }
  }
  # Now we divide by n
  pmse.cv = pmse.lambda/dim(X)[1]
  return(pmse.cv)
}
k = 10
lambda.max = 1e5
n.lambdas = 100
lambda.v = exp(seq(0,log(lambda.max+1),length=n.lambdas))-1
pmse.cv = pmse_kcv(X,Y,lambda.v,k)
posmin = which.min(pmse.cv)
lambda.cv = lambda.v[posmin]# 5.428
cv.ridge.1$lambda.1se
cv.ridge.2$lambda.1se # 4.7815 more close to lambda.cv
cv.ridge.1$lambda.min
cv.ridge.2$lambda.min

