# Function computing the degrees of freedom for the vector of lambdas.
df.fun = function(X, lambda.v){
  n.lambdas = length(lambda.v)
  XtX <- t(X)%*%X # X'X
  XtX.vaps <- eigen(XtX,symmetric = TRUE, only.values = TRUE)$values
  df.lambda = matrix(0,nrow=n.lambdas, ncol=1)
  for (l in 1:n.lambdas){ 
    lambda = lambda.v[l]
    df.lambda[l,] = sum(XtX.vaps/(XtX.vaps+lambda))
  } 
  return(df.lambda)
}


prostate <- read.table("prostate_data.txt", header=TRUE, row.names = 1)

# 1. RIDGE REGRESSION PENALIZATION PARAMETER

# Inputs
training.sample = which(prostate$train==TRUE)
validation.sample = which(prostate$train==FALSE)
lambda.max = 1e5
n.lambdas = 100
lambda.v = exp(seq(0,log(lambda.max+1),length=n.lambdas))-1  # candidate lambdas

# Define training and validation sets
Y = scale(prostate$lpsa[training.sample], center=TRUE, scale=FALSE)
X = scale( as.matrix(prostate[training.sample,1:8]), center=TRUE, scale=TRUE)
Yv = scale(prostate$lpsa[validation.sample], center=TRUE, scale=FALSE)
Xv = scale( as.matrix(prostate[validation.sample,1:8]), center=TRUE, scale=TRUE)

pmse_val <- function(X, Y, Xv, Yv, lambda.v) {
  nv <- dim(Xv)[1]
  p <- dim(X)[2]
  
  # Calculating the coefficients beta for each lambda
  pmse.lambda = matrix(0,nrow=length(lambda.v), ncol=1)
  for (l in 1:length(lambda.v)){ 
    H.lambda.aux = t(solve(t(X)%*%X + lambda.v[l]*diag(1,p))) %*% t(X) # (X'X + lId)^{-1} X'
    beta.lambda = H.lambda.aux %*% Y
    
    # Now we can calculate the PMSE for lambda
    error.lambda = Yv - Xv %*% beta.lambda
    pmse.lambda[l,] = (1/nv) * t(error.lambda) %*% error.lambda
  } 
  return(pmse.lambda)
}

# Function to plot PMSE vs DF.
plot.pmse.df = function(df.lambda, pmse.lambda, title, bool_log){
  if (bool_log == TRUE){
    df.lambda = log(1+lambda.v)-1
  }
  plot(range(df.lambda), range(pmse.lambda),type="n",xlab="df(lambda)",ylab="PMSE(lambda)",
       main=title)
  lines(df.lambda, pmse.lambda, col=4)
  points(df.lambda, pmse.lambda, pch=19,cex=1,col=4)
  posmin = which.min(pmse.lambda)
  points(df.lambda[posmin], pmse.lambda[posmin], pch=19,cex=1,col=2)
  abline(v=df.lambda[posmin],col=4,lty=3,lwd=2)
  
  lambda.text = paste("l = ", round(lambda.v[posmin],2),sep="")
  legend('topright', legend=c(lambda.text), 
         pch = c(16), col = c(2))
}

df.lambda = df.fun(X, lambda.v)
pmse.lambda = pmse_val(X, Y, Xv, Yv, lambda.v)
plot.pmse.df(df.lambda, pmse.lambda, "PMSE(lambda) in validation set", FALSE)

# 2.
pmse_kCV = function(X, Y, lambda.v, k){
  n.lambdas = length(lambda.v)
  p = dim(X)[2]
  
  # Split the data into k-folds
  set.seed(123) # for reproducibility
  folds = sample(rep(1:k, length.out = dim(X)[1])) # assign a set to each obs.
  pmse.lambda = matrix(0,nrow=n.lambdas, ncol=1)
  for (i in 1:k) {
    # Get the training and validation sets for this fold
    x_train <- X[folds != i,]
    y_train <- Y[folds != i]
    x_val <- X[folds == i,]
    y_val <- Y[folds == i]
    nv = dim(x_val)[1]

    # Now we fit the model for each lambda
    for (l in 1:n.lambdas){ 
      beta.lambda = t(solve(t(x_train)%*%x_train + lambda.v[l]*diag(1,p))) %*% t(x_train) %*% y_train
      
      # Now we can calculate the PMSE for lambda
      error.lambda = y_val - x_val %*% beta.lambda
      pmse.lambda[l] = pmse.lambda[l] + t(error.lambda)%*%error.lambda
    } 
  }
  # Now we divide by n
  pmse.cv = pmse.lambda/dim(X)[1]
  return(pmse.cv)
}

k = 10
pmse.cv = pmse_kCV(X, Y, lambda.v, k)

# Plot
title = paste('Evolution of PMSE as a function of df [', k, sep="")
title = paste(title, '-fold CV]', sep='')
plot.pmse.df(df.lambda, pmse.cv, title, FALSE)



# 3.
# Let's train the regression to find the optimal lambda.
# (1) Validation 
pmse.val = pmse_val(X, Y, Xv, Yv, lambda.v)
posmin = which.min(pmse.val)
lambda.val = lambda.v[posmin]
df.lambda.val = df.lambda[posmin]
# (2) 5-fold cross-validation
pmse.5 = pmse_kCV(X, Y, lambda.v, 5)
posmin = which.min(pmse.5)
lambda.5 = lambda.v[posmin]
df.lambda.5 = df.lambda[posmin]
# (2) 10-fold cross-validation
pmse.10 = pmse_kCV(X, Y, lambda.v, 10)
posmin = which.min(pmse.10)
lambda.10 = lambda.v[posmin]
df.lambda.10 = df.lambda[posmin]

# We can compare with leave-one-out CV easily with our function taking k=n
pmse.loo = pmse_kCV(X, Y, lambda.v, dim(X)[1])
posmin = which.min(pmse.loo)
lambda.loo = lambda.v[posmin]
df.lambda.loo = df.lambda[posmin]

# And compute generalized-cross-validation as follows:
n = dim(X)[1]
p = dim(X)[2]
pmse.gcv = numeric(n.lambdas)
for (l in 1:n.lambdas){
  lambda = lambda.v[l]
  H.lambda = X %*% t(solve(t(X)%*%X + lambda*diag(1,p))) %*% t(X) # X (X'X + lId)^{-1} X'
  nu = sum(diag(H.lambda))
  pmse.gcv[l] = sum(((Y-H.lambda%*% Y)/(1-nu/n))^2 )/n
}
posmin = which.min(pmse.gcv)
lambda.gcv = lambda.v[posmin]
df.lambda.gcv = df.lambda[posmin]


# NOW WE CAN COMPARE THE RESULTS:
torange = c(pmse.val, pmse.5, pmse.10, pmse.loo, pmse.gcv)
plot(range(df.lambda), c(max(torange),min(torange)),type="n",xlab="df(lambda)",ylab="PMSE(lambda)",
     main='lambda choice for different methods')

lines(df.lambda, pmse.val, col=4)
points(df.lambda, pmse.val, pch=19,cex=.7,col=4)
abline(v=df.lambda.val,col=4,lty=3,lwd=2)

lines(df.lambda, pmse.5, col=1)
points(df.lambda, pmse.5, pch=19,cex=.7,col=1)
abline(v=df.lambda.5,col=1,lty=3,lwd=2)

lines(df.lambda, pmse.10, col=3)
points(df.lambda, pmse.10, pch=19,cex=.7,col=3)
abline(v=df.lambda.10,col=3,lty=3,lwd=2)

lines(df.lambda, pmse.loo, col=6)
points(df.lambda, pmse.loo, pch=19,cex=.7,col=6)
abline(v=df.lambda.loo,col=6,lty=3,lwd=2)

lines(df.lambda, pmse.gcv, col=7)
points(df.lambda, pmse.gcv, pch=19,cex=.7,col=7)
abline(v=df.lambda.gcv,col=7,lty=3,lwd=2)

legend.text = c(paste("validation set => l=", round(lambda.val,2)),
                paste("5-fold CV => l=", round(lambda.5,2)),
                paste("10-fold CV => l=", round(lambda.10,2)),
                paste("LOO-CV => l=", round(lambda.loo,2)),
                paste("GCV => l=", round(lambda.gcv,2)))

legend('topright', legend=legend.text, 
       pch = c(16,16,16,16,16), col = c(4,1,3,6,7))


# NOW IMPLEMENTATION IN THE BOSTON HOUSE DATAFRAME
library(Hmisc)

load("boston.Rdata")

# Select only the variables of interest to the model, in the order desired:
boston = boston.c[c("CRIM","ZN","INDUS","CHAS","NOX","RM","AGE","DIS",
                  "RAD","TAX","PTRATIO","B","LSTAT","MEDV")]

# Add description to favor interpretation of the model:
labs = list(
  'per capita crime rate by town',
  'proportion of residential land zoned for lots over 25,000 sq.ft.',
  'proportion of non-retail business acres per town',
  'Charles River dummy variable (= 1 if tract bounds river; 0 otherwise)',
  'nitric oxides concentration (parts per 10 million)',
  'average number of rooms per dwelling',
  'proportion of owner-occupied units built prior to 1940',
  'weighted distances to five Boston employment centres',
  'index of accessibility to radial highways',
  'full-value property-tax rate per $10,000',
  'pupil-teacher ratio by town',
  '1000(Bk - 0.63)^2 where Bk is the proportion of blacks by town',
  '% lower status of the population',
  'Median value of owner-occupied homes in $1000’s'
)

label(boston, which=NULL) = labs

# Add a description to favor interpretation of the model:
names(boston.c)  # names of the variables
describe(boston,1) # name, label and summary statistics
nobs = dim(boston)[1] # 506 observations

Y = scale(boston$MEDV, center=TRUE, scale=FALSE)

# Since CHAS is a binary variable, we don't have to include it in the scale. 
X1 = scale(boston[,c(-4,-14)], center=TRUE, scale=TRUE)
X = cbind(X1, as.numeric(boston$CHAS)-1) #Include it again 
colnames(X) = c(colnames(X1), "CHAS")
p = dim(X)[2] # 13 variables

# Lambdas to try:
lambda.max <- 1e5
n.lambdas <- 100
lambda.v <- exp(seq(0,log(lambda.max+1),length=n.lambdas))-1  # candidate lambdas

# Now calculate degrees of freedom:
XtX <- t(X)%*%X # X'X
XtX.vaps <- eigen(XtX,symmetric = TRUE, only.values = TRUE)$values
df.lambda = matrix(0,nrow=n.lambdas, ncol=1)
for (l in 1:n.lambdas){ 
  lambda = lambda.v[l]
  df.lambda[l,] = sum(XtX.vaps/(XtX.vaps+lambda))
} 

# (0) Validation set

# Split the dataset into validation set and train set:
ntrain <- round(0.7 * nobs)
indtrain <- sample(1:nobs, ntrain, replace = FALSE)
Xtrain <- X[indtrain, ]
Ytrain <- Y[indtrain,]
Xval <- X[-indtrain,]
Yval <- Y[-indtrain,]

pmse_val <- pmse_val(Xtrain,Ytrain, Xval, Yval, lambda.v)
posmin = which.min(pmse_val)
lambda.valid = lambda.v[posmin]
df.lambda.valid = df.lambda[posmin]

# (1) LOO cross-validation
pmse.loo = pmse_kCV(X, Y, lambda.v, nobs)
posmin = which.min(pmse.loo)
lambda.loo = lambda.v[posmin]
df.lambda.loo = df.lambda[posmin]

# (2) 5-fold cross-validation
pmse.5 = pmse_kCV(X, Y, lambda.v, 5)
posmin = which.min(pmse.5)
lambda.5 = lambda.v[posmin]
df.lambda.5 = df.lambda[posmin]

# (3) 10-fold cross-validation
pmse.10 = pmse_kCV(X, Y, lambda.v, 10)
posmin = which.min(pmse.10)
lambda.10 = lambda.v[posmin]
df.lambda.10 = df.lambda[posmin]

# We plot to compare the results:
torange = c(pmse.5, pmse.10, pmse.loo)
plot(range(df.lambda), c(max(torange),min(torange)),type="n",xlab="df(lambda)",ylab="PMSE(lambda)",
     main='lambda choice for different methods')


lines(df.lambda, pmse_val, col=4)
points(df.lambda, pmse_val, pch=19,cex=.7,col=4)
abline(v=df.lambda.valid ,col=4,lty=3,lwd=2)


lines(df.lambda, pmse.5, col=1)
points(df.lambda, pmse.5, pch=19,cex=.7,col=1)
abline(v=df.lambda.5,col=1,lty=3,lwd=2)

lines(df.lambda, pmse.10, col=3)
points(df.lambda, pmse.10, pch=19,cex=.7,col=3)
abline(v=df.lambda.10,col=3,lty=3,lwd=2)

lines(df.lambda, pmse.loo, col=6)
points(df.lambda, pmse.loo, pch=19,cex=.7,col=6)
abline(v=df.lambda.loo,col=6,lty=3,lwd=2)

legend.text = c(paste("validation set => l=", round(lambda.valid,2)),
                paste("5-fold CV => l=", round(lambda.5,2)),
                paste("10-fold CV => l=", round(lambda.10,2)),
                paste("LOO-CV => l=", round(lambda.loo,2))
                )


legend('top', legend=legend.text, 
       pch = c(16,16,16,16), col = c(4,1,3,6))


# The 5- and 10-fold CV provide the same minimum lambda, at around 5.8
# corresponding to 12.5 degrees of freedom.
# We will use lambda so that we have 12 degrees of freedom:
lambda.12 <- approx(x=df.lambda, y=lambda.v, xout=12)$y # 13.47
beta.12 = t(solve(t(X)%*%X + lambda*diag(1,p))) %*% t(X) %*% Y
beta.12.sorted = sort(beta.12, decreasing = TRUE)
print(beta.12)




  






