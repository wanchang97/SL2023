setwd("c:/UPC master course/Statistical Learning/2-Regression Estimators Regression ridge Lasso/1-Ridge Regression Assignment")

library(MASS)
data(Boston)
help(Boston)
dim(Boston)
################################################################################
# 1. Data preprocession
################################################################################
train.sample <- Boston[1:200,]
validation.sample <- Boston[201:400,]
test.sample <- Boston[401:dim(Boston)[1],]
Y <- scale(train.sample$medv, center=TRUE, scale=FALSE)
X <- scale(as.matrix(train.sample[,1:13]), center=TRUE, scale=TRUE)
n <- dim(X)[1]
p <- dim(X)[2]

XtX <- t(X)%*%X 
d2 <- eigen(XtX,symmetric = TRUE, only.values = TRUE)$values

(cond.number <- max(d2)/min(d2))
lambda.max <- 1e5
n.lambdas <- 25
lambda.v <- exp(seq(0,log(lambda.max+1),length=n.lambdas))-1
################################################################################
# estimated coefficients path
################################################################################
beta.path <- matrix(0,nrow=n.lambdas, ncol=p)
diag.H.lambda <- matrix(0,nrow=n.lambdas, ncol=n)
for (l in 1:n.lambdas){ 
  lambda <- lambda.v[l]
  H.lambda.aux <- solve(XtX + lambda*diag(1,p)) %*% t(X) 
  beta.path[l,] <-  H.lambda.aux %*% Y
  H.lambda <- X %*% H.lambda.aux 
  diag.H.lambda[l,] <- diag(H.lambda)
} 
plot(c(-1,log(1+lambda.v[n.lambdas])), range(beta.path),type="n",
     xlab="log(1+lambda)",ylab="coefficients")
abline(h=0,lty=2)
for(j in 1:p){
  lines(log(1+lambda.v),beta.path[,j],col=4)
  points(log(1+lambda.v),beta.path[,j],pch=19,cex=.7,col=4)
}
text(0*(1:p), beta.path[1,],names(Boston)[1:p],pos=2)

################################################################################
# effective degrees of freedom
################################################################################
df.v <- numeric(n.lambdas)
for (l in 1:n.lambdas){
  lambda <- lambda.v[l]
  df.v[l] <- sum(d2/(d2+lambda)) 
}
plot(log(1+lambda.v),df.v)
points(0*df.v,df.v,col=2,pch=19)
text(0*df.v,df.v,round(df.v,2),col=2,pch=19,pos=4)


# linear interpolation to obtain the values lambda.vv
# such that the corresponding df are 0,1,...,8 (approx.)
lambda.vv <- approx(x=df.v,y=lambda.v,xout=0:8)$y
lambda.vv[1] <- lambda.v[n.lambdas]
df.vv <- numeric(length(lambda.vv))
for (l in 1:length(lambda.vv)){
  lambda <- lambda.vv[l]
  df.vv[l] <- sum(d2/(d2+lambda)) 
}
print(df.vv)

# another way to compute df's
trace.H.lambda <- apply(diag.H.lambda,1,sum)
print(trace.H.lambda - df.v)

# estimated coefficients path against effective degrees of freedom
plot(c(0,p+1), range(beta.path),type="n",xlab="df(lambda)",ylab="coefficients")
abline(h=0,lty=2)
for(j in 1:p){
  lines(df.v,beta.path[,j],col=4)
  points(df.v,beta.path[,j],pch=19,cex=.7,col=4)
}
text(p+0*(1:p), beta.path[1,],names(prostate)[1:p],pos=4)


################################################################################
# choosing lambda based on the Cross Validation method on the validation sample
################################################################################
lambda_val <- function(X.train,Y.train,X.val,Y.val,lambda.v){
  ##############################################################################
#  Input:
#    X.train, Y.train  -> training samples(already scaled)
#    X.val,Y.val       -> validating samples
#    lambda.v          -> contains all lambda to be chosen for the optimization step
#  Output
#   lambda.val         -> the optimal lambda deducted from the Cross Validation
  ##############################################################################
  
  ##############################################################################
  # Step 1 calculate the beta in the train sample
  ##############################################################################
  n.train <- dim(X.train)[1]
  p <- dim(X.train)[2] # the number of parameters are the same for all sample
  XtX.train <- t(X.train)%*%X.train
  d2.train <- eigen(XtX.train,symmetric= TRUE, only.values= TRUE) $values
  (cond.number <- max(d2.train)/min(d2.train))  # for the OLS regression estimator
  #### Set up the optimization range and resolution
  #lambda.max <- 1e5
  #n.lambdas <- 25
  #lambda.v <- exp(seq(0,log(lambda.max+1),length=n.lambdas))-1
  n.lambdas = length(lambda.v)
  ### Calculate the beta for all lambda and save in beta.path
  beta.path <- matrix(0,nrow=n.lambdas,ncol=p)
  diag.W.lambda <- matrix(0,nrow=n.lambdas, ncol=n.train)
  for (l in 1:n.lambdas){
    lambda <- lambda.v[l]
    W.lambda.aux <- solve(XtX.train + lambda*diag(1,p)) %*% t(X.train) # matrix multiplication
    beta.path[l,] <- W.lambda.aux %*%Y.train
    W.lambda <- X.train %*% W.lambda.aux
    diag.W.lambda[l,] <- diag(W.lambda)}
  # Plot beta to check
  if (1) {
    plot(c(-1,log(1+lambda.v[n.lambdas])),range(beta.path),type="n",xlab = "log(1+lambda)",ylab = "coefficients")
    abline(h=0,lty=2)
    for(j in 1:p){
      lines(log(1+lambda.v),beta.path[,j],col=4)
      points(log(1+lambda.v),beta.path[,j],pch=19,cex=.7,col=4)
    }  
    #text(0*(1:p),beta.path[1,],names(prostate)[1:p],pos=2)  
    
    title("beta v.s lambda") 
  }
  ##############################################################################
  # Step 2 calculate the MSPE in the validation sample
  ##############################################################################
  n.val <- dim(X.val)[1]
  XtX.val <- t(X.val)%*%X.val
  d2.val <- eigen(XtX.val,symmetric= TRUE, only.values= TRUE) $values
  MSPE.val <- numeric(n.lambdas)
  for (l in 1:n.lambdas){
    for (i in 1:n.val){
      Xi <- X.val[i,];Yi <- Y.val[i]
      hat.Y.evali <- Xi %*% beta.path[l,]
      MSPE.val[l] <- MSPE.val[l] + (hat.Y.evali-Yi)^2
    }
    MSPE.val[l]<- MSPE.val[l]/n.val}
  ##############################################################################
  # Step 3 choose the optimal lambda which minimizing the MSPE
  ##############################################################################
  lambda.opt <- lambda.v[which.min(MSPE.val)]
  return (lambda.opt)
}

#### Call the function lambda.val
library(MASS)
data(Boston)
help(Boston)
dim(Boston)
train.sample <- Boston[1:200,]
validation.sample <- Boston[201:400,]
test.sample <- Boston[401:dim(Boston)[1],]
Y.train <- scale(train.sample$medv, center=TRUE, scale=FALSE)
X.train <- scale(as.matrix(train.sample[,1:13]), center=TRUE, scale=TRUE)
Y.val <- scale(validation.sample$medv, center=TRUE, scale=FALSE)
X.val <- scale(as.matrix(validation.sample[,1:13]), center=TRUE, scale=TRUE)
lambda.max <- 1e5
n.lambdas <- 25
lambda.v <- exp(seq(0,log(lambda.max+1),length=n.lambdas))-1
lambda.CV.val <- lambda_val(X.train,Y.train,X.val,Y.val,lambda.v)
### Fit the model on the test sample
Y_Ridge_Regression <- function(X.test,Y.test,lambda.opt){
  n.test <- dim(X.test)[1]
  p <- dim(X.test)[2] 
  Y_estimator <- X.test %*% solve(t(X.test) %*% X.test+lambda.opt*diag(p)) %*% t(X.test) %*% Y.test
  return(Y_estimator) 
}
Y.test <- scale(test.sample$medv, center=TRUE, scale=FALSE)
X.test <- scale(as.matrix(test.sample[,1:13]), center=TRUE, scale=TRUE)
Y.test.estimator <- Y_Ridge_Regression(X.test,Y.test,lambda.CV.val)
W = X.train%*% solve(t(X.train%*%X.train+lambda.CV.val*diag(p)))%*%t(X.train)%*%Y.train

###############################################################################
# function to get the effective degree of freedom from lambda
###############################################################################
df_ridge <- function(lambda,X.train){
  XtX <- t(X.train)%*%X.train
  d2 <- eigen(XtX,symmetric= TRUE, only.values= TRUE) $values
  df <- sum(d2/(d2+lambda))
  return (df)
}
###############################################################################
# function leave one out cross validation
###############################################################################
MSPE_CV_one <- function(X,Y,lambda.v){
  n.lambdas <- length(lambda.v)
  n <- dim(X)[1]
  p <- dim(X)[2]
  MSPE.CV.one <- numeric(n.lambdas) # n.lambdas 0s initialize MSPE.CV
  for (l in 1:n.lambdas){
    for (i in 1:n){
      m.Y.i <- 0 #? 
      X.i <- X[-i,]; Y.i <- Y[-i]-m.Y.i# remove the xi,yi from data
      Xi <- X[i,];Yi <- Y[i]
      beta.i <- solve(t(X.i)%*%X.i+lambda.v[l]*diag(1,p)) %*% t(X.i) %*% Y.i
      hat.Yi <- Xi %*% beta.i
      MSPE.CV.one[l] <-MSPE.CV.one[l] + (hat.Yi-Yi)^2
    }
    MSPE.CV.one[l]<- MSPE.CV.one[l]/n  
  }
  return (MSPE.CV.one)
}

###############################################################################
# function leave one out cross validation a more computational effective way
###############################################################################
MSPE_CV_one_efficient <- function(X,Y,lambda.v){
  n.lambdas <- length(lambda.v)
  n <- dim(X)[1]
  p <- dim(X)[2]
  MSPE.CV.one <- numeric(n.lambdas)
  ## calculate beta.path and diag.W.lambda
  beta.path <- matrix(0,nrow=n.lambdas,ncol=p)
  diag.W.lambda <- matrix(0,nrow=n.lambdas, ncol=n)
  for (l in 1:n.lambdas){
    lambda <- lambda.v[l]
    W.lambda.aux <- solve(XtX + lambda*diag(1,p)) %*% t(X) # matrix multiplication
    beta.path[l,] <- W.lambda.aux %*%Y
    W.lambda <- X %*% W.lambda.aux
    diag.W.lambda[l,] <- diag(W.lambda)}
  
  for (l in 1:n.lambdas){
    lambda <- lambda.v[l]
    hat.Y <- X %*% beta.path[l,]
    #nu <- sum(diag.W.lambda[l,])
    for (i in 1:n){
      MSPE.CV.one[l] <-MSPE.CV.one[l] + ((Y[i]-hat.Y[i])/(1-diag.W.lambda[l,i]))^2
    }
    MSPE.CV.one[l] <- MSPE.CV.one[l]/n
  }
}



###############################################################################
# function k-fold cross validation
###############################################################################
set.seed(777)
MSPE_CV_k <- function(X,Y,lambda.v,k){
  n.lambdas <- length(lambda.v)
  n <- dim(X)[1]
  p <- dim(X)[2]
  MSPE.CV.k <- numeric(n.lambdas) # n.lambdas 0s initialize MSPE.CV
  for (l in 1:n.lambdas){
    for (i in 1:k){
      # split data into training and validation
      picked <- sample(seq_len(n),size=k)
      # m.Y.k <- 0 #? 
      X.k <- X[-picked,]; Y.k <- Y[-picked]# validation sample
      Xk <- X[picked,];Yk <- Y[picked] # training sample
      beta.k <- solve(t(X.k)%*%X.k+lambda.v[l]*diag(1,p)) %*% t(X.k) %*% Y.k
      hat.Yk <- Xk %*% beta.k
      MSPE.CV.k[l] <-MSPE.CV.k[l] + (hat.Yk-Yk)^2
    }
    MSPE.CV.k[l]<- MSPE.CV.k[l]/k  
  }
  return (MSPE.CV.k)
}

###############################################################################
# calculated the coefficients
beta_path <- function(X,Y,lambda.v){
  n.lambdas <- length(lambda.v)
  n <- dim(X)[1]
  p <- dim(X)[2]
  beta.path <- matrix(0,nrow=n.lambdas,ncol=p)
  diag.W.lambda <- matrix(0,nrow=n.lambdas, ncol=n)
  for (l in 1:n.lambdas){
    lambda <- lambda.v[l]
    W.lambda.aux <- solve(XtX + lambda*diag(1,p)) %*% t(X) # matrix multiplication
    beta.path[l,] <- W.lambda.aux %*%Y
    W.lambda <- X %*% W.lambda.aux
    diag.W.lambda[l,] <- diag(W.lambda)}
  return(list(beta.path,diag.W.lambda))
}
###############################################################################


###############################################################################
# function generalized cross validation
###############################################################################
MSPE_GCV <- function(X,Y,lambda.v){
  n.lambdas <- length(lambda.v)
  n <- dim(X)[1]
  p <- dim(X)[2]
  MSPE.GCV <- numeric(n.lambdas)
  ## calculate beta.path and diag.W.lambda
  beta.path <- matrix(0,nrow=n.lambdas,ncol=p)
  diag.W.lambda <- matrix(0,nrow=n.lambdas, ncol=n)
  for (l in 1:n.lambdas){
    lambda <- lambda.v[l]
    W.lambda.aux <- solve(XtX + lambda*diag(1,p)) %*% t(X) # matrix multiplication
    beta.path[l,] <- W.lambda.aux %*%Y
    W.lambda <- X %*% W.lambda.aux
    diag.W.lambda[l,] <- diag(W.lambda)}

  for (l in 1:n.lambdas){
    lambda <- lambda.v[l]
    hat.Y <- X %*% beta.path[l,]
    nu <- sum(diag.W.lambda[l,])
    MSPE.GCV[l] <- sum( ((Y-hat.Y)/(1-nu/n))^2 )/n}
  return (MSPE.GCV)}