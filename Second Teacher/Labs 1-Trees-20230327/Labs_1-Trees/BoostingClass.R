
# Classification Boosting
##################################################

if(!require("ISLR")) install.packages("ISLR")
library("ISLR")

if(!require("gbm")) install.packages("gbm")
library("gbm")

if(!require("pROC")) install.packages("pROC")
library("pROC")


data("Carseats")


Carseats$High=ifelse(Carseats$Sales<=8,0,1)

Carseats <- Carseats[,-1]


set.seed(2)
pt <- 1/2
train <- sample(1:nrow(Carseats),pt*nrow(Carseats))
Carseats.test=Carseats[-train,]
High.test=Carseats$High[-train]



set.seed(1)


boost.Carseats=gbm(High~.,data=Carseats[train,],distribution="adaboost",
                   n.trees=5000,interaction.depth=4)

summary(boost.Carseats)

par(mfrow=c(1,2))
plot(boost.Carseats,i="Price")
plot(boost.Carseats,i="ShelveLoc")
par(mfrow=c(1,1))

pretty.gbm.tree(boost.Carseats,i.tree = 1)


predict.Carseats=predict(boost.Carseats,newdata=Carseats[-train,],
                         n.trees=5000, type = "response")


list.pred <- data.frame("Actual" = High.test, 
           "PredictedProbability" = predict.Carseats)

head(list.pred)



auc <- roc(list.pred[,1], list.pred[,2])
plot.roc(auc)
print(auc$auc)

######

boost.Carseats=gbm(High~.,data=Carseats[train,],distribution="adaboost",
                   n.trees=2000,interaction.depth=4, cv.folds = 3)

summary(boost.Carseats)

gbm.perf(boost.Carseats)



# Regression Boosting

library(MASS)

set.seed(1234)
pt <- 1/2
train <- sample(1:nrow(Boston),pt*nrow(Boston))
boston.test <- Boston$medv[-train]


set.seed(1)
boost.boston=gbm(medv~.,data=Boston[train,],distribution="gaussian",
                 n.trees=5000,interaction.depth=4)

summary(boost.boston)
par(mfrow=c(1,2))
plot(boost.boston,i="rm")
plot(boost.boston,i="lstat")
par(mfrow=c(1,1))

yhat.boost=predict(boost.boston,newdata=Boston[-train,],n.trees=5000)
mean((yhat.boost-boston.test)^2)

# new model
boost.boston=gbm(medv~.,data=Boston[train,],distribution="gaussian",
                 n.trees=5000,interaction.depth=4,shrinkage=0.2,verbose=F)
yhat.boost=predict(boost.boston,newdata=Boston[-train,],n.trees=5000)
mean((yhat.boost-boston.test)^2)

