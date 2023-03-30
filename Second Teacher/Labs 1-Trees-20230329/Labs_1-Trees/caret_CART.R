## ----setup, include=FALSE------------------------------------------------
options(width=100) 
if(!require("knitr")) install.packages("knitr")
library("knitr")
#getOption("width")
knitr::opts_chunk$set(comment=NA,echo = TRUE, cache=TRUE)

## ----s1, include=TRUE, message=FALSE, warning=FALSE----------------------
if(!require("caret")) install.packages("caret")
if(!require("mlbench")) install.packages("mlbench")
library("caret")
library("mlbench")

## ----s1.1----------------------------------------------------------------
data(Sonar)
names(Sonar)
set.seed(1234) # Control of data generation
inTrain <- createDataPartition(y=Sonar$Class, p=.75, list=FALSE)
str(inTrain)
training <- Sonar[inTrain,]
testing <- Sonar[-inTrain,]
nrow(training)

## ----s2------------------------------------------------------------------
CART1Model <- train (Class ~ ., 
                   data=training, 
                   method="rpart1SE",
                   preProc=c("center","scale"))
CART1Model

## ------------------------------------------------------------------------
ctrl <- trainControl(method = "repeatedcv", repeats=3)
CART1Model3x10cv <- train (Class ~ ., 
                         data=training, 
                         method="rpart1SE",
                         trControl=ctrl,
                         preProc=c("center","scale"))

CART1Model3x10cv

## ------------------------------------------------------------------------
ctrl <- trainControl(method = "repeatedcv", repeats=3, classProbs=TRUE,
summaryFunction=twoClassSummary)
CART1Model3x10cv <- train (Class ~ ., 
                         data=training, 
                         method="rpart1SE", 
                         trControl=ctrl, 
                         metric="ROC", 
                         preProc=c("center","scale"))

CART1Model3x10cv

## ------------------------------------------------------------------------
CART2Fit3x10cv <- train (Class ~ ., 
                       data=training, 
                       method="rpart", 
                       trControl=ctrl, 
                       metric="ROC", 
                       preProc=c("center","scale"))
CART2Fit3x10cv
plot(CART2Fit3x10cv)

CART2Fit3x10cv <- train (Class ~ ., 
                       data=training, 
                       method="rpart", 
                       trControl=ctrl, 
                       metric="ROC",  
                       tuneLength=10,
                       preProc=c("center","scale"))
CART2Fit3x10cv
plot(CART2Fit3x10cv)

## ---- eval= TRUE---------------------------------------------------------
CART2Probs <- predict(CART2Fit3x10cv, newdata = testing, type = "prob")
CART2Classes <- predict(CART2Fit3x10cv, newdata = testing, type = "raw")
confusionMatrix(data=CART2Classes,testing$Class)

## ---- eval= TRUE---------------------------------------------------------
resamps=resamples(list(CART2=CART2Fit3x10cv,CART1=CART1Model3x10cv))
summary(resamps)
xyplot(resamps,what="BlandAltman")
diffs<-diff(resamps)
summary(diffs)

