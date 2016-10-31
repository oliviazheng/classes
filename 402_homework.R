rm(list = ls(all = TRUE))
setwd('C:\\Users\\ohz\\Dropbox\\402 Regression')
diabetes<-read.csv('diabetes.csv',header = TRUE)
library(boot)
library(cubature)
library(np)
library(MASS)
require(mgcv)
require(graphics)
patient<-diabetes$patient
age<-diabetes$age
base.deficit<-diabetes$base.deficit
c.peptide<-diabetes$c.peptide

###########
#Problem 1#
###########
addfit<-gam(c.peptide~s(age)+s(base.deficit),data=diabetes)
plot(addfit,scale=0,se=2,shade=TRUE,resid=TRUE,pages=1)

x1.points <- seq(-29,0,length.out=100)
x2.points <- seq(.9,16,length.out=100)
x12grid <- expand.grid(base.deficit=x1.points,age=x2.points)
y.out <- matrix(0,100,100)
y.out <- predict(addfit,newdata=x12grid)
require(lattice)
wireframe(y.out~x12grid$base.deficit*x12grid$age,scales=list(arrows=FALSE),
  xlab="base.deficit",ylab="age",zlab="c.peptide")

###########
#Problem 2#
###########
df5<-data.frame(age=rep(5,200),base.deficit=seq(from=-29,to=0,length=1000))
df10<-data.frame(age=rep(10,200),base.deficit=seq(from=-29,to=0,length=1000))
df12<-data.frame(age=rep(12,200),base.deficit=seq(from=-29,to=0,length=1000))

plot(seq(from=-29,to=0,length=1000),predict(addfit,df5),xlab="Base.Deficit",ylab="Predicted c.peptide",pch=".")
points(seq(from=-29,to=0,length=1000),predict(addfit,df10),pch=".",col="blue")
points(seq(from=-29,to=0,length=1000),predict(addfit,df12),pch=".",col="red")

###########
#Problem 3#
###########
np<-npreg(c.peptide~base.deficit+age,data=diabetes)
plot(np,phi=50)

y.out2 <- matrix(0,100,100)
y.out2 <- predict(np,newdata=x12grid)
wireframe(y.out2~x12grid$base.deficit*x12grid$age,scales=list(arrows=FALSE),
  xlab="base.deficit",ylab="age",zlab="c.peptide")

plot(seq(from=-29,to=0,length=1000),predict(np,newdata=df5),pch=".")
points(seq(from=-29,to=0,length=1000),predict(np,newdata=df10),col="blue",pch=".")
points(seq(from=-29,to=0,length=1000),predict(np,newdata=df12),col="red",pch=".")

###########
#Problem 4#
###########
nfolds<-5
case.folds<-sample(rep(1:nfolds,length.out=nrow(diabetes)))
fold.mses<-matrix(0,nrow=nfolds,ncol=2)
colnames(fold.mses)<-c("additive","kernel")
for(fold in 1:nfolds){
  train<-diabetes[case.folds!=fold,]
  test<-diabetes[case.folds==fold,]
  additive<-gam(c.peptide~s(age)+s(base.deficit),data=train)
  kernel<-npreg(c.peptide~base.deficit+age,data=train)
  add.preds<-predict(additive,newdata=test)
  kern.preds<-predict(kernel,newdata=test)
  fold.mses[fold,1]<-mean((test$c.peptide-add.preds)^2)
  fold.mses[fold,2]<-mean((test$c.peptide-kern.preds)^2)
}
colMeans(fold.mses)
#additive wins with no interaction

###########
#Problem 5#
###########
p<-c("base.deficit","age")

percent.predictionsA=function(data=diabetes,predictors=p){
  change=rep(NA,length(predictors))
  names(change)=predictors
  baseline<-colMeans(data[,predictors])
  baseline<-as.data.frame(t(baseline))
  baseline.pred<-predict(addfit,newdata=baseline)
  for(p in predictors){
    tmp=baseline
    tmp[,p]=baseline[,p]+1
    prd=predict(addfit,newdata=tmp)
    change[p]=abs(prd-baseline.pred)
  }
  return(change)
}
percent.predictionsA()

percent.predictionsB=function(data=diabetes,predictors=p){
  change=rep(NA,length(predictors))
  names(change)=predictors
  baseline<-colMeans(data[,predictors])
  baseline<-as.data.frame(t(baseline))
  baseline.pred<-predict(addfit,newdata=baseline)
  for(p in predictors){
    tmp=baseline
    tmp[,p]=baseline[,p]+.1*sd(data[,p])
    prd=predict(addfit,newdata=tmp)
    change[p]=abs(prd-baseline.pred)
  }
  return(change)
}

percent.predictionsB()

resample <- function(x) {
  sample(x,size=length(x),replace=TRUE)
}

func<-function(subset1) {
  fitA<-gam(c.peptide~s(age1)+s(base.deficit),data=diabetes,subset=subset1)
  fitB<-gam(c.peptide~s(age)+s(bd10),data=diabetes,subset=subset1)
  return(list(fitA=fitA,fitB=fitB))
}

resample.dia <- function() {
  sample.rows <- resample(1:nrow(diabetes))
  return(sample.rows)
}

evalA<-seq(from=.9,to=15.6,length.out=43)
evalB<-seq(from=-29,to=-.2,length.out=43)
evaluation.points <- data.frame(age=evalA,base.deficit=evalB)

eval<- function(add) {
  return(predict(add,newdata=evaluation.points))
}

test <- func(1:nrow(diabetes))
main.curveA <- eval(test$fitA)
main.curveB <- eval(test$fitB)

se <- function(B,alpha) {
  tbootA <- replicate(B,eval(func(resample.dia())$fitA))
  tbootB<-replicate(B,eval(func(resample.dia())$fitB))
  diff<-abs(tbootA-tbootB)
  seA<-apply(tbootA,2,sd)
  seB<-apply(tbootB,2,sd)
  seD<-apply(diff,2,sd)
  
  low.quantilesA <- apply(tbootA,1,quantile,probs=alpha/2)
  high.quantilesA <- apply(tbootA,1,quantile,probs=1-alpha/2)
  low.cisA <- 2*main.curveA - high.quantilesA
  high.cisA <- 2*main.curveA - low.quantilesA
  cisA <- rbind(mean(low.cisA),mean(high.cisA))

  low.quantilesB <- apply(tbootB,1,quantile,probs=alpha/2)
  high.quantilesB <- apply(tbootB,1,quantile,probs=1-alpha/2)
  low.cisB <- 2*main.curveB - high.quantilesB
  high.cisB <- 2*main.curveB - low.quantilesB
  cisB <- rbind(mean(low.cisB),mean(high.cisB))
  
  sedf<-data.frame(seA=seA,seB=seB,seD=seD)
  #return(sedf)
  return(list(cisA=cisA,cisB=cisB))
}
colMeans(se(50,.1))

#Confidence Intervals
se(1000,.1)
##$cisA
##         [,1]
##[1,] 3.823553
##[2,] 6.183828
##
##$cisB
##         [,1]
##[1,] 3.989556
##[2,] 5.467770
