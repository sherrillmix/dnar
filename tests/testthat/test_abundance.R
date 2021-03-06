context("Test abundance functions")
test_that("Test shannon",{
  expect_equal(shannon(0),0) #a bit undefined
  expect_equal(shannon(1),0)
  expect_equal(shannon(c(1,1)),-log(1/2))
  expect_equal(shannon(c(1,1),2),-log2(1/2))
  expect_equal(shannon(rep(1,100)),-log(1/100))
  expect_equal(shannon(rep(1,100),10),-log10(1/100))
  expect_equal(shannon(1:3),sum(1:3*-log(1:3/6)/6))
  expect_error(shannon(c(rep(1,100),-1)),'positive')
  expect_equal(shannon(c(c(.25,.25),standardize=FALSE)),-.5*log(1/4))
})

test_that("Test rao",{
  expect_equal(rao(1,matrix(0)),0)
  expect_equal(rao(rep(1,2),matrix(c(0,1,1,0),nrow=2)),.5)
  expect_equal(rao(rep(1,2),matrix(c(0,1,1,0),nrow=2)),.5)
  expect_equal(rao(rep(1,2),matrix(c(0,2,2,0),nrow=2)),1)
  expect_error(rao(rep(1,2),as.matrix(dist(1:3))),'match')
  expect_error(rao(rep(1,3),as.matrix(dist(1:2))),'match')
  expect_equal(rao(rep(1,3),as.matrix(dist(1:3))),8/9)
  expect_equal(rao(rep(1,4),as.matrix(dist(1:4))),20/16)
  expect_equal(rao(1:3,as.matrix(dist(1:3))),28/36)
})

test_that("Test jensenShannon",{
  expect_equal(jensenShannon(1,1),0)
  expect_equal(jensenShannon(rep(1,20),rep(1,20)),0)
  expect_equal(jensenShannon(0:2,2:0),2/3)
  expect_equal(jensenShannon(0:1,1:0),1)
  expect_equal(jensenShannon(0:2,2:0,3),1+1/3*log(1/3,3)+2/3*log(2/3,3))
  expect_error(jensenShannon(rep(1,20),rep(1,21)),'[Ll]ength')
  counts1<-sample(1:100)
  counts2<-sample(1:100)
  expect_equal(jensenShannon(counts1,counts2),jensenShannon(counts2,counts1))
})

test_that("Test kullback",{
  expect_equal(kullback(1,1),0)
  expect_equal(kullback(0:1,1:0),NaN)
  expect_equal(kullback(c(1,1),0:1),Inf)
  expect_equal(kullback(0:1,c(1,1)),NaN)
  expect_equal(kullback(c(.5,.5),c(.2,.8),standardize=FALSE), sum(.5*log2(.5/c(.2,.8))))
  expect_equal(kullback(c(20,20),c(20,80),standardize=TRUE), sum(.5*log2(.5/c(.2,.8))))
  expect_equal(kullback(c(.5,.5),c(.2,.8)), sum(.5*log2(.5/c(.2,.8))))
  expect_equal(kullback(c(.5,.5),c(.2,.8),base=3), sum(.5*log(.5/c(.2,.8),3)))
  expect_equal(kullback(c(1,1,1),c(.1,.1,.8),base=3), sum(1/3*log(1/3/c(.1,.1,.8),3)))
  expect_equal(kullback(c(.1,.1,.8),c(1,1,1)), sum(c(.1,.1,.8)*log2(c(.1,.1,.8)*3)))
})

test_that("Test chao",{
  expect_equal(chao(-1),0)
  expect_equal(chao(0),0)
  expect_equal(chao(1),1)
  expect_equal(chao(rep(1,2)),3)
  expect_equal(chao(rep(1,3)),6)
  expect_equal(chao(rep(1,4)),10)
  expect_equal(chao(c(rep(1,2),rep(2,1))),3.5)
  expect_equal(chao(c(rep(1,2),rep(2,1),rep(100,100))),103.5)
  expect_equal(chao(c(rep(1,4),rep(100,100))),110)
  expect_equal(chao(c(rep(1,4),rep(2,4),rep(100,100))),108+4*3/2/5)
})

test_that("Test rarefy",{
  expect_equal(sum(rarefy(rep(1,100),0,reps=1000)),0)
  expect_equal(unique(unlist(rarefy(rep(1,100),1,reps=1000))),1)
  expect_equal(unique(unlist(rarefy(rep(1,100),100,reps=1000))),100)
  expect_equal(unique(unlist(rarefy(100,100,reps=1000))),1)
  expect_equal(unique(unlist(rarefy(rep(9,10),90,reps=1000,minObs=10))),0)
  expect_equal(unique(unlist(rarefy(rep(9,10),90,reps=1000,minObs=9))),10)
  out<-data.frame('50%'=1:2,'2.5%'=1:2,'97.5%'=1:2,row.names=c(1,100),check.names=FALSE)
  expect_equal(rarefy(c(50,50),c(1,100),reps=1000,probs=c(.5,.025,.975)),out)
  expect_equal(rarefy(c(50,50),2,reps=1000,chao=TRUE)[,'100%'],3)
  expect_gt(rarefy(1:100,20,reps=1000,chao=TRUE)[,'100%'],100)
})

test_that("Test rareEquation",{
  expect_equal(rareEquation(rep(1,100),0:100),structure(0:100,.Names=0:100))
  expect_equal(rareEquation(c(50,50),51),c('51'=2))
  expect_equal(rareEquation(2,3),c('3'=NaN))
  expect_equal(as.vector(round(rareEquation(1:20,40))),as.vector(round(rarefy(1:20,40,reps=1000)[,'50%'])))
  #example calc from http://ww2.tnstate.edu/ganter/B412%20ExtraRarefaction.html
  expect_equal(round(rareEquation(c(42,23,16,14,6,5),25),2),c('25'=5.53))
  expect_lt(abs(unlist(rarefy(c(50,50),3,reps=1000,statFunc=mean))-rareEquation(c(50,50),3)),.1)
  expect_lt(abs(unlist(rarefy(rep(50,50),40,reps=1000,statFunc=mean))-rareEquation(rep(50,50),40)),.3)
  expect_lt(abs(unlist(rarefy(1:10,10,reps=1000,statFunc=mean))-rareEquation(1:10,10)),.3)
  expect_lt(abs(unlist(rarefy(1:10,10,reps=1000,statFunc=mean,minObs=3))-rareEquation(1:10,10,minObs=3)),.2)
  expect_lt(abs(unlist(rarefy(1:100,200,reps=1000,statFunc=mean,minObs=5))-rareEquation(1:100,200,minObs=5)),.5)
  expect_equal(rareEquation(1:10,200),c('200'=NaN))
  expect_equal(rareEquation(1:10,-1),c('-1'=NaN))
  expect_equal(rareEquation(1:100,1:99,minObs=100),structure(rep(0,99),.Names=1:99))
  expect_equal(rareEquation(1:100,sum(1:100),minObs=100),structure(1,.Names=sum(1:100)))
})

test_that("Test pRare",{
  expect_equal(sum(sapply(0:10,pRare,10,10,10)),1)
  expect_equal(sum(sapply(0:11,pRare,11,12,13)),1)
  expect_equal(pRare(21,20,10,10),0)
  expect_equal(pRare(0,20,10,10),0)
  expect_lt(pRare(1,20,10,10),pRare(2,20,10,10))
  expect_lt(pRare(1,20,10,10),pRare(2,20,10,10))
  expect_lt(pRare(2,21,10,10),pRare(2,20,10,10))
  expect_equal(pRare(1,2,2,2),1/3)
  expect_equal(pRare(2,2,2,2),2/3)
  expect_equal(pRare(1,2,2,3),3*2/15)
  expect_equal(pRare(2,2,2,3),1-3*2/15)
  expect_equal(pRare(1,2,3,3),3*3/36)
  expect_equal(pRare(2,2,3,3),1-3*3/36)
  expect_equal(pRare(3,3,3,3),3^3/84)
})

test_that("Test chooseAtLeastOneFromEach",{
  expect_equal(chooseAtLeastOneFromEach(10,10,1),1)
  expect_equal(chooseAtLeastOneFromEach(10,10,2),2^10)
  expect_equal(chooseAtLeastOneFromEach(10,10,3),3^10)
  expect_equal(chooseAtLeastOneFromEach(11,10,3),3^10*20/2)
  expect_equal(chooseAtLeastOneFromEach(11,10,2),2^10*10/2)
  expect_equal(chooseAtLeastOneFromEach(11,10,4),4^10*30/2)
  expect_equal(chooseAtLeastOneFromEach(3,2,3),3^2*4/2)
  expect_equal(chooseAtLeastOneFromEach(3,2,2),4)
  expect_equal(chooseAtLeastOneFromEach(4,2,4),choose(8,4)-2)
  expect_equal(chooseAtLeastOneFromEach(3,2,4),choose(8,3)-choose(4,3)*2)
  expect_equal(chooseAtLeastOneFromEach(30,10,3),1)
  expect_equal(chooseAtLeastOneFromEach(29,10,3),30)
  expect_equal(chooseAtLeastOneFromEach(9,10,2),0)
  expect_equal(chooseAtLeastOneFromEach(-1,10,2),0)
  expect_equal(chooseAtLeastOneFromEach(10,10,0),0)
  expect_equal(chooseAtLeastOneFromEach(10,1,0),0)
  expect_equal(chooseAtLeastOneFromEach(10,0,1),0)
  expect_equal(chooseAtLeastOneFromEach(10,-1,-1),0)
})

test_that("Test unifracMatrix",{
  x<-matrix(c('a','b','c','a','b','cd'),ncol=3,byrow=TRUE)
  y<-matrix(c('a','b','c','a','c','e'),ncol=3,byrow=TRUE)
  z<-matrix(c('x','b','c'),ncol=3,byrow=TRUE)
  NAs<-matrix(c(NA,NA,'c'),ncol=3,byrow=TRUE)
  expect_equal(unifracMatrix(x,x,TRUE),0)
  expect_equal(unifracMatrix(y,y,FALSE),0)
  expect_equal(unifracMatrix(x,z,TRUE),1)
  expect_equal(unifracMatrix(y,z,FALSE),1)
  expect_equal(unifracMatrix(x,y,FALSE),3/6)
  expect_equal(unifracMatrix(x,y[c(1,1),],FALSE),1/4)
  expect_equal(unifracMatrix(x[rep(1:2,20:21),],y[rep(1:2,90:91),],FALSE),3/6)
  expect_equal(unifracMatrix(x,rbind(y,z),FALSE),6/9)
  expect_equal(unifracMatrix(x[1,,drop=FALSE],y[2,,drop=FALSE],TRUE),4/6)
  expect_equal(unifracMatrix(x,y,TRUE),.5*4/6)
  expect_equal(unifracMatrix(x,y[c(1,1),],TRUE),.5*2/6)
  expect_equal(unifracMatrix(x,y[c(1,1,1),],TRUE),.5*2/6)
  expect_equal(unifracMatrix(x,y[c(1,1,2),],TRUE),((2/3-1/2+1/3+1/2)+(1-2/3+1/3))/6)
  expect_equal(unifracMatrix(y,x,TRUE),unifracMatrix(x,y,TRUE))
  expect_equal(unifracMatrix(y,x),unifracMatrix(x,y))
  expect_equal(unifracMatrix(x,y),unifracMatrix(x,y,checkUpstream=FALSE))
  expect_equal(unifracMatrix(x[1,,drop=FALSE],z,checkUpstream=FALSE),1/3)
  expect_equal(unifracMatrix(x[1,,drop=FALSE],z),1)
  expect_equal(unifracMatrix(x,z),unifracMatrix(x,NAs))
  expect_equal(unifracMatrix(x,z,checkUpstream=TRUE),unifracMatrix(x,NAs,checkUpstream=TRUE))
  expect_error(unifracMatrix(x,NULL),'list')
  expect_error(unifracMatrix(list(x),y),'list')
  expect_equal(unifracMatrix(list(x)),matrix(0))
  expect_equal(unifracMatrix(list(x,x,x,x)),matrix(0,ncol=4,nrow=4))
  xyDist<-unifracMatrix(x,y,weighted=TRUE)
  outMat<-matrix(rep(c(0,xyDist,0,xyDist,xyDist,0,xyDist,0),each=2),nrow=4)
  expect_equal(unifracMatrix(list(x,x,y,y),weighted=TRUE),outMat)
  xyDist<-unifracMatrix(x,rbind(y,x),weighted=FALSE)
  outMat<-matrix(c(0,0,xyDist,0,0,xyDist,xyDist,xyDist,0),nrow=3)
  expect_equal(unifracMatrix(list(x,x,rbind(y,x)),weighted=FALSE),outMat)
  outMat<-matrix(c(0,1,1,0),nrow=2)
  expect_equal(unifracMatrix(list(x,NAs)),outMat)
  expect_equal(dim(unifracMatrix(rep(list(x,y),20))),c(40,40))
  expect_message(unifracMatrix(list(x,NAs),vocal=TRUE),'1')
  expect_message(unifracMatrix(list(x,NAs),vocal=TRUE),'2')
  expect_equal(unifracMatrix(list(x,y,z),weighted=FALSE),unifracMatrix(list(x,y,z),weighted=FALSE,mc.cores=2))
  out<-unifracMatrix(list(x,y,z,y,z,x),weighted=FALSE)
  expect_equal(out[upper.tri(out)],t(out)[upper.tri(out)])
  expect_equal(out[lower.tri(out)],t(out)[lower.tri(out)])
})


test_that("Test cumpaste",{
  expect_equal(cumpaste(c('AA','BBB','C')),c('AA','AA BBB','AA BBB C'))
  expect_equal(cumpaste(c(NA,'BBB','C')),c('NA','NA BBB','NA BBB C'))
  expect_equal(cumpaste(c('','BBB','C'),sep='><'),c('','><BBB','><BBB><C'))
})


test_that("Test rarefyCounts",{
  expect_equal(replicate(10,sum(rarefyCounts(1:100,123))),rep(123,10))
  expect_equal(replicate(10,length(rarefyCounts(c(1,1,1,10000),1))),rep(4,10))
  expect_equal(replicate(10,sum(rarefyCounts(c(1,1,1,1),3)>0)),rep(3,10))
  expect_equal(rarefyCounts(100:1,sum(100:1)),100:1)
  expect_equal(rarefyCounts(1:200,sum(1:200)),1:200)
  expect_equal(rarefyCounts(1:200,0),rep(0,200))
  expect_equal(rarefyCounts(c('a'=1,'b'=2,'zz'=3),0),c('a'=0,'b'=0,'zz'=0))
  expect_error(rarefyCounts(1:200,-1),'invalid')
  expect_error(rarefyCounts(1:2,10),'larger')
})

test_that("Test gini",{
  #http://shlegeris.com/gini
  expect_equal(gini(1:10),.3)
  expect_equal(gini(20:1),gini(1:20))
  expect_equal(gini(0:1),.5)
  expect_equal(gini(c(rep(0,99),1)),.99)
  expect_equal(gini(rep(1,99)),0)
})
