context("Helper functions")

test_that("Test isError",{
	expect_equal(isError(c(as.list(1:100),list(simpleError('Test')),as.list(letters))),rep(c(FALSE,TRUE,FALSE),c(100,1,26)))
})

test_that("Test index matrix",{
	expect_equal(indexMatrix(1,1,matrix(1)), 1)
	expect_equal(indexMatrix(2,1,matrix(1:2,nrow=2)), 2)
	expect_equal(indexMatrix(2,1,matrix(1:100,nrow=20)), 2)
	expect_error(indexMatrix(NA,1,matrix(1:2,nrow=2)), "NA")
	expect_error(indexMatrix(1,NA,matrix(1:2,nrow=2)), "NA")
	expect_error(indexMatrix(100,1,matrix(1:2,nrow=2)), "rows outside")
	expect_error(indexMatrix(1,100,matrix(1:2,nrow=2)), "cols outside")
	expect_error(indexMatrix(1:2,1,matrix(1:2,nrow=2)), "different")
	expect_equal(indexMatrix(c(1:3,20),c(1:3,1),matrix(1:100,nrow=20,byrow=TRUE)), c(1,7,13,96))
	expect_equal(indexMatrix(c(1:3,20),c(1:3,1),matrix(1:100,nrow=20,byrow=TRUE),returnIndex=TRUE), c(1,22,43,20))
	expect_equal(indexMatrix(letters[c(1:3,20)],LETTERS[c(1:3,1)],matrix(1:100,nrow=20,byrow=TRUE,dimnames=list(letters[1:20],LETTERS[1:5]))), c(1,7,13,96))
})

test_that("Test adjustWindow",{
	expect_error(adjustWindow(),NA)
	noColumns<-is.na(as.integer(Sys.getenv('COLUMNS')))
	expect_equal(c(ifelse(noColumns,NA,options('width')$width),1)[1],as.integer(Sys.getenv('COLUMNS')))
})

test_that("Test object.sizes",{
	test<-list(xxx=1:10, yyy=1:20)
	testEnv<-as.environment(test)
	first<-object.sizes(testEnv)
	test<-c(test,'abc'=list(rnorm(1000)))
	testEnv2<-as.environment(test)
	second<-object.sizes(testEnv2)
	#expect_equal(sort(c(names(first),'abc')),sort(names(second)))
})

test_that("Test stopError",{
	expect_error(stopError('Test'),'Test')
	expect_error(stopError('Test ',sum(1:3),' test'),'Test 6 test')
})

test_that("Test conservativeBoundary",{
	expect_equal(conservativeBoundary(c(-1,1),0),0)
	expect_equal(conservativeBoundary(c(-1,1),.5),.5)
	expect_equal(conservativeBoundary(c(-1,1),2),1)
	expect_equal(conservativeBoundary(c(-1,1),-1.1),-1)
	expect_equal(conservativeBoundary(c(-Inf,1),0),0)
	expect_equal(conservativeBoundary(c(-1,Inf),0),0)
	expect_equal(conservativeBoundary(c(10,Inf),0),10)
})

test_that("Test lagNA",{
	expect_equal(lagNA(1:10,1),c(2:10,NA))
	expect_equal(lagNA(1:10,5),c(6:10,rep(NA,5)))
	expect_equal(lagNA(1:10,5,-99),c(6:10,rep(-99,5)))
	expect_equal(lagNA(1:10,-1),c(NA,1:9))
	expect_equal(lagNA(1:10,-5),c(rep(NA,5),1:5))
	expect_equal(lagNA(1:10,-5,-99),c(rep(-99,5),1:5))
	expect_equal(lagNA(1:20,-50),rep(NA,20))
	expect_equal(lagNA(1:20,-50),rep(NA,20))
	expect_equal(lagNA(1:20,-20),rep(NA,20))
	expect_equal(lagNA(1:20,20),rep(NA,20))
	expect_equal(lagNA(1:20,-19),rep(c(NA,1),c(19,1)))
	expect_equal(lagNA(1:20,19),rep(c(20,NA),c(1,19)))
})

test_that("Test movingStat",{
	expect_equal(movingStat(1:10,min,2),c(1,1,1,2:8))
	expect_equal(movingStat(10:1,min,2),rev(c(1,1,1,2:8)))
	expect_equal(movingStat(20:1,max,5),c(rep(20,6),19:6))
	expect_equal(movingStat(20,max,5),c(20))
	expect_equal(movingStat(1:100,mean,0),1:100)
})

test_that("Test allSameLength",{
	expect_equal(allSameLength(1:10,2:11,letters[1:10]),TRUE)
	expect_equal(allSameLength(1:10,2:11,letters[1:11]),FALSE)
	expect_equal(allSameLength(1:11,2:11,letters[1:11]),FALSE)
	expect_equal(allSameLength(as.list(1:10),1:10),TRUE)
	test<-lapply(1:100,function(x)rnorm(100))
	expect_equal(do.call(allSameLength,test),TRUE)
	expect_equal(do.call(allSameLength,c(test,list(1:99))),FALSE)
})

test_that("Test anyArgsNA",{
	expect_equal(anyArgsNA(1:10,2:11,letters[1:10]),FALSE)
	expect_equal(anyArgsNA(1:100),FALSE)
	expect_equal(anyArgsNA(list(1:100,1:200,1:500)),FALSE)
	expect_equal(anyArgsNA(list(1:100,1:200,1:500,NA)),FALSE)
	expect_equal(anyArgsNA(list(1:100,1:200,1:500,NA),recursive=TRUE),TRUE)
	expect_equal(anyArgsNA(1:10,2:11,list(1:100,1:200,1:500,NA),recursive=TRUE),TRUE)
	expect_equal(anyArgsNA(1:10,2:11,list(1:100,1:200,1:500,list(list(NA))),recursive=TRUE),TRUE)
	expect_equal(do.call(anyArgsNA,as.list(1:1000)),FALSE)
	expect_equal(do.call(anyArgsNA,c(as.list(1:1000),list(NA))),TRUE)
})


test_that("Test cacheOperation",{
	cache<-tempfile()
	expect_equal(cacheOperation(cache,mean,1:10),mean(1:10))
	expect_message(cacheOperation(cache,mean,1:10,VOCAL=TRUE),'[Cc]ache')
	expect_message(cacheOperation(cache,mean,1:11,VOCAL=FALSE,OVERWRITE=TRUE),NA)
	expect_message(cacheOperation(cache,mean,1:10,VOCAL=FALSE,OVERWRITE=TRUE),NA)
	expect_equal(cacheOperation(cache,mean,1:10),mean(1:10))
	tmp<-new.env()
	load(cache,tmp)
	expect_equal(get('out',tmp),mean(1:10))
	expect_error(cacheOperation(cache,mean,1:11),'match')
	expect_error(cacheOperation(cache,median,1:10),'match')
	expect_error(cacheOperation(cache,median,1:11),'match')
	expect_equal(cacheOperation(cache,median,1:11,OVERWRITE=TRUE),median(1:11))
	expect_equal(cacheOperation(cache,mean,x=1:10,OVERWRITE=TRUE,EXCLUDE='x'),mean(1:10))
	expect_equal(cacheOperation(cache,mean,x=1:20,EXCLUDE='x'),mean(1:10)) #incorrect answer but expected when the md5 check is excluded
})

test_that("Test cv.glm.par",{
	test<-data.frame(x=1:10,y=1:10*4)
	expect_equal(cv.glm.par(glm(y~x,data=test),test)$pred,test$y) 
	expect_output(cv.glm.par(glm(y~x,data=test),test,vocal=TRUE),'\\.')
	expect_output(cv.glm.par(glm(y~x,data=test),test,vocal=FALSE),NA)
	test<-data.frame(x=1:11,y=c(1:10*4,1000))
	expect_equal(cv.glm.par(glm(y~x,data=test),test)$pred[11],11*4)
})


test_that("Test cleanMclapply",{
	expect_equal(cleanMclapply(1:10,2,function(x)x^2),as.list((1:10)^2))
	expect_equal(cleanMclapply(1:100,2,function(x)log(x^2)),lapply(1:100,function(x)log(x^2)))
	y<-10
	env<-environment()
	expect_error(cleanMclapply(1:10,2,function(x)x^2+y),"[Pp]roblem")
	expect_equal(cleanMclapply(1:10,2,function(x,y)x^2+y,y=y,envir=env),as.list((1:10)^2+10))
	expect_equal(cleanMclapply(1:10,2,function(x)x^2+y-z,extraCode='y<-10;z<-5'),as.list((1:10)^2+10-5))
	expect_message(cleanMclapply(1:10,2,function(x)x^2),'Logs')
	expect_message(cleanMclapply(1:10,2,function(x)x^2,VOCAL=FALSE),NA)
})

test_that("Test orderIn",{
	expect_equal(orderIn(1:10,1:10),1:10)
	expect_equal(orderIn(1:10,10:1),10:1)
	expect_equal(orderIn(1:10,10:1,decreasing=TRUE),1:10)
	expect_equal(orderIn(1:10,10:1,orderFunc=rank),10:1)
	expect_equal(orderIn(1:5,c(5,3,1,2,4)),c(5,3,1,2,4))
	expect_equal(orderIn(1:5,c(5,3,1,2,4,1:5,1:5)),c(5,3,1,2,4))
	expect_equal(orderIn(c('z','a','t','c'),letters),c(2,4,3,1))
	expect_error(orderIn(1:5,c(1),strict=TRUE),'not in')
	expect_equal(orderIn(1:5,c(1)),1:5)
})

test_that("Test fillList",{
	expect_equal(fillList(list(c('a'=1,'b'=2),c('c'=3))),list(c('a'=1,'b'=2,'c'=NA),c('c'=3,'a'=NA,'b'=NA)))
	expect_equal(fillList(list(c('a'=1,'b'=2),c('c'=3)),fill=-99),list(c('a'=1,'b'=2,'c'=-99),c('c'=3,'a'=-99,'b'=-99)))
	expect_equal(fillList(list(c('a'=1,'b'=2),c('c'=3)),namesList='a'),list(c('a'=1),c('a'=NA,'b'=99)[1]))
})

test_that("Test convertLineToUser",{
	plotFile<-tempfile()
	pdf(plotFile)
	plot(1:10,xlim=c(1,10),ylim=c(5,15),xaxs='i',yaxs='i')
	expect_equal(convertLineToUser(0,1),5)
	expect_equal(convertLineToUser(0,2),1)
	expect_equal(convertLineToUser(0,3),15)
	expect_equal(convertLineToUser(0,4),10)
	expect_gt(convertLineToUser(1,4),convertLineToUser(0,4))
	expect_lt(convertLineToUser(-1,4),convertLineToUser(0,4))
	expect_lt(convertLineToUser(-1,4),par('usr')[2])
	expect_lt(convertLineToUser(1,2),convertLineToUser(0,2))
	expect_gt(convertLineToUser(-1,2),convertLineToUser(0,2))
	expect_gt(convertLineToUser(-1,2),par('usr')[1])
	expect_error(convertLineToUser(0,5),'axis')
	dev.off()
})

test_that("Test convertUserToLine",{
	plotFile<-tempfile()
	pdf(plotFile)
	plot(1:10,xlim=c(1,10),ylim=c(5,15),xaxs='i',yaxs='i')
	expect_equal(convertUserToLine(5,1),0)
	expect_equal(convertUserToLine(1,2),0)
	expect_equal(convertUserToLine(15,3),0)
	expect_equal(convertUserToLine(10,4),0)
	expect_gt(convertUserToLine(1,4),convertUserToLine(0,4))
	expect_lt(convertUserToLine(-1,4),convertUserToLine(0,4))
	expect_equal(convertUserToLine(par('usr')[2],4),0)
	expect_equal(convertUserToLine(par('usr')[4],3),0)
	expect_equal(convertUserToLine(par('usr')[1],2),0)
	expect_equal(convertUserToLine(par('usr')[3],1),0)
	expect_lt(convertUserToLine(1,2),convertLineToUser(0,2))
	expect_gt(convertUserToLine(-1,2),convertUserToLine(0,2))
	expect_error(convertUserToLine(0,5),'axis')
	expect_equal(convertUserToLine(convertLineToUser(-10:10,1),1),-10:10)
	expect_equal(convertUserToLine(convertLineToUser(-10:10,2),2),-10:10)
	expect_equal(convertUserToLine(convertLineToUser(-10:10,3),3),-10:10)
	expect_equal(convertUserToLine(convertLineToUser(-10:10,4),4),-10:10)
	dev.off()
})

test_that("Test arrow",{
	plotFile<-tempfile()
	pdf(plotFile)
	plot(1:10,xlim=c(1,10),ylim=c(5,15),xaxs='i',yaxs='i')
	coords<-arrow(1,5,2,shaft=.2,point=.4)
	expect_equal(length(unique(na.omit(coords$y))),5)
	expect_equal(sort(unique(round(na.omit(coords$y)-2,10))),c(-.4,-.2,0,.2,.4))
	expect_equal(length(unique(na.omit(coords$x))),3)
	expect_equal(max(coords$x,na.rm=TRUE),5)
	expect_equal(min(coords$x,na.rm=TRUE),1)
	expect_true(is.list(arrow(1,5,2,concat=FALSE)))
	expect_true(is.list(arrow(1:10,11,2,concat=FALSE)))
	expect_equal(length(arrow(1:10,11,2,concat=FALSE)),10)
	expect_equal(sum(is.na(arrow(1:10,20,2)$x)),10)
	expect_equal(sum(is.na(arrow(1:10,20,2)$y)),10)
	expect_error(arrow(5.1,5,2),'border')
	dev.off()
})

test_that("Test stackRegions",{
	expect_equal(stackRegions(1:10,1:10),rep(1,10))
	expect_equal(stackRegions(1:10,1:10+.1),rep(1,10))
	expect_equal(stackRegions(1:10,1:10+100),1:10)
	expect_equal(stackRegions(1:10,1:10+4.9),rep(1:5,2))
	expect_error(stackRegions(1:10,1:11),'length')
	expect_error(stackRegions(1:10,c(-1,2:10)),'less')
})

test_that("Test wilsonInt",{
	#values calculated from http://epitools.ausvet.com.au/content.php?page=CIProportion
	expect_equal(wilsonInt(1,9),c(.0178762131,.4041500268))
	expect_equal(wilsonInt(90,10),c(0.8256343385,0.9447708629))
	expect_equal(wilsonInt(90:89,10:9),wilsonInt(90,10)) #just takes the first one
	expect_equal(wilsonInt(1,9),rev(1-wilsonInt(9,1))) 
	expect_equal(wilsonInt(900,9),rev(1-wilsonInt(9,900))) 
	expect_error(wilsonInt(900,-1),'less')
	expect_error(wilsonInt(-1,10),'less')
	expect_error(wilsonInt(-1,-1),'less')
})

test_that('Test mostAbundant',{
	expect_equal(mostAbundant(c(1:10,10)),'10')
	expect_equal(mostAbundant(c('c',letters)),'c')
	expect_true(mostAbundant(c(1:10,10,1)) %in% c('10','1'))
	expect_equal(mostAbundant(c('c',rep(letters,10))),'c')
})

test_that('Test table2vector',{
	expect_equal(table2vector(table(c(3,3,4,rep(1:4,20)))),c('1'=20,'2'=20,'3'=22,'4'=21))
	expect_equal(table2vector(table(rep('a',1000))),c('a'=1000))
	expect_error(table2vector(table(1:10,1:10)),'dimension')
	expect_error(table2vector(table(1:10,rep(1,10))),'dimension')
})

test_that('Test table2vector',{
	expect_equal(escapeRegexBracketChars(c(']','\\','-','A','b')),c('\\]','\\\\','\\-','A','b'))
	expect_equal(escapeRegexBracketChars(c('a','b','c'),c('a','c')),c('\\a','b','\\c'))
})

test_that("Test fillDown",{
  expect_equal(fillDown(c(1:5,NA,NA,6,NA,7)),c(1:5,5,5,6,6,7))
  expect_equal(fillDown(c('a','c','d',' ',NA,'')),c('a','c','d',' ',' ',' '))
  expect_equal(fillDown(c('a','c','d',' ',NA,''),c(' ',NA)),c('a','c','d','d','d',''))
  expect_equal(fillDown(c('a','c','d','',NA,''),c('')),c('a','c','d','d',NA,NA))
  expect_error(fillDown(c(NA,'c','d','')),'[Ff]irst')
})
