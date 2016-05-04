
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


