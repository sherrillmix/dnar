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

