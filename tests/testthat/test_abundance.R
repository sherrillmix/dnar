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

