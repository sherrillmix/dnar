
context("Simple DNA string manipulation")
test_that("Test reverseString",{
	expect_that(reverseString("1234\nabc_"), equals("_cba\n4321"))
	expect_that(reverseString(reverseString("1234\nabc_")), equals("1234\nabc_"))
	testStrings<-replicate(20,paste(sample(c(letters,LETTERS),1000,TRUE),collapse=''))
	revTestStrings<-sapply(strsplit(testStrings,''),function(x)paste(rev(x),collapse=''))
	expect_that(reverseString(testStrings), equals(revTestStrings))
})


test_that("Test compliment",{
	expect_that(complimentDna(c("ACTG",'A')), equals(c("TGAC",'T')))
	expect_that(complimentDna(complimentDna(c("AATTGGCCA",'','A'))), equals(c("AATTGGCCA",'','A')))
})

test_that("Test reverse compliment",{
	expect_that(revComp(c("ACTG",'A')), equals(c("CAGT",'T')))
	expect_that(revComp(revComp(c("AATTGGCCA",'','A'))), equals(c("AATTGGCCA",'','A')))
})
