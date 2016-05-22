context("Test readFile functions")

test_that("Test read.fastq",{
})

test_that("Test read.fa",{
	dat<-'>read1\nAATTGG\n>read2\nGGTTCC'
	out<-data.frame('name'=c('read1','read2'),'seq'=c('AATTGG','GGTTCC'),stringsAsFactors=FALSE)
	expect_equal(read.fa(textConnection(dat)),out)
})
