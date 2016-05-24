context("Test readFile functions")

test_that("Test read.fastq",{
})

test_that("Test read.fa",{
	expect_equal(read.fa(textConnection('')),NULL)
	dat<-'>read1\nAATTGG\n>read2\nGGT\nTCC\nGGG'
	out<-data.frame('name'=c('read1','read2'),'seq'=c('AATTGG','GGTTCCGGG'),stringsAsFactors=FALSE)
	expect_equal(read.fa(textConnection(dat)),out)
	tmp<-tempfile()
	writeLines(dat,tmp)
	expect_equal(read.fa(tmp),out)
	tmpGz<-tempfile(fileext='.gz')
	gzHandle<-gzfile(tmpGz)
	writeLines(dat,gzHandle)
	close(gzHandle)
	expect_equal(read.fa(tmpGz),out)
	dat<-'>read1\nAATTGG\n>read2\nGGTTCCGGG'
	expect_equal(read.fa(textConnection(dat)),read.fa(textConnection(dat),assumeSingleLine=TRUE))
	dat<-'>read1\nAATTGG      \n>read2\nGGTTCCGGG         '
	expect_equal(read.fa(textConnection(dat)),out)
	dat<-'read1\nAATTGG\n>read2\nGGTTCCGGG'
	#error if missing >
	expect_error(read.fa(textConnection(dat),assumeSingleLine=TRUE),'read')
	out2<-data.frame('name'=c('read2'),'seq'=c('GGTTCCGGG'),stringsAsFactors=FALSE)
	#throw out first
	expect_equal(read.fa(textConnection(dat)),out2)
})
