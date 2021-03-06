context("454 functions")

test_that("Test seq2flow",{
	expect_equal(seq2flow('TTTTCCCCTCCCCCC'),c('T'=4,'A'=0,'C'=4,'G'=0,'T'=1,'A'=0,'C'=6))
	expect_equal(seq2flow('T',outputLength=4),c('T'=1,'A'=0,'C'=0,'G'=0))
	expect_equal(seq2flow('C',outputLength=1),c('T'=0,'A'=0,'C'=1)) #just a minimum length not maximum
	expect_equal(seq2flow('CCCT',outputLength=1),c('T'=0,'A'=0,'C'=3,'G'=0,'T'=1))
	expect_error(seq2flow('TTTTCCCZCTCCCCCC'),'character')
	expect_equal(seq2flow('TTTTCCCCTCCCCCC',c('T','C','Z','!','#')),c('T'=4,'C'=4,'Z'=0,'!'=0,'#'=0,'T'=1,'C'=6))
	expect_equal(seq2flow('C','C'),c('C'=1))
	expect_error(seq2flow('C',c('C','C')),'[Rr]epeat')
})

test_that("Test indexFlow",{
	expect_equal(indexFlow(seq2flow('TACG'),1:4),1:4)
	expect_equal(indexFlow(seq2flow('GTACG'),1:5),3+1:5)
	expect_equal(indexFlow(seq2flow('C'),1),3)
	expect_equal(indexFlow(seq2flow('C','C'),1),1)
	expect_equal(indexFlow(seq2flow('TTTGCCCCAA'),1:13),rep(c(1,4,7,10,Inf),c(3,1,4,2,3)))
	expect_error(indexFlow(seq2flow('TTTGCCCCAA'),-1),'Negative')
	expect_error(indexFlow(seq2flow('TTTGCCCCAA'),0),'Negative')
})

test_that("Test flow2seq",{
	expect_equal(flow2seq(c('A'=2,'Z'=5,'_'=3)),'TTAAAAACCC') #doesnt use names
	expect_equal(flow2seq(c('A'=2,'Z'=5,'_'=3),flowOrder=1:3),'1122222333')
	expect_equal(flow2seq(c('A'=2,'Z'=5,'_'=3),flowOrder=rep('A',3)),'AAAAAAAAAA')
	dna<-paste(sample(c('A','C','G','T'),1000,TRUE),collapse='')
	expect_equal(flow2seq(seq2flow(dna)),dna)
	expect_equal(flow2seq(seq2flow(dna,c('A','C','G','T','Z')),c('A','C','G','T','Z')),dna)
})
