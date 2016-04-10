context("Helper functions")
test_that("Test index matrix",{
	expect_that(indexMatrix(1,1,matrix(1)), equals(1))
	expect_that(indexMatrix(2,1,matrix(1:2,nrow=2)), equals(2))
	expect_that(indexMatrix(2,1,matrix(1:100,nrow=20)), equals(2))
	expect_that(indexMatrix(100,1,matrix(1:2,nrow=2)), throws_error("rows outside"))
	expect_that(indexMatrix(1,100,matrix(1:2,nrow=2)), throws_error("cols outside"))
	expect_that(indexMatrix(1:2,1,matrix(1:2,nrow=2)), throws_error("different"))
	expect_that(indexMatrix(c(1:3,20),c(1:3,1),matrix(1:100,nrow=20,byrow=TRUE)), equals(c(1,7,13,96)))
	expect_that(indexMatrix(c(1:3,20),c(1:3,1),matrix(1:100,nrow=20,byrow=TRUE),returnIndex=TRUE), equals(c(1,22,43,20)))
})


context("Simple DNA string manipulation")
test_that("Test ambiguous2regex",{
	expect_that(ambiguous2regex(c("AATTCCGG",'GN@!','RYHB')), equals(c("AATTCCGG",'G[ACGT]@!','[AG][CT][ACT][CGT]')))
	expect_that(grep(ambiguous2regex('ANRT'),c('ATGT','ACAT','AGTT')), equals(1:2))
})

test_that("Test bases2ambiguous",{
	expect_that(bases2ambiguous(c("AATTCCTG",'AAATCTAG','AGCTCACG','AAGTCCGT')), equals(c("ARNTCHNK")))
	expect_that(bases2ambiguous(c("AATTCCTG",'AAATCTAG','AGCTCACG','AAGT1CGT')), throws_error("[ACTG]{4}"))
	expect_that(bases2ambiguous(c("AATTCCGG",'AAACCGG','AACTCCGG','AAGTCCGG')), throws_error("length"))
	expect_that(bases2ambiguous(c("AATTCCGGTTATT")), equals(c("AATTCCGGTTATT")))
	seqs<-c("ATTTCCGG",'AGACCGGT','AACTCCGG','AAGTCCGG')
	expect_that(grep(ambiguous2regex(bases2ambiguous(seqs)),seqs),equals(1:4))
})

test_that("Test expandAmbiguous",{
	expect_that(expandAmbiguous(c("AATTCCTG",'AAATCTAG','AGCTCACG','AAGTCCGT12!@')), equals(list("AATTCCTG",'AAATCTAG','AGCTCACG','AAGTCCGT12!@')))
	expect_that(lapply(expandAmbiguous(c("AAGN",'ARV')),sort), equals(list(c("AAGA","AAGC","AAGG","AAGT"),c("AAA","AAC","AAG","AGA","AGC","AGG"))))
	expect_that(sapply(expandAmbiguous(c("AAN",'AAANN','NANAN',"HBVD","RYMKSWHBVD")),length), equals(c(4,16,64,81,2^6*3^4)))
})

test_that("Test countNmers",{
	expect_that(countNmers(c('AATT'),2),equals(c('AA'=1,'AT'=1,'TT'=1)))
	expect_that(countNmers(c('AATTT'),2),equals(c('AA'=1,'AT'=1,'TT'=2)))
	expect_that(countNmers(c('GTTTT'),2),equals(c('GT'=1,'TT'=3)))
	expect_that(countNmers(c('GTTTT'),3),equals(c('GTT'=1,'TTT'=2)))
	expect_that(countNmers(c('GTTTT'),100),gives_warning('shorter'))
})

test_that("Test reverseString",{
	expect_that(reverseString("1234\nabc_"), equals("_cba\n4321"))
	expect_that(reverseString(reverseString("1234\nabc_")), equals("1234\nabc_"))
	testStrings<-replicate(20,paste(sample(c(letters,LETTERS),1000,TRUE),collapse=''))
	revTestStrings<-sapply(strsplit(testStrings,''),function(x)paste(rev(x),collapse=''))
	expect_that(reverseString(testStrings), equals(revTestStrings))
	expect_that(reverseString("(0)1[23]4"), equals("4[32]1(0)"))
	expect_that(reverseString("(0)1[23]4",brackets=FALSE), equals("4]32[1)0("))
})


test_that("Test compliment",{
	expect_that(complimentDna(c("ACTG",'A')), equals(c("TGAC",'T')))
	expect_that(complimentDna(c(rep(paste(rep(c('AA','TTT','C','GGGG'),100),collapse=''),100),'A')), equals(c(rep(paste(rep(c('TT','AAA','G','CCCC'),100),collapse=''),100),'T')))
	expect_that(complimentDna(complimentDna(c("AATTGGCCA",'','A'))), equals(c("AATTGGCCA",'','A')))
	expect_that(complimentDna(c("XZ!@")), equals(c("XZ!@")))
	expect_that(complimentDna(c("CAGT","MNK")), equals(c("GTCA","KNM")))
	expect_that(complimentDna(c("CAGT","MNK"),ambigs=FALSE), equals(c("GTCA","MNK")))
})

test_that("Test reverse compliment",{
	expect_that(revComp(c("ACTG",'A')), equals(c("CAGT",'T')))
	expect_that(revComp(revComp(c("AATTGGCCA",'','A'))), equals(c("AATTGGCCA",'','A')))
	expect_that(revComp(c("AATTGGCCA",'','A')), equals(complimentDna(reverseString(c("AATTGGCCA",'','A')))))
	expect_that(revComp(c("A[CT]G",'HDMNK')), equals(c("C[AG]T",'MNKHD')))
})
