context("Helper functions")
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



context("Simple DNA string manipulation")
test_that("Test ambiguous2regex",{
	expect_equal(ambiguous2regex(c("AATTCCGG",'GN@!','RYHB','')), c("AATTCCGG",'G[ACGT]@!','[AG][CT][ACT][CGT]',''))
	expect_equal(grep(ambiguous2regex('ANRT'),c('ATGT','ACAT','AGTT')), 1:2)
})

test_that("Test bases2ambiguous",{
	expect_equal(bases2ambiguous(c("AATTCCTG",'AAATCTAG','AGCTCACG','AAGTCCGT')), c("ARNTCHNK"))
	expect_error(bases2ambiguous(c("AATTCCTG",'AAATCTAG','AGCTCACG','AAGT1CGT')), "[ACTG]{4}")
	expect_error(bases2ambiguous(c("AATTCCGG",'AAACCGG','AACTCCGG','AAGTCCGG')), "length")
	expect_equal(bases2ambiguous(c("AATTCCGGTTATT")), c("AATTCCGGTTATT"))
	seqs<-c("ATTTCCGG",'AGACCGGT','AACTCCGG','AAGTCCGG')
	expect_equal(grep(ambiguous2regex(bases2ambiguous(seqs)),seqs),1:4)
})

test_that("Test expandAmbiguous",{
	expect_equal(expandAmbiguous(c("AATTCCTG",'AAATCTAG','AGCTCACG','AAGTCCGT12!@')), list("AATTCCTG",'AAATCTAG','AGCTCACG','AAGTCCGT12!@'))
	expect_equal(lapply(expandAmbiguous(c("AAGN",'ARV','')),sort), list(c("AAGA","AAGC","AAGG","AAGT"),c("AAA","AAC","AAG","AGA","AGC","AGG"),''))
	expect_equal(sapply(expandAmbiguous(c("AAN",'AAANN','NANAN',"HBVD","RYMKSWHBVD")),length), c(4,16,64,81,2^6*3^4))
})

test_that("Test countNmers",{
	expect_equal(countNmers(c('AATT'),2),c('AA'=1,'AT'=1,'TT'=1))
	expect_equal(countNmers(c('AATTT'),2),c('AA'=1,'AT'=1,'TT'=2))
	expect_equal(countNmers(c('GTTTT'),2),c('GT'=1,'TT'=3))
	expect_equal(countNmers(c('GTTTT'),3),c('GTT'=1,'TTT'=2))
	expect_warning(countNmers(c('GTTTT'),100),'shorter')
})

test_that("Test gcPercent",{
	expect_equal(gcPercent(c('AATT','AATTGGCC','AGAGAG','GGCCCC','ATGCGC')),c(0,.5,.5,1,4/6))
	expect_equal(gcPercent(c('AATT','AATTGGCC','AGAGAG','GGCCCC','ATGCGC'),'A'),c(.5,.25,.5,0,1/6))
	expect_equal(gcPercent(c('AATT','AATTGGCC','AGAGAG','GGCCCC','ATGCGC'),c('A','T')),c(1,.5,.5,0,2/6))
	expect_equal(gcPercent(paste(rep('ATGC',1000),collapse='')),.5)
	expect_equal(gcPercent(''),NaN)
})

test_that("Test toggleCase",{
	expect_equal(toggleCase(c("1234\nabc_","AbCdE",'12345!@','')), c("1234\nABC_","aBcDe",'12345!@',''))
	expect_equal(toggleCase(paste(letters,LETTERS,collapse='\n',sep='_@!$%^')), paste(LETTERS,letters,collapse='\n',sep='_@!$%^'))
})

test_that("Test highlightString",{
	expect_equal(highlightString('AGA',c('TTTAGATTTAGAT','','!@#@$AGATTT123','TAIQWEUQOWIEUO','AGAGA')), c('TTTagaTTTagaT','','!@#@$agaTTT123','TAIQWEUQOWIEUO','agaGA'))
	expect_equal(highlightString('aga',tolower(c('TTTAGATTTAGAT','','!@#@$AGATTT123','TAIQWEUQOWIEUO','AGAGA'))), toggleCase(c('TTTagaTTTagaT','','!@#@$agaTTT123','TAIQWEUQOWIEUO','agaGA')))
	expect_equal(highlightString('AG1A',c('TTTAG1ATTTAG1AT')), c('TTTag1aTTTag1aT'))
})

test_that("Test highlightDifferences",{
	expect_equal(highlightDifferences('AGA1123','AGA1123'), 'AGA1123')
	expect_equal(highlightDifferences('AGAAAT','AGACAT'), 'AGAaAT')
	expect_equal(highlightDifferences('agaaat','agacat'), 'agaAat')
	expect_equal(highlightDifferences('',''), '')
	expect_error(highlightDifferences("ACA","ACAA"), "length")
})

test_that("Test seqSplit",{
	expect_equal(seqSplit('ATA','ACA','A'), matrix(c('A','T','A','A','C','A','A','.','.'),nrow=3,byrow=TRUE))
	expect_equal(seqSplit('ATA',c('ACA','A'),fill='Z'), matrix(c('A','T','A','A','C','A','A','Z','Z'),nrow=3,byrow=TRUE))
	expect_equal(dim(seqSplit('','','')), c(3,0))
	expect_equal(dim(seqSplit(replicate(100,paste(sample(letters,1000,TRUE),collapse='')))), c(100,1000))
	expect_error(seqSplit("ACA","ACAA",fill=NULL), "length")
})

test_that("Test dna2codons",{
	expect_equal(dna2codons(c('GTTGAA','ACGTTT123','GTAAA')),list(c('GTT','GAA'),c('ACG','TTT','123'),c('GTA')))
	expect_equal(dna2codons(c('GTTGAA','ACGTTT123','GTAAA'),frame=0:2),list(c('GTT','GAA'),c('CGT','TT1'),c('AAA')))
	expect_equal(dna2codons(c('GTTGAA','ACGTTT123','GTAAA'),frame=2),list(c('TGA'),c('GTT','T12'),c('AAA')))
	expect_warning(dna2codons(c('GT')),'shorter')
	expect_warning(dna2codons(c('GTT','GT')),'shorter')
	expect_warning(dna2codons(c('GTT'),frame=1),'shorter')
	expect_warning(dna2codons(c('GTTA'),frame=2),'shorter')
})

test_that("Test codons2aa",{
	expect_equal(codons2aa(c('AUG','GTT','ACG','GTA')),c('M','V','T','V'))
	expect_equal(codons2aa(c('AUG','GTT','A1G','GTA'),naReplace='@',warn=FALSE),c('M','V','@','V'))
	expect_equal(codons2aa(c('AUG','GTT','A1G','GTA'),'abbr',naReplace='@',warn=FALSE),c('Met','Val','@','Val'))
	expect_warning(codons2aa(c('123')),'Unknown')
})

test_that("Test dna2aa",{
	expect_warning(dna2aa(c('AA1')),'Unknown')
	expect_equal(dna2aa(c('AA1'),warn=FALSE,naReplace='@'),'@')
	expect_equal(dna2aa(c('ATGAA1'),warn=FALSE,naReplace='z'),'Mz')
	expect_equal(dna2aa(c('ATGAA')),'M')
	expect_equal(dna2aa(c('ATGAAT','ATGGATGATTAGGAG')),c('MN','MDDXE'))
})

test_that("Test aa2codons",{
	expect_error(aa2codons(c('1','G')),'Unknown')
	expect_equal(lapply(aa2codons(c('M','A','A','X')),function(x)sort(strsplit(gsub('[()]','',x),'\\|')[[1]])),list(c('AUG'),c('GCA','GCC','GCG','GCU'),c('GCA','GCC','GCG','GCU'),c('UAA','UAG','UGA'))) 
	expect_equal(lapply(aa2codons(c('M','A','A','X'),regex=FALSE),sort),list(c('AUG'),c('GCA','GCC','GCG','GCU'),c('GCA','GCC','GCG','GCU'),c('UAA','UAG','UGA'))) 
	expect_equal(lapply(aa2codons(c('M','A','A','X')),function(x)sort(strsplit(gsub('[()]','',x),'\\|')[[1]])),lapply(aa2codons(c('M','A','A','X'),regex=FALSE),sort)) 
	expect_equal(lapply(aa2codons(c('Met','Ala','Ala','Xxx',"Glu"),type='abbr',regex=FALSE),sort),lapply(aa2codons(c('M','A','A','X','E'),regex=FALSE),sort)) 
	expect_equal(lapply(aa2codons(c('Methionine','Alanine','Alanine'),type='name',regex=FALSE),sort),lapply(aa2codons(c('M','A','A'),regex=FALSE),sort)) 
})


test_that("Test aa2dna",{
	expect_error(aa2dna(c('1','G')),'Unknown')
	expect_equal(grep(aa2dna(c('MAAX')),c('ATGGCTGCTTAG','ATGGCTGCTTAA','GTGGCTGCTTAA','ATGGCTGCTTAC','TTTTTATGGCTGCTTAGTTTT')),c(1,2,5)) 
	expect_equal(sapply(strsplit(aa2dna(c('MAAX','MMFFYYS',paste(sample(dnar::aminoAcids$code,50,TRUE),collapse=''))),'[()]+'),length),c(4,7,50)+1)
})

test_that("Test gap2NoGap",{
	expect_equal(gap2NoGap("AA--CC--TT",1:10), c(1,2,2,2,3,4,4,4,5,6))
	expect_equal(gap2NoGap("AA**CC..TT",1:10), c(1,2,2,2,3,4,4,4,5,6))
	expect_equal(is.na(gap2NoGap("AA**CC..TT",11)), TRUE)
	expect_equal(gap2NoGap(sprintf("A%sT",paste(rep('-',1000),collapse='')),1:1002), c(rep(1,1001),2))
})

test_that("Test noGap2Gap",{
	expect_equal(noGap2Gap("AA--CC--TT",1:6), c(1,2,5,6,9,10))
	expect_equal(noGap2Gap("AA**CC..TT",1:6), c(1,2,5,6,9,10))
	expect_equal(is.na(noGap2Gap("AA**CC..TT",7)), TRUE)
	seq<-sprintf("A%sT-A--T---T",paste(rep('-',1000),collapse=''))
	expect_equal(noGap2Gap(seq,1:2), c(1,1002))
	expect_equal(gap2NoGap(seq,noGap2Gap(seq,1:5)), 1:5)
})

test_that("Test binary2range",{
	expect_equal(binary2range(c(TRUE,TRUE,FALSE,FALSE,FALSE)), data.frame('start'=1,'end'=2))
	expect_equal(binary2range(c(FALSE,FALSE,FALSE,TRUE,TRUE)), data.frame('start'=4,'end'=5))
	expect_equal(binary2range(c(TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE)), data.frame('start'=c(1,3,8),'end'=c(1,4,9)))
	expect_equal(binary2range(c(rep(TRUE,100),rep(FALSE,200),rep(TRUE,300),FALSE,rep(TRUE,20))), data.frame('start'=c(1,301,602),'end'=c(100,600,621)))
	expect_equal(binary2range(rep(FALSE,100)), data.frame('start'=NA,'end'=NA)[0,])
})

test_that("Test index2range",{
	expect_equal(index2range(c(1:100,200:250)), data.frame('start'=c(1,200),'end'=c(100,250)))
	expect_equal(index2range(c(1:100,2:110)), data.frame('start'=c(1),'end'=c(110)))
	seq<-c(rep(TRUE,100),rep(FALSE,200),rep(TRUE,300),FALSE,rep(TRUE,20))
	expect_equal(index2range(which(seq)), data.frame('start'=c(1,301,602),'end'=c(100,600,621)))
	expect_equal(index2range(which(seq)), binary2range(seq))
	expect_equal(index2range(c()), data.frame('start'=NA,'end'=NA)[0,])
})

test_that("Test reverseString",{
	expect_equal(reverseString("1234\nabc_"), "_cba\n4321")
	expect_equal(reverseString(reverseString("1234\nabc_")), "1234\nabc_")
	testStrings<-replicate(20,paste(sample(c(letters,LETTERS),1000,TRUE),collapse=''))
	revTestStrings<-sapply(strsplit(testStrings,''),function(x)paste(rev(x),collapse=''))
	expect_equal(reverseString(testStrings), revTestStrings)
	expect_equal(reverseString("(0)1[23]4"), "4[32]1(0)")
	expect_equal(reverseString("(0)1[23]4",brackets=FALSE), "4]32[1)0(")
})

test_that("Test compliment",{
	expect_equal(complimentDna(c("ACTG",'A')), c("TGAC",'T'))
	expect_equal(complimentDna(c(rep(paste(rep(c('AA','TTT','C','GGGG'),100),collapse=''),100),'A')), c(rep(paste(rep(c('TT','AAA','G','CCCC'),100),collapse=''),100),'T'))
	expect_equal(complimentDna(complimentDna(c("AATTGGCCA",'','A'))), c("AATTGGCCA",'','A'))
	expect_equal(complimentDna(c("XZ!@")), c("XZ!@"))
	expect_equal(complimentDna(c("CAGT","MNK")), c("GTCA","KNM"))
	expect_equal(complimentDna(c("CAGT","MNK"),ambigs=FALSE), c("GTCA","MNK"))
})

test_that("Test reverse compliment",{
	expect_equal(revComp(c("ACTG",'A')), c("CAGT",'T'))
	expect_equal(revComp(revComp(c("AATTGGCCA",'','A'))), c("AATTGGCCA",'','A'))
	expect_equal(revComp(c("AATTGGCCA",'','A')), complimentDna(reverseString(c("AATTGGCCA",'','A'))))
	expect_equal(revComp(c("A[CT]G",'HDMNK')), c("C[AG]T",'MNKHD'))
})

test_that("Test degap",{
	expect_equal(degap(c("...A---CT--G!@#$%^&(),,...",'A...')), c("ACTG!@#$%^&(),,",'A'))
	expect_equal(degap(c("^^^A--CT^^G...",'A%%%'),c('^','%')), c("A--CTG...",'A'))
	expect_equal(degap(c("^^^A--CT^^G...",'A%%%'),c('^')), c("A--CTG...",'A%%%'))
	expect_equal(degap(c("A--CT\\\\^^G",'A%%%'),c('\\')), c("A--CT^^G",'A%%%'))
	expect_equal(degap(c("A--CT[][]G",'A%%%'),c('[',']')), c("A--CTG",'A%%%'))
})
