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
	expect_that(ambiguous2regex(c("AATTCCGG",'GN@!','RYHB','')), equals(c("AATTCCGG",'G[ACGT]@!','[AG][CT][ACT][CGT]','')))
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
	expect_that(lapply(expandAmbiguous(c("AAGN",'ARV','')),sort), equals(list(c("AAGA","AAGC","AAGG","AAGT"),c("AAA","AAC","AAG","AGA","AGC","AGG"),'')))
	expect_that(sapply(expandAmbiguous(c("AAN",'AAANN','NANAN',"HBVD","RYMKSWHBVD")),length), equals(c(4,16,64,81,2^6*3^4)))
})

test_that("Test countNmers",{
	expect_that(countNmers(c('AATT'),2),equals(c('AA'=1,'AT'=1,'TT'=1)))
	expect_that(countNmers(c('AATTT'),2),equals(c('AA'=1,'AT'=1,'TT'=2)))
	expect_that(countNmers(c('GTTTT'),2),equals(c('GT'=1,'TT'=3)))
	expect_that(countNmers(c('GTTTT'),3),equals(c('GTT'=1,'TTT'=2)))
	expect_that(countNmers(c('GTTTT'),100),gives_warning('shorter'))
})

test_that("Test gcPercent",{
	expect_that(gcPercent(c('AATT','AATTGGCC','AGAGAG','GGCCCC','ATGCGC')),equals(c(0,.5,.5,1,4/6)))
	expect_that(gcPercent(c('AATT','AATTGGCC','AGAGAG','GGCCCC','ATGCGC'),'A'),equals(c(.5,.25,.5,0,1/6)))
	expect_that(gcPercent(c('AATT','AATTGGCC','AGAGAG','GGCCCC','ATGCGC'),c('A','T')),equals(c(1,.5,.5,0,2/6)))
	expect_that(gcPercent(paste(rep('ATGC',1000),collapse='')),equals(.5))
	expect_that(gcPercent(''),equals(NaN))
})

test_that("Test toggleCase",{
	expect_that(toggleCase(c("1234\nabc_","AbCdE",'12345!@','')), equals(c("1234\nABC_","aBcDe",'12345!@','')))
	expect_that(toggleCase(paste(letters,LETTERS,collapse='\n',sep='_@!$%^')), equals(paste(LETTERS,letters,collapse='\n',sep='_@!$%^')))
})

test_that("Test highlightString",{
	expect_that(highlightString('AGA',c('TTTAGATTTAGAT','','!@#@$AGATTT123','TAIQWEUQOWIEUO','AGAGA')), equals(c('TTTagaTTTagaT','','!@#@$agaTTT123','TAIQWEUQOWIEUO','agaGA')))
	expect_that(highlightString('aga',tolower(c('TTTAGATTTAGAT','','!@#@$AGATTT123','TAIQWEUQOWIEUO','AGAGA'))), equals(toggleCase(c('TTTagaTTTagaT','','!@#@$agaTTT123','TAIQWEUQOWIEUO','agaGA'))))
	expect_that(highlightString('AG1A',c('TTTAG1ATTTAG1AT')), equals(c('TTTag1aTTTag1aT')))
})

test_that("Test highlightDifferences",{
	expect_that(highlightDifferences('AGA1123','AGA1123'), equals('AGA1123'))
	expect_that(highlightDifferences('AGAAAT','AGACAT'), equals('AGAaAT'))
	expect_that(highlightDifferences('agaaat','agacat'), equals('agaAat'))
	expect_that(highlightDifferences('',''), equals(''))
	expect_that(highlightDifferences("ACA","ACAA"), throws_error("length"))
})

test_that("Test seqSplit",{
	expect_that(seqSplit('ATA','ACA','A'), equals(matrix(c('A','T','A','A','C','A','A','.','.'),nrow=3,byrow=TRUE)))
	expect_that(seqSplit('ATA',c('ACA','A'),fill='Z'), equals(matrix(c('A','T','A','A','C','A','A','Z','Z'),nrow=3,byrow=TRUE)))
	expect_that(dim(seqSplit('','','')), equals(c(3,0)))
	expect_that(dim(seqSplit(replicate(100,paste(sample(letters,1000,TRUE),collapse='')))), equals(c(100,1000)))
	expect_that(seqSplit("ACA","ACAA",fill=NULL), throws_error("length"))
})

test_that("Test dna2codons",{
	expect_that(dna2codons(c('GTTGAA','ACGTTT123','GTAAA')),equals(list(c('GTT','GAA'),c('ACG','TTT','123'),c('GTA'))))
	expect_that(dna2codons(c('GTTGAA','ACGTTT123','GTAAA'),frame=0:2),equals(list(c('GTT','GAA'),c('CGT','TT1'),c('AAA'))))
	expect_that(dna2codons(c('GTTGAA','ACGTTT123','GTAAA'),frame=2),equals(list(c('TGA'),c('GTT','T12'),c('AAA'))))
	expect_that(dna2codons(c('GT')),gives_warning('shorter'))
	expect_that(dna2codons(c('GTT','GT')),gives_warning('shorter'))
	expect_that(dna2codons(c('GTT'),frame=1),gives_warning('shorter'))
	expect_that(dna2codons(c('GTTA'),frame=2),gives_warning('shorter'))
})

test_that("Test codons2aa",{
	expect_that(codons2aa(c('AUG','GTT','ACG','GTA')),equals(c('M','V','T','V')))
	expect_that(codons2aa(c('AUG','GTT','A1G','GTA'),naReplace='@',warn=FALSE),equals(c('M','V','@','V')))
	expect_that(codons2aa(c('AUG','GTT','A1G','GTA'),'abbr',naReplace='@',warn=FALSE),equals(c('Met','Val','@','Val')))
	expect_that(codons2aa(c('123')),gives_warning('Unknown'))
})

test_that("Test dna2aa",{
	expect_that(dna2aa(c('AA1')),gives_warning('Unknown'))
	expect_that(dna2aa(c('AA1'),warn=FALSE,naReplace='@'),equals('@'))
	expect_that(dna2aa(c('ATGAA1'),warn=FALSE,naReplace='z'),equals('Mz'))
	expect_that(dna2aa(c('ATGAA')),equals('M'))
	expect_that(dna2aa(c('ATGAAT','ATGGATGATTAGGAG')),equals(c('MN','MDDXE')))
})

test_that("Test aa2codons",{
	expect_that(aa2codons(c('1','G')),throws_error('Unknown'))
	expect_that(lapply(aa2codons(c('M','A','A','X')),function(x)sort(strsplit(gsub('[()]','',x),'\\|')[[1]])),equals(list(c('AUG'),c('GCA','GCC','GCG','GCU'),c('GCA','GCC','GCG','GCU'),c('UAA','UAG','UGA')))) 
	expect_that(lapply(aa2codons(c('M','A','A','X'),regex=FALSE),sort),equals(list(c('AUG'),c('GCA','GCC','GCG','GCU'),c('GCA','GCC','GCG','GCU'),c('UAA','UAG','UGA')))) 
	expect_that(lapply(aa2codons(c('M','A','A','X')),function(x)sort(strsplit(gsub('[()]','',x),'\\|')[[1]])),equals(lapply(aa2codons(c('M','A','A','X'),regex=FALSE),sort))) 
	expect_that(lapply(aa2codons(c('Met','Ala','Ala','Xxx',"Glu"),type='abbr',regex=FALSE),sort),equals(lapply(aa2codons(c('M','A','A','X','E'),regex=FALSE),sort))) 
	expect_that(lapply(aa2codons(c('Methionine','Alanine','Alanine'),type='name',regex=FALSE),sort),equals(lapply(aa2codons(c('M','A','A'),regex=FALSE),sort))) 
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
