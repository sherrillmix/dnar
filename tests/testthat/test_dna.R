context("Test dna.R functions")
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

test_that("Test blatFindGaps",{
	expect_equal(blatFindGaps(c('1,100,200','1','1,10,20,30'),c('4,200,300','100','1,1001,2001,3001'),c('99,50,100','10000','3,4,5,6')),
		list(
			data.frame(qStartAfter=c(99,149),qGaps=c(0,50),tStartAfter=c(102,249),tGaps=c(97,50)),
			data.frame(qStartAfter=-99,qGaps=-99,tStartAfter=-99,tGaps=-99)[0,],
			data.frame(qStartAfter=c(3,13,24),qGaps=c(6,6,5),tStartAfter=c(3,1004,2005),tGaps=c(997,996,995))
		)
	)
	expect_error(blatFindGaps(c('1,100,200','1','1,10,20,30'),c('4,200,300','100','1,1001,2001,3001'),c('99,50,100','10000')),'Length')
	expect_error(blatFindGaps(c('1,100,200','1'),c('4,200,300','100','1,1001,2001,3001'),c('99,50,100','10000')),'Length')
	expect_error(blatFindGaps(c('1'),c('4,200,300','100'),c('99,50,100','10000')),'Length')
	expect_error(blatFindGaps(c('1,100','1'),c('4,200,300','100'),c('99,50,100','10000')),'length')
	expect_error(blatFindGaps(c('1,100,200','1'),c('4,200,300','100'),c('99,50','10000')),'length') 
	expect_error(blatFindGaps(c('1,100,200','1'),c('4,200,300','100'),c('150,99,50','10000')),'negative') 
	expect_error(blatFindGaps(c('1,100,300','1'),c('4,200,300','100'),c('50,101,50','10000')),'negative') 
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


test_that("Test checkOverlap",{
	queries<-data.frame('chr'=letters[rep(1:3,each=2)],'start'=c(1,1000,-1,1000,-10000,1),'end'=c(1010,2000,50,1050,10000,1))
	targets<-data.frame('chr'=letters[rep(1:2,each=2)],'start'=c(1,1111,1,1000),'end'=c(1000,9999,2,1000),'name'=letters[1:4])
	out<-checkOverlap(queries$chr,queries$start,queries$end,targets$chr,targets$start,targets$end,targets$name)
	expect_equal(out,c('a','a|b','c','d','',''))
	out<-checkOverlap(queries$chr,queries$start,queries$end,targets$chr,targets$start,targets$end,targets$name,sep=' ')
	expect_equal(out,c('a','a b','c','d','',''))
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

test_that("Test complement",{
	expect_equal(complementDna(c("ACTG",'A')), c("TGAC",'T'))
	expect_equal(complementDna(c(rep(paste(rep(c('AA','TTT','C','GGGG'),100),collapse=''),100),'A')), c(rep(paste(rep(c('TT','AAA','G','CCCC'),100),collapse=''),100),'T'))
	expect_equal(complementDna(complementDna(c("AATTGGCCA",'','A'))), c("AATTGGCCA",'','A'))
	expect_equal(complementDna(c("XZ!@")), c("XZ!@"))
	expect_equal(complementDna(c("CAGT","MNK")), c("GTCA","KNM"))
	expect_equal(complementDna(c("CAGT","MNK"),ambigs=FALSE), c("GTCA","MNK"))
})

test_that("Test reverse complement",{
	expect_equal(revComp(c("ACTG",'A')), c("CAGT",'T'))
	expect_equal(revComp(revComp(c("AATTGGCCA",'','A'))), c("AATTGGCCA",'','A'))
	expect_equal(revComp(c("AATTGGCCA",'','A')), complementDna(reverseString(c("AATTGGCCA",'','A'))))
	expect_equal(revComp(c("A[CT]G",'HDMNK')), c("C[AG]T",'MNKHD'))
})

test_that("Test degap",{
	expect_equal(degap(c("...A---CT--G!@#$%^&(),,...",'A...')), c("ACTG!@#$%^&(),,",'A'))
	expect_equal(degap(c("^^^A--CT^^G...",'A%%%'),c('^','%')), c("A--CTG...",'A'))
	expect_equal(degap(c("^^^A--CT^^G...",'A%%%'),c('^')), c("A--CTG...",'A%%%'))
	expect_equal(degap(c("A--CT\\\\^^G",'A%%%'),c('\\')), c("A--CT^^G",'A%%%'))
	expect_equal(degap(c("A--CT[][]G",'A%%%'),c('[',']')), c("A--CTG",'A%%%'))
})

test_that("Test parseRegion",{
	expect_equal(parseRegion('chr1:1234-2345'), data.frame('chr'='chr1','start'=1234,'end'=2345,'strand'='*',stringsAsFactors=FALSE))
	expect_equal(parseRegion(c('chr1:1234-2345','chrX:435-567-','chr!@#@$!:435-567+')), data.frame('chr'=c('chr1','chrX','chr!@#@$!'),'start'=c(1234,435,435),'end'=c(2345,567,567),'strand'=c('*','-','+'),stringsAsFactors=FALSE))
	expect_warning(parseRegion(c('chr1:1234-2345','chrX:435-567-','chrY:43a5-567+')), 'NA')
	expect_error(parseRegion(c('asdfea','asd')), 'asdfea, asd.*parsed')
})

test_that("Test pasteRegion",{
	expect_equal(pasteRegion('chr1',1234,2345), 'chr1:1234-2345')
	expect_equal(pasteRegion(c('chr1','chrX'),c(1234,1e9+1),c(2345,1e11+10)), c('chr1:1234-2345','chrX:1000000001-100000000010'))
	expect_equal(pasteRegion(c('chr1','chrX','chrZaasd'),c(1234,1e9+1,10),c(2345,1e11+10,20),c('-','*','+')), c('chr1:1234-2345-','chrX:1000000001-100000000010*','chrZaasd:10-20+'))
	expect_error(pasteRegion(c('chr1','chrX'),c(1234,1e9+1,10),c(2345,1e11+10,20),c('-','*','+')), 'length')
})

test_that("Test samFlag",{
	expect_equal(samFlag(c(1+2^(1:8),2^(1:8)),'paired'), rep(c(TRUE,FALSE),each=8))
	expect_equal(samFlag(c(1+2^(1:8),2^(1:8)),1), rep(c(TRUE,FALSE),each=8))
	expect_equal(samFlag(as.character(c(1+2^(1:8),2^(1:8))),'paired'), rep(c(TRUE,FALSE),each=8))
	expect_equal(samFlag(c(64+2^(1:4)),'first'), rep(TRUE,4))
	expect_equal(samFlag(c(64+2^(0:4)),c('paired','first')), c(TRUE,rep(FALSE,4)))
	expect_equal(samFlag(c(64+2^(0:4)+2^(1:5)),c('properPair','unmapped','first')), c(FALSE,TRUE,rep(FALSE,3)))
	expect_error(samFlag(1:10,'notARealFlag'),'Unknown')
})

test_that("Test trimNs",{
	expect_equal(trimNs(c('AAA','TTAT'),3),c('AAA','TTAT'))
	expect_equal(trimNs(c('AAA','TTAT'),4),c('','TTAT'))
	expect_equal(trimNs(c('AAA','TTAT'),c(3,4)),c('','TTAT'))
	expect_equal(trimNs(c('AAA','TTAT'),c(4,3)),c('','TTAT'))
	expect_equal(trimNs(c('AAA','TTAT'),c(3,3)),c('AAA','TTAT'))
	expect_equal(trimNs(c('A','NTN','NNACACANN',''),c(0,0)),c('A','NTN','NNACACANN',''))
	expect_equal(trimNs(c('A','NTN','NNACACANN'),c(1,0)),c('A','TN','ACACANN'))
	expect_equal(trimNs(c('A','NTN','NNACACANN'),c(0,1)),c('A','NT','NNACACA'))
	expect_equal(trimNs(c('A','NTN','NANACNNACANN'),2),c('','','ACNNACA'))
	expect_equal(trimNs(c('A','NTNNNNNNNNNNNN','NNNAANNNNNNNNNNNNNNAANNNNNNN'),2),c('','','AANNNNNNNNNNNNNNAA'))
	expect_equal(trimNs(c('ZZZNNNANNNZZ!Z',''),2,c('Z','!')),c('NNNANNN',''))
	expect_equal(trimNs(c('[-^ANNA-N^^^]]'),2,c('[',']','^','-')),c('ANNA'))
})

test_that("Test_cigarToBlock",{
	expect_equal(cigarToBlock('100M',2),data.frame(qStarts='1',tStarts='2',sizes='100',stringsAsFactors=FALSE))
	expect_equal(cigarToBlock(
			c('100M','10H10M1000H','10M10I20M10D2M','10S10M10H','100M10000N100M'),
			1:5
		),data.frame(
			qStarts=c('1','1','1,21,41','11','1,101'),
			tStarts=c('1','2','3,13,43','4','5,10105'),
			sizes=c('100','10','10,20,2','10','100,100'),
		stringsAsFactors=FALSE)
	)
	out<-data.frame(
		start=c(1,2,3,13,43,4,5,10105),
		end=c(100,11,12,32,44,13,104,10204),
		id=c(1,2,3,3,3,4,5,5),
		qStart=c(1,1,1,21,41,11,1,101),
		qEnd=c(100,10,10,40,42,20,100,200)
	)
	expect_equal(cigarToBlock(
			c('100M','10H10M1000H','10M10I20M10D2M','10S10M10H','100M10000N100M'),
			1:5,
			startEnds=TRUE
	),out)
	expect_error(cigarToBlock('100M',1:2),'length')
	expect_error(cigarToBlock(c('100M','10M'),1),'length')
	expect_error(cigarToBlock(c('100M','10M1Z'),1:2),'operation')
})

test_that("Test blat2exons",{
	out<-data.frame(
		chrom=rep(c('chr1','chr2'),c(5,3)),
		name=rep(c('read1','read2','read1'),c(2,3,3)),
		exonName=sprintf('%s_ex%d',rep(c('read1','read2','read1'),c(2,3,3)),c(1:3,2:1,1:3)),
		start=c(1,100,50,150,300,1000,2000,3000),
		end=c(20,299,99,199,399,1998,2998,4233),
		strand=rep(c('+','-','+'),c(2,3,3)),
		stringsAsFactors=FALSE
	)
	blat<-data.frame(
		'chr'= c('chr1','chr1','chr2'),
		'qName' = c('read1','read2','read1'),
		'starts' = c('1,100','50,150,300','1000,2000,3000'),
		'blocks' = c('20,200','50,50,100','999,999,1234'),
		'strand' = c('+','-','+'),
		stringsAsFactors=FALSE)
	extraCols<-data.frame('x'=1:3,'y'=2:4)
	extraColsOut<-extraCols[rep(1:3,c(2,3,3)),]
	extraColsOut<-extraCols[rep(1:3,c(2,3,3)),]
	rownames(extraColsOut)<-1:nrow(extraColsOut)
	extraSplits<-data.frame('xx'=c('1,2','4,4,4','5,6,7'),'yy'=c('a,b','aa,bb,cc','zz,yy,xx'),stringsAsFactors=FALSE)
	extraSplitsOut<-data.frame('xx'=as.character(c(1,2,4,4,4,5,6,7)),'yy'=c('a','b','aa','bb','cc','zz','yy','xx'),stringsAsFactors=FALSE)
	expect_equal(blat2exons(blat$chr,blat$qName,blat$starts,blat$blocks,blat$strand),out)
	expect_equal(blat2exons(blat$chr[1:2],blat$qName[1:2],blat$starts[1:2],blat$blocks[1:2],blat$strand[1:2]),out[1:5,])
	expect_equal(blat2exons(blat$chr,blat$qName,blat$starts,blat$blocks,blat$strand,extraCols=extraCols),cbind(out,extraColsOut))
	expect_equal(blat2exons(blat$chr,blat$qName,blat$starts,blat$blocks,blat$strand,extraSplits=extraSplits),cbind(out,extraSplitsOut))
	expect_equal(blat2exons(blat$chr,blat$qName,blat$starts,blat$blocks,blat$strand,extraSplits=extraSplits,extraCols=extraCols),cbind(out,extraColsOut,extraSplitsOut))
	expect_error(blat2exons(blat$chr,blat$qName,blat$starts,blat$blocks,c(blat$strand,'+')),'length')
	expect_error(blat2exons(blat$chr,c(blat$qName,'asd'),blat$starts,blat$blocks,c(blat$strand)),'length')

})

test_that("Test findIntrons",{
	expect_equal(
		findIntrons(c('1,100','50,150,300','1000,2000,3000'),c('20,299','99,199,399','1998,2998,4223')),
		data.frame('startAfter'=c(21,100,200,1999,2999)-1,'length'=c(99,149,299,1999,2999)-c(21,100,200,1999,2999)+1,'name'=c('1_in1','2_in1','2_in2','3_in1','3_in2'),stringsAsFactors=FALSE)
	)
	expect_equal(
		findIntrons(c('1,100,300','50,150,300','1000,2000,3000'),c('20,299,350','99,199,399','1998,2998,4223')),
		data.frame('startAfter'=c(21,300,100,200,1999,2999)-1,'length'=c(99,298,149,299,1999,2999)-c(21,299,100,200,1999,2999)+1,'name'=c('1_in1','1_in2','2_in1','2_in2','3_in1','3_in2'),stringsAsFactors=FALSE)
	)
	expect_equal(
		findIntrons(c('1,100,300','50,150,300','1000,2000,3000','1'),c('20,299,350','99,199,399','1998,2998,4223','10000'),names=letters[1:4]),
		data.frame('startAfter'=c(21,300,100,200,1999,2999)-1,'length'=c(99,298,149,299,1999,2999)-c(21,299,100,200,1999,2999)+1,'name'=c('a_in1','a_in2','b_in1','b_in2','c_in1','c_in2'),stringsAsFactors=FALSE)
	)
	expect_equal(
		findIntrons(c('1,100','50,150,300','1000,2000,3000'),c('20,299','99,199,399','1998,2998,4223')),
		findIntrons(list(c(1,100),c(50,150,300),c(1000,2000,3000)),list(c(20,299),c(99,199,399),c(1998,2998,4223)))
	)
	expect_equal(
		findIntrons(c('1000,2000,3000'),c('1998,2998,4223')),
		findIntrons(c(1000,2000,3000),c(1998,2998,4223))
	)
	expect_error(findIntrons(c('1,100,300','50,150,300','1000,2000,3000'),c('20,299,350','99,199,399','1998,2998,4223','1')),'length')
	expect_error(findIntrons(c('1,100,300','50,150,300','1000,2000,3000','1'),c('20,299,350','99,199,399','1998,2998,4223')),'length')
	expect_error(findIntrons(c('1,100,300,10000','50,150,300','1000,2000,3000'),c('20,299,350','99,199,399','1998,2998,4223')),'length')
	expect_error(findIntrons(c('1,100,300','50,150,300','1000,2000,3000'),c('20,299,350,10000','99,199,399','1998,2998,4223')),'length')
	expect_error(findIntrons(c('1,100,299','50,150,300','1000,2000,3000'),c('20,299,350','99,199,399','1998,2998,4223')),'negative')
})

test_that("Test pwm",{
	expect_equal(pwm(c('ACA','ACA','ACT')), matrix(c(1,0,0,0,0,1,0,0,2/3,0,0,1/3),nrow=4,ncol=3,dimnames=list(c('A','C','G','T'))))
	expect_equal(pwm(c('ACAT','ACAC','ACTG','ACTz')), matrix(c(1,0,0,0,0,1,0,0,.5,0,0,.5,0,1/3,1/3,1/3),nrow=4,ncol=4,dimnames=list(c('A','C','G','T'))))
	expect_equal(pwm(c('ACAT','ACAC','ACTG','ACTz'),chars=c('A','C','G')), matrix(c(1,0,0,0,1,0,1,0,0,0,1/2,1/2),nrow=3,ncol=4,dimnames=list(c('A','C','G'))))
	expect_equal(pwm(c('ACAT','ACAC','ACTG','ACTz'),chars=c('A','C','G'),priors=c('A'=1,'C'=1,'G'=1)), matrix(c(5/7,1/7,1/7,1/7,5/7,1/7,3/5,1/5,1/5,1/5,2/5,2/5),nrow=3,ncol=4,dimnames=list(c('A','C','G'))))
	expect_error(pwm(c('ACAT','ACAC','ACTG','ACT')),'length')
	expect_error(pwm(c('ACAT','ACAC','ACTG','ACTT'),chars=c('A','C'),priors=c('A'=1,'T'=1)),'match')
})

test_that("Test scoreFromPWM",{
	pwm<-matrix(c(1,0,0,0,0,1,0,0,2/3,0,0,1/3),nrow=4,ncol=3,dimnames=list(c('A','C','G','T')))
	expect_equal(scoreFromPWM(c('ACA','ACT','ACA','AGT','GCT'),pwm),log(c(2/3,1/3,2/3,0,0)))
	expect_error(scoreFromPWM(c('ACA','ACT','ACA','AGT','GCTA'),pwm),'length')
	expect_error(scoreFromPWM(c('ACA','ACT','ACA','AGT','GCT'),cbind(pwm,c(0,1,0,0))),'length')
	expect_warning(scoreFromPWM(c('ACA','ACT','ACA','AGT','GCZ'),pwm),'base')
	seqs<-rep(paste(rep('ACGT',100),collapse=''),200)
	pwm<-pwm(seqs)
	expect_equal(scoreFromPWM(seqs,pwm),log(rep(1,200)))
	seqs<-replicate(200,paste(sample(LETTERS,500,TRUE),collapse=''))
	pwm<-pwm(seqs,LETTERS)
	expect_true(all(scoreFromPWM(seqs,pwm)>log(1/26)*500))
})


