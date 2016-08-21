context("Test readFile functions")

test_that("Test read.fastq",{
  fastq<-c(
    "@seq1",
    "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT",
    "+",
    "!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65",
    "@seq2",
    "GACCGGAACT",
    "+",
    "!''*(!!!*!"
  )
  out<-data.frame('name'=c('seq1','seq2'),'seq'=c('GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT','GACCGGAACT'),'qual'=c("!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65","!''*(!!!*!"),stringsAsFactors=FALSE)
  expect_equal(read.fastq(textConnection(fastq)),out)
  tmpFile<-tempfile()
  writeLines(fastq,tmpFile)
  expect_equal(read.fastq(tmpFile),out)
  shortq<-fastq
  shortq[8]<-substring(shortq[8],1,nchar(shortq[8])-1)
  expect_error(read.fastq(textConnection(shortq)),'length')
  shortq<-fastq
  shortq[6]<-substring(shortq[6],1,nchar(shortq[6])-1)
  expect_error(read.fastq(textConnection(shortq)),'length')
  expect_true(all(grepl('^[0-9 ]+',read.fastq(textConnection(fastq),convert=TRUE)$qual)))
  expect_equal(read.fastq(textConnection('@1\nAA\n+\nI@'),convert=TRUE)$qual,'40 31')
  tmpGz<-tempfile(fileext='.gz')
  gzHandle<-gzfile(tmpGz)
  writeLines(fastq,gzHandle)
  close(gzHandle)
  expect_equal(read.fastq(tmpGz),out)
  expect_lt(file.size(tmpGz),file.size(tmpFile))
  #replace first G with Z
  expect_warning(read.fastq(textConnection(sub('G','Z',fastq))),'character')
})

test_that("Test intsToQual",{
  expect_equal(intsToQual(2:5),'#$%&')
  expect_equal(intsToQual(2:4,34),'$%&')
})

test_that("Test qualToInts",{
  expect_equal(qualToInts('#$%&'),list(2:5))
  expect_equal(qualToInts('$%&',34),list(2:4))
  expect_equal(qualToInts(c(intsToQual(1:40),intsToQual(20:10))),list(1:40,20:10))
  expect_equal(qualToInts(c(intsToQual(1:40,10),intsToQual(20:10,10)),10),list(1:40,20:10))
})

test_that("Test read.phd",{ tmpFile<-tempfile()
  phd<-'Some\nextra lines\nBEGIN_DNA\nc 9 6\nt 9 18\nc 10 26\nc 19 38\nA 40 39\nC 15 45\nEND_DNA\nmore\nextra lines'
  expect_equal(read.phd(textConnection(phd))[['seq']],'ctccAC')
  expect_equal(read.phd(textConnection(phd))[['qual']],'9 9 10 19 40 15') 
  expect_equal(read.phd(textConnection(phd),19)[['seq']],'cA')
  expect_equal(read.phd(textConnection(phd),19)[['qual']],'19 40') 
  expect_equal(read.phd(textConnection(phd),20)[['seq']],'A')
  expect_equal(read.phd(textConnection(phd),20)[['qual']],'40') 
  writeLines(phd,tmpFile)
  expect_equal(read.phd(textConnection(phd)),read.phd(tmpFile)) 
})

test_that("Test readFaDir",{
  tmpDir<-tempdir()
  files<-file.path(tmpDir,sprintf('%d.fa',1:9))
  seqs<-lapply(files,function(x)generateFasta(sample(1:100,1),10:1000,c('A','C','T','G','\n')))
  out<-do.call(rbind,seqs)
  out$seq<-gsub('\n','',out$seq)
  out$name<-sub('^>','',out$name)
  out$file<-rep(basename(files),sapply(seqs,nrow))
  mapply(function(x,file)write.fa(x$name,x$seq,file),seqs,files)
  read<-readFaDir(tmpDir)
  expect_equal(read,out)
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
  expect_equal(read.fa(textConnection('')),NULL)
  expect_equal(read.fa(textConnection(rep('\n\n',10))),NULL)
  expect_equal(read.fa(textConnection(rep('AA\nCC\n',10))),NULL)
})

test_that("Test generateFasta",{
  fasta<-generateFasta(1000,90:110,c('A','Z'))
  expect_equal(nrow(fasta),1000)
  expect_true(all(grepl('^>',fasta$name)))
  expect_true(!any(grepl('[^AZ]',fasta$seq)))
  expect_true(all(nchar(fasta$seq) %in% c(90:110)))
  expect_equal(length(unique(fasta$name)),1000)
  fastq<-generateFasta(100,90:110,generateQuals=TRUE,qualRange=2:39)
  expect_equal(nrow(fastq),100)
  expect_equal(nchar(fastq$seq),nchar(fastq$qual))
  expect_true(all(unlist(qualToInts(fastq$qual)) %in% 2:39))
})

test_that("Test write.fa",{
  fa<-generateFasta(100,10:300)
  fa$noGreater<-sub('^>','',fa$name)
  tmpFile<-tempfile()
  write.fa(fa$noGreater,fa$seq,tmpFile)
  expect_equal(read.fa(tmpFile)$seq,fa$seq)
  expect_equal(read.fa(tmpFile)$name,fa$noGreater)
  fa<-generateFasta(100,10:300,c('A','C','T','G','\n'))
  fa$noGreater<-sub('^>','',fa$name)
  tmpGz<-tempfile(fileext='.gz')
  write.fa(fa$name,fa$seq,tmpGz)
  expect_equal(read.fa(tmpGz)$seq,gsub('\n','',fa$seq))
  expect_equal(read.fa(tmpGz)$name,fa$noGreater)
  expect_lt(file.size(tmpGz),file.size(tmpFile))
})

test_that("Test fillZeros",{
  expect_equal(fillZeros(c(1:4,7,9),data.frame('x'=1:6,'aa'=letters[1:6],stringsAsFactors=FALSE)),list('pos'=1:9,'data'=data.frame('x'=c(1:4,0,0,5,0,6),'aa'=c('a','b','c','d',0,0,'e',0,'f'),stringsAsFactors=FALSE)))
  expect_equal(fillZeros(c(1:4,7,9),data.frame('x'=letters[1:6]),'a')[['data']],data.frame('x'=c('a','b','c','d','a','a','e','a','f')))
  expect_equal(fillZeros(c(13,sample(14:999,40),1000),data.frame('x'=1:42))[['pos']],13:1000)
  expect_error(fillZeros(c(1:4,7,9),data.frame('x'=1:7)),'length')
  expect_error(fillZeros(c(1:4,7,9,10),data.frame('x'=1:6)),'length')
  expect_equal(fillZeros(-c(1:4,7,9),data.frame('x'=1:6)),list('pos'=-9:-1,'data'=data.frame('x'=rev(c(1:4,0,0,5,0,6)))))
})

test_that("Test fillCover",{
  cover<-data.frame('pos'=c(2,10,11,14),'counts1'=1:4,'counts2'=2:5,'xx'=1)
  out<-data.frame('pos'=2:14,'counts1'=rep(c(1,0,2,3,0,4),c(1,7,1,1,2,1)),'counts2'=rep(c(2,0,3,4,0,5),c(1,7,1,1,2,1)),'xx'=1)
  expect_equal(fillCover(cover),out)
  expect_error(fillCover(cbind(cover,'yy'=2:5)),'nonunique')
})

test_that("Test readBlast",{
  blastString<-c(
    "read1	target1	34.56	100	2	1	101	200	400	500	1e-11	96",
    "read2	target2	100	1000	200	10	102	1100	401	2000	1e-21	123"
  )
  out<-data.frame('qName'=c('read1','read2'),'tName'=c('target1','target2'),'percID'=c(34.56,100),'alignLength'=c(100,1000),'mismatch'=c(2,200),'nGap'=c(1,10),'qStart'=101:102,'qEnd'=c(200,1100),'tStart'=400:401,'tEnd'=c(500,2000),'eValue'=c(1e-11,1e-21),'bit'=c(96,123),score=c(97,790),stringsAsFactors=FALSE)
  expect_equal(read.blast(textConnection(blastString)),out)
})

test_that("Test write.fastq",{
  out<-data.frame('name'=c('seq1','seq2'),'seq'=c('GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT','GACCGGAACT'),'qual'=c("!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65","!''*(!!!*!"),stringsAsFactors=FALSE)
  tmpFile<-tempfile()
  write.fastq(out$name,out$seq,out$qual,tmpFile)
  expect_equal(read.fastq(tmpFile),out)
  tmpGz<-tempfile(fileext='.gz')
  write.fastq(out$name,out$seq,out$qual,tmpGz)
  expect_equal(read.fastq(tmpGz),out)
  expect_lt(file.size(tmpGz),file.size(tmpFile))
  expect_error(write.fastq(out$name[-1],out$seq,out$qual,tmpFile),'length')
  expect_error(write.fastq(out$name,out$seq[-1],out$qual,tmpFile),'length')
  expect_error(write.fastq(out$name,out$seq,out$qual[-1],tmpFile),'length')
  fastq<-generateFasta(100,10:300,c('A','C','T','G'),generateQuals=TRUE)
  write.fastq(fastq$name,fastq$seq,fastq$qual,tmpFile)
  expect_equal(read.fastq(tmpFile),fastq)
})
