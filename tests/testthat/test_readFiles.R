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
})


test_that("Test readFaDir",{
  tmpDir<-tempdir()
  files<-file.path(tmpDir,sprintf('%d.fa',1:9))
  seqs<-lapply(files,function(x)generateFakeFasta(sample(1:100,1),10:1000,c('A','C','T','G','\n')))
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

test_that("Test generateFakeFasta",{
  fasta<-generateFakeFasta(1000,90:110,c('A','Z'))
  expect_equal(nrow(fasta),1000)
  expect_true(all(grepl('^>',fasta$name)))
  expect_true(!any(grepl('[^AZ]',fasta$seq)))
  expect_true(all(nchar(fasta$seq) %in% c(90:110)))
  expect_equal(length(unique(fasta$name)),1000)
})

test_that("Test write.fa",{
  fa<-generateFakeFasta(100,10:300)
  fa$noGreater<-sub('^>','',fa$name)
  tmp<-tempfile()
  write.fa(fa$noGreater,fa$seq,tmp)
  expect_equal(read.fa(tmp)$seq,fa$seq)
  expect_equal(read.fa(tmp)$name,fa$noGreater)
  fa<-generateFakeFasta(100,10:300,c('A','C','T','G','\n'))
  fa$noGreater<-sub('^>','',fa$name)
  tmpGz<-tempfile(fileext='.gz')
  write.fa(fa$name,fa$seq,tmpGz)
  expect_equal(read.fa(tmpGz)$seq,gsub('\n','',fa$seq))
  expect_equal(read.fa(tmpGz)$name,fa$noGreater)
})
