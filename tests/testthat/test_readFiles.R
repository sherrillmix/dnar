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
})
