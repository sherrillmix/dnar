% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readFiles.R
\name{readFaDir}
\alias{readFaDir}
\title{Read in a bunch of fasta files in a target directory}
\usage{
readFaDir(
  dir = ".",
  suffix = "\\\\.(fn?a|fasta)$",
  recursive = FALSE,
  vocal = FALSE,
  ...
)
}
\arguments{
\item{dir}{target directory}

\item{suffix}{regex to select file names}

\item{recursive}{if TRUE recurse through data directory}

\item{vocal}{if TRUE display status message for each file loading}

\item{...}{additional arguments to read.fa}
}
\value{
data.frame with columns name, seq and file
}
\description{
Read in a bunch of fasta files in a target directory
}
\examples{
tmpDir<-tempdir()
files<-file.path(tmpDir,sprintf('\%d.fa',1:9))
for(ii in files){
  tmp<-generateFasta(5,10:20,c('A','C','T','G'))
  write.fa(tmp$name,tmp$seq,ii)
}
seqs<-readFaDir(tmpDir)
}
