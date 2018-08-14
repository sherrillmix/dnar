#' Converts ambiguous dna to an appropriate regular expression
#'
#' @param dna input sequence with ambiguous bases
#' @export
#' @return regular expression to match the ambiguous dna
#' @references \url{https://en.wikipedia.org/wiki/Nucleic_acid_notation}
#' @examples
#' dnaRegex<-ambiguous2regex('ANRT')
#' grepl(dnaRegex,c('ATGT','ACAT','AGTT'))
ambiguous2regex<-function(dna){
  dna<-toupper(dna)
  for (i in names(dnar::ambiguousBaseCodes)){
    dna<-gsub(i,paste('[',dnar::ambiguousBaseCodes[i],']',sep=''),dna)
  }
  return(dna)
}

#' Convert set of single bases to ambiguous code
#'
#' @param bases vector of character sequences all the same number of characters
#' @export
#' @return single character string with ambiguous nucleotide codes representing any differences between sequences 
#' @references \url{https://en.wikipedia.org/wiki/Nucleic_acid_notation}
#' @examples
#' seqs<-c('ACTAGG','ACGTGG','ACCTGG','ACATGG')
#' ambig<-bases2ambiguous(seqs)
#' pattern<-ambiguous2regex(ambig)
#' grepl(pattern,seqs)
bases2ambiguous<-function(bases){
  if(any(grepl('[^ACTG]',bases)))stop(simpleError('Non ACTG characters found in bases2ambiguous'))
  bases<-sort(unique(bases))
  nBases<-nchar(bases[1])
  if(any(nchar(bases)!=nBases))stop(simpleError('Different length strings found in bases2ambiguous'))
  if(length(bases)==1)return(bases)
  if(nBases>1)return(paste(sapply(1:nBases,function(x)bases2ambiguous(substring(bases,x,x))),collapse=''))
  else return(dnar::reverseAmbiguous[paste(bases,collapse='')])
}


#' Convert ambiguous dna to all possible sequences
#'
#' @param dna vector dna containing ambiguous bases
#' @export
#' @return list with each entry containing all combinations of ambiguous bases for that entry of the dna vector
#' @references \url{https://en.wikipedia.org/wiki/Nucleic_acid_notation}
#' @examples
#' expandAmbiguous(c('AACACANAT',"ACCHB"))
expandAmbiguous<-function(dna){
  dna<-toupper(dna)
  ambigRegex<-sprintf('[%s]',paste(names(dnar::ambiguousBaseCodes),collapse=''))
  out<-lapply(dna,function(x){
    pos<-regexpr(ambigRegex,x)  
    if(pos!=-1){
      replaces<-strsplit(dnar::ambiguousBaseCodes[substring(x,pos,pos)],'')[[1]]
      new<-rep(x,length(replaces))
      for(i in 1:length(replaces))substring(new[i],pos,pos)<-replaces[i]
      out<-unlist(expandAmbiguous(new))
    }else{out<-x}
    return(out)
  })
  return(out)
}


#' Count nMers in a string
#'
#' @param string string to be searched for nmers
#' @param k length of nmers
#' @export
#' @return a sorted named vector giving the identity and counts for nMers
#' @examples
#' countNmers('AATTAATT',2)
#' countNmers('AATTAATT',3)
countNmers<-function(string,k=10){
  if(nchar(string)<k){
    warning('string shorter than k')
    return(0)
  }
  indices<-1:(nchar(string)-k+1)
  substrings<-substring(string,indices,indices+k-1)
  return(table2vector(sort(table(substrings))))
}

#' Calculate GC percentage for strings
#'
#' @param seqs sequences to get GC content of
#' @param chars characters to count (G,C by default)
#' @export
#' @return proportion of GC in sequence
#' @examples
#' gcPercent(c('AAATA','AATTC','GCGCTT'))
#' gcPercent(c('12345','22234','01011101'),'1')
gcPercent<-function(seqs,chars=c('C','G')){
  regex<-sprintf('[%s]+',paste(chars,collapse=''))
  gcs<-nchar(gsub(regex,'',seqs,perl=TRUE))
  return(1-gcs/nchar(seqs))
}


#' Swap case of a string
#'
#' @param string String to have case toggled
#' @export
#' @return String with lowercase and uppercase swapped
#' @examples
#' toggleCase(c('AGTCcccT','AbCdEfG'))
toggleCase<-function(string){
  chartr(paste(letters,LETTERS,collapse='',sep=''),paste(LETTERS,letters,collapse='',sep=''),string)
}

#' Find pattern in strings and flip case
#'
#' Use the case of a string to highlight a pattern. For example, highlighting AC in TTTTGTACAAATC as TTTTGTacAAATC. Note that it will not find overlapping patterns e.g. looking for AGA in AGAGA will only find the first AGA. Non-letter characters can be included in the pattern but their case will remain unchanged.
#'
#' @param pattern Pattern to look for in strings
#' @param strings Strings in which to look for pattern
#' @export
#' @return Strings with case of any occurences of pattern swapped
#' @examples
#' highlightString('AGA',c('TTAGATTAGA','AGAGA'))
#' highlightString('aga',tolower(c('TTAGATTAGA','AGAGA')))
highlightString<-function(pattern,strings){
  locs<-gregexpr(pattern,strings)
  nLetter<-nchar(pattern)
  strings<-mapply(function(x,y){y<-y[y!=-1];for(i in y)substring(x,i,i+nLetter-1)<-toggleCase(substring(x,i,i+nLetter-1));return(x)},strings,locs,USE.NAMES=FALSE)
  return(strings)
}

#' Switch the case of positions in seq1 that mismatch seq2
#'
#' @param seq1 sequence to highlight differences in 
#' @param seq2 sequence to compare to
#' @export
#' @return seq1 with bases that differ from seq2 having their case toggled (e.g. upper to lower case)
#' highlightDifferences('aaccat','aacgat')
#' highlightDifferences('TATATTCTTT','TATATCCTTT')
highlightDifferences<-function(seq1,seq2){
  if(nchar(seq1)!=nchar(seq2))stop(simpleError('Highlighting differences in different length sequences not supported'))
  seqMat<-seqSplit(seq1,seq2)
  diffs<-which(apply(seqMat,2,function(x)x[1]!=x[2]))
  for(ii in diffs)substring(seq1,ii,ii)<-toggleCase(substring(seq1,ii,ii))
  return(seq1)
}

#' Convenience function for splitting a bunch of sequences into a matrix
#'
#' @param ... Various sequences to split into a matrix
#' @param fill A character to pad ends of sequences. An error is generated if any sequence differs in length if NULL
#' @export
#' @return A matrix of single characters each row corresponding to a read
#' @examples
#' seqSplit('ACACA','ACA')
#' seqSplit('ACACA','ACA',fill='-')
#' seqSplit(c('ACACA','ACA'),'TCACA')
seqSplit<-function(...,fill='.'){
  seqs<-c(...)
  seqN<-nchar(seqs)
  maxN<-max(seqN)
  if(is.null(fill)&&any(seqN!=maxN))stop(simpleError('All sequences not same length'))
  seqs<-paste(seqs,sapply(maxN-seqN,function(x)paste(rep(fill,x),collapse='')),sep='')
  return(do.call(rbind,strsplit(seqs,'')))
}

#' Convert DNA string into seperate codons
#'
#' @param dna Vector of DNA strings
#' @param frame Starting frame for codons (0=start on first base, 1=on second, 2=on third)
#' @export
#' @return List with a vector of 3 base codons for each input DNA string
#' @examples
#' dna2codons('ACATTTGGG')
#' dna2codons('ACATTTGGG',1)
#' dna2codons(c('ACATTTGGG','AAATTAGGC'))
#' dna2codons(c('ACATTTGGG','AAATTAGGC'),c(1,2))
dna2codons<-function(dna,frame=0){
  frame<-frame+1
  out<-mapply(function(seq,start){
    if(nchar(seq)-start+1<3){
      warning('DNA shorter than 3 bases in dna2codons')
      return(NULL)
    }
    starts<-seq(start,nchar(seq)-2,3)
    return(substring(seq,starts,starts+2))
  },dna,frame,SIMPLIFY=FALSE,USE.NAMES=FALSE)
  return(out)
}

#' Convert codon to amino acid
#'
#' @param codons Vector of 3 base codons
#' @param type Amino acid info to return; code for single letter, name for full name, or abbr for 3-letter abbreviation
#' @param naReplace Replace NAs with this string. Anything other than a single character will result in a protein longer than the number of codons
#' @param warn Warn if any unknown codons are found. Usually due to ACTG characters
#' @export
#' @return Vector of amino acids
#' @examples
#' codons2aa(c('AUG','GTT','AGG','GTA'))
#' codons2aa(dna2codons('ATGTTTATTGGG')[[1]])
codons2aa<-function(codons,type='code',naReplace='z',warn=TRUE){
  if(!type %in% c('code','name','abbr'))stop(simpleError('Invalid amino acid type'))
  codons<-gsub('T','U',toupper(codons))
  aas<-dnar::aminoAcids[,type,drop=FALSE][codons,1]
  if(warn&any(is.na(aas)))warning('Unknown codons in codons2aa')
  aas[is.na(aas)]<-naReplace
  return(aas)
}

#' Convert dna/rna to amino acids
#'
#' @param dna A string of DNA/RNA
#' @param frame Starting frame (0=start on first base, 1=on second, 2=on third)
#' @param ... Additional arguments to codons2aa
#' @export
#' @return A string of amino acids
#' @examples
#' dna2aa(c('ATGAGATGCAGTTAA','AATTTA'))
dna2aa<-function(dna,frame=0,...){
  codons<-dna2codons(dna,frame)  
  output<-sapply(codons,function(xx)paste(codons2aa(xx,...),collapse=''))
  return(output)
}

#' Find the DNA to code for a given amino acid
#'
#' @param aas single amino acids to convert to dna
#' @param type code or abbr or name
#' @param regex if true return a (X|Y) regex of codons else return vector
#' @export
#' @return A list of vectors of codons if regex=FALSE else a regular expression matching the possible codons
#' @examples
#' aa2codons(c('A','T','X'))
#' aa2codons(c('A','T','X'),regex=FALSE)
#' aa2codons(c('Glu','Ala','Xxx'),type='abbr')
aa2codons<-function(aas,type='code',regex=TRUE){
  selectors<-lapply(aas,function(aa)dnar::aminoAcids[,type]==aa)
  aaFound<-sapply(selectors,any)
  if(any(!aaFound))stopError('Unknown amino acid ',paste(aas[!aaFound],collapse=', '))
  codons<-lapply(selectors,function(selector)dnar::aminoAcids[selector,'codon'])
  if(regex)codons<-sapply(codons,function(codon)sprintf('(%s)',paste(codon,collapse='|')))
  return(codons)
}

#' Convert amino acids to dna regex
#'
#' @param aas A vector of amino acids to turn to dna regex
#' @export
#' @return A vector with a regular expression for each element of aas
#' @examples
#' aa2dna(c('MATX','MGCATKRVX'))
aa2dna<-function(aas){
  splitAas<-strsplit(aas,'')
  out<-sapply(splitAas,function(x)paste(aa2codons(x),collapse=''))
  out<-gsub('U','T',out)
  return(out)
}

# Check overlap between two sets of coordinates
#
# @param chrs Chromosomes of query
# @param starts Start coordinates of query
# @param ends End coordinates of query
# @param tChrs Chromsomes of target
# @param tStarts Start coordinates of target
# @param tEnds End coordinates of target
# @param tNames Target names 
# @param sep Matches are pasted together separated by this
# @export
# @return '|' separated vector of tNames within overlap or '' if no overlapping target
# @examples
# queries<-data.frame(
#   'chr'=letters[rep(1:3,each=2)],
#   'start'=c(1,1000,-1,1000,-10000,1),
#   'end'=c(1010,2000,50,1050,10000,1)
# )
# targets<-data.frame('chr'=letters[rep(1:2,each=2)],
#   'start'=c(1,1111,1,1000),
#   'end'=c(1000,9999,2,1000),
#   'name'=letters[1:4]
# )
# checkOverlap(
#   queries$chr,queries$start,queries$end,
#   targets$chr,targets$start,targets$end,targets$name
# )
#checkOverlap<-function(chrs,starts,ends,tChrs,tStarts,tEnds,tNames,sep='|'){
#  queries<-GenomicRanges::GRanges(seqnames=chrs,ranges=IRanges::IRanges(starts,end=ends))
#  targets<-GenomicRanges::GRanges(seqnames=tChrs,ranges=IRanges::IRanges(tStarts,end=tEnds))
#  overlaps<-GenomicRanges::findOverlaps(queries,targets)
#  pasted<-tapply(tNames[overlaps@to],overlaps@from,paste,collapse=sep)
#  overlapNames<-pasted[as.character(1:length(chrs))]
#  overlapNames[is.na(overlapNames)]<-''
#  overlapNames<-as.vector(unname(overlapNames))
#  return(overlapNames)
#}

#' Convert gapped coordinates to what the coordinates would be without gaps
#'
#' @param gapSeq the reference sequence containing gaps
#' @param coords coordinates on the gapped gapSeq to be converted into equivalent nongap coordinatess
#' @export
#' @return equivalent gapped coordinates
#' @examples
#' gap2NoGap('AA-TT-GG',7)
#' gap2NoGap('AAA----TT-GG',1:10)
gap2NoGap<-function(gapSeq,coords){
  gapSeqSplit<-strsplit(gapSeq,'')[[1]]
  nonDash<-!gapSeqSplit %in% c('*','.','-')
  newCoords<-cumsum(nonDash)
  coords[coords<1|coords>length(newCoords)]<-NA
  return(newCoords[coords])
}

#' Convert ungapped coordinates to what the coordinates would be with gaps
#'
#' @param gapSeq the reference sequence containing gaps
#' @param coords coordinates on the ungapped gapSeq to be converted into equivalent gap coordinates
#' @export
#' @return equivalent ungapped coordinates
#' @examples
#' noGap2Gap('AA-TT-GG',5)
#' noGap2Gap('AAA----TT-GG',1:7)
noGap2Gap<-function(gapSeq,coords){
  gapSeqSplit<-strsplit(gapSeq,'')[[1]]
  nonDash<-which(!gapSeqSplit %in% c('.','*','-'))
  coords[coords<1|coords>length(nonDash)]<-NA
  return(nonDash[coords])
}

#' Find covered ranges in binary data
#'
#' @param index logical vector containing TRUE in the regions of interest
#' @export
#' @return data.frame with rows for each contiguous region of interest with columns start and end of regions
#' @examples
#' binary2range(c(FALSE,TRUE,FALSE,TRUE,TRUE,TRUE))
binary2range<-function(index){
  return(index2range(which(index)))
}

#' Find covered ranges with numerical index data
#'
#' @param index numeric indices indicating positions of interest
#' @export
#' @return data.frame with rows for each contiguous region of interest with columns start and end of regions
#' @examples
#' index2range(c(1:100,300:450))
#' index2range(c(1:100,2:110))
index2range<-function(index){
  if(length(index)==0)return(data.frame('start'=NA,'end'=NA)[0,])
  index<-sort(unique(index))
  diffs<-c(diff(index),1)
  ends<-c(which(diffs>1),length(index))
  starts<-c(1,ends[-length(ends)]+1)
  return(data.frame('start'=index[starts],'end'=index[ends]))
}

#' Reverse strings
#'
#' @param strings vector of strings to be reversed
#' @param brackets if TRUE then reverse brackets and parenthesis, [ goes to ]. keeps brackets in order after reversing
#' @export
#' @return vector with the reversed strings
#' @examples
#' reverseString(c("ABCDEFG","Reverse"))
reverseString<-function(strings,brackets=TRUE){
  output<-sapply(strings,function(x)intToUtf8(rev(utf8ToInt(x))),USE.NAMES=FALSE) #http://stackoverflow.com/questions/13612967/how-to-reverse-a-string-in-r
  #slower
  #output<-sapply(strsplit(strings,''),function(x)paste(rev(x),collapse=''))  
  if(brackets)output<-chartr('[]()','][)(',output)
  return(output)
}

#' Complement DNA 
#'
#' @param dnas vector of sequences
#' @param ambigs if TRUE complement ambiguous bases
#' @export
#' @return vector with the DNA sequences complemented
#' @examples
#' complementDna(c('CTAG','ATCCAC'))
#' complementDna(c('CT[AC]G','ATNRY'))
complementDna<-function(dnas,ambigs=TRUE){
  finds<-'TGAC'
  replaces<-'ACTG'
  #deal with ambiguous
  if(ambigs){
    sortAmbig<-sapply(lapply(strsplit(dnar::ambiguousBaseCodes,''),sort),paste,collapse='',USE.NAMES=FALSE)
    revAmbig<-sapply(strsplit(complementDna(sortAmbig,ambigs=FALSE),''),function(x)paste(sort(x),collapse=''),USE.NAMES=FALSE)
    ambigComp<-names(sortAmbig)[sapply(revAmbig,function(x)which(x==sortAmbig))]
    finds<-sprintf('%s%s',finds,paste(names(sortAmbig),collapse=''))
    replaces<-sprintf('%s%s',replaces,paste(ambigComp,collapse=''))
  }
  return(chartr(finds,replaces,dnas))
}
#' Reverse complement dna
#'
#' @param dnas vector of sequences
#' @export
#' @return vector of reverse complemented dna sequences
#' @examples
#' revComp(c('CTAG','ATCCAC'))
#' revComp(c('CT[AC]G','ATNRY'))
revComp<-function(dnas){
  return(complementDna(reverseString(dnas),TRUE))
}

#' Remove gap characters from DNA sequences
#'
#' @param seq vector of gapped DNA sequences
#' @param gaps vector of characters to be considered gaps
#' @export
#' @return vector of DNA sequences with gaps removed
#' @examples
#' degap(c('...ACTATATA----ACATG--G..','ATTAT--T'))
#' degap(c('1000111011'),'0')
degap<-function(seq,gaps=c('*','-','.')){
  escapedGaps<-escapeRegexBracketChars(gaps)
  gsub(sprintf('[%s]+',paste(escapedGaps,collapse='')),'',seq,perl=TRUE)
}

#' Parse a region string into chr, start, end and strand
#'
#' @param reg vector of region strings in the format "chrX:123545-123324" or "chr12:1234-1236+"
#' @export
#' @return data.frame with one row per reg string and columns chr, start, end and strand
#' @examples
#' parseRegion(c('chr1:12345-23456','chrX:2222:3333'))
#' parseRegion(c('chr1:12345-23456+','chrX:2222:3333-','chrX:2222:3333*'))
parseRegion<-function(reg){
  strand<-ifelse(substring(reg,nchar(reg)) %in% c('*','-','+'),substring(reg,nchar(reg)),'*')
  reg<-sub('[-+*]$','',reg)
  splits<-strsplit(reg,'[:-]')
  if(any(sapply(splits,length)!=3))stop(simpleError(sprintf('Region %s not parsed',paste(reg[sapply(splits,length)!=3],collapse=', '))))
  out<-data.frame('chr'=sapply(splits,'[[',1),stringsAsFactors=FALSE)
  out[,c('start','end')]<-as.numeric(do.call(rbind,lapply(splits,'[',2:3)))
  out$strand<-strand
  return(out)
}

#' Make region from chr start end
#'
#' @param chrs vector of chromosomes
#' @param starts vector of starts
#' @param ends vector of ends
#' @param strands optional vector of strands
#' @export
#' @return vector of region strings
#' @examples
#' pasteRegion(c('chr1','chr2'),c(100,1235),c(200,2346))
#' pasteRegion(c('chr1','chr2'),c(100,1235),c(200,2346),c('*','-'))
pasteRegion<-function(chrs,starts,ends,strands=''){
  sprintf('%s:%s-%s%s',chrs,sub('^[ \t\r\n]+','',format(starts,scientific=FALSE)),sub('^[ \t\r\n]+','',format(ends,scientific=FALSE)),strands)
}

#' Test sam flags for values
#'
#' @param flags vector of integer flags from sam
#' @param test either character vector of flag short names or integers
#' @export
#' @return logical vector with element for each element in flags with TRUE if flag meets all tests and FALSE otherwise
#' @examples
#' print(samFlags)
#' samFlag(1:10,'paired')
#' samFlag(1:10,'unmapped')
samFlag<-function(flags,test='paired'){
  test<-unique(test)
  if(!is.numeric(test)){
    if(!all(test %in% dnar::samFlags$short))stop(simpleError(sprintf('Unknown flag please select from %s',paste(dnar::samFlags$short,collapse=', '))))
    test<-dnar::samFlags[dnar::samFlags$short %in% test,'bit']
  }else{
    test<-as.integer(test)
  }
  testInt<-0
  for(i in test)testInt<-bitops::bitOr(testInt,i)
  return(bitops::bitAnd(flags,testInt)==testInt)
}


#' Trim ambiguous sequence from ends of reads
#'
#' Trim sequences to the first occurrence of a given length of unambiguous bases, e.g. to trim noisy starts and ends from Sanger sequences. Note that sequences shorter than nonNStretch will be completely deleted since they do not contain any stretch of unambiguous bases to pass the threshold.
#'
#' @param seqs sequences to be trimmed
#' @param nonNStretch a one or two element (specifying start, end) to delete until finding a stretch of nonNs greater this length (0 means no trimming)
#' @param nChars a vector of characters to count as Ns
#' @export
#' @return Vector sequences with ends trimmed
#' @examples
#' trimNs('NNNANNACTNGAGANNNNAAN',3)
trimNs<-function(seqs,nonNStretch=c(10,10),nChars=c("N")){
  if(length(nonNStretch)==1)nonNStretch<-rep(nonNStretch,2)
  regex<-sprintf('[^%s]{%d,}',paste(escapeRegexBracketChars(nChars),collapse=''),nonNStretch)
  if(nonNStretch[1]>0){
    starts<-regexpr(regex[1],seqs,perl=TRUE)
    seqs[starts==-1]<-''
    seqs<-substring(seqs,starts)
  }
  if(nonNStretch[2]>0){
    regexResults<-gregexpr(regex[2],seqs,perl=TRUE)
    starts<-sapply(regexResults,utils::tail,1)
    lengths<-sapply(regexResults,function(x)utils::tail(attr(x,'match.length'),1))
    seqs[starts==-1]<-''
    seqs<-substring(seqs,1,starts+lengths-1)
  }
  return(seqs)
}

#' Convert cigar and starts to qStarts, tStarts, blockSizes as in blat
#'
#' @param cigars vector of SAM cigar strings
#' @param starts vector of starting positions in target
#' @param startEnds if FALSE single line with comma separated starts, if TRUE data.frame with single start, end and id column
#' @export
#' @return dataframe with qStarts,tStarts,sizes if !startEnds or dataframe with starts, ends and ids if startEnds
#' @examples
#' cigars<-c('100M','10H10M1000H','10M10I20M10D2M','10S10M10H','100M10000N100M')
#' starts<-1:5
#' cigarToBlock(cigars,starts)
#' cigarToBlock(cigars,starts,startEnds=TRUE)
cigarToBlock<-function(cigars,starts,startEnds=FALSE){
  #M=match, I=insertion in query, D=deletion in query, N="intron" deletion in query, S=soft clipping (clip sequence), H=hard clipping (sequence was already clipped)
  nAligns<-length(starts)
  if(nAligns!=length(cigars))stop(simpleError('Cigars and starts not same length'))
  supportedOps<-c('M','I','D','N','S','H')
  if(any(grep(sprintf('[^0-9%s]',paste(supportedOps,collapse='')),cigars)))stop(simpleError(sprintf('Unknown operation. Only %s cigar operations supported',paste(supportedOps,collapse=''))))
  tPos<-starts
  qPos<-rep(1,nAligns)
  stillWorking<-rep(TRUE,nAligns)
  qStarts<-tStarts<-rep('',nAligns)
  blockSizes<-tStarts<-rep('',nAligns)
  if(startEnds){
    startEndsOut<-data.frame('start'=-1,'end'=-1,'id'=-1,'qStart'=-1,'qEnd'=-1)[0,]
    ids<-1:length(starts)
  }
  N_DUMMY_GAPS<-1000000
  dummyGaps<-paste(rep('-',N_DUMMY_GAPS),collapse='')
  while(any(stillWorking)){
    for(i in supportedOps){
      regex<-sprintf('^([0-9]+)%s',i)
      matches<-grep(regex,cigars[stillWorking])
      if(!any(matches))next()
      num<-as.numeric(sub(sprintf('%s.*',regex),'\\1',cigars[stillWorking][matches]))
      cigars[stillWorking][matches]<-sub(regex,'',cigars[stillWorking][matches])
      if(i %in% c('M')){
        if(any(is.na(num)))browser()
        tStarts[stillWorking][matches]<-sprintf('%s,%d',tStarts[stillWorking][matches],tPos[stillWorking][matches])
        qStarts[stillWorking][matches]<-sprintf('%s,%d',qStarts[stillWorking][matches],qPos[stillWorking][matches])
        blockSizes[stillWorking][matches]<-sprintf('%s,%d',blockSizes[stillWorking][matches],num)
        if(startEnds){
          startEndsOut<-rbind(startEndsOut,data.frame('start'=tPos[stillWorking][matches],'end'=tPos[stillWorking][matches]+num-1,'id'=ids[stillWorking][matches],'qStart'=qPos[stillWorking][matches],'qEnd'=qPos[stillWorking][matches]+num-1))
        }
      }
      if(i %in% c('M','D','N')){
        tPos[stillWorking][matches]<-num+tPos[stillWorking][matches]
      }
      if(i %in% c('M','I','S')){
        qPos[stillWorking][matches]<-num+qPos[stillWorking][matches]
      }
      stillWorking[stillWorking]<-nchar(cigars[stillWorking])>0
    }
  }
  qStarts<-sub('^,','',qStarts)
  tStarts<-sub('^,','',tStarts)
  blockSizes<-sub('^,','',blockSizes)
  if(startEnds){
    startEndsOut<-startEndsOut[order(startEndsOut$id,startEndsOut$start),]
    rownames(startEndsOut)<-NULL
    return(startEndsOut)
  }
  else return(data.frame('qStarts'=qStarts,'tStarts'=tStarts,'sizes'=blockSizes,stringsAsFactors=FALSE))
}

#' Convert blocks from blat into alignment
#'
#' @param seqs vector of sequences
#' @param tSeqs vector of target sequences or a single target
#' @param qStarts comma separated starts of query matches e.g. from blat or cigarToBlock (1 based)
#' @param tStarts comma separated starts of target matches e.g. from blat or cigarToBlock (1 based)
#' @param sizes comma separated lengths of matches e.g. from blat or cigarToBlock
#' @export
#' @return data.frame with sequences aligned in columns qSeq and tSeq
#' @examples
#' blockToAlign('ABEC','BCD','2,4','1,2','1,1')
#' blockToAlign(c('ABC','CDE'),c('BCD','ABCDE'),c('1,3','1'),c('1,2','3'),c('1,1','3'))
#' blocks<-cigarToBlock(c('3M','2M1I1M'),1:2)
#' blockToAlign(
#'   c('ATC','ATGC'),
#'   c('ATCT','GATCT'),
#'   blocks$qStarts,blocks$tStarts,blocks$sizes
#' )
blockToAlign<-function(seqs,tSeqs,qStarts,tStarts,sizes){
  nSeqs<-length(seqs)
  if(length(tSeqs)==1)tSeqs<-rep(tSeqs,length(seqs))
  if(length(tSeqs)!=nSeqs)stop(simpleError('Target seqs not same length as seqs and not a single sequence'))
  if(length(qStarts)!=nSeqs||length(tStarts)!=nSeqs||length(sizes)!=nSeqs)stop(simpleError('All arguments to blockToAlign not same length'))
  tPieces<-blat2exons(1:length(seqs),1:length(seqs),tStarts,sizes)
  tPieces$seq<-substring(tSeqs[tPieces$name],tPieces$start,tPieces$end)
  qPieces<-blat2exons(1:length(seqs),1:length(seqs),qStarts,sizes)
  qPieces$seq<-substring(seqs[qPieces$name],qPieces$start,qPieces$end)
  if(any(tPieces$exonName!=qPieces$exonName))stop(simpleError('Problem constructing fragments'))
  gaps<-blatFindGaps(qStarts,tStarts,sizes)
  gaps<-mapply(function(gap,tSeq,qSeq){
    gap<-as.data.frame(gap,stringsAsFactors=FALSE)
    #depending on substring(seq,x,x-1) returning ''
    if(nrow(gap)==0)return(gap)
    gap$tSeq<-substring(tSeq,gap$tStartAfter+1,gap$tStartAfter+gap$tGaps)
    gap$qSeq<-substring(qSeq,gap$qStartAfter+1,gap$qStartAfter+gap$qGaps)
    dummy<-paste(rep('-',max(nchar(c(gap$tSeq,gap$qSeq)))),collapse='')
    gap[,c('tSeq','qSeq')]<-t(apply(gap[,c('tSeq','qSeq')],1,function(x)paste(x,substring(dummy,1,max(nchar(x))-nchar(x)),sep='')))
    return(gap)
  },gaps,tSeqs,seqs,SIMPLIFY=FALSE)
  qSeqs<-sapply(1:length(seqs),function(ii){
    x<-qPieces[qPieces$name==ii,]
    thisGaps<-gaps[[ii]]
    return(paste(c(x$seq,thisGaps$qSeq)[order(c(x$start,thisGaps$qStartAfter+.5))],collapse=''))
  })
  tSeqs<-sapply(1:length(seqs),function(ii){
    x<-tPieces[tPieces$name==ii,]
    thisGaps<-gaps[[ii]]
    return(paste(c(x$seq,thisGaps$tSeq)[order(c(x$start,thisGaps$tStartAfter+.5))],collapse=''))
  })
  return(data.frame('qSeq'=qSeqs,'tSeq'=tSeqs,stringsAsFactors=FALSE))
}

#' Take the output from blat and make continous reads out of it
#'
#' @param seqs comma separated target sequences from blat
#' @param starts comma separated starting location for each piece of sequence 0 indexed
#' @param fillStarts start base for each output sequence 0 indexed
#' @param fillEnds total length for each output sequence
#' @param gaps list of matrices with columns tGaps and qGaps
#' @export
#' @return vector of sequences
cutBlatReads<-function(seqs,starts,fillStarts=NULL,fillEnds=NULL,gaps=NULL){
  if(length(seqs)!=length(starts))stop(simpleError('Length of seqs and starts not equal'))
  if(!is.null(fillEnds)){
    if(length(seqs)!=length(fillEnds))stop(simpleError('Length of seqs and fillEnds not equal'))
    if(any(fillEnds<1))stop(simpleError('Total length must be positive'))
  }
  if(!is.null(gaps)&&length(gaps)!=length(starts))stop(simpleError('Gaps not the same length as seqs'))
  if(!is.null(fillStarts)&&length(seqs)!=length(fillStarts))stop(simpleError('Length of seqs and fillStarts not equal'))
  lengthMax<-1000000
  dots<-rep('.',lengthMax)
  seqList<-strsplit(seqs,',')
  startList<-lapply(strsplit(starts,','),as.numeric)
  if(any(sapply(seqList,length)!=sapply(startList,length)))stop(simpleError('Number of sequences and starts do not match'))
  if(is.null(fillEnds)){
    fillEnds<-mapply(function(x,y)x[length(x)]+nchar(y[length(y)])-1,startList,seqList)
  }
  if(is.null(fillStarts)){
    fillStarts<-sapply(startList,min)
  }
  output<-mapply(function(seqs,starts,fillStart,fillEnd){
    if(any(starts<fillStart))browser()#stop(simpleError('fillStart less than starts'))
    numChars<-nchar(seqs)
    starts<-starts-fillStart
    ends<-starts+numChars-1
    if(any(ends>fillEnd))stop(simpleError('Sequence goes past fillEnd'))
    goodPos<-unlist(mapply(function(x,y)x:y,starts,ends))+1
    if(any(table(goodPos)>1))stop(simpleError('Base covered more than once'))
    seqLength<-diff(c(fillStart,fillEnd))
    if(seqLength>lengthMax){
      warning('Sequence longer than ',lengthMax,' bases. Returning NA')
      return(NA)
    }
    seqSplit<-dots[1:seqLength]
    seqSplit[goodPos]<-unlist(strsplit(seqs,''))
    return(paste(seqSplit,collapse=''))
  },seqList,startList,fillStarts,fillEnds)
  return(output)
}

#' Find gaps in blat coordinates
#'
#' @param qStarts comma separated starts for query seq
#' @param tStarts comma separated starts for target seq
#' @param blockSizes comma separated block sizes
#' @export
#' @return data.frame giving qStartAfter, qGaps, tGapStartAfter and tGaps
#' @examples
#' blatFindGaps('1,100,200','4,200,300','99,50,100')
blatFindGaps<-function(qStarts,tStarts,blockSizes){
  if(length(qStarts)!=length(tStarts)||length(blockSizes)!=length(qStarts))stop(simpleError('Lengths of starts and blocksizes not equal'))
  qStartList<-lapply(strsplit(qStarts,','),as.numeric)
  tStartList<-lapply(strsplit(tStarts,','),as.numeric)
  blockSizeList<-lapply(strsplit(blockSizes,','),as.numeric)
  blockNums<-sapply(qStartList,length)
  if(any(blockNums!=sapply(tStartList,length))||any(blockNums!=sapply(blockSizeList,length)))stop(simpleError('Comma separated lists of starts and blocksizes not equal length'))

  gaps<-mapply(function(qStarts,tStarts,blockSizes){
    tEnds<-tStarts+blockSizes-1
    qEnds<-qStarts+blockSizes-1
    qGaps<-findIntrons(list(qStarts),list(qEnds))
    tGaps<-findIntrons(list(tStarts),list(tEnds))
    if(nrow(qGaps)!=nrow(tGaps))stop(simpleError('Mismatch in target and query introns'))
    out<-data.frame('qStartAfter'=qGaps$startAfter,'qGaps'=qGaps$length,'tStartAfter'=tGaps$startAfter,'tGaps'=tGaps$length)
    return(out)
  },qStartList,tStartList,blockSizeList,SIMPLIFY=FALSE)

  return(gaps)
}


#' Take coordinates of blat matches and split into a single line for each exon
#'
#' @param chroms Chromosome or other identifier
#' @param names Name of gene of other container of exons
#' @param starts Start coordinates of exons
#' @param ends End coordinates of exons if lengths is FALSE or length of exon if lengths is TRUE
#' @param strands Strand (for numbering exons in reverse on - strand)
#' @param lengths logical whether ends are end coordinates or lengths
#' @param extraCols a dataframe of extra columns (1 per batch of starts) to be added to the output
#' @param extraSplits a dataframe of extra comma-separated values (1 string of comma separated values per batch of starts, 1 value per start-stop pair) to be added to the output
#' @param prefix prefix to be added to exon names
#' @param adjustStart add adjustStart to starts (good for 0 index start, 1 index ends of UCSC)
#' @export
#' @return data.frame with a row for each exon or piece of alignment and columns chrom, name, exonName, start, end and strand
#' @examples
#' blat2exons(
#'   c('chr1','chr1','chr2'),
#'   c('read1','read2','read1'),
#'   c('1,100','50,150,300','1000,2000,3000'),
#'   c('20,200','50,50,100','999,999,1234'),
#'   c('+','-','+')
#' )
blat2exons<-function(chroms,names,starts,ends,strands=rep('+',length(names)),lengths=TRUE,extraCols=NULL,extraSplits=NULL,prefix='ex',adjustStart=0){
  if(any(c(length(chroms),length(names),length(strands),length(starts))!=length(ends)))stop(simpleError('Different lengths for chrom, strand, starts, lengths'))
  startsList<-strsplit(starts,',')
  if(adjustStart!=0)startsList<-lapply(startsList,function(x)as.numeric(x)+adjustStart)
  exonCounts<-exonCountsStrand<-sapply(startsList,length)
  endsList<-strsplit(ends,',')
  if(any(exonCounts!=sapply(endsList,length)))stop(simpleError("Split starts and ends not equal length"))
  if(lengths)endsList<-mapply(function(x,y)as.numeric(x)+as.numeric(y)-1,endsList,startsList,SIMPLIFY=FALSE)
  #make sure extraCols and extraSplits are dataframes
  if(is.vector(extraSplits))extraSplits<-data.frame(extraSplits,stringsAsFactors=FALSE)
  if(is.vector(extraCols))extraCols<-data.frame(extraCols,stringsAsFactors=FALSE)
  endsUnlist<-as.numeric(unlist(endsList))
  startsUnlist<-as.numeric(unlist(startsList))
  exonCountsStrand[strands=='-']<- -exonCountsStrand[strands=='-']
  exonNums<-unlist(lapply(exonCountsStrand,function(x){if(x<0){-x:1}else{1:x}}))
  output<-data.frame('chrom'=rep(chroms,exonCounts),'name'=rep(names,exonCounts),'exonName'=sprintf('%s_%s%s',rep(names,exonCounts),prefix,exonNums),'start'=startsUnlist,'end'=endsUnlist,'strand'=rep(strands,exonCounts),stringsAsFactors=FALSE)
  if(!is.null(colnames(extraCols))){
    for(i in colnames(extraCols)){
      output[,i]<-rep(extraCols[,i],exonCounts)
    }
  }
  if(!is.null(colnames(extraSplits))){
    for(i in colnames(extraSplits)){
      thisData<-unlist(strsplit(extraSplits[,i],','))
      if(length(thisData)!=nrow(output))message('Problem splitting extraSplits ',i)
      output[,i]<-thisData
    }
  }
  return(output)
}

#' Convert set of starts and ends of exons to introns
#'
#' @param starts starts for exons given as a comma separated string for each "gene" or list of numeric vectors
#' @param ends ends for exons given as a comma separated string for each "gene" or list of numeric vectors
#' @param names vector of each "gene". Output data.frame will contain a column name with name_in + intron number
#' @param additionalColumns a data.frame with a row for each starts and columns giving additional information
#' @export
#' @return A data.frame with columns startAfter giving the base before the intron, length giving the length of intron and name of intron where name is the input names + "_in" + number. Note that startAfter and length are used to allow zero length introns
#' @examples
#' findIntrons(c('1,100','50,150,8000'),c('50,150','148,185,9000'))
findIntrons<-function(starts,ends,names=1:length(startList),additionalColumns=NULL){
  if(length(starts)!=length(ends))stop(simpleError('Starts and ends not same length'))
  if(!is.null(additionalColumns)&&length(ends)!=nrow(additionalColumns))stop(simpleError('additionalColumns not same length as starts and ends'))
  if(is.character(starts)){
    startList<-lapply(strsplit(starts,','),as.numeric)
    endList<-lapply(strsplit(ends,','),as.numeric)
  }else if(is.numeric(starts)){
    startList<-list(starts)
    endList<-list(ends)
  }else{
    startList<-starts
    endList<-ends
  }
  inEnds<-lapply(startList,function(x)x[-1]-1)
  inStarts<-lapply(endList,function(x)x[-length(x)]+1)
  if(any(sapply(inEnds,length)!=sapply(inStarts,length)))stop(simpleError('Intron ends and starts not equal length'))
  problemIndex<-mapply(function(x,y)which(x<y),inEnds,inStarts,SIMPLIFY=FALSE)
  problems<-which(sapply(problemIndex,length)>0)
  for(problem in problems){
    for(thisIndex in problemIndex[[problem]]){
      #check if these are actually nogap neighbor exons
      if(inStarts[[problem]][thisIndex]-1!=inEnds[[problem]][thisIndex]){
        stop(simpleError('Introns have negative length'))
      }
    }
    #if these weren't nogap neighbor exons we would have errored out so we can delete them
    #inStarts[[problem]]<-inStarts[[problem]][-problemIndex[[problem]]]
    #inEnds[[problem]]<-inEnds[[problem]][-problemIndex[[problem]]]
  }
  inCount<-sapply(inEnds,length)
  if(any(inCount!=sapply(inStarts,length)))stop(simpleError('Intron ends and starts not equal length after removing neighbors'))
  introns<-data.frame('startAfter'=unlist(inStarts)-1,'length'=unlist(inEnds)-unlist(inStarts)+1,stringsAsFactors=FALSE)
  introns$name<-paste(rep(names,inCount),unlist(lapply(inCount,function(x)if(x==0)return(NULL) else return(1:x))),sep='_in')
  if(!is.null(additionalColumns))introns<-cbind(introns,additionalColumns[rep(1:nrow(additionalColumns),inCount),])
  return(introns)
}

#' Generate a PWM for a set of sequences
#'
#' @param seqs vector of strings, sequences to form a position weight matrix from (all same length)
#' @param chars allowed characters
#' @param priors additional counts to add to each column with (default is a lazy way to get a named 0 vector)
#' @export
#' @return position weight matrix with length(chars) rows and nchar(seqs) columns
#' @examples
#' pwm(c('ACTG','ACTT','ACTT','ACGT'))
#' pwm(c('ACTG','ACTT','ACTT','ACGT'),priors=c('C'=1,'G'=1,'T'=1,'A'=1))
pwm<-function(seqs,chars=c('C','G','T','A'),priors=table(chars)-1){
  nBases<-nchar(seqs[1])
  if(any(nchar(seqs)!=nBases))stop(simpleError('All sequences not same length for PWM'))
  if(any(sort(names(priors))!=sort(chars)))stop(simpleError('Priors and possible characters do not match'))
  seqMat<-do.call(rbind,strsplit(seqs,''))
  out<-apply(seqMat,2,function(x){
    x<-x[x %in% chars]
    baseTable<-(table(c(chars,x))-1)
    baseTable<-baseTable+priors[names(baseTable)]
    return(baseTable/(length(x)+sum(priors)))
  })
  return(out)
}

#' Calculate PWM scores for sequences
#'
#' @param seqs vector of strings, sequences to score
#' @param pwm position weight matrix with uniqueBases x seqLength dimensions
#' @export
#' @return vector with a score for each sequence
#' @examples
#' seqs<-c('ACA','ACA','ACA','ACT','ACT','ACC')
#' seqPwm<-pwm(seqs)
#' scoreFromPWM(seqs,seqPwm)
scoreFromPWM<-function(seqs,pwm){
  bases<-rownames(pwm)  
  nBases<-length(bases)
  seqLength<-nchar(seqs[1])
  if(any(nchar(seqs)!=seqLength))stop(simpleError('All sequences not same length for PWM'))
  if(seqLength!=ncol(pwm))stop(simpleError('Sequences and PWM not same length'))
  badSeqs<-grepl(sprintf('[^%s]',paste(bases,collapse='')),seqs)
  if(any(badSeqs))warning('Unknown bases for PWM scoring')
  seqs<-seqs[!badSeqs]

  nSeqs<-length(seqs)
  scores<-rep(NA,length(seqs))
  seqMat<-do.call(rbind,strsplit(seqs[!badSeqs],''))
  ids<-1:nBases
  names(ids)<-bases
  #assuming as.vector takes items columnwise
  probs<-matrix(indexMatrix(ids[as.vector(seqMat)],rep(1:seqLength,each=nSeqs),pwm),nrow=nSeqs)
  scores[!badSeqs]<-apply(log(probs),1,sum)
  return(scores)
}

#' Calculate sequence for PWM
#'
#' @param pwm a matrix with rows unique bases and columns positions
#' @param cutoff consider any base with proportion greater than this cutoff present
#' @param singleCutoff only label base with a single bases if the proportion is greater than this cutoff
#' @export
#' @return a single character string containing the conensus sequence using ambiguous base codes
#' @seealso \code{\link{pwm}}, \code{\link{ambiguous2regex}}
#' @examples
#' seqs<-c('ACA','ACA','ACA','ACT','ACT','ACC')
#' seqPwm<-pwm(seqs)
#' pwmToSeq(seqPwm)
pwmToSeq<-function(pwm,cutoff=.2,singleCutoff=cutoff){
  isPresent<-apply(pwm,2,function(xx)if(sum(xx>cutoff)>1)xx>=cutoff else xx>=singleCutoff)
  bases<-matrix(rownames(pwm),nrow=nrow(pwm),ncol=ncol(pwm))
  bases[!isPresent]<-NA
  out<-paste(apply(bases,2,function(xx){
    if(all(is.na(xx)))return('N')
    else bases2ambiguous(xx[!is.na(xx)])
  }),collapse='')
  return(out)
}

