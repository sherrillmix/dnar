#' Read fastq file
#'
#' @param fileName name of fastq file
#' @param convert if TRUE convert condensed quals to space separated numeric quals
#' @param baseQual integer to subtract from quality characters to convert into quality space
#' @export
#' @return dataframe with name, seq, qual
#' @examples
#' fastq<-c(
#'   "@SEQ_ID",
#'   "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT",
#'   "+",
#'	  "!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65"
#' )
#' read.fastq(textConnection(fastq))
read.fastq<-function(fileName,convert=FALSE,baseQual=33){
	#assuming no comments and seq and qual on a single line each
	#assuming any line starting with @ 2 lines later by + is the block and no extra chars (who designed this format?)
	x<-readLines(fileName)
	plusLines<-grep('^\\+',x)
	atLines<-grep('^@',x)
	#make sure matching + and @
	plusLines<-plusLines[plusLines %in% (atLines+2)]
	atLines<-atLines[atLines %in% (plusLines-2)]
	if(any(grep('[^ACTGN]',x[atLines+1])))warning('Non ATCGN chars found in sequence')
	if(length(plusLines)!=length(atLines))stop(simpleError('Problem finding @ + lines'))
	output<-data.frame('name'=sub('^@','',x[atLines]), 'seq'=x[atLines+1], 'qual'=x[atLines+3],stringsAsFactors=FALSE)
	if(any(nchar(output$seq)!=nchar(output$qual)))stop(simpleError('Sequence and qual lengths do not match'))
	if(convert)output$qual<-unlist(lapply(x[atLines + 3],function(x)paste(as.integer(charToRaw(x))-baseQual,collapse=' ')))
	return(output)
}

#' Read a sanger phred .phd file
#'
#' @param fileName name of file
#' @param trimQual trim bases with quality lower than trimQual from start and end
#' @export
#' @return vector of sequence and space seperated qualities
read.phd<-function(fileName,trimQual=-Inf){
	tmp<-readLines(fileName)
	skip<-grep('BEGIN_DNA',tmp)[1]
	lastLine<-grep('END_DNA',tmp)[1]
	thisData<-do.call(rbind,strsplit(tmp[(skip+1):(lastLine-1)],'[ \t]'))
	lims<-range(which(as.numeric(thisData[,2])>trimQual))
	return(c('seq'=paste(thisData[lims[1]:lims[2],1],collapse=''),'qual'=paste(thisData[lims[1]:lims[2],2],collapse=' ')))
}

#' Read in a bunch of fasta files in a target directory
#'
#' @param dir target directory
#' @param suffix regex to select file names
#' @param recursive if TRUE recurse through data directory
#' @param vocal if TRUE display status message for each file loading
#' @param ... additional arguments to read.fa
#' @export
#' @return data.frame with columns name, seq and file
#' @examples
#' tmpDir<-tempdir()
#' files<-file.path(tmpDir,sprintf('%d.fa',1:9))
#' for(ii in files){
#'	  tmp<-generateFakeFasta(5,10:20,c('A','C','T','G'))
#'	  write.fa(tmp$name,tmp$seq,ii)
#' }
#' seqs<-readFaDir(tmpDir)
readFaDir<-function(dir='.',suffix='\\.(fn?a|fasta)$',recursive=FALSE,vocal=FALSE,...){
	faFiles<-list.files(dir,suffix,recursive=recursive)
	if(length(faFiles)<1)stop(simpleError('No fa files found'))
	for(ii in faFiles){
		if(vocal)message('Working on ',ii)
		tmp<-read.fa(file.path(dir,ii),...)
		tmp$file<-ii
		if(exists('allFa'))allFa<-rbind(allFa,tmp)
		else allFa<-tmp
	}
	return(allFa)
}

#' Read a fasta file 
#' 
#' @param fileName name of file
#' @param assumeSingleLine if TRUE don't process sequence lines. just assume they're one line per sequence
#' @param ... additional arguments to readLines
#' @export
#' @return  data.frame with columns name and seq
#' @examples
#' file<-tempfile()
#' x<-generateFakeFasta(100,bases=c('A','C','T','G','\n','-'))
#' write.fa(x$name,x$seq,file)
#' y<-read.fa(file)
read.fa<-function(fileName,assumeSingleLine=FALSE,...){
	x<-readLines(fileName,warn=FALSE,...)
	if(length(x)==0|sum(nchar(x))==0)return(NULL)
	if(assumeSingleLine){
		output<-data.frame('name'=x[seq(1,length(x),2)],'seq'=x[seq(2,length(x),2)],stringsAsFactors=FALSE)
		if(any(grepl('^>',output$seq,perl=TRUE)|!grepl('^>',output$name,perl=TRUE)))stop(simpleError('Problem reading single line fasta'))
		output$name<-substring(output$name,2)
	}else{
		x<-x[!grepl('^[;#]',x,perl=TRUE)&x!='']
		nameLines<-grep('^>',x,perl=TRUE)
		thisNames<-sub('^>','',x[nameLines],perl=TRUE)
		seqs<-apply(cbind(nameLines+1,c(nameLines[-1]-1,length(x))),1,function(coords){
			if(coords[1]<=coords[2])return(paste(x[coords[1]:coords[2]],collapse=''))
			else return('')
		})
		output<-data.frame('name'=thisNames,'seq'=seqs,stringsAsFactors=FALSE)
	}
	output$seq<-sub(' +$','',output$seq,perl=TRUE)
	return(output)
}



#' Generate fake name and sequence data
#'
#' @param nSeq number of sequences to generate
#' @param nChar possible lengths of sequences
#' @param bases possible bases
#' @export
#' @return data.frame with nSeq rows and columns name and seq
generateFakeFasta<-function(nSeq=10000,nChar=100:1000,bases=c('A','C','T','G','-','N')){
	nChars<-sample(nChar,nSeq,TRUE)
	names<-sprintf('>%d_%s',1:nSeq,replicate(nSeq,paste(sample(c(letters,LETTERS),30,TRUE),collapse='')))
	seqs<-sapply(nChars,function(x)paste(sample(bases,x,TRUE),collapse=''))
	return(data.frame('name'=names,'seq'=seqs))
}


#' Call samtools view on a sam/bam file
#'
#' @param fileName vector of sam/bam files to read
#' @param samArgs 1 element character vector of args to pass to samtools
#' @param samtoolsBinary location of samtools
#' @param vocal if TRUE print status messages
#' @param samCommand which samtools tool to run
#' @param ... additional arguments for read.sam
#' @export
#' @return data.frame with columns as read.sam
samView<-function(fileName,samArgs='',...,samtoolsBinary='samtools',vocal=FALSE,samCommand='view'){
	if(length(fileName)>1){ #recurse
		allOut<-do.call(rbind,lapply(fileName,function(x){
			out<-samView(x,samArgs,...,samtoolsBinary=samtoolsBinary)
			out$file<-x
			return(out)
		}))
		return(allOut)
	}else{ #do samtools on single file
		cmd<-sprintf('%s %s %s %s',samtoolsBinary,samCommand,fileName,samArgs[1])
		if(vocal)message(cmd)
		samOut<-textConnection(system(cmd,intern=TRUE))
		if(vocal)message('Parsing')
		if(length(readLines(samOut,n=1))<1){
			out<-NULL
		}else{
			if(samCommand=='view')out<-read.sam(samOut,skips=0,...)
			else if(samCommand=='depth'||grepl('bam2depth',samtoolsBinary))out<-utils::read.table(samOut,sep='\t',stringsAsFactors=FALSE,...)
			else out<-readLines(samOut)
		}
		close(samOut)
		return(out)
	}
}


#' Read a sam file
#'
#' @param fileName name of file
#' @param nrows number of rows to return. -1 for all rows
#' @param skips number of lines to skip (if negative find @ headers automtically)
#' @param condense throw out a bunch of columns for smaller file size?
#' @export
#' @return dataframe with columns qName, flag, tName, pos, cigar, seq
read.sam<-function(fileName,nrows=-1,skips=-1,condense=TRUE){
	colNames<-c('qName','flag','tName','pos','mapq','cigar','mrnm','mpos','isize','seq','qual','tags')
	if(condense)colClasses<-c('character','numeric','character','numeric','null','character','null','null','null','character','null')
	else colClasses<-c('character','numeric','character',rep('numeric',2),rep('character',2),rep('numeric',2),rep('character',2))
	
	if(skips<0){
		testLines<-readLines(fileName,n=1000)
		atLines<-which(grepl('^@',testLines))
		if(any(atLines))skips<-max(atLines)
		else skips=-1
		message('Found ',skips,' line header')
	}

	#-15 to exclude tags that can be variable length and tab seperated
	#will need to modify if we want tags
	x<-scan(fileName, what = list('character'='','numeric'=1,'null'=NULL)[colClasses[-12]], fill=TRUE,skip=skips,nmax=nrows,flush=TRUE)
	#x<-read.table(fileName,skip=skips,sep="\t",stringsAsFactors=FALSE,colClasses=colClasses,nrows=nrows,col.names=colNames)
	isNull<-sapply(x,is.null)
	x<-do.call(data.frame,c(x[!isNull],stringsAsFactors=FALSE))
	colnames(x)<-colNames[-12][!isNull]
	return(x)
}


#' Fill zeros in cover data
#'
#' @param cover output from pullRegion with missing zero positions
#' @param posCol name of position column
#' @param countCols vector of column names for columns containing counts
#' @keywords internal 
#' @return filled in data.frame
fillZeros<-function(cover,posCol='pos',countCols=colnames(cover)[grep('counts',colnames(cover))]){
	repeatedCols<-!colnames(cover) %in% c(posCol,countCols)
	if(any(apply(cover[,repeatedCols,drop=FALSE],2,function(x)length(unique(x)))>1))stop(simpleError('Found nonunique extra columns in fillZeros'))
	cover<-cover[order(cover[,posCol]),]
	diffs<-diff(cover[,posCol])
	missingZeros<-which(diffs>1)
	if(!any(missingZeros))return(cover)
	missingPos<-unlist(mapply(function(start,end)start:end,cover[missingZeros,posCol]+1,cover[missingZeros+1,posCol]-1,SIMPLIFY=FALSE))
	filler<-cover[rep(1,length(missingPos)),]
	filler[,countCols]<-0
	filler[,posCol]<-missingPos
	out<-rbind(cover,filler)
	out<-out[order(out[,posCol]),]
	return(out)
}


#' Read a tab-delimited blast file
#'
#' @param fileName name of file
#' @param skips number of lines to skip (0 for a normal blast file)
#' @param nrows number of rows to read in (-1 for all)
#' @param calcScore if TRUE, calculate score column
#' @export
#' @return dataframe of blat data
read.blast<-function(fileName,skips=0,nrows=-1,calcScore=TRUE){
	x<-utils::read.table(fileName,skip=skips,sep="\t",stringsAsFactors=FALSE,colClasses=c(rep('character',2),rep('numeric',10)),nrows=nrows)
	colnames(x)<-c('qName','tName','percID','alignLength','mismatch','nGap','qStart','qEnd','tStart','tEnd','eValue','hspBit')
	#Score equation from blat's webpage
	if(calcScore)x$score<-x$alignLength-x$mismatch-x$nGap
	return(x)
}

#' Read a large blat file piece by piece 
#'
#' @section Warning:
#' Assumes blat files is sorted by qName
#'
#' @param fileName name of file
#' @param nHeader number of lines to skip at start of file
#' @param nrows number of lines ot read in per shot
#' @param isGz if TRUE then file is assumed to be gzipped
#' @param filterFunc a function taking a data.frame of blat data and returning a TRUE/FALSE vector specifying rows to keep 
#' @param ... additional arguments to read.blat
#' @export
#' @return data.frame of blat data
readLargeBlat<-function(fileName,nHeader=5,nrows=1e6,isGz=grepl('.gz$',fileName),filterFunc=function(x)rep(TRUE,nrow(x)),...){
	#R complains about seek and gzipped files. should work but...
	#openFile<-gzfile(fileName,'r')
	#so we'll do the stupid way if it's a gzipped file
	if(isGz){
		openFile<-pipe(sprintf('zcat %s',fileName),'r')
	}else{
		openFile<-file(fileName,'r')
	}
	on.exit(close(openFile))
	#burn off the header
	readLines(openFile,n=nHeader)
	leftOverData<-read.blat(openFile,skips=0,nrows=1,...)
	out<-leftOverData[0,]
	counter<-1
	finished<-FALSE
	while(!finished){
		message('Working on ',counter*nrows)
		counter<-counter+1
		thisData<-read.blat(openFile,skips=0,nrows=nrows,...)
		finished<-nrow(thisData)<1
		thisData<-rbind(leftOverData,thisData)
		if(finished){
			#we're done
			leftOverData<-thisData[0,]
		}else{
			#we're not sure if the last entry is complete
			lastSelector<-thisData$qName==utils::tail(thisData$qName,1)
			leftOverData<-thisData[lastSelector,]
			thisData<-thisData[!lastSelector,]
		}
		#if all reads were same qName then we could have empty (otherwise shouldn't occur)
		if(nrow(thisData)>0){
			selector<-filterFunc(thisData)
			if(any(is.na(selector)))stop(simpleError('NA returned from selection function'))
			message('Filtering ',sum(!selector),' of ',length(selector),' matches for ',length(unique(thisData$qName)),' reads')
			out<-rbind(out,thisData[selector,])
		}
	}
	return(out)
}

#' Read a blat file
#'
#' @param fileName name of file
#' @param skips number of lines to skip (5 for a normal blat file)
#' @param nCols 21 (no alignments) or 23 columns (alignments)
#' @param calcScore calculate score column?
#' @param fixStarts convert start to 1-based?
#' @param filterFunction if not NULL then pass the output data.frame to this function and filter by logical value returned
#' @param ... additional arguments to read.table
#' @export
#' @return dataframe of blat data
read.blat<-function(fileName,skips=5,calcScore=TRUE,fixStarts=TRUE,nCols=21,filterFunction=NULL,...){
	#if we don't read in all lines at once it becomes a pain to deal with open/closed connections when we want to peak 
	
	#test for empty file
	#if(length(allLines)<1)return(NULL) 
	#testDf<-read.table(textConnection(allLines[1]),sep="\t",stringsAsFactors=FALSE,...)
	if(!nCols %in% c(21,23))stop(simpleError('Please select 21 or 23 columns'))
	colNames<-c('match','mismatch','repmatch','ns','qGaps','qGapBases','tGaps','tGapBases','strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blocks','blockSizes','qStarts','tStarts')
	colClasses<-c(rep('numeric',8),rep('character',2),rep('numeric',3),'character',rep('numeric',4),rep('character',3))
	if(nCols==23){
		colNames<-c(colNames,'qAlign','tAlign')
		colClasses<-c(colClasses,rep('character',2))
	}
	#textConnection is too slow to be useful here
	x<-utils::read.table(fileName,sep="\t",stringsAsFactors=FALSE,skip=skips,colClasses=colClasses,col.names=colNames,...)

	#Score equation from blat's webpage
	if(calcScore)x$score<-x$match-x$mismatch-x$qGaps-x$tGaps
	if(fixStarts&nrow(x)>0){
		#blat uses 0-index starts and 1-index ends
		#put starts in 1-index
		x$tStart<-x$tStart+1
		x$qStart<-x$qStart+1
		x$tStarts<-sapply(strsplit(x$tStarts,','),function(x)paste(as.numeric(x)+1,collapse=','))
		x$qStarts<-sapply(strsplit(x$qStarts,','),function(x)paste(as.numeric(x)+1,collapse=','))
		x$qStartsBak<-x$qStarts
		negSelect<-x$strand=='-'
		#qSize-(blockSizes-1)-(qStarts-1)
		if(any(negSelect)){
			x$qStarts[negSelect]<-mapply(function(qStart,blockSize,qSize){
				paste(as.numeric(qSize)-as.numeric(qStart)-as.numeric(blockSize)+2,collapse=',')
			},strsplit(x$qStarts[negSelect],','),strsplit(x$blockSizes[negSelect],','),x$qSize[negSelect])
		}
	}
	if(!is.null(filterFunction)){
		x<-x[filterFunction(x),]
	}
	return(x)
}

#' Read a UCSC wiggle file
#'
#' @param fileName wiggle file to read
#' @export
#' @return data.frame with columns start, value, end and chr
read.wiggle<-function(fileName){
	x<-readLines(fileName)
	trackLines<-grep('chrom=',x)
	if(length(trackLines)==0)stop(simpleError('No track lines found in wiggle file'))
	out<-do.call(rbind,mapply(function(startLine,endLine){
		regex<-'^.*chrom=([^ \t]+).*$'
		if(!grepl(regex,x[startLine]))stop(simpleError("Can't assign chromosome in wiggle file"))
		chr<-sub(regex,'\\1',x[startLine])
		regex<-'^.*span=([^ \t]+).*$'
		if(grepl(regex,x[startLine]))span<-as.numeric(sub(regex,'\\1',x[startLine]))
		else span<-1
		out<-data.frame(do.call(rbind,strsplit(x[(startLine+1):endLine],'[\t ]')),stringsAsFactors=FALSE)
		colnames(out)<-c('start','value')
		out$start<-as.numeric(out$start)
		out$value<-as.numeric(out$value)
		out$end<-out$start+span-1
		out$chr<-chr
		return(out)
	},trackLines,c(trackLines[-1]-1,length(x)),SIMPLIFY=FALSE))
	return(out)
}

#' Read an ace file
#'
#' @section Warning:
#' Only reads 1 contig
#'
#' @param aceFile string of file name or file handle
#' @param dropMosaik Remove extraneous line from Mosaik labelled either .MosaikAnchor.C1 or MosaikReference
#' @param checkSnps Find SNPs from pyroBayes
#' @param vocal if TRUE print status messages
#' @export
#' @return list of reference sequence in [[1]], aligned reads in [[2]] (note reads are not globally aligned still need to use start coordinate to place globally), snps in [[3]]
read.ace<-function(aceFile,dropMosaik=TRUE,checkSnps=TRUE,vocal=TRUE){
	#debug<-FALSE
	#if(!exists(x)|!debug)x<-readLines(aceFile)
	x<-readLines(aceFile)
	asLines<-grep('^AS [0-9]+ [0-9]+ *$',x,perl=TRUE)
	if(length(asLines)!=1)stop(simpleError('Incorrect number of AS lines found'))
	asLine<-strsplit(x[asLines],' ')[[1]]
	numContigs<-asLine[2]
	numReads<-asLine[3]
	if(numContigs!=1)stop(simpleError('Sorry this function only handles 1 contig'))
	if(vocal)message('Expecting ',numReads,' reads')
	coLines<-grep('^CO [^ ]+ [0-9]+ [0-9]+ [0-9]+ [UC] *$',x,perl=TRUE)
	if(length(coLines)!=1)stop(simpleError('Incorrect number of CO lines found'))
	coLine<-strsplit(x[coLines],' ')[[1]]
	contigName<-coLine[2]
	contigBaseNum<-coLine[3]
	contigReadNum<-coLine[4]
	if(vocal)message('Contig ',contigName,' has ',contigBaseNum,' bases and ',contigReadNum,' reads')
	bqLines<-grep('^BQ *$',x)
	if(length(bqLines)!=1)stop(simpleError('Incorrect number of BQ lines found'))
	contigSeq<-paste(x[(coLines+1):(bqLines-1)],sep='',collapse='')
	if(nchar(contigSeq)!=contigBaseNum)stopError('Found ',nchar(contigSeq),' bases in contig but was expecting ',contigBaseNum,' bases')	
	#skipping BQ
	afLines<-grep('^AF [^ ]+ [UC] [0-9]+ *$',x)
	if(length(afLines)!=contigReadNum)stopError(length(afLines),' AF lines found but was expecting ',contigReadNum)
	afLine<-strsplit(x[afLines],' ')
	reads<-as.data.frame(do.call(rbind,lapply(afLine,function(x)return(x[c(2,3,4)]))),stringsAsFactors=FALSE)
	colnames(reads)<-c('name','dir','start')
	reads$start<-as.numeric(reads$start)
	#skipping BS
	rdLines<-grep('^RD [^ ]+ [0-9]+ [0-9]+ [0-9]+ *$',x)
	if(length(rdLines)!=contigReadNum)stopError(length(afLines),' RD lines found but was expecting ',contigReadNum)
	qaLines<-grep('^QA [0-9]+ [0-9]+ [0-9]+ [0-9]+ *$',x)
	if(length(qaLines)!=contigReadNum)stopError(length(qaLines),' QA lines found but was expecting ',contigReadNum)
	#Assuming RD lines are followed immediately by QA
	readSeqs<-apply(cbind(rdLines,qaLines),1,function(x,y){return(paste(y[(x[1]+1):(x[2]-1)],sep='',collapse=''))},x)
	reads$seq<-readSeqs
	qaLine<-strsplit(x[qaLines],' ')
	qa<-as.data.frame(do.call(rbind,lapply(qaLine,function(x)return(x[c(2,3,4,5)]))),stringsAsFactors=FALSE)
	colnames(qa)<-c('qStart','qEnd','aStart','aEnd')
	reads<-cbind(reads,qa)	
	reads$length<-nchar(reads$seq)
	ctLines<-grep('^CT *\\{ *$',x)
	if(length(ctLines)>0&&checkSnps){
		ctLine<-x[ctLines+1]
		snpCT<-grep('PB++',ctLine)
		psnpLine<-x[ctLines+2][snpCT]
		ctLine<-ctLine[snpCT]
		if(length(snpCT)!=length(ctLine)|length(grep('pSnp=',psnpLine))!=length(psnpLine))stop(simpleError('psnp lines not the same length as PB++ CT lines'))
		if(vocal)message('Found ',length(ctLine),' CT lines of which ',length(snpCT),' were from PB++')
		ct<-strsplit(ctLine,' ')
		notes<-as.data.frame(do.call(rbind,lapply(ct,function(x)return(x[c(4,5,1)]))),stringsAsFactors=FALSE)
		colnames(notes)<-c('start','stop','contig')
		notes$psnp<-sub('pSnp=','',psnpLine)
	}else{
		if(vocal)message('Not looking for snp notes')
		notes<-c()
	}
	if(dropMosaik&(reads[1,'name']=='.MosaikAnchor.C1'|reads[1,'name']=='.MosaikReference')){
		reads<-reads[-1,]
		if(vocal)message('Removing .MosaikAnchor.C1')
	}
	return(list('contig'=contigSeq,'reads'=reads,'notes'=notes))
}

#' Read a bed file
#'
#' @param fileName .bed file to be read in
#' @param startAddOne if TRUE add 1 to starts to adjust for 0-based start, 1-based ends in UCSC formats
#' @export
#' @return list with a dataframe with columns chr, start and end for each track
read.bed<-function(fileName,startAddOne=FALSE){
	x<-readLines(fileName)
	tracks<-grep('track',x)
	if(length(tracks)==0){
		tracks<-c(0)
		trackNames<-'main'
	} else trackNames<-gsub('.*name=([^ ]+).*','\\1',x[tracks])
	message('Found ',length(tracks),' tracks')
	tracks<-c(tracks,length(x)+1)
	output<-list()
	for(i in 1:(length(tracks)-1)){
		output[[trackNames[i]]]<-data.frame(do.call(rbind,strsplit(x[(tracks[i]+1):(tracks[i+1]-1)],'\t')),stringsAsFactors=FALSE)
		thisNames<-c('chr','start','end')
		if(ncol(output[[trackNames[i]]])==4)thisNames<-c(thisNames,'name')
		if(ncol(output[[trackNames[i]]])==12)thisNames<-c(thisNames,'name','score','strand','thickStart','thickEnd','rgb','blockCount','blockSizes','blockStarts')
		colnames(output[[trackNames[i]]])<-thisNames
		output[[trackNames[i]]][,'start']<-as.numeric(output[[trackNames[i]]][,'start'])
		output[[trackNames[i]]][,'end']<-as.numeric(output[[trackNames[i]]][,'end'])
		if(startAddOne)output[[trackNames[i]]][,'start']<-output[[trackNames[i]]][,'start']+1
	}
	return(output)
}

#' Pull counts from a bam file for target region
#'
#' @param reg region in the format "chrX:123545-123324"
#' @param files bam files to pull the counts from
#' @param bam2depthBinary bam2depth executable file
#' @param fillMissingZeros if TRUE fill in uncovered positions with zeros otherwise may be left out of output by bam2depth
#' @export
#' @return data.frame with columns chr, pos and counts 
pullRegion<-function(reg,files,bam2depthBinary='./bam2depth',fillMissingZeros=TRUE){
	region<-parseRegion(reg)
	region$start<-region$start+1 #bam2depth using ucsc 0-start, 1-ends
	samArg<-sprintf('-r %s',reg)
	fileArg<-paste(files,collapse=' ')
	cover<-samView(fileArg,samArgs=samArg,samCommand='',samtoolsBinary=bam2depthBinary,colClasses=c('character',rep('numeric',length(files)+1)))
	if(is.null(cover))cover<-do.call(data.frame,c(list('XXX'),as.list(-(1:(length(files)+1)))))[0,]
	countCols<-sprintf('counts%d',1:(ncol(cover)-2))
	colnames(cover)<-c('chr','pos',countCols)
	if(fillMissingZeros){
		#deal with missing start or ends
		filler<-cover[c(1,1),]
		filler$chr<-region$chr
		filler[,countCols]<-0
		filler$pos<-unlist(region[,c('start','end')])
		if(!region$start %in% cover$pos)cover<-rbind(filler[1,],cover)
		if(!region$end %in% cover$pos)cover<-rbind(cover,filler[2,])
		cover<-fillZeros(cover)
	}
	return(cover)
}

#' Writes sequences to a fasta file
#' @param names names for the > lines
#' @param dna the sequences
#' @param fileName file to write to
#' @param isGz if TRUE write to a gzipped file
#' @export
#' @return NULL
write.fa<-function(names,dna,fileName,isGz=grepl('.gz$',fileName)){
	if(any(grepl('^[^>]',names)))names<-paste('>',names,sep='')
	dna<-sub(' +$','',dna,perl=TRUE)
	output<-paste(names,dna,sep="\n")
	if(isGz)fileName<-gzfile(fileName)
	writeLines(output,sep="\n",con=fileName)
	if(isGz)close(fileName)
	return(NULL)
}

#' Writes a fastq file
#'
#' @param names sequence names for @/+ line
#' @param seqs the sequences strings
#' @param quals the qualities strings
#' @param fileName file to write to
#' @param isGz if TRUE gz compress output
#' @export
#' @return NULL
write.fastq<-function(names,seqs,quals,fileName,isGz=grepl('.gz$',fileName)){
	names1<-paste('@',names,sep='')
	names2<-paste('+',names,sep='')
	output<-paste(names1,seqs,names2,quals,sep="\n")
	if(isGz)fileName<-gzfile(fileName)
	writeLines(output,sep="\n",con=fileName)
	if(isGz)close(fileName)
	return(NULL)
}


