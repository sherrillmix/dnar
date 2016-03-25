#read fastq file
#fileName: name of fastq file
#convert: convert condensed quals to numeric quals?
#returns: dataframe with name, seq, qual
read.fastq<-function(fileName,convert=TRUE){
	#assuming no comments and seq and qual on a single line each
	#assuming any line starting with @ 2 lines later by + is the block and no extra chars (who designed this format?)
	#as.integer(charToRaw())
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
	if(convert)output$qual<-unlist(lapply(x[atLines + 3],function(x)paste(as.integer(charToRaw(x))-33,collapse=' ')))
	return(output)
}

#read a sanger phred .phd file
#fileName:name of file
#trimEnds:trim off low quality bases and quals at start and end?
#trimQual:trim bases with quality lower than trimQual from start and end
#returns: vector of sequence and space seperated qualities
read.phd<-function(fileName,trimEnds=TRUE,trimQual=30){
	tmp<-readLines(fileName)
	skip<-grep('BEGIN_DNA',tmp)[1]
	lastLine<-grep('END_DNA',tmp)[1]
	thisData<-do.call(rbind,strsplit(tmp[(skip+1):(lastLine-1)],'[ \t]'))
	if(trimEnds) lims<-range(which(as.numeric(thisData[,2])>trimQual))
	else lims<-c(1,nrow(thisData))
	return(c('seq'=paste(thisData[lims[1]:lims[2],1],collapse=''),'qual'=paste(thisData[lims[1]:lims[2],2],collapse=' ')))
}

#read in a bunch of fasta files in a target directory
#dir: target directory
#suffix: regex to select file names
#recursive: recurse through data directory?
#vocal: status message for each file loading
#...: extra arguments to read.fa
readFaDir<-function(dir='.',suffix='\\.(fn?a|fasta)$',recursive=FALSE,vocal=FALSE,...){
	faFiles<-list.files(dir,suffix,recursive=recursive)
	if(length(faFiles)<1)stop(simpleError('No fa files found'))
	for(i in faFiles){
		if(vocal)message('Working on ',i)
		tmp<-read.fa(sprintf('%s/%s',dir,i),...)
		tmp$file<-i
		if(exists('allFa'))allFa<-rbind(allFa,tmp)
		else allFa<-tmp
	}
	return(allFa)
}


#read a fasta file
#fileName:name of file
#longNameTrim:trim off anything after a space for name column (preserve original in long name)
#assumeSingleLine:don't process sequence lines. just assume they're one line per sequence
#returns: dataframe with columns name (name between > and the first ' '), seq (sequence), and longName (the whole > line)
read.fa<-function(fileName,longNameTrim=TRUE,assumeSingleLine=FALSE){
	x<-readLines(fileName)
	if(assumeSingleLine){
		output<-data.frame('longName'=x[seq(1,length(x),2)],'seq'=x[seq(2,length(x),2)],stringsAsFactors=FALSE)
		if(any(grepl('^>',output$seq,perl=TRUE)|!grepl('^>',output$longName,perl=TRUE)))stop(simpleError('Problem reading single line fasta'))
		output$longName<-substring(output$longName,2)
	}else{
		selector<-grep('^[^>].* .*[^ ]$',x,perl=TRUE)
		x[selector]<-paste(x[selector],' ',sep='')
		if(length(x)==0)return(NULL)
		x<-x[!grepl('^[#;]',x,perl=TRUE)&x!='']
		y<-paste(x,collapse="\n")
		splits<-strsplit(y,'>',fixed=TRUE)[[1]][-1]
		splits2<-strsplit(splits,"\n",fixed=TRUE)
		output<-lapply(splits2,function(x){return(c(x[1],paste(x[-1],collapse='')))})
		output<-as.data.frame(do.call(rbind,output),stringsAsFactors=FALSE)
		colnames(output)<-c('longName','seq')
	}
	if(longNameTrim){
		output$name<-unlist(lapply(strsplit(output$longName,' ',fixed=TRUE),function(x)x[1]))
		output<-output[,3:1]
	}else{
		colnames(output)<-c('name','seq')	
	}
	output$seq<-gsub(' +$','',output$seq,perl=TRUE)
	return(output)
}
#alternative version of the above (a bit quicker)
read.fa2<-function(fileName=NULL,longNameTrim=TRUE,x=NULL,...){
	if(is.null(fileName)&is.null(x))stop(simpleError('Please specify fileName or string vector x'))
	if(is.null(x))x<-readLines(fileName,warn=FALSE,...)
	if(length(x)==0)return(NULL)
	x<-x[!grepl('^[;#]',x,perl=TRUE)&x!='']
	nameLines<-grep('^>',x,perl=TRUE)
	thisNames<-sub('^>','',x[nameLines],perl=TRUE)
	hasSpaces<-any(grep(' [^ ]',x[-nameLines],perl=TRUE))
	seqs<-apply(cbind(nameLines+1,c(nameLines[-1]-1,length(x))),1,function(coords){
		if(coords[1]<=coords[2])return(paste(x[coords[1]:coords[2]],collapse=''))
		else return('')
	})
	#seqs<-gsub('  +',' ',seqs,perl=TRUE)
	seqs<-sub(' +$','',seqs,perl=TRUE)

	output<-data.frame('longName'=thisNames,'seq'=seqs,stringsAsFactors=FALSE)
	if(longNameTrim){
		output$name<-unlist(lapply(strsplit(output$longName,' ',fixed=TRUE),function(x)x[1]))
		output<-output[,3:1]
	}else{
		colnames(output)<-c('name','seq')	
	}
	return(output)
}

#call samtools view on a sam/bam file
#fileName: sam/bam file to read
#samArgs: 1 element character vector of args to pass to samtools
#samtoolsBinary: location of samtools
#vocal: print status messages?
#samCommand: which samtools tool to run
#...: additional arguments for read.sam
#samtoolsBinary: location of samtools
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
			else if(samCommand=='depth'||grepl('bam2depth',samtoolsBinary))out<-read.table(samOut,sep='\t',stringsAsFactors=FALSE,...)
			else out<-readLines(samOut)
		}
		close(samOut)
		return(out)
	}
}


#read a sam file
#fileName: name of file
#nrows: number of rows to return
#skips: number of lines to skip (if negative find @ headers automtically)
#condense: throw out a bunch of columns for smaller file size?
#returns: dataframe with columns 
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


#cover: output from pullRegion with missing zero positions
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


#read a tab-delimited blast file
#fileName: name of file
#skips: number of lines to skip (0 for a normal blast file)
#nrows: number of rows to read in (-1 for all)
#calcScore: calculate score column?
#returns: dataframe of blat data
readBlast<-function(fileName,skips=0,nrows=-1,calcScore=TRUE){
	x<-read.table(fileName,skip=skips,sep="\t",stringsAsFactors=FALSE,colClasses=c(rep('character',2),rep('numeric',10)),nrows=nrows)
	colnames(x)<-c('qName','tName','percID','alignLength','mismatch','nGap','qStart','qEnd','tStart','tEnd','eValue','hspBit')
	#Score equation from blat's webpage
	if(calcScore)x$score<-x$alignLength-x$mismatch-x$nGap
	return(x)
}

readBlast2<-function(fileName,excludeUnculture=TRUE){
	x<-readLines(fileName)
	queries<-grep('Query=',x)
	queryData<-gsub('Query= *','',x[queries])
	datas<-grep('\\|.*[0-9]+  +[0-9]',x)
	if(excludeUnculture)datas<-datas[!grepl('(Uncultured)|(Unidentified)',x[datas])]
	queryVec<-queryData[sapply(datas,function(x)max(which(queries<x)))]
	dataData<-gsub('  +','\t',x[datas])
	dataMat<-do.call(rbind,strsplit(dataData,'\t'))
	out<-data.frame(queryVec,dataMat,stringsAsFactors=FALSE)
	colnames(out)<-c('query','gNum','descr','score','eval')
	return(out)
}


#read a large blat file piece by piece NOTE: assumes blat files is sorted by qName
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
	leftOverData<-readBlat(openFile,skips=0,nrows=1,...)
	out<-leftOverData[0,]
	counter<-1
	finished<-FALSE
	while(!finished){
		message('Working on ',counter*nrows)
		counter<-counter+1
		thisData<-readBlat(openFile,skips=0,nrows=nrows,...)
		finished<-nrow(thisData)<1
		thisData<-rbind(leftOverData,thisData)
		if(finished){
			#we're done
			leftOverData<-thisData[0,]
		}else{
			#we're not sure if the last entry is complete
			lastSelector<-thisData$qName==tail(thisData$qName,1)
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

#read a blat file
#fileName: name of file
#skips: number of lines to skip (5 for a normal blat file)
#nCols: 21 (no alignments) or 23 columns (alignments)
#calcScore: calculate score column?
#fixStarts: convert start to 1-based?
#filterFunction: if not NULL then pass the output data.frame to this function and filter by logical value returned
#...: additional arguments to read.table
#returns: dataframe of blat data
readBlat<-function(fileName,skips=5,calcScore=TRUE,fixStarts=TRUE,nCols=21,filterFunction=NULL,...){
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
	x<-read.table(fileName,sep="\t",stringsAsFactors=FALSE,skip=skips,colClasses=colClasses,col.names=colNames,...)

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

#fileName: wiggle file to read
readWiggle<-function(fileName){
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

#read an ace file
#LIMITATION: only reads 1 contig
#aceFile: string of file name or file handle
#dropMosaik: Remove extraneous line from Mosaik labelled either .MosaikAnchor.C1 or MosaikReference
#checkSnps: Find SNPs from pyroBayes
#returns: list of reference sequence in [[1]], aligned reads in [[2]] (note reads are not globally aligned still need to use start coordinate to place globally), snps in [[3]]
readAce<-function(aceFile,dropMosaik=TRUE,checkSnps=TRUE,vocal=TRUE){
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
	if(nchar(contigSeq)!=contigBaseNum)stopError('Found ',char(contigSeq),' bases in contig but was expecting ',contigBaseNum,' bases')	
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
	#readDir<-lapply(afLine,function(x)return(x[3]))
	#readStart<-lapply(afLine,function(x)return(x[4]))
	#reads<-data.frame('name'=readNames,'dir'=readDir,'start'=readStart)
}

#take the output from an ace file and fill in the starts and ends of sequences with gaps to make a global alignment
#seqs: sequences from ace file
#starts: starting location for each sequence 
#low: start point for global alignment
#high: end point for global alignment
#lengths: lengths of seqs (probably can go with default in 99% of cases)
#filter: filter out reads falling outside range?
#returns: list of aligned seqs in [[1]], logical vector of whether read fell within cut region in [[2]]
cutReads<-function(seqs,starts,low=min(starts),high=max(starts+nchar(seqs)-1),lengths=nchar(seqs),filter=TRUE){
	#make a string of dots for cutting
	dots<-paste(rep('.',max(high-low+1)),collapse='')
	debug<-TRUE
	goodReads<-findReads(low,starts,lengths,high)
	if(!any(goodReads)){
		if(filter)return(FALSE)
		else return(rep(dots,length(seqs)))
	}
	thisSeqs<-seqs[goodReads];	starts<-starts[goodReads];	lengths<-lengths[goodReads]
	cutlow<-low-starts+1
	startDash<-0-cutlow+1
	startDash[startDash<0]<-0
	cutlow[cutlow<1]<-1
	cuthigh<-high-starts+1
	endDash<-high-starts-lengths+1
	endDash[endDash<0]<-0
	cuts<-substr(thisSeqs,cutlow,cuthigh)
	predots<-substring(dots,1,startDash)
	postdots<-substring(dots,1,endDash)
	cuts<-paste(predots,cuts,postdots,sep='')
	if(filter){
		return(list(cuts,goodReads))
	}else{
		out<-rep(NA,length(seqs))
		out[goodReads]<-cuts
		out[!goodReads]<-dots
		return(out)
	}
}


parseGff<-function(gffFile,individuals=NULL,contig=individuals[1]){
	message('Using ',contig,' as main contig')
	gff<-tryCatch(read.table(gffFile,stringsAsFactors=FALSE),error=function(e)return(data.frame('contig'=1,'program'=1,'type'=1,'start'=1,'stop'=1,'psnp'=1,'dummy1'=1,'dummy2'=1,'extra'=1)[0,]))
	colnames(gff)<-c('contig','program','type','start','stop','psnp','dummy1','dummy2','extra')
	gff<-gff[,!colnames(gff) %in% c('dummy1','dummy2')]
	if(nrow(gff)>0){
	details<-strsplit(gff$extra,'\\;')
	details<-do.call(rbind,details)
	gff$alleles<-gsub('^alleles=','',details[,1])
	gff$genos<-gsub('^individualGenotypes=','',details[,2])
	gff$genoProbs<-gsub('^individualGenotypeProbabilities=','',details[,3])
	gff$alleleCounts<-gsub('^individualAlleleCounts=','',details[,4])
	if(!is.null(individuals)){
		alleleCounts<-strsplit(gff$alleleCounts,',')
		alleleData<-lapply(alleleCounts,function(x){
			output<-c()
			for(i in individuals){
				countLabel<-paste('coverage_',i,sep='');alleleLabel<-paste('allele_',i,sep='');maxAlleleLabel<-paste('max_',i,sep='');pLabel<-paste('pMax_',i,sep='')
				thisIndex<-grep(i,x)
				if(length(thisIndex)==1){
					splits<-strsplit(sub('^[^:]+:','',x[thisIndex]),'&')[[1]]	
					alleles<-gsub('^([ACTG-]+)\\|.*$','\\1',splits)
					counts<-as.numeric(gsub('^[ACTG-]+\\|([0-9]+).*$','\\1',splits))
					if(!is.null(contig) & i!=contig){
						contigAllele<-output[paste('max_',contig,sep='')]
						diffCounts<-counts[alleles!=contigAllele]
						diffAlleles<-alleles[alleles!=contigAllele]
					}else{
						diffCounts<-counts
						diffAlleles<-alleles
					}
					maxCount<-max(diffCounts)
					output[countLabel] <- sumCount <- sum(counts)
					output[maxAlleleLabel]<-ifelse(maxCount==0,'',diffAlleles[diffCounts==maxCount][1])
					output[alleleLabel]<-paste(alleles,counts,sep='=',collapse='|')
					output[pLabel]<-maxCount/sumCount
				}else{
					if(length(thisIndex)>1)warning('Found multiple matches for allele counts for individual ',i)
					output[countLabel]<-0
					output[maxAlleleLabel]<-''
					output[alleleLabel]<-''
					output[pLabel]<-NA
				}
			}
			return(output)
		})
		alleleData<-do.call(rbind,alleleData)
		gff[,colnames(alleleData)]<-alleleData
		for(i in grep('coverage_',colnames(gff)))gff[,i]<-as.numeric(gff[,i])
	}
	}
	return(gff)
}


#names: names of reads e.g. >DRDR12A125
#samples: sample ID for each above read
#barcodes: list of vectors of barcodes or barcode/primers with each entry indexed by sample
#outdir: fileBinirectory to put seperate sffs
#sffDir: directory to look for sffs
#baseName: prepend to sample names 
#sffBinDir: location of sfffile binary
makeSeperateSffs<-function(names,samples,barcodes,outDir,sffDir,baseName='reads',sffBinDir='',vocal=TRUE){
	if(length(names)!=length(samples))stop(simpleError('Length of names and sample assignments differ'))
	if(sffBinDir!='')sffBinDir<-sprintf('%s/',sffBinDir)
	#apparently sfffile can't find sff files in a directory even though it says it can so we'll just feed it every sff
	sffFiles<-paste(list.files(sffDir,'\\.sff$',full.names=TRUE),collapse=' ')

	#can't let sfffile break this up for us because some people have overlapping barcodes among samples
	if(!all(unique(samples) %in% names(barcodes)))stop(simpleError('Please give barcodes for every sample'))
	for(i in unique(samples)){
		#outFile<-sprintf('%s/%s%s.sff',outDir,baseName,i)
		#don't need to include sample since sfffile will do it automatically from the barcode file
		outFile<-sprintf('%s/%s.sff',outDir,baseName)
		finalOutFile<-sprintf('%s/%s.%s.sff',outDir,baseName,i)
		capsOutFile<-sprintf('%s/%s.%s.sff',outDir,baseName,toupper(i))
		nameFile<-tempfile()
		selector<-samples==i
		writeLines(names[selector],nameFile)
		barcodeFile<-tempfile()
		thisBarcodes<-barcodes[[i]]
		barLines<-paste(sprintf('mid = "%s", "%s",2;',i,thisBarcodes),collapse='\n')
		writeLines(c('GSMIDs','{',barLines,'}'),barcodeFile)
		cmd<-sprintf('%ssfffile -i %s -o %s -mcf %s -s %s',sffBinDir,nameFile,outFile,barcodeFile,sffFiles)
		if(vocal)message(cmd)
		returnCode<-system(cmd,intern=TRUE)
		if(length(grep('reads written into the SFF file',returnCode))==0)stop(simpleError('Some error in sfffile'))
		if(length(grep('reads written into the SFF file',returnCode))!=1)stop(simpleError('Too many sff files made'))
		thisRegex<-sprintf(' *%s: *([0-9]+) reads written into the SFF file\\.',toupper(i))
		numberWritten<-returnCode[grep(thisRegex,returnCode)]
		numberWritten<-gsub(thisRegex,'\\1',numberWritten)
		if(sum(selector)!=numberWritten)browser()#stop(simpleError('Reads written do not match selected reads'))
		cmd<-sprintf('%ssffinfo -a %s',sffBinDir,capsOutFile)
		if(vocal)message(' ',cmd)
		if(length(system(cmd,intern=TRUE))!=sum(selector))stop(simpleError('sffinfo gives different number of reads from selected'))
		file.rename(capsOutFile,finalOutFile)
		unlink(barcodeFile)
		unlink(nameFile)
	}
	return(TRUE)
}

#reads bed file
#returns list with a dataframe (columns chr,start,end) for each track
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




#reg: region in the format "chrX:123545-123324"
#files: bam files to pull the counts from
#bam2depthBinary: bam2depth executable file
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
		filler<-cover[rep(1,2),]
		filler$chr<-region$chr
		filler[,countCols]<-0
		filler$pos<-unlist(region[,c('start','end')])
		if(!region$start %in% cover$pos)cover<-rbind(filler[1,],cover)
		if(!region$end %in% cover$pos)cover<-rbind(cover,filler[2,])
		cover<-fillZeros(cover)
	}
	return(cover)
}

#writes a fasta file
#names: the > line
#dna: the sequences
#fileName: file to write to
#addBracket: add > to start of names?
write.fa<-function(names,dna,fileName,addBracket=FALSE,isGz=grepl('.gz$',fileName)){
	if(addBracket|any(grep('^[^>]',names)))names<-paste('>',names,sep='')
	dna<-sub(' +$','',dna,perl=TRUE)
	output<-paste(names,dna,sep="\n")
	if(isGz)fileName<-gzfile(fileName)
	writeLines(output,sep="\n",con=fileName)
	if(isGz)close(fileName)
}

#writes a fastq file
#names: sequence names for @/+ line
#seqs: the sequences
#quals: the qualities
#fileName: file to write to
write.fastq<-function(names,seqs,quals,fileName,isGz=grepl('.gz$',fileName)){
	names1<-paste('@',names,sep='')
	names2<-paste('+',names,sep='')
	output<-paste(names1,seqs,names2,quals,sep="\n")
	if(isGz)fileName<-gzfile(fileName)
	writeLines(output,sep="\n",con=fileName)
	if(isGz)close(fileName)
}


