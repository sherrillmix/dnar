nimblegenPrimers<-c('CTCGAGAATTCTGGATCCTC','GAGGATCCAGAATTCTCGAGTT')
primers454<-c("GCCTCCCTCGCGCCATCAG","GCCTTGCCAGCCCGCTCAG")

#convenience function for stop(simpleError())
stopError<-function(...){
	stop(simpleError(paste(...,sep='')))
}


#pastes characters into string
#chars: an array of characters
#returns: string
c2s<-function (chars){
	return(paste(chars, collapse = ""))
}

#breaks string into characters
#string: a character string
#returns: an array of characters
s2c<-function (string){
	if (is.character(string) && length(string) == 1) {
		return(strsplit(string,'')[[1]])
	} else {
		warning("Wrong argument type in s2c(), NA returned")
		return(NA)
	}
}

checkOverlap<-function(starts,ends,tStarts,tEnds,tNames,allCover=FALSE,allCoverFuzz=0){
	overlapNames<-apply(cbind(starts,ends),1,function(x,tStarts,tEnds,tNames){
		if(!allCover)thisNames<-tNames[x[1]<=tEnds&x[2]>=tStarts]
		else thisNames<-tNames[tStarts-allCoverFuzz<=x[1]&tEnds+allCoverFuzz>=x[2]]
		if(length(thisNames)==0)return('')
		else return(paste(thisNames,collapse='|'))
	},tStarts,tEnds,tNames)
	return(overlapNames)
}


matchBlocks<-function(starts,ends,tStarts,tEnds){
	if(length(starts)!=length(ends))stop(simpleError('starts and ends not same length'))
	if(length(tStarts)!=length(tEnds))stop(simpleError('tStarts and tEnds not same length'))
	bases<-as.vector(unlist(apply(cbind(starts,ends),1,function(x)x[1]:x[2])))
	tBases<-as.vector(unlist(apply(cbind(tStarts,tEnds),1,function(x)x[1]:x[2])))
	tBases<-sort(unique(tBases))
	index<-tBases %in% bases
	if(any(index))	tOverlap<-mean(index[min(which(index)):max(which(index))])
	else tOverlap<-0
	bases<-sort(unique(bases))
	index<-bases %in% tBases
	overlap<-mean(index)
	return(c(overlap,tOverlap))
}


checkCover<-function(tStart,tEnd,starts,ends){
	if(length(starts)!=length(ends))stop(simpleError('Starts and ends not same length'))
	goodStarts<-findReads(tStart,starts,ends-starts+1,tEnd)
	starts<-starts[goodStarts]
	ends<-ends[goodStarts]
	output<-rep(0,tEnd-tStart+1)
	names(output)<-tStart:tEnd
	tmp<-apply(cbind(starts,ends),1,function(x,start,end){selector<-x[1]:x[2];selector<-paste(selector[selector<=end&selector>=start]);output[selector]<<-output[selector]+1},tStart,tEnd)
	return(output)
}

findReads<-function(low,starts,lengths,high=low){
	if(length(starts)!=length(lengths))stopError('Length of starts and lengths not equal')
	output<-rep(TRUE,length(starts))
	output[starts>high]<-FALSE
	output[starts+lengths<low]<-FALSE
	return(output)
}


#changed argument order may need adjusted
checkSO<-'~/scripts/R/c/checkCover.so'
loader<-try(dyn.load(checkSO),TRUE)
if (any(grep("Error",loader))){
	checkCoverage<-function(starts,lengths,totalNumBases=max(starts+lengths),range=FALSE,coverMin=0){
		if(length(starts)!=length(lengths))stop(simpleError('Starts and lengths not same length'))
		if(any(starts+lengths-1>totalNumBases))stop(simpleError('totalNumBases < starts + lengths'))
		if(!range){
			output<-rep(0,totalNumBases)
			counter<-1
			tmp<-apply(cbind(starts,lengths),1,function(x){if(counter%%1000==0)message(counter);counter<<-counter+1;output[x[1]:(x[1]+x[2]-1)]<<-output[x[1]:(x[1]+x[2]-1)]+1})
		}else{
			holder<-c()
			counter<-1
			tmp<-apply(cbind(starts,lengths),1,function(x){
				if(counter%%1000==0)message(counter)
				counter<<-counter+1
				indices<-paste(x[1]:(x[1]+x[2]-1))
				set<-indices %in% names(holder)
				holder[indices[set]]<<-holder[indices[set]]+1
				holder[indices[!set]]<<-1
			})
			holder<-as.numeric(names(holder)[holder>coverMin])
			output<-index2range(holder)
		}
		return(output)
	}
}else{
	checkCoverage <- function(starts, lengths, outLength=max(starts+lengths-1)) {
		if(length(starts)!=length(lengths))stop(simpleError('starts and lengths different lengths'))
		if(any(starts+lengths-1>outLength))stop(simpleError('starts+length-1 greater than outLength'))
		ans<-.C('checkCover',as.integer(rep(0,outLength)),as.integer(starts),as.integer(starts+lengths-1),as.integer(length(starts)))	
		return(ans[[1]])
	}
}

gap2NoGap<-function(gapSeq,coords){
	gapSeqSplit<-strsplit(gapSeq,'')[[1]]
	nonDash<-!gapSeqSplit %in% c('*','.','-')
	newCoords<-cumsum(nonDash)
	return(newCoords[coords])
}

noGap2Gap<-function(gapSeq,coords){
	gapSeqSplit<-strsplit(gapSeq,'')[[1]]
	nonDash<-which(!gapSeqSplit %in% c('.','*','-'))
	return(nonDash[coords])
}

binary2range<-function(index){
	return(index2range(which(index)))
}

index2range<-function(index){
	index<-sort(unique(index))
	diffs<-c(diff(index),1)
	ends<-c(which(diffs>1),length(index))
	starts<-c(1,ends[-length(ends)]+1)
	return(data.frame('start'=index[starts],'end'=index[ends]))
}

#reverse compliment dna
#dna: sequence
#returns: reverse complimented dna
revComp<-function(dna){
	revdna<-c2s(rev(s2c(dna)))
	revdna<-chartr('TGAC][','ACTG[]',revdna)
	return(revdna)
}

read.fastq<-function(fileName,convert=TRUE){
	#assuming no comments and seq and qual on a single line each
	#assuming any line starting with @ 2 lines later by + is the block and no extra chars (who designed this format?)
	#as.integer(charToRaw())
	x<-readLines(fileName)
	plusLines<-grep('^\\+',x)
	atLines<-grep('^@',x)
	plusLines<-plusLines[plusLines %in% (atLines+2)]
	atLines<-atLines[atLines %in% (plusLines-2)]
	if(any(grep('[^ACTGN]',x[atLines+1])))warning('Non ATCGN chars found in sequence')
	if(length(plusLines)!=length(atLines))stopError('Problem finding @ + lines')
	output<-data.frame('name'=sub('^@','',x[atLines]), 'seq'=x[atLines+1], 'qual'=x[atLines+3],stringsAsFactors=FALSE)
	if(any(nchar(output$seq)!=nchar(output$qual)))stopError('Sequence and qual lengths do not match')
	if(convert)output$qual<-unlist(lapply(x[atLines + 3],function(x)paste(as.integer(charToRaw(x))-33,collapse=' ')))
	return(output)
}

#converts ambigous dna to an appropriate regular expression
#input sequence with ambigous bases
#returns: regular expression
ambigous2regex<-function(dna){
	dna<-toupper(dna)
	ambig<-c(
		'R'='AG',
		'Y'='CT',
		'M'='AC',
		'K'='GT',
		'S'='CG',
		'W'='AT',
		'H'='ACT',
		'B'='CGT',
		'V'='ACG',
		'D'='AGT',
		'N'='ACGT'
	)
	for (i in names(ambig)){
		dna<-gsub(i,paste('[',ambig[i],']',sep=''),dna)
	}
	return(dna)
}
#read a fasta file
#fileName:name of file
#returns: dataframe with columns name (name between > and the first ' '), seq (sequence), and longName (the whole > line)
read.fa<-function(fileName,longNameTrim=TRUE){
	x<-readLines(fileName)
	if(length(x)==0)return(NULL)
	x<-x[!1:length(x) %in% c(grep('^#',x,perl=TRUE),grep('^$',x,perl=TRUE))]
	y<-paste(x,collapse="\n")
	splits<-strsplit(y,'>',fixed=TRUE)[[1]][-1]
	splits2<-strsplit(splits,"\n",fixed=TRUE)
	output<-lapply(splits2,function(x){return(c(x[1],paste(x[-1],collapse='')))})
	output<-as.data.frame(do.call(rbind,output),stringsAsFactors=FALSE)
	colnames(output)<-c('longName','seq')
	if(longNameTrim){
		output$name<-unlist(lapply(strsplit(output$longName,' ',fixed=TRUE),function(x)x[1]))
		output<-output[,3:1]
	}else{
		colnames(output)<-c('name','seq')	
	}
	return(output)
}
#writes a fasta file
#names: the > line
#dna: the sequences
#fileName: file to write to
#addBracket: add > to start of names?
write.fa<-function(names,dna,fileName,addBracket=FALSE){
	if(addBracket|any(grep('^[^>]',names)))names<-paste('>',names,sep='')
	output<-paste(names,dna,sep="\n")
	writeLines(output,sep="\n",con=fileName)
}

readBlat<-function(fileName){
	x<-read.table(fileName,skip=5,sep="\t",stringsAsFactors=FALSE,colClasses=c(rep('numeric',8),rep('character',2),rep('numeric',3),'character',rep('numeric',4),rep('character',3)))
	colnames(x)<-c('match','mismatch','repmatch','ns','qGaps','qGapBases','tGaps','tGapBases','strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blocks','blockSizes','qStarts','tStarts')
	#Score equation from blat's webpage
	x$score<-x$match-x$mismatch-x$qGaps-x$tGaps
	return(x)
}

trimEnd<-function(seqs,revCompPrimer,trimmed=rep(FALSE,length(seqs)),minSubstring=8){
	if(length(seqs)!=length(trimmed)) stop(simpleError('Already trimmed vector and seqs vector not same length'))
	trimSeq<-seqs
	for(i in nchar(revCompPrimer):minSubstring){
			  thisRegex<-paste(substr(revCompPrimer,1,i),'$',sep='')
			  selector<-!trimmed & 1:length(seqs) %in% grep(thisRegex,seqs)
			  trimSeq[selector]<-gsub(thisRegex,'',trimSeq[selector])
			  trimmed<-trimmed|selector
			  message('Found ',sum(selector),' seqs matching ',thisRegex)
	}
	message('Trimmed ',sum(selector),' reads out of ',length(seqs))
	return(list(trimSeq,trimmed))
}

parseAce<-function(aceFile,dropMosaik=TRUE,checkSnps=TRUE){
	#debug<-FALSE
	#if(!exists(x)|!debug)x<-readLines(aceFile)
	x<-readLines(aceFile)
	asLines<-grep('^AS [0-9]+ [0-9]+ *$',x,perl=TRUE)
	if(length(asLines)!=1)stopError('Incorrect number of AS lines found')
	asLine<-strsplit(x[asLines],' ')[[1]]
	numContigs<-asLine[2]
	numReads<-asLine[3]
	if(numContigs!=1)stopError('Sorry this function only handles 1 contig')
	message('Expecting ',numReads,' reads')
	coLines<-grep('^CO [^ ]+ [0-9]+ [0-9]+ [0-9]+ [UC] *$',x,perl=TRUE)
	if(length(coLines)!=1)stopError('Incorrect number of CO lines found')
	coLine<-strsplit(x[coLines],' ')[[1]]
	contigName<-coLine[2]
	contigBaseNum<-coLine[3]
	contigReadNum<-coLine[4]
	message('Contig ',contigName,' has ',contigBaseNum,' bases and ',contigReadNum,' reads')
	bqLines<-grep('^BQ *$',x)
	if(length(bqLines)!=1)stopError('Incorrect number of BQ lines found')
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
		if(length(snpCT)!=length(ctLine)|length(grep('pSnp=',psnpLine))!=length(psnpLine))stopError('psnp lines not the same length as PB++ CT lines')
		message('Found ',length(ctLine),' CT lines of which ',length(snpCT),' were from PB++')
		ct<-strsplit(ctLine,' ')
		notes<-as.data.frame(do.call(rbind,lapply(ct,function(x)return(x[c(4,5,1)]))),stringsAsFactors=FALSE)
		colnames(notes)<-c('start','stop','contig')
		notes$psnp<-sub('pSnp=','',psnpLine)
	}else{
		message('Not looking for snp notes')
		notes<-c()
	}
	if(dropMosaik&(reads[1,'name']=='.MosaikAnchor.C1'|reads[1,'name']=='.MosaikReference')){
		reads<-reads[-1,]
		message('Removing .MosaikAnchor.C1')
	}
	return(list('contig'=contigSeq,'reads'=reads,'notes'=notes))
	#readDir<-lapply(afLine,function(x)return(x[3]))
	#readStart<-lapply(afLine,function(x)return(x[4]))
	#reads<-data.frame('name'=readNames,'dir'=readDir,'start'=readStart)
}

cutReads<-function(seqs,starts,low,high,lengths=nchar(seqs)){
	debug<-TRUE
	goodReads<-findReads(low,starts,lengths,high)
	if(!any(goodReads))return(FALSE)
	thisSeqs<-seqs[goodReads];	starts<-starts[goodReads];	lengths<-lengths[goodReads]
	cutlow<-low-starts+1
	startDash<-0-cutlow+1
	startDash[startDash<0]<-0
	cutlow[cutlow<1]<-1
	cuthigh<-high-starts+1
	endDash<-high-starts-lengths+1
	endDash[endDash<0]<-0
	cuts<-substr(thisSeqs,cutlow,cuthigh)
	#make a string of dots for cutting
	dots<-paste(rep('.',100000),collapse='')
	predots<-substring(dots,1,startDash)
	postdots<-substring(dots,1,endDash)
	cuts<-paste(predots,cuts,postdots,sep='')
	return(list(cuts,goodReads))
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

seq2flow<-function(seq,flowOrder=c('T','A','C','G'),outputLength=NULL){
	seqSplit<-strsplit(seq,'')[[1]]
	dif<-seqSplit!=c(seqSplit[-1],'DUMMY')
	strLengths<-diff(c(0,which(dif)))
	chars<-seqSplit[dif]
	nextChars<-c(chars[-1],NA)
	numFlows<-length(flowOrder)
	distMat<-matrix(NA,nrow=numFlows,ncol=numFlows,dimnames=list(flowOrder,flowOrder))
	distMat[,1]<-length(flowOrder):1
	for(i in 2:numFlows)distMat[,i]<-c(distMat[numFlows,i-1],distMat[-numFlows,i-1])
	dists<-rep(NA,length(chars))
	for(i in flowOrder){
		dists[chars==i&is.na(nextChars)]<-NA
		dists[chars==i&!is.na(nextChars)]<-distMat[i,nextChars[chars==i&!is.na(nextChars)]]
	}
	dists<-c(which(flowOrder==chars[1]),dists[-length(dists)])
	if(is.null(outputLength))outputLength<-sum(dists)
	output<-rep(0,outputLength)
	output[cumsum(dists)]<-strLengths
	names(output)<-rep(flowOrder,length.out=length(output))
	return(output)
}

flow2seq<-function(flow,flowOrder=c('T','A','C','G')){
	chars<-rep(flowOrder,length.out=length(flow))
	output<-paste(rep(chars,round(flow)),collapse='')
	return(output)
}

#reads bed file
#returns list with a dataframe (columns chr,start,end) for each track
read.bed<-function(fileName){
	x<-readLines(fileName)
	tracks<-grep('track',x)
	if(length(tracks)==0)tracks<-c(0)
	message('Found ',length(tracks),' tracks')
	trackNames<-gsub('.*name=([^ ]+).*','\\1',x[tracks])
	tracks<-c(tracks,length(x)+1)
	output<-list()
	for(i in 1:(length(tracks)-1)){
		output[[trackNames[i]]]<-data.frame(do.call(rbind,strsplit(x[(tracks[i]+1):(tracks[i+1]-1)],'\t')),stringsAsFactors=FALSE)
		colnames(output[[trackNames[i]]])<-c('chr','start','end')
		output[[trackNames[i]]][,'start']<-as.numeric(output[[trackNames[i]]][,'start'])
		output[[trackNames[i]]][,'end']<-as.numeric(output[[trackNames[i]]][,'end'])
	}
	return(output)
}


