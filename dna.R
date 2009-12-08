nimblegenPrimers<-c('CTCGAGAATTCTGGATCCTC','GAGGATCCAGAATTCTCGAGTT')
primers454<-c("GCCTCCCTCGCGCCATCAG","GCCTTGCCAGCCCGCTCAG")
primerTitanium<-c('CCATCTCATCCCTGCGTGTCTCCGACTCAG','CCTATCCCCTGTGTGCCTTGGCAGTCTCAG')

#convenience function for stop(simpleError())
stopError<-function(...){
	stop(simpleError(paste(...,sep='')))
}

#convenience function for selecting elements from a matrix
indexMatrix<-function(x,y,mat){
	if(!is.integer(x)){tmp<-1:nrow(mat);names(tmp)<-rownames(mat);x<-tmp[x]}
	if(!is.integer(y)){tmp<-1:ncol(mat);names(tmp)<-colnames(mat);y<-tmp[y]}
	if(length(x)!=length(y)|max(x)>nrow(mat)|max(y)>ncol(mat))stop(simpleError("Dimensions don't match up"))
	index<-(y-1)*nrow(mat)+x
	return(mat[index])
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

#calculate shannon diversity of a vector
shannon<-function(x,base=exp(1)){
	return(-sum(x/sum(x)*log(x/sum(x),base)))
}




#calculate heights for a 'weblogo' like plot from a matrix of rows acgt, cols base position
#http://en.wikipedia.org/wiki/Weblogo
calcWebLogo<-function(baseMat,num=rep(9999,ncol(baseMat))){
	baseMat<-apply(baseMat,1,function(x)x/sum(x))
	shannon<-apply(rbind(baseMat,num),1,function(x){x[1:4]*(2-(shannon(x[1:4],base=2)+3/2/ln(2)/x[5]))})

########WORK HERE

}

#calculates chao diversity index
#counts: a vector of counts with one entry per "species"
#returns: chao index
chao<-function(counts){
	return(length(counts)+sum(counts==1)*(sum(counts==1)-1)/2/(sum(counts==2)+1))
}

#calculate bootstrapped rarefactions 
#species: ids of species
#counts: corresponding counts of species
#samples: vector of numbers of draws to sample at
#reps: how many random samples to take at each step
#quants: quantiles to return
#chaoAdjust: calculate chao-predicted species number on each random draw
#returns: list containing dataframe of calculated quantiles and vector of samples argument
#chao<-tapply(rep(ace[[2]]$otu2,ace[[2]]$num),rep(ace[[2]]$ampPat2,ace[[2]]$num),function(x){y<-table(x);return(length(y)+sum(y==1)*(sum(y==1)-1)/2/(sum(y==2)+1))})
rarefy<-function(species,counts=rep(1,length(species)),samples=seq(10,sum(counts),10),reps=10000,quants=c(.5,.025,.975),chaoAdjust=FALSE,debug=FALSE,replaceSpecies=FALSE,minCount=0){
	if(length(species)!=length(counts))stop(simpleError('Length of species and counts not equal'))
	if(minCount>0){
		speciesCounts<-tapply(counts,species,sum)
		badSpecies<-names(speciesCounts)[speciesCounts<minCount]
		message('Bad species:',paste(badSpecies,collapse=','))
		counts<-counts[!species %in% badSpecies]
		species<-species[!species %in% badSpecies]
	}
	species<-rep(species,counts)
	if(debug)message('Number of samples: ',length(samples),' Last sample:',samples[length(samples)])
	output<-lapply(samples,function(sample,reps,species,debug){
		if(debug)message('Sample ',sample,' started')
		numSpecies<-sapply(1:reps,function(rep,species,sample,chaoAdjust){
			thisSpecies<-sample(species,sample,replace=replaceSpecies)
			if(chaoAdjust){
				return(chao(table(thisSpecies)))	
			} else return(length(unique(thisSpecies)))
		},species,sample,chaoAdjust)	
		estimate<-quantile(numSpecies,quants)
		return(estimate)
	},reps,species,debug)
	output<-do.call(rbind,output)
	return(list(output,samples))
}

#calculate rarefaction using formula
#sample: vector of numbers of individuals per species
#step: size of sampling steps for rarefaction
quickRare<-function(sample,step=10){
	sampleSize<-20;
	steps<-unique(c(seq(step,sum(sample),step),sum(sample)))
	output<-sapply(steps,function(x)rareEquation(sample,x))
	return(cbind(output,steps))

}

rareEquation<-function(sample,sampleSize){
	#numbers too big
	#output2<-length(sample)-choose(sum(sample),sampleSize)^-1*sum(choose(sum(sample)-sample,sampleSize))
	#message(output2)
	#no way to log sum 
	#logSum<-log(sum(choose(sum(sample)-sample,sampleSize)))
	#output<-length(sample) - exp(- lchoose(sum(sample),sampleSize) + logSum)
	output<-sum(1-exp(lchoose(sum(sample)-sample,sampleSize)-lchoose(sum(sample),sampleSize)))
	
	if(is.na(output)|is.infinite(output))browser()
	return(output)
}


#convert dna string into seperate codons
#dna: single string of DNA
#frame: starting frame (0=start on first base, 1=on second, 2=on third)
dna2codons<-function(dna,frame=0){
	if(nchar(dna)<3){
		warning('DNA less than 3 bases long')
		return(NULL)
	}
	frame<-frame+1
	starts<-seq(frame,nchar(dna)-2,3)
	##############WORK HERE
	return(substring(dna,starts,starts+2))
}

#convert codon to amino acid (using aminoAcids table below)
#codons: vector of 3 base codons
#type: amino acid info to return; code for single letter, name for full name, or abbr for 3-letter abbreviation
codon2aa<-function(codons,type='code'){
	if(!type %in% c('code','name','abbr'))stop(simpleError('Invalid amino acid type'))
	codons<-gsub('T','U',toupper(codons))
	return(aminoAcids[,type,drop=FALSE][codons,1])
}

#convert dna/rna to amino acids
#dna: a string of DNA/RNA
#frame: starting frame (0=start on first base, 1=on second, 2=on third)
#debug: print debug info?
dna2aa<-function(dna,frame=0,debug=FALSE){
	codons<-dna2codons(dna,frame)	
	if(debug)print(codons)
	output<-paste(codon2aa(codons),collapse='')
	return(output)
}


#find a single codon at a given position in dna
#dna: a string of DNA/RNA
#start: start coordinate of exon
#frame: starting frame (0=start on first base, 1=on second, 2=on third)
#strand: strand of dna (i.e. revcomp the dna first if '-')
dnaPos2aa<-function(dna,pos,start=1,end=nchar(dna),frame=0,strand='+',refStart=1,debug=FALSE){
	if(length(pos)==1&length(start)>1)pos<-rep(pos,length(start))
	pos[pos>end|pos<start|frame==-1]<--99999999
	start<-start-refStart+1
	end<-end-refStart+1
	pos<-pos-refStart+1
	#adjust for frame
	end[strand=='-']<-end[strand=='-']-(3-frame[strand=='-'])%%3
	start[strand!='-']<-start[strand!='-']+(3-frame[strand!='-'])%%3
	#adjust end and start for outisde dna 
	#print(c(start,end))
	start[strand=='-'&start<1]<-1
	end[strand!='-'&end>nchar(dna)]<-nchar(dna[strand!='-'&end>nchar(dna)])
	start[strand!='-'&start<1]<-1+(3-(1-start[strand!='-'&start<1]))%%3
	selector<-strand=='-'&end>nchar(dna)
	end[selector]<-nchar(dna[selector])-(3-(end[selector]-nchar(dna[selector])))%%3

	dna<-substring(dna,start,end)
	if(any(nchar(dna)!=end-start+1&nchar(dna)!=0))stop(simpleError("Missing DNA"))
	dna[strand=='-']<-revComp(dna[strand=='-'])
	pos[strand=='-']<-end[strand=='-']-pos[strand=='-']+1
	pos[strand!='-']<-pos[strand!='-']-start[strand!='-']+1
	startPos<-floor((pos-1)/3)*3+1
	if(debug)message(paste('StartPos: ',startPos,' EndPos:',startPos+2,' Cut:',substring(dna,startPos,startPos+2),sep='',collapse='\n'))
	if(debug)message(paste('DNA: ',dna,' pos: ',pos,' start: ',start,' end: ',end,' strand-:',strand=='-',sep='',collapse='\n'))
	return(codon2aa(substring(dna,startPos,startPos+2)))
		
}


#check overlap between two sets of coordinates
#starts: start coordinates
#ends: end coordinates
#tStarts: start coordinates of target
#tEnds: end coordinates of target
#tNames: target names 
#allCover: if TRUE entire start and end of base must fall within targets +- allCoverFuzz
#allCoverFuzz: extra overlap to consider on each end of target when using allCover
#returns: | seperated vector of tNames with overlap ('' if no overlapping target)
checkOverlap<-function(starts,ends,tStarts,tEnds,tNames,allCover=FALSE,allCoverFuzz=0,sep='|'){
	overlapNames<-apply(cbind(as.numeric(starts),as.numeric(ends)),1,function(x,tStarts,tEnds,tNames){
		if(!allCover)thisNames<-tNames[x[1]<=tEnds&x[2]>=tStarts]
		else thisNames<-tNames[tStarts-allCoverFuzz<=x[1]&tEnds+allCoverFuzz>=x[2]]
		if(length(thisNames)==0)return('')
		else return(paste(unique(thisNames),collapse=sep))
	},as.numeric(tStarts),as.numeric(tEnds),tNames)
	return(overlapNames)
}

checkOverlapMulti<-function(starts,ends,tStarts,tEnds,tNames,qChrom,tChrom,vocal=FALSE,...){
	results<-rep(NA,length(qChrom))
	for(i in unique(qChrom)){
		if(vocal)message('Working on ',i)
		selector<-qChrom==i
		tSelector<-tChrom==i
		if(any(selector)&any(tSelector))results[selector]<-checkOverlap(starts[selector],ends[selector],tStarts[tSelector],tEnds[tSelector],tNames[tSelector],...)
	}
	return(results)
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
#alternative check coverage (better for sparse coverage in high numbers of bases)
#starts:starts of coverage ranges
#ends:ends of coverage ranges
checkCoverage2<-function(starts,ends){
	cover<-unlist(apply(cbind(starts,ends),1,function(x)x[1]:x[2]))
	output<-table(cover)
	output<-data.frame('pos'=as.numeric(names(output)),'cover'=as.vector(output),stringsAsFactors=FALSE)
	return(output)
}

#starts:starts of coverage ranges
#ends:ends of coverage ranges
startStop2Range<-function(starts,stops){
	cover<-unique(unlist(mapply(function(x,y)x:y,starts,stops)))
	ranges<-index2range(cover)
	return(ranges)
}

#read in a bunch of fasta files in a target directory
#dir: target directory
#suffix: regex to select file names
readFaDir<-function(dir='.',suffix='\\.fn?a$'){
	faFiles<-list.files(dir,suffix)
	for(i in faFiles){
		tmp<-read.fa(sprintf('%s/%s',dir,i))
		tmp$file<-i
		if(exists('allFa'))allFa<-rbind(allFa,tmp)
		else allFa<-tmp
	}
	return(allFa)
}
	

#make a .bedgraph file for use on UCSC browser
#fileName: name of file
#chroms: chromosomes of coverage
#starts: start coordinates of coverage
#ends: end coordinates of coverage
#values: values for each start-end range (or NULL to calculate counts) (if not NULL should be at most one value per individual base)
#header: string to add to .bedgraph header e.g. name="" description="" visibility=full 
#vocal: turn on various progress messages
#autoScale: add autoscale=off to header and calculate maximum (or Ymax) for entire data
#proportion: divide all values by this value (e.g. number of reads in sample)
#yMin: minimum value for y scale if autoScale
#yMax: maximum value for y scale if autoScale (if NULL find max in data)
#side effect: writes .bedgraph to fileName
#returns: Calculated dataframe of chrom, start, end, value
makeBedGraph<-function(fileName,chroms,starts,ends,values=NULL,header='',vocal=FALSE,autoScale=TRUE,proportion=NULL,yMin=0,yMax=NULL){
	if(length(chroms)!=length(starts)|length(chroms)!=length(ends))stop(simpleError("chroms, starts, ends not equal length"))
	if(is.null(values)){
		if(vocal)message('Calculating counts')
		uniqueChroms<-unique(chroms)
		for(chrom in uniqueChroms){
			if(vocal)message('Calculating cover')
			cover<-checkCoverage2(starts[chroms==chrom],ends[chroms==chrom])
			if(vocal)message('Calculating ranges')
			ranges<-tapply(cover$pos,cover$cover,index2range)
			rangeNums<-lapply(ranges,nrow)
			tmp<-do.call(rbind,ranges)
			colnames(tmp)<-c('starts','ends')
			tmp$values<-as.numeric(rep(names(ranges),rangeNums))
			tmp$chroms<-chrom
			if(chrom==uniqueChroms[1])data<-tmp
			else data<-rbind(data,tmp)
		}
	}else{
		data<-cbind(chroms,starts,ends,values)
	}
	data<-data[,c('chroms','starts','ends','values')]
	colnames(data)<-c('chrom','start','end','value')
	if(!is.null(proportion))data$value<-data$value/proportion
	#make sure we don't get any 6.231e-09
	if(is.null(yMax))yMax<-max(data$value)
	data$value<-formatC(data$value,format='fg')
	if(vocal)message('Outputting to ',fileName)
	if(!grepl('^track type=',header))header<-sprintf('track type=bedGraph %s',header)
	if(autoScale)header<-sprintf('%s viewLimits=%s:%s autoScale=off',header,formatC(yMin,format='fg'),formatC(yMax,format='fg'))
	output<-c(header,paste(data$chrom,format(data$start,scientific=FALSE),format(data$end,scientific=FALSE),data$value,sep='\t'))
	writeLines(output,fileName)
	return(data)
}

#write psl file for ucsc genome browser
#blat: dataframe containing columns 
#file: file to be writted
#header: header line for psl
write.psl<-function(blat,file,header=''){
	if(!grepl('^track ',header))header<-sprintf('track %s',header)
	pslCols<-c(	 'match', 'mismatch', 'repmatch', 'ns', 'qGaps', 'qGapBases', 'tGaps', 'tGapBases', 'strand', 'qName', 'qSize', 'qStartBak', 'qEnd', 'tName', 'tSize', 'tStartBak', 'tEnd', 'blocks', 'blockSizes', 'qStarts', 'tStarts')
	pslSelector<-pslCols %in% colnames(blat)
	if(any(!pslSelector))stop(simpleError('Missing columns: ',paste(pslCols[!pslSelector],collapse=', ')))
	writeLines(header,file)
	write.table(blat[,pslCols],file,append=TRUE,quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
}


#convert gapped coordinates to what the coordinates would be without gaps
#gapSeq: the sequence containing gaps
#coords: coordinates on the gapped gapSeq to be converted into equivalent nongap coordinatess
#returns: equivalent gapped coordinates
gap2NoGap<-function(gapSeq,coords){
	gapSeqSplit<-strsplit(gapSeq,'')[[1]]
	nonDash<-!gapSeqSplit %in% c('*','.','-')
	newCoords<-cumsum(nonDash)
	coords[coords<1|coords>length(newCoords)]<-NA
	return(newCoords[coords])
}

#convert ungapped coordinates to what the coordinates would be with gaps
#gapSeq: the sequence containing gaps
#coords: coordinates on the ungapped gapSeq to be converted into equivalent gap coordinatess
#returns: equivalent ungapped coordinates
noGap2Gap<-function(gapSeq,coords){
	gapSeqSplit<-strsplit(gapSeq,'')[[1]]
	nonDash<-which(!gapSeqSplit %in% c('.','*','-'))
	coords[coords<1|coords>length(nonDash)]<-NA
	return(nonDash[coords])
}

#find covered ranges in binary data
#index:logical index
#returns: dataframe with start and end of ranges
binary2range<-function(index){
	return(index2range(which(index)))
}

#find ranges with at least one coverage in numerical index data
#index: numeric indices
#returns: dataframe with start and end of ranges
index2range<-function(index){
	index<-sort(unique(index))
	diffs<-c(diff(index),1)
	ends<-c(which(diffs>1),length(index))
	starts<-c(1,ends[-length(ends)]+1)
	return(data.frame('start'=index[starts],'end'=index[ends]))
}

#reverse strings
#strings: vector of strings to be revered
reverseString<-function(strings){
	output<-sapply(strsplit(strings,''),function(x)paste(rev(x),collapse=''))	
	return(output)
}
#compliment dna 
#dnas: vector of sequences
complimentDna<-function(dnas,brackets=TRUE){
	finds<-'TGAC'
	replaces<-'ACTG'
	if(brackets){finds<-paste(finds,'[]',sep='');replaces<-paste(replaces,'][',sep='')}
	return(chartr(finds,replaces,dnas))
}

#reverse compliment dna
#dnas: vector of sequences
#returns: reverse complimented dna
revComp<-function(dnas){
	return(complimentDna(reverseString(dnas),TRUE))
}

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
	plusLines<-plusLines[plusLines %in% (atLines+2)]
	atLines<-atLines[atLines %in% (plusLines-2)]
	if(any(grep('[^ACTGN]',x[atLines+1])))warning('Non ATCGN chars found in sequence')
	if(length(plusLines)!=length(atLines))stopError('Problem finding @ + lines')
	output<-data.frame('name'=sub('^@','',x[atLines]), 'seq'=x[atLines+1], 'qual'=x[atLines+3],stringsAsFactors=FALSE)
	if(any(nchar(output$seq)!=nchar(output$qual)))stopError('Sequence and qual lengths do not match')
	if(convert)output$qual<-unlist(lapply(x[atLines + 3],function(x)paste(as.integer(charToRaw(x))-33,collapse=' ')))
	return(output)
}

#read solexa fastq file
#fileName: name of solexa fastq file
#convert: convert condensed quals to numeric quals?
#returns: dataframe with solexa info seq and quals
read.solexa<-function(fileName,convert=TRUE){
	x<-readLines(fileName)
	output<-do.call(rbind,strsplit(x,':'))
	output<-cbind(output,do.call(rbind,strsplit(output[,5],'[#/]')))
	output<-as.data.frame(output[,-5])
	colnames(output)<-c('instrument','lane','tile','x','seq','rawQual','y','barcode','pair')
	if(convert)output$qual<-unlist(lapply(output$rawQual,function(x)paste(as.integer(charToRaw(x))-33,collapse=' ')))
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




#read a fasta file
#fileName:name of file
#returns: dataframe with columns name (name between > and the first ' '), seq (sequence), and longName (the whole > line)
read.fa<-function(fileName,longNameTrim=TRUE){
	x<-readLines(fileName)
	selector<-grep('^[^>].* .*[^ ]$',x,perl=TRUE)
	x[selector]<-paste(x[selector],' ',sep='')
	if(length(x)==0)return(NULL)
	x<-x[!grepl('^#',x,perl=TRUE)&x!='']
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
	output$seq<-gsub(' +$','',output$seq)
	return(output)
}
#alternative version of the above (a bit quicker)
read.fa2<-function(fileName,longNameTrim=TRUE){
	x<-readLines(fileName,warn=FALSE)
	if(length(x)==0)return(NULL)
	x<-x[!grepl('^#',x,perl=TRUE)&x!='']
	nameLines<-grep('^>',x)
	thisNames<-sub('^>','',x[nameLines])
	seqs<-apply(cbind(nameLines+1,c(nameLines[-1]-1,length(x))),1,function(coords){
		if(coords[1]<=coords[2])return(paste(x[coords[1]:coords[2]],collapse=''))
		else return('')
	})

	output<-data.frame('longName'=thisNames,'seq'=seqs,stringsAsFactors=FALSE)
	if(longNameTrim){
		output$name<-unlist(lapply(strsplit(output$longName,' ',fixed=TRUE),function(x)x[1]))
		output<-output[,3:1]
	}else{
		colnames(output)<-c('name','seq')	
	}
	return(output)
}

#return order of one vector in another
#query: values to be sorted in target order
#target: order for query to be sorted into
#strict: if true error if query not in target, if false append unknown queries to end of target
#...: arguments for order
orderIn<-function(query,target,strict=FALSE,...){
	if(any(!query %in% target)){
		if(strict){
			stop(simpleError('Query not in target'))
		}else{
			target<-c(target,unique(query[!query %in% target]))
		}
	}
	newOrder<-order(sapply(query,function(x)min(which(x==target))),...)
	return(newOrder)
}



#trim leading and trailing space characters
#x: vector of strings
trim<-function(x){
	sub('\\s+$','',sub('^\\s+','',x),perl=TRUE)
}

#loop through list, get unique names and make sure every element has those names
#x: list to loop through
#namesList: names to make sure every element has (and delete extras)
#fill: value to insert in missing elements
fillList<-function(x,namesList=unique(unlist(lapply(x,names))),fill=NA){
	output<-lapply(x,function(x){
		x<-x[names(x) %in% namesList]
		x[namesList[!namesList %in% names(x)]]<-fill
		return(x)
	})	
	return(output)
}

#parse a line of >ASDASD extra args test=1 test2=1 test3=asdasd asdasd 
#kind of slow
#assumes anything after = goes with the name= so put any extra arguments first
#nameLine: vector of lines
parseEqualLines<-function(nameLine,firstDel=TRUE){
	nameSplit<-strsplit(nameLine,' ',fixed=TRUE)
	output<-lapply(nameSplit,function(x){
		equals<-grep('=',x)	
		equalNames<-sub('=.*$','',x[equals])
		ends<-c(equals[-1]-1,length(x))
		values<-mapply(function(start,end){
			sub('^[^=]*=','',paste(x[start:end],collapse=' '))
		},equals,ends)
		names(values)<-equalNames
		return(values)
	})
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


#start a gfServer on an open port and return port
#nibDir: directory containing nib files to align against
#options: options to gfServer
#gfServer: path to gfServer program
#startPort: begin looking for open (no other R called gfServer using) ports on startPort and add 1 until finding one
#nibSuffix: nib files end in nibSuffix
#wait: number of 5 second intervals to wait for gfServer to finish starting before giving up
startBlat<-function(nibDir,options='',gfServer='gfServer',startPort=37900,nibSuffix='.nib',wait=120){
	port<-startPort
	counter<-1
	while(checkBlat(port)){
		port<-port+1
		if(counter>30)stop(simpleError('All 30 checked ports taken'))
		counter<-counter+1
	}
	system(sprintf('%s start localhost %d %s %s/*%s',gfServer,port,options,nibDir,nibSuffix),wait=FALSE)
	counter<-1
	while(system(sprintf('%s status localhost %d >/dev/null',gfServer,port),ignore.stderr=TRUE)!=0){
		Sys.sleep(5)
		if(counter>wait){
			killBlat(port)
			stopError('gfServer took longer than ',wait*5,' seconds to start')
		}
		counter<-counter+1
	}
	return(port)
}

#check if blat is running on port port
#port: port to check if blat is running on
checkBlat<-function(port){
	return(system(sprintf('pgrep -f "gfServer *start *localhost *%d">/dev/null',port))==0)
}

#run blat on port port
#faFile: path to fasta file for alignment
#port: port to run blat on
#gfClient: path to gfClient program
#gfClientOptions: arguments for gfClient
#...: arguments for startBlat
runBlat<-function(faFile,gfClientOptions='',outFile=gsub('\\.fn?a$','.blat',faFile),gfClient='gfClient',...){
	port<-startBlat(...)
	
	system(sprintf('%s localhost %d / %s %s %s',gfClient,port,faFile,gfClientOptions,outFile))

	killBlat(port)
}

#kill blat running on port port
#port: pkill blat running on port port
killBlat<-function(port){
	if(!checkBlat(port)){
		stopError('Blat does not appear to be running on port ',port)
	}
	code<-system(sprintf('pkill -f "gfServer *start *localhost *%d">/dev/null',port))
	if(checkBlat(port)){
		stopError('Could not kill blat on port ',port)
	}	
}


#read a blat file
#fileName: name of file
#skips: number of lines to skip (5 for a normal blat file)
#nrows: number of rows to read in (-1 for all)
#calcScore: calculate score column?
#returns: dataframe of blat data
readBlat<-function(fileName,skips=5,nrows=-1,calcScore=TRUE){
	x<-read.table(fileName,skip=skips,sep="\t",stringsAsFactors=FALSE,colClasses=c(rep('numeric',8),rep('character',2),rep('numeric',3),'character',rep('numeric',4),rep('character',3)),nrows=nrows)
	colnames(x)<-c('match','mismatch','repmatch','ns','qGaps','qGapBases','tGaps','tGapBases','strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blocks','blockSizes','qStarts','tStarts')
	#Score equation from blat's webpage
	if(calcScore)x$score<-x$match-x$mismatch-x$qGaps-x$tGaps
	return(x)
}

trimEnd<-function(seqs,revCompPrimer,trimmed=rep(FALSE,length(seqs)),minSubstring=8){
	if(length(seqs)!=length(trimmed)) stop(simpleError('Already trimmed vector and seqs vector not same length'))
	trimSeq<-seqs
	trimmedBak<-trimmed
	for(i in nchar(revCompPrimer):minSubstring){
			  thisRegex<-paste(substr(revCompPrimer,1,i),'$',sep='')
			  selector<-!trimmed & 1:length(seqs) %in% grep(thisRegex,seqs)
			  trimSeq[selector]<-gsub(thisRegex,'',trimSeq[selector])
			  trimmed<-trimmed|selector
			  message('Found ',sum(selector),' seqs matching ',thisRegex)
	}
	message('Trimmed ',sum(trimmed&!trimmedBak),' reads out of ',sum(!trimmedBak))
	return(list(trimSeq,trimmed))
}

#read an ace file
#LIMITATION: only reads 1 contig
#aceFile: string of file name or file handle
#dropMosaik: Remove extraneous line from Mosaik labelled either .MosaikAnchor.C1 or MosaikReference
#checkSnps: Find SNPs from pyroBayes
#returns: list of reference sequence in [[1]], aligned reads in [[2]] (note reads are not globally aligned still need to use start coordinate to place globally), snps in [[3]]
parseAce<-function(aceFile,dropMosaik=TRUE,checkSnps=TRUE,vocal=TRUE){
	#debug<-FALSE
	#if(!exists(x)|!debug)x<-readLines(aceFile)
	x<-readLines(aceFile)
	asLines<-grep('^AS [0-9]+ [0-9]+ *$',x,perl=TRUE)
	if(length(asLines)!=1)stopError('Incorrect number of AS lines found')
	asLine<-strsplit(x[asLines],' ')[[1]]
	numContigs<-asLine[2]
	numReads<-asLine[3]
	if(numContigs!=1)stopError('Sorry this function only handles 1 contig')
	if(vocal)message('Expecting ',numReads,' reads')
	coLines<-grep('^CO [^ ]+ [0-9]+ [0-9]+ [0-9]+ [UC] *$',x,perl=TRUE)
	if(length(coLines)!=1)stopError('Incorrect number of CO lines found')
	coLine<-strsplit(x[coLines],' ')[[1]]
	contigName<-coLine[2]
	contigBaseNum<-coLine[3]
	contigReadNum<-coLine[4]
	if(vocal)message('Contig ',contigName,' has ',contigBaseNum,' bases and ',contigReadNum,' reads')
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
#returns: list of aligned seqs in [[1]], logical vector of whether read fell within cut region in [[2]]
cutReads<-function(seqs,starts,low=min(starts),high=max(starts+nchar(seqs)-1),lengths=nchar(seqs)){
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

#take the output from blat and make continous reads out of it
#seqs: comma seperated target sequences from blat
#starts: comma seperated starting location for each piece of sequence 0 indexed
#fillStarts: start base for each output sequence 0 indexed
#fillEnds: total length for each output sequence
#gaps: list of matrices with columns tGaps and qGaps
#returns: dataframe of sequence
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

#qStarts: comma seperated starts for query seq
#tStarts: comma seperated starts for target seq
#blockSizes: comma seperated block sizes
blatFindGaps<-function(qStarts,tStarts,blockSizes){
	if(length(qStarts)!=length(tStarts)||length(blockSizes)!=length(qStarts))stop(simpleError('Lengths of starts and blocksizes not equal'))
	qStartList<-lapply(strsplit(qStarts,','),as.numeric)
	tStartList<-lapply(strsplit(tStarts,','),as.numeric)
	blockSizeList<-lapply(strsplit(blockSizes,','),as.numeric)
	blockNums<-sapply(qStartList,length)
	if(any(blockNums!=sapply(tStartList,length))||any(blockNums!=sapply(blockSizeList,length)))stop(simpleError('Comma seperated lists of starts and blocksizes not equal length'))
	gaps<-mapply(function(qStarts,tStarts,blockSizes){
		tEnds<-tStarts+blockSizes-1
		qEnds<-qStarts+blockSizes-1

		tGapLengths<-tStarts[-1]-tEnds[-length(tEnds)]-1
		qGapLengths<-qStarts[-1]-qEnds[-length(qEnds)]-1
		return(cbind('qGaps'=qGapLengths,'tGaps'=tGapLengths,'qGapStartAfter'=qEnds[-length(qEnds)],'tGapStartAfter'=tEnds[-length(tEnds)]))
	},qStartList,tStartList,blockSizeList)

	return(gaps)
}



#take coordinates of blat matches and split into a single for each exon
#chroms: Chromosome or other identifier
#names: Name of gene of other container of exons
#starts: Start coordinates of exons
#ends: End coordinates of exons if lengths is FALSE or length of exon if lengths is TRUE
#strands: Strand (for numbering exons in reverse on - strand)
#lengths: logical whether ends are end coordinates or lengths
#extraCols: a dataframe of extra columns (1 per batch of starts) to be added to the output
#extraSplits: a dataframe of extra comma-separated values (1 string of comma separated values per batch of starts, 1 value per start-stop pair) to be added to the output
#introns: also output introns (the spaces between exons)
#adjustStart: add adjustStart to starts (good for 0 index start, 1 index ends of UCSC)
#example: with(blat,blat2exons(tName,qName,tStarts,blockSizes,strand))
blat2exons<-function(chroms,names,starts,ends,strands=rep('+',length(names)),lengths=TRUE,extraCols=NULL,extraSplits=NULL,introns=FALSE,prefix='ex',adjustStart=0){
	if(any(c(length(chroms),length(names),length(strands),length(starts))!=length(ends)))stop(simpleError('Different lengths for chrom, strand, starts, lengths'))
	startsList<-strsplit(starts,',')
	if(adjustStart!=0)startsList<-lapply(startsList,function(x)as.numeric(x)+1)
	exonCounts<-exonCountsStrand<-sapply(startsList,length)
	endsList<-strsplit(ends,',')
	if(lengths)endsList<-mapply(function(x,y)as.numeric(x)+as.numeric(y)-1,endsList,startsList)
	if(any(exonCounts!=sapply(endsList,length)))stop(simpleError("Starts and ends not equal"))
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
	if(introns){
		inEnds<-sapply(startsList,function(x)as.numeric(x[-1])-1)	
		inStarts<-sapply(endsList,function(x)as.numeric(x[-length(x)])+1)	
		if(any(sapply(inEnds,length)!=sapply(inStarts,length)))stop(simpleError('Intron ends and starts not equal length'))
		problemIndex<-mapply(function(x,y)which(x<y),inEnds,inStarts)
		problems<-which(sapply(problemIndex,length)>0)
		for(problem in problems){
			for(thisIndex in problemIndex[[problem]]){
				#check if these are actually nogap neighbor exons
				if(inStarts[[problem]][thisIndex]-1!=inEnds[[problem]][thisIndex]){
					browser()
					stop(simpleError('Introns have negative length'))
				}
			}
			#if these weren't nogap neighbor exons we would have errored out so we can delete them
			inStarts[[problem]]<-inStarts[[problem]][-problemIndex[[problem]]]
			inEnds[[problem]]<-inEnds[[problem]][-problemIndex[[problem]]]
		}
		inCount<-sapply(inEnds,length)
		if(any(inCount!=sapply(inStarts,length)))stop(simpleError('Intron ends and starts not equal length after removing neighbors'))
		selector<-inCount>0
		introns<-blat2exons(chroms[selector],names[selector],sapply(inStarts[selector],paste,collapse=','),sapply(inEnds[selector],paste,collapse=','),strands[selector],FALSE,extraCols[selector,,drop=FALSE],NULL,FALSE,'in',adjustStart=0)
		introns$isIntron<-TRUE
		output$isIntron<-FALSE
		output<-rbind(output,introns)
	}
#this is too slow
#	thisExons<-apply(cbind(names,chroms,strands,starts,ends),1,function(x){
#		exonStarts<-as.numeric(strsplit(x['starts'],',')[[1]])
#		exonEnds<-as.numeric(strsplit(x['ends'],',')[[1]])
#		if(length(exonStarts)!=length(exonEnds))stop(simpleError('Exon starts and ends not equal length'))
#		if(lengths)exonEnds<-exonStarts+exonEnds-1
#		exonNum<-1:length(exonStarts)
#		if(x['strands']=='-')exonNum<-rev(exonNum)
#		thisNames<-paste(x['names'],'_ex',1:length(exonStarts),sep='')
#		return(data.frame('chr'=x['chroms'],'name'=thisNames,'start'=exonStarts,'end'=exonEnds,stringsAsFactors=FALSE,row.names=thisNames))
#	})
#	output<-do.call(rbind,thisExons)
	return(output)
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

#flow: vector of flowgram (e.g. produced by seq2flow)
#coords: bp of desired flow position
#return: flow number for each coord
indexFlow<-function(flow,coords){
	indices<-cumsum(flow)	
	output<-unlist(lapply(coords,function(x,y)return(min(which(y>=x))),indices))
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
		if(ncol(output[[trackNames[i]]])>3)thisNames<-c(thisNames,'name')
		colnames(output[[trackNames[i]]])<-thisNames
		output[[trackNames[i]]][,'start']<-as.numeric(output[[trackNames[i]]][,'start'])
		output[[trackNames[i]]][,'end']<-as.numeric(output[[trackNames[i]]][,'end'])
	}
	return(output)
}


degap<-function(seqs,extraChars=''){
	return(gsub(sprintf('[.*%s-]+',extraChars),'',seqs,perl=TRUE))
}



#blat: dataframe from readBlat
#ambigousThreshold: Throw out reads with no match better than ambigousThreshold
#matchThreshold: Throw reads with more than ambigousNumThreshold matches above matchThreshold
#ambigousNumThreshold: Throw reads with more than ambigousNumThreshold matches above matchThreshold
#debug: Display debug messages
trimBlat<-function(blat,ambigousThreshold,matchThreshold,ambigousNumThreshold=1,debug=FALSE){
		goodHits<-tapply(blat$score,blat$qName,function(x,y)sum(x>=max(x)*y),ambigousThreshold)
		selector<-goodHits>ambigousNumThreshold
		message(sum(selector),' reads have at least ',ambigousNumThreshold,' matches more than max(match)*',ambigousThreshold,'. Discarding.')
		selector<-!blat$qName %in% names(goodHits)[selector]
		if(debug)print(t(t(tapply(blat[!selector,'qName'],blat[!selector,'file'],function(x)length(unique(x))))))
		blat<-blat[selector,]
		goodHits<-tapply(blat$score,blat$qName,function(x)max(x))
		goodHits<-data.frame('qName'=names(goodHits),'maximumScore'=goodHits,stringsAsFactors=FALSE)
		numCheck<-nrow(blat)
		blat<-merge(blat,goodHits,all.x=TRUE)
		if(nrow(blat)!=numCheck|any(is.na(blat$maximumScore)))stop(simpleError('Problem finding max score'))
		numCheck<-length(unique(blat$qName))
		blat<-blat[blat$score==blat$maximumScore,]
		if(ambigousNumThreshold==1&nrow(blat)!=numCheck)stop(simpleError('Problem selecting max score read'))
		if(ambigousNumThreshold!=1)message('Returning multiple matches')
		selector<-apply(blat[,c('match','qSize')],1,function(x,y)x['match']<(x['qSize'])*y,matchThreshold)
		message(sum(selector),' reads have a match less than qSize*',matchThreshold,'. Discarding.')
		if(debug){print(t(t(tapply(blat[selector,'qName'],blat[selector,'file'],function(x)length(unique(x))))));browser()}
		blat<-blat[!selector,]
		return(blat)
}

#blat: dataframe from readBlat
#ambigousThreshold: Throw out reads with no match better than ambigousThreshold
#matchThreshold: Throw reads with more than ambigousNumThreshold matches above matchThreshold
#ambigousNumThreshold: Throw reads with more than ambigousNumThreshold matches above matchThreshold
#debug: Display debug messages
#returnSelector: if TRUE return logical selection matrix
trimBlat2<-function(blat,ambigousThreshold,matchThreshold,ambigousNumThreshold=1,debug=FALSE,returnSelector=FALSE){
		selectTable<-data.frame('BADAmbig'=rep(NA,nrow(blat)),'BADNotMax'=rep(NA,nrow(blat)),'BADMatch'=rep(NA,nrow(blat)))
		goodHits<-tapply(blat$score,blat$qName,function(x,y)sum(x>=max(x)*y),ambigousThreshold)
		selector<-goodHits>ambigousNumThreshold
		message(sum(selector),' reads have at least ',ambigousNumThreshold,' matches more than max(match)*',ambigousThreshold,'. Discarding.')
		selector<-!blat$qName %in% names(goodHits)[selector]
		if(debug)print(t(t(tapply(blat[!selector,'qName'],blat[!selector,'file'],function(x)length(unique(x))))))
		blat<-blat[selector,]
		selectTable$BADAmbig<-!selector

		goodHits<-tapply(blat$score,blat$qName,function(x)max(x))
		blat$maximumScore<-goodHits[blat$qName]
		numCheck<-length(unique(blat$qName))
		selector<-blat$score==blat$maximumScore
		blat<-blat[selector,]
		selectTable[!selectTable$BADAmbig,'BADNotMax']<-!selector
		if(ambigousNumThreshold==1&nrow(blat)!=numCheck)stop(simpleError('Problem selecting max score read'))
		if(ambigousNumThreshold!=1)message('Returning multiple matches')
		selector<-blat[,'match']<blat[,'qSize']*matchThreshold
		message(sum(selector),' reads have a match less than qSize*',matchThreshold,'. Discarding.')
		if(debug){print(t(t(tapply(blat[selector,'qName'],blat[selector,'file'],function(x)length(unique(x))))));browser()}
		blat<-blat[!selector,]
		selectTable[!selectTable$BADAmbig,'BADMatch'][!selectTable$BADNotMax[!selectTable$BADAmbig]]<-selector

		if(returnSelector) return(selectTable)
		else return(blat)
}



#function to generate a smoothing window
#sigma: distribution parameter
#numPoints: number of points for window
#author:Kyle Bittinger
#usage:window <- gaussWindow(0.5, 50);smoothData <- convolve(data, window, type="open")
gaussWindow <- function(sigma, numPoints) {
	# formula from http://en.wikipedia.org/wiki/Window_function
	N <- numPoints - 1
	xs <- 0:N

	# modified form of gaussian equation
	numer <- xs - N/2
	denom <- sigma * N/2
	ys <- exp(- 0.5 * ((numer / denom) ^ 2))

	return(ys)
}

#function to smooth data
#data:vector of data to smooth
#window: smoothiing window froom gaussWindow
#author:Kyle Bittinger
#usage:smoothData <- smooth(data, window)
smooth <- function(data, window,truncateWindow=TRUE) {
	smoothData <- convolve(data, window, type="open")

	# smoothData contains length(window)-1 extra points
	N <- length(window)
	trimSize <- N %/% 2 # integer division
	extraTrim <- N %% 2 # modulus

	# trim extra point from beginning of data if necessary
	firstIndex <- trimSize + extraTrim
	lastIndex  <- length(smoothData) - trimSize
	numObs<-length(data)

	windowSums<-rep(sum(window),numObs)
	if(truncateWindow){
		windowSums[1:trimSize]<-windowSums[1:trimSize]-rev(cumsum(window[1:trimSize]))
		reverseIndices<-numObs-(0:(trimSize-2+extraTrim))
		windowSums[reverseIndices]<-windowSums[reverseIndices]-rev(cumsum(rev(window)[1:(trimSize-1+extraTrim)]))
	}
	return(smoothData[firstIndex:lastIndex]/windowSums)
}











aminoAcids<-data.frame('codon'=c('UUU','UUC','UCU','UCC','UAU','UAC','UGU','UGC','UUA','UCA','UAA','UGA','UUG','UCG','UAG','UGG','CUU','CUC','CCU','CCC','CAU','CAC','CGU','CGC','CUA','CUG','CCA','CCG','CAA','CAG','CGA','CGG','AUU','AUC','ACU','ACC','AAU','AAC','AGU','AGC','AUA','ACA','AAA','AGA','AUG','ACG','AAG','AGG','GUU','GUC','GCU','GCC','GAU','GAC','GGU','GGC','GUA','GUG','GCA','GCG','GAA','GAG','GGA','GGG'),'abbr'=c('Phe','Phe','Ser','Ser','Tyr','Tyr','Cys','Cys','Leu','Ser','Ochre','Opal','Leu','Ser','Amber','Trp','Leu','Leu','Pro','Pro','His','His','Arg','Arg','Leu','Leu','Pro','Pro','Gln','Gln','Arg','Arg','Ile','Ile','Thr','Thr','Asn','Asn','Ser','Ser','Ile','Thr','Lys','Arg','Met','Thr','Lys','Arg','Val','Val','Ala','Ala','Asp','Asp','Gly','Gly','Val','Val','Ala','Ala','Glu','Glu','Gly','Gly'),'code'=c('F','F','S','S','Y','Y','C','C','L','S','X','X','L','S','X','W','L','L','P','P','H','H','R','R','L','L','P','P','Q','Q','R','R','I','I','T','T','N','N','S','S','I','T','K','R','M','T','K','R','V','V','A','A','D','D','G','G','V','V','A','A','E','E','G','G'),'name'=c('Phenylalanine','Phenylalanine','Serine','Serine','Tyrosine','Tyrosine','Cysteine','Cysteine','Leucine','Serine','Stop','Stop','Leucine','Serine','Stop','Tryptophan','Leucine','Leucine','Proline','Proline','Histidine','Histidine','Arginine','Arginine','Leucine','Leucine','Proline','Proline','Glutamine','Glutamine','Arginine','Arginine','Isoleucine','Isoleucine','Threonine','Threonine','Asparagine','Asparagine','Serine','Serine','Isoleucine','Threonine','Lysine','Arginine','Methionine','Threonine','Lysine','Arginine','Valine','Valine','Alanine','Alanine','Aspartic acid','Aspartic acid','Glycine','Glycine','Valine','Valine','Alanine','Alanine','Glutamic acid','Glutamic acid','Glycine','Glycine'),row.names=c('UUU','UUC','UCU','UCC','UAU','UAC','UGU','UGC','UUA','UCA','UAA','UGA','UUG','UCG','UAG','UGG','CUU','CUC','CCU','CCC','CAU','CAC','CGU','CGC','CUA','CUG','CCA','CCG','CAA','CAG','CGA','CGG','AUU','AUC','ACU','ACC','AAU','AAC','AGU','AGC','AUA','ACA','AAA','AGA','AUG','ACG','AAG','AGG','GUU','GUC','GCU','GCC','GAU','GAC','GGU','GGC','GUA','GUG','GCA','GCG','GAA','GAG','GGA','GGG'),stringsAsFactors=FALSE)
