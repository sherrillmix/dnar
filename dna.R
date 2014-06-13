nimblegenPrimers<-c('CTCGAGAATTCTGGATCCTC','GAGGATCCAGAATTCTCGAGTT')
primers454<-c("GCCTCCCTCGCGCCATCAG","GCCTTGCCAGCCCGCTCAG")
primerTitanium<-c('CCATCTCATCCCTGCGTGTCTCCGACTCAG','CCTATCCCCTGTGTGCCTTGGCAGTCTCAG')
#convenience function to resize R console window
adjustWindow<-function()options(width=as.integer(Sys.getenv('COLUMNS')))
#convenience function to list objects by size
object.sizes<-function(env=.GlobalEnv)sort(sapply(ls(env=env),function(x)object.size(get(x))),decreasing=TRUE)

ambigousBaseCodes<-c(
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

#get the conservative edge of the odds ratio from fisher's exact test
#x: a 2x2 numeric matrix to feed to fisher.test
conservativeOddsRatio<-function(x){
	conf.int<-fisher.test(x)$conf.int
	out<-ifelse(all(conf.int>1),conf.int[1],ifelse(all(conf.int<1),conf.int[2],1))
	return(out)
}

#convenience function for lagging with NA
lagNA<-function(x,lag=1,fill=NA){
	out<-x	
	if(lag>0)out<-c(out[-(1:lag)],rep(fill,lag))
	if(lag<0)out<-c(rep(fill,abs(lag)),out[-(length(x)-(1:lag)+1)])
	return(out)
}


#convenience function for stop(simpleError())
stopError<-function(...){
	stop(simpleError(paste(...,sep='')))
}

#calculating moving average/max/min etc
#vec: to average over
#statFunc: function to use
#spacer: how many +- to average
movingStat<-function(vec,statFunc=max,spacer=2){
	n<-length(vec)
	sapply(1:n,function(x)statFunc(vec[max(1,x-spacer):min(n,x+spacer)]))
}

#convenience function to check all args are same length
#...: arguments to check the lengths of 
allSameLength<-function(...){
	args<-list(...)
	ns<-sapply(args,length)
	return(all(ns==ns[1]))
}

#convenience function to check all args don't have NAs
#...: arguments to check the lengths of 
anyNa<-function(...){
	args<-list(...)
	hasNa<-sapply(args,function(x)any(is.na(x)))
	return(any(hasNa))
}

#cache an operation to disk or load if present
#cacheFile: file location to save data to 
#operation: a function taking ... arguments and returning object to be stored
#...: arguments for operation function
#OVERWRITE: If FALSE throw an error if hash of ... changes from cached values. If TRUE redo operation and overwrite cache without asking
#VOCAL: If TRUE report on status of caching
#EXCLUDE: Vector of names of arguments to exclude from md5 digest comparison (for very large arguments)
cacheOperation<-function(cacheFile,operation,...,OVERWRITE=FALSE,VOCAL=TRUE,EXCLUDE=NULL){
	#avoid evaluation EXCLUDEd args until necessary since they're probably big
	unevalArgs<-match.call(expand.dots=FALSE)$'...'
	varSelector<-if(is.null(names(unevalArgs)))rep(TRUE,length(unevalArgs)) else !names(unevalArgs) %in% EXCLUDE
	allArgs<-lapply(unevalArgs[varSelector],eval)
	if(!require(digest))stop(simpleError('caching requires digest'))
	#try to prevent function scope from changing things
	md5<-digest(lapply(allArgs,function(x)if(is.function(x))deparse(x)else x))
	if(file.exists(cacheFile)){
		tmp<-new.env()
		load(cacheFile,env=tmp)
		if(with(tmp,md5)==md5){
			if(VOCAL)message('Cache ',cacheFile,' does exist. Loading data')
			return(with(tmp,out))
		}else{
			if(!OVERWRITE)stop(simpleError(sprintf('Input does not match cache for %s. Please delete cache file',cacheFile)))
			if(VOCAL)message('Cache hash does not match args. Rerunning operation')
		}
	}

	if(VOCAL)message('Cache ',cacheFile,' does not exist. Running operation')
	#make sure we have all the args if we excluded some
	if(!is.null(EXCLUDE))allArgs<-list(...)
	out<-do.call(operation,allArgs)
	save(md5,out,file=cacheFile)
	return(out)
}


#string: string to be hashed
#hashSize: size of hashed strings
#everyBase: generate a hash every X strings
#start: output labels start from e.g. label starting from 1000
hashString<-function(string,hashSize,everyBase=1,start=1){
	cuts<-seq(1,nchar(string)-hashSize+1,everyBase)
	hashes<-substring(string,cuts,cuts+hashSize-1)
	return(data.frame('forw'=hashes,'revComp'=revComp(hashes),'start'=cuts,'end'=cuts+hashSize-1,stringsAsFactors=FALSE))
}

#string: string to have case toggled
toggleCase<-function(string){
	chartr(paste(letters,LETTERS,collapse='',sep=''),paste(LETTERS,letters,collapse='',sep=''),string)
}

#pattern: look for pattern in string
#strings: look for pattern in strings
highlightString<-function(pattern,strings){
	locs<-gregexpr(pattern,strings)
	nLetter<-nchar(pattern)
	strings<-mapply(function(x,y){y<-y[y!=-1];for(i in y)substring(x,i,i+nLetter-1)<-toggleCase(substring(x,i,i+nLetter-1));return(x)},strings,locs,USE.NAMES=FALSE)
	return(strings)
}

#convenience function for binding a bunch of sequences together
#...: various sequences to split into a matrix
#fill: fill sequence to pad ends of sequences
seqSplit<-function(...,fill=NULL){
	seqs<-c(...)
	seqN<-nchar(seqs)
	maxN<-max(seqN)
	if(is.null(fill)&&any(seqN!=maxN))stop(simpleError('All sequences not same length'))
	else seqs<-paste(seqs,sapply(maxN-seqN,function(x)paste(rep(fill,x),collapse='')),sep='')
	return(do.call(rbind,strsplit(seqs,'')))
}


#convenience function for selecting elements from a matrix
indexMatrix<-function(x,y,mat,returnIndex=FALSE){
	mat<-as.matrix(mat)
	if(!is.integer(x)){tmp<-1:nrow(mat);names(tmp)<-rownames(mat);x<-tmp[x]}
	if(!is.integer(y)){tmp<-1:ncol(mat);names(tmp)<-colnames(mat);y<-tmp[y]}
	if(length(x)!=length(y)|max(x)>nrow(mat)|max(y)>ncol(mat))stop(simpleError("Dimensions don't match up"))
	index<-(y-1)*nrow(mat)+x
	if(returnIndex)return(index)
	else return(mat[index])
}

#convenience function for picking first most abundant of a set of values
#values: a vector of values
#returns: first most abundant value as a string
mostAbundant<-function(values){
	tmp<-table(values)
	return(names(tmp)[tmp==max(tmp)][1])
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
shannon<-function(x,base=exp(1),standardize=TRUE){
   if(any(x<0))stop(simpleError('Please give positive values of x'))
   x<-x[x>0]
   if(standardize)prop<-x/sum(x)
   else prop<-x
   return(-sum(prop*log(prop,base)))
}


#calculate rao diversity of a vector
rao<-function(x,dist){
	zeros<-x==0
	x<-x[!zeros]
	dist<-dist[!zeros,!zeros]
	props<-x/sum(x)
	propTable<-outer(props,props)
	if(any(dim(propTable)!=dim(dist)))stop(simpleError('Distance and proportion vectors do not match'))
	return(sum(propTable*dist))
}

#jensen shannon divergence between two probability distributions
jensenShannon<-function(x,y,base=2){
	propX<-x/sum(x)
	propY<-y/sum(y)
	shannonCalc<-shannon((propX+propY)/2,base,FALSE)-shannon(propX,base,FALSE)/2-shannon(propY,base,FALSE)/2
	kullbackCalc<-.5*kullback(propX,(propY+propX)/2,base,FALSE)+.5*kullback(propY,(propY+propX)/2,base,FALSE)
	if(round(shannonCalc,5)!=round(kullbackCalc,5))stop(simpleError('Problem calculating shannon jensen'))
	return(kullbackCalc)
}

#Kullback–Leibler divergence between two probability distributions
kullback<-function(x,y,base=2,standardize=TRUE){
	selector<-y>0&x>0
	x<-x[selector];y<-y[selector]
	if(standardize){
		propX<-x/sum(x)
		propY<-y/sum(y)
	}else{
		propX<-x
		propY<-y
	}
	sum(propX*log(propX/propY,base))
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
	counts<-counts[counts>0]
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
#maxN: maximum sample size tested
quickRare<-function(sample,step=10,maxN=sum(sample)){
	sampleSize<-20;
	steps<-unique(c(seq(step,maxN,step),maxN))
	output<-sapply(steps,function(x)rareEquation(sample,x))
	return(data.frame('rare'=output,'sampleN'=steps))
}

#speciesCounts: vector of counts for each "species" e.g. c(10,100,5)
#sampleSize: single value of number of draws from sample 
rareEquation<-function(speciesCounts,sampleSize){
	#numbers too big
	#output2<-length(sample)-choose(sum(sample),sampleSize)^-1*sum(choose(sum(sample)-sample,sampleSize))
	#message(output2)
	#no way to log sum 
	#logSum<-log(sum(choose(sum(sample)-sample,sampleSize)))
	#output<-length(sample) - exp(- lchoose(sum(sample),sampleSize) + logSum)
	#zeros can take computational time
	speciesCounts<-speciesCounts[speciesCounts>0]
	output<-sum(1-exp(lchoose(sum(speciesCounts)-speciesCounts,sampleSize)-lchoose(sum(speciesCounts),sampleSize)))
	if(is.na(output)||is.infinite(output))browser()
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
#naReplace: replace unknown codons with
codon2aa<-function(codons,type='code',naReplace='z',warn=TRUE){
	if(!type %in% c('code','name','abbr'))stop(simpleError('Invalid amino acid type'))
	codons<-gsub('T','U',toupper(codons))
	aas<-aminoAcids[,type,drop=FALSE][codons,1]
	if(warn&any(is.na(aas)))warning('Unknown codons')
	aas[is.na(aas)]<-naReplace
	return(aas)
}

#convert dna/rna to amino acids
#dna: a string of DNA/RNA
#frame: starting frame (0=start on first base, 1=on second, 2=on third)
#debug: print debug info?
dna2aa<-function(dna,frame=0,debug=FALSE,...){
	codons<-dna2codons(dna,frame)	
	if(debug)print(codons)
	output<-paste(codon2aa(codons,...),collapse='')
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

#find closest region for query
#qStart: query start coordinate
#qEnd: query end coordinate
#tStarts: query start coordinates
#tEnds: query end coordinates
closestRegion<-function(qStart,qEnd,tStarts,tEnds){
	distRight<-tStarts-qEnd
	distLeft<-qStart-tEnds
	isRight<-distRight>0
	isLeft<-distLeft>0
	!isRight&!isLeft
	dists<-ifelse(!isRight&!isLeft,0,ifelse(isLeft,distLeft,distRight)) #0 if overlap
	#should probably deal with best overlap better
	return(which.min(dists))
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

checkOverlapIRange<-function(starts,ends,names,tStarts,tEnds,tNames,qChrom,tChrom,vocal=FALSE,...){
	library(IRanges)
	queries<-RangedData(IRanges(starts,end=ends),space=qChrom,names=names)
	targets<-RangedData(IRanges(tStarts,end=tEnds),space=tChrom,names=tNames)
	out<-rdapply(RDApplyParams(queries,function(x){
		thisChr<-unique(space(x))
		message('Working on ',thisChr)
		if(length(thisChr)>1)stop(simpleError('Problem with rdapply'))
		tTree<-IntervalTree(ranges(targets)[[thisChr]])
		tData<-targets[thisChr]
		overlaps<-as.data.frame(slot(findOverlaps(ranges(x)[[1]],tTree),'matchMatrix'))
		xRanges<-ranges(x)[[1]]
		tRanges<-ranges(tData)[[1]]
#		browser()
		#zz<-apply(overlaps,1,function(y){cat('.');matchStats(xRanges[y[1],],tRanges[y[2],])})
		out<-by(overlaps,overlaps$query,function(thisQ){
			query<-xRanges[thisQ$query[1],]#ranges(x[thisQ$query[1],])[[1]]
			matches<-tData[thisQ$subject,]
			mRanges<-tRanges[thisQ$subject]
			if(nrow(matches)==0)return(NULL)
			targetNames<-unique(matches$names)
			stats<-do.call(rbind,lapply(targetNames,function(tName){
				target<-mRanges[matches$names==tName,]
				matchStats(query,target)
			}))
			stats$name<-targetNames
			scores<-stats$qUncover+stats$tUncover
			stats<-stats[scores==min(scores),]
			return(stats)
		})
		return(out)
	}))
#		out<-lapply(unique(x$names),function(thisQ){
#	cat('.')
#			query<-ranges(x[x$names==thisQ,])[[1]]
#			overlap<-slot(findOverlaps(query,tTree),'matchMatrix')
#			matches<-tData[overlap[,'subject'],]
#			if(nrow(matches)==0)return(NULL)
#			targetNames<-unique(matches$names)
#			stats<-do.call(rbind,lapply(targetNames,function(tName){
#				target<-ranges(matches[matches$names==tName,])[[1]]
#				print(system.time(out<-matchStats(query,target)))
#				return(out)
#			}))
#			stats$name<-targetNames
#			scores<-stats$qUncover+stats$tUncover
#			stats<-stats[scores==min(scores),]
#			return(stats)
#		})
#		browser()
#		return(out)
#	}))
	return(out)
}

matchStats<-function(qRanges,tRanges){
	qStart<-min(start(qRanges))
	qEnd<-max(end(qRanges))
	tRanges<-restrict(tRanges,qStart,qEnd)
	tGaps<-gaps(tRanges,start=qStart,end=qEnd)
	qGaps<-gaps(qRanges,start=qStart,end=qEnd)
	uncoveredQuery<-sum(width(safeIntersect(tGaps,qRanges)))
	uncoveredTarget<-sum(width(safeIntersect(qGaps,tRanges)))
	coveredBases<-sum(width(safeIntersect(qRanges,tRanges)))
	coveredGaps<-sum(width(safeIntersect(qGaps,tGaps)))
	output<-data.frame('qUncover'=uncoveredQuery,'tUncover'=uncoveredTarget,'coverMatch'=coveredBases,'gapMatch'=coveredGaps)
	return(output)
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

safeIntersect<-function(x,y){
	if(any(width(x)>0)&&any(width(y)>0)) return(intersect(x,y))
	else return(IRanges(0,width=0))
}


findBestBlockMatch<-function(starts,ends,tStarts,tEnds,tNames,vocal=FALSE,returnMismatch=FALSE,sep='|'){
	if(length(starts)!=length(ends))stop(simpleError('starts and ends not same length'))
	if(length(tStarts)!=length(tEnds))stop(simpleError('tStarts and tEnds not same length'))
	if(length(tStarts)!=length(tNames))stop(simpleError('tStarts and tNames not same length'))
	targets<-RangedData(ranges = IRanges(tStarts,end=tEnds), space=tNames)
	query<-IRanges(starts,end=ends)
	matchData<-do.call(rbind,lapply(ranges(targets),function(x,y)matchStats(y,x),query))
	matchData$score<-matchData$coverMatch+matchData$gapMatch
	matchData<-matchData[matchData$score==max(matchData$score),]
	if(nrow(matchData)!=1&&sum(duplicated(matchData))!=nrow(matchData)-1){message('Multiple equal matches');browser()}
	output<-paste(rownames(matchData),collapse=sep)
	if(returnMismatch)output<-list(output,matchData)
	return(output)
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
	tmp<-apply(cbind(starts,ends),1,function(x,start,end){
		selector<-x[1]:x[2]
		selector<-paste(selector[selector<=end&selector>=start])
		output[selector]<<-output[selector]+1}
	,tStart,tEnd)
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
	checkCoverage <- function(starts, lengths=NULL,ends=starts+lengths-1,outLength=max(ends),trimToStart=FALSE) {
		if(is.null(starts))return(NULL)
		if(length(starts)!=length(ends))stop(simpleError('starts and ends different lengths'))
		if(any(ends>outLength))stop(simpleError('starts+length-1 greater than outLength'))
		if(any(ends<starts))stop(simpleError('ends less than starts'))
		if(trimToStart){
			firstNonZero<-min(starts)
			outLength<-outLength-firstNonZero+1
			starts<-starts-firstNonZero+1
			ends<-ends-firstNonZero+1
		}
		ans<-.C('checkCover',as.integer(rep(0,outLength)),as.integer(starts),as.integer(ends),as.integer(length(starts)))	
		if(trimToStart){
			names(ans[[1]])<-format((1:outLength)+firstNonZero-1) #need format or 7000000 becomes 7e7
		}
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
	data<-data[order(data$start),]
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
reverseString<-function(strings,faster=TRUE){
	if(faster) output<-sapply(strings,function(x)intToUtf8(rev(utf8ToInt(x)))) #http://stackoverflow.com/questions/13612967/how-to-reverse-a-string-in-r
	else output<-sapply(strsplit(strings,''),function(x)paste(rev(x),collapse=''))	
	return(output)
}
#compliment dna 
#dnas: vector of sequences
complimentDna<-function(dnas,brackets=TRUE,ambigs=TRUE){
	finds<-'TGAC'
	replaces<-'ACTG'
	#deal with ambigous
	if(ambigs){
		sortAmbig<-sapply(lapply(strsplit(ambigousBaseCodes,''),sort),paste,collapse='')
		revAmbig<-sapply(strsplit(complimentDna(sortAmbig,ambigs=FALSE),''),function(x)paste(sort(x),collapse=''))
		ambigComp<-names(sortAmbig)[sapply(revAmbig,function(x)which(x==sortAmbig))]
		finds<-sprintf('%s%s',finds,paste(names(sortAmbig),collapse=''))
		replaces<-sprintf('%s%s',replaces,paste(ambigComp,collapse=''))
	}

	if(brackets){finds<-paste(finds,'[]',sep='');replaces<-paste(replaces,'][',sep='')}
	return(chartr(finds,replaces,dnas))
}

#reverse compliment dna
#dnas: vector of sequences
#returns: reverse complimented dna
revComp<-function(dnas){
	return(complimentDna(reverseString(dnas),TRUE))
}

#find a rough maximum independent set
#mat: logical matrix of connection or not connected
#vocal: show progress of matrix and output set?
mis<-function(mat,vocal=FALSE){
	#very basic error checks
	if(is.null(colnames(mat)))colnames(mat)<-1:ncol(mat)
	if(ncol(mat)!=nrow(mat))stop(simpleError('Square matrix only for MIS'))
	if(any(!mat %in% c(0,1,TRUE,FALSE)))stop(simpleError('Only binary/logical matrices for MIS'))

	output<-c()
	while(nrow(mat)>0){
		if(vocal){message('Matrix:');print(mat)}
		connects<-apply(mat,1,sum)
		if(any(connects==1)){
			#anything only connected to itself can go
			remove<-connects==1
			output<-c(output,colnames(mat)[remove])
			mat<-mat[!remove,!remove,drop=FALSE]
		}else{
			#find first node with minimum connections
			remove<-which.min(connects)
			output<-c(output,colnames(mat)[remove])
			mat<-mat[!mat[,remove],!mat[,remove],drop=FALSE]
		}
		if(vocal){message('Current Set:');print(output)}
	}
	return(output)
}

#brute force graph coloring
graphColor<-function(mat,reps=3){
	N<-nrow(mat)
	if(N==1)return(1)
	if(ncol(mat)!=N)stop(simpleError('Square matrix only for coloring'))
	if(all(mat))return(1:N)
	#mask out diagonal
	mat<-mat&!diag(N)
	degree<-apply(mat,1,sum)
	if(all(degree==0))return(rep(1,N))
	possibleOut<-list(rep(NA,reps))
	for(i in 1:reps){
		colors<-lapply(1:N,function(x)1:c(max(degree)+1))
		colorLength<-sapply(colors,length)
		colorSelector<-colorLength>1
		while(any(colorSelector)){
			current<-which(colorSelector)[order(colorLength[colorSelector],-degree[colorSelector],runif(sum(colorSelector)),decreasing=TRUE)][1]
			colors[[current]]<-min(colors[[current]])
			neighbors<-which(mat[current,]&colorSelector)
			sapply(neighbors,function(x,y)colors[[x]]<<-colors[[x]][colors[[x]]!=y],colors[[current]])
			colorLength<-sapply(colors,length)
			colorSelector<-colorLength>1
		}
		possibleOut[[i]]<-unlist(colors)
	}
	best<-which.min(sapply(possibleOut,max))
	return(possibleOut[[best]])
}

#read DIMACS file
read.col<-function(file){
	x<-readLines(file)
	y<-x[grep('^e',x)]
	y<-gsub('^e ','',y)
	z<-do.call(rbind,lapply(strsplit(y,'\\s'),as.numeric))
	info<-strsplit(x[grep('^p',x)],'\\s')[[1]]
	nodes<-as.numeric(info[3])
	edges<-as.numeric(info[4])
	if(any(z>nodes))stop(simpleError('Nodes not numbered correctly'))
	if(nrow(z)!=edges)stop(simpleError('Edges not numbered correctly'))
	output<-matrix(FALSE,ncol=nodes,nrow=nodes)
	apply(z,1,function(x){
		output[x[1],x[2]]<<-TRUE
		output[x[2],x[1]]<<-TRUE
	})
	return(output)
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
	#make sure matching + and @
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
#limit: Only read in limit lines
#vocal: Echo progress messages
#returns: dataframe with solexa info seq and quals
read.solexa<-function(fileName,convert=TRUE,limit=-1,vocal=FALSE){
	x<-readLines(fileName,n=limit)
	if(vocal)message('Done reading in')
	output<-do.call(rbind,strsplit(x,':'))
	rm(x)
	if(vocal)message('Do.call done')
	output<-cbind(output,do.call(rbind,strsplit(output[,5],'[#/]')))
	if(vocal)message('2nd do.call done')
	output<-as.data.frame(output[,-5],stringsAsFactors=FALSE)
	colnames(output)<-c('instrument','lane','tile','x','seq','rawQual','y','barcode','pair')
	if(vocal)message('Start convert')
	if(convert)output$qual<-sapply(output$rawQual,function(x)paste(as.integer(charToRaw(x))-64,collapse=' '))
	if(vocal)message('Convert done')
	output$name<-paste(output$instrument,output$lane,output$tile,output$x,output$y,output$barcode,output$pair,sep='_')
	return(output)
}

#converts ambigous dna to an appropriate regular expression
#input sequence with ambigous bases
#returns: regular expression
ambigous2regex<-function(dna){
	dna<-toupper(dna)
	for (i in names(ambigousBaseCodes)){
		dna<-gsub(i,paste('[',ambigousBaseCodes[i],']',sep=''),dna)
	}
	return(dna)
}
#convert set of single bases to ambigous code
#bases: vector of single character base strings
bases2ambigous<-function(bases){
	bases<-sort(unique(bases))
	nBases<-nchar(bases[1])
	if(any(nchar(bases)!=nBases))stop(simpleError('Convert bases to ambigous requires same length sequences'))
	if(length(bases)==1)return(bases)
	if(nBases>1)return(paste(sapply(1:nBases,function(x)bases2ambigous(substring(bases,x,x))),collapse=''))
	else return(names(ambigousBaseCodes)[paste(bases,collapse='')==ambigousBaseCodes])
}



#convert ambigous dna to all possible sequences
#dna: vector dna containing ambigous bases
#unlist: return a unlisted vector instead of a list
#returns: list with each entry containing all combinations for that entry of the vector
expandAmbigous<-function(dna,delist=FALSE){
	dna<-toupper(dna)
	ambigRegex<-sprintf('[%s]',paste(names(ambigousBaseCodes),collapse=''))
	out<-lapply(dna,function(x){
		pos<-regexpr(ambigRegex,x)	
		if(pos!=-1){
			replaces<-strsplit(ambigousBaseCodes[substring(x,pos,pos)],'')[[1]]
			new<-rep(x,length(replaces))
			for(i in 1:length(replaces))substring(new[i],pos,pos)<-replaces[i]
			out<-expandAmbigous(new,TRUE)
		}else{out<-x}
		return(out)
	})
	if(delist)out<-unlist(out)
	return(out)
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
read.fa2<-function(fileName,longNameTrim=TRUE,...){
	x<-readLines(fileName,warn=FALSE,...)
	if(length(x)==0)return(NULL)
	x<-x[!grepl('^[;#]',x,perl=TRUE)&x!='']
	nameLines<-grep('^>',x)
	thisNames<-sub('^>','',x[nameLines])
	if(any(grep(' [^ ]',x[-nameLines][3],perl=TRUE)))hasSpaces<-TRUE
	else hasSpaces<-FALSE
	seqs<-apply(cbind(nameLines+1,c(nameLines[-1]-1,length(x))),1,function(coords){
		if(coords[1]<=coords[2])return(paste(x[coords[1]:coords[2]],collapse=ifelse(hasSpaces,' ','')))
		else return('')
	})
	seqs<-gsub('  +',' ',seqs,perl=TRUE)
	seqs<-sub(' $','',seqs,perl=TRUE)

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

#reg: region in the format "chrX:123545-123324"
#files: bam files to pull the counts from
#bam2depthBinary: bam2depth executable file
pullRegion<-function(reg,files,bam2depthBinary='./bam2depth'){
	samArg<-sprintf('-r %s',reg)
	fileArg<-paste(files,collapse=' ')
	cover<-samView(fileArg,samArgs=samArg,samCommand='',samtoolsBinary=bam2depthBinary,colClasses=c('character',rep('numeric',length(files)+1)))
	if(is.null(cover))cover<-do.call(data.frame,c(list('XXX'),as.list(-(1:(length(files)+1)))))[0,]
	colnames(cover)<-c('chr','pos',sprintf('counts%d',1:(ncol(cover)-2)))
	return(cover)
}

#reg: region in the format "chrX:123545-123324"
parseRegion<-function(reg){
	splits<-strsplit(reg,'[:-]')
	if(any(sapply(splits,length)!=3))stop(simpleError('Region not parsed'))
	out<-data.frame('chr'=sapply(splits,'[[',1),stringsAsFactors=FALSE)
	out[,c('start','end')]<-as.numeric(do.call(rbind,lapply(splits,'[',2:3)))
	return(out)
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

#convert cigar and starts to qStarts, tStarts, blockSizes as in blat
#cigars: vector of SAM cigar strings
#starts: vector of starting positions in target
#seqs: vector of query sequences (only if alignments desired)
#tSeq: single target sequence (only if alignments desired)
#returns: dataframe with qStarts,tStarts,sizes if !startEnds or dataframe with starts, ends and ids if startEnds or dataframe with pairwise query alignments, qAlign, and target alignments, tAlign if provided seqs and tSeq
cigarToBlock<-function(cigars,starts,startEnds=FALSE,seqs=NULL,tSeq=NULL){
	#M=match, I=insertion in query, D=deletion in query, N="intron" deletion in query, S=soft clipping (clip sequence), H=hard clipping (sequence was already clipped)
	nAligns<-length(starts)
	if(nAligns!=length(cigars))stop(simpleError('Cigars and starts not same length'))
	supportedOps<-c('M','I','D','N','S','H')
	if(any(grep(sprintf('[^0-9%s]',paste(supportedOps,collapse='')),cigars)))stop(simpleError(sprintf('Only %s cigar operations supported',paste(supportedOps,collapse=''))))
	tPos<-starts
	qPos<-rep(1,nAligns)
	stillWorking<-rep(TRUE,nAligns)
	qStarts<-tStarts<-rep('',nAligns)
	blockSizes<-tStarts<-rep('',nAligns)
	isAlign<-FALSE
	if(!is.null(seqs)&&!is.null(tSeq)){
		if(length(seqs)!=nAligns)stop(simpleError('Cigars and seqs not same length'))
		tSeq<-tSeq[1]
		isAlign<-TRUE
		aligns<-data.frame('qAlign'=rep('',nAligns),'tAlign'=rep('',nAligns),stringsAsFactors=FALSE)
	}
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
				if(isAlign)aligns[stillWorking,][matches,c('tAlign')]<-sprintf('%s%s',aligns[stillWorking,][matches,c('tAlign')],substring(tSeq,tPos[stillWorking][matches],tPos[stillWorking][matches]+num-1))
				tPos[stillWorking][matches]<-num+tPos[stillWorking][matches]
			}
			if(i %in% c('M','I','S')){
				if(isAlign)aligns[stillWorking,][matches,c('qAlign')]<-sprintf('%s%s',aligns[stillWorking,][matches,c('qAlign')],substring(seqs[stillWorking][matches],qPos[stillWorking][matches],qPos[stillWorking][matches]+num-1))
				qPos[stillWorking][matches]<-num+qPos[stillWorking][matches]
			}
			if(i!='M'&&isAlign){
				if(isAlign)if(any(num>N_DUMMY_GAPS))stop(simpleError('Gap larger than ',N_DUMMY_GAPS))
				if(i %in% c('D','N'))aligns[stillWorking,][matches,c('qAlign')]<-sprintf('%s%s',aligns[stillWorking,][matches,c('qAlign')],substring(dummyGaps,1,num))
				if(i %in% c('I','S'))aligns[stillWorking,][matches,c('tAlign')]<-sprintf('%s%s',aligns[stillWorking,][matches,c('tAlign')],substring(dummyGaps,1,num))
			}
			stillWorking[stillWorking]<-nchar(cigars[stillWorking])>0
		}
	}
	qStarts<-sub('^,','',qStarts)
	tStarts<-sub('^,','',tStarts)
	blockSizes<-sub('^,','',blockSizes)
	if(isAlign)return(aligns)
	if(startEnds)return(startEndsOut)
	else return(data.frame('qStarts'=qStarts,'tStarts'=tStarts,'sizes'=blockSizes,stringsAsFactors=FALSE))
}

#seqs: vector of sequences
#tSeqs: vector of target sequences or a single target
#qStarts: comma separated starts of query matches e.g. from blat or cigarToBlock (1 based)
#tStarts: comma separated starts of target matches e.g. from blat or cigarToBlock (1 based)
#sizes: comma separated lengths of matches e.g. from blat or cigarToBlock
blockToAlign<-function(seqs,tSeqs,qStarts,tStarts,sizes){
	nSeqs<-length(seqs)
	if(length(tSeqs)==1)tSeqs<-rep(tSeqs,length(seqs))
	if(length(tSeqs)!=length(seqs))stop(simpleError('Target seqs not same length as seqs and not a single sequence'))
	if(!all.equal(length(seqs),length(qStarts),length(tStarts),length(sizes)))stop(simpleError('All arguments not same length'))
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
		gap$tSeq<-substring(tSeq,gap$tGapStartAfter+1,gap$tGapStartAfter+gap$tGaps)
		gap$qSeq<-substring(qSeq,gap$qGapStartAfter+1,gap$qGapStartAfter+gap$qGaps)
		dummy<-paste(rep('-',max(nchar(c(gap$tSeq,gap$qSeq)))),collapse='')
		gap[,c('tSeq','qSeq')]<-t(apply(gap[,c('tSeq','qSeq')],1,function(x)paste(x,substring(dummy,1,max(nchar(x))-nchar(x)),sep='')))
		return(gap)
	},gaps,tSeqs,seqs,SIMPLIFY=FALSE)
	qSeqs<-sapply(1:length(seqs),function(ii){
		x<-qPieces[qPieces$name==ii,]
		thisGaps<-gaps[[ii]]
		return(paste(c(x$seq,thisGaps$qSeq)[order(c(x$start,thisGaps$qGapStartAfter+.5))],collapse=''))
	})
	tSeqs<-sapply(1:length(seqs),function(ii){
		x<-tPieces[tPieces$name==ii,]
		thisGaps<-gaps[[ii]]
		return(paste(c(x$seq,thisGaps$tSeq)[order(c(x$start,thisGaps$tGapStartAfter+.5))],collapse=''))
	})
	return(data.frame('qSeq'=qSeqs,'tSeq'=tSeqs,stringsAsFactors=FALSE))
}

#return order of one vector in another
#query: values to be sorted in target order
#target: order for query to be sorted into
#strict: if true error if query not in target, if false append unknown queries to end of target
#...: arguments for order
orderIn<-function(query,target,strict=FALSE,orderFunc=order,...){
	if(any(!query %in% target)){
		if(strict){
			stop(simpleError('Query not in target'))
		}else{
			target<-c(target,unique(query[!query %in% target]))
		}
	}
	newOrder<-orderFunc(sapply(query,function(x)min(which(x==target))),...)
	return(newOrder)
}

#flags: vector of integer flags from sam
#test: either character vector of flag short names or integers
samFlag<-function(flags,test='paired'){
	if(!require(bitops))stop(simpleError('Sam flag requires bitops package'))
	test<-unique(test)
	if(!is.numeric(test)){
		if(!all(test %in% samFlags$short))stop(simpleError(sprintf('Unknown flag please select from %s',paste(samFlags$short,collapse=', '))))
		test<-samFlags[samFlags$short %in% test,'bit']
	}else{
		test<-as.integer(test)
	}
	testInt<-0
	for(i in test)testInt<-bitOr(testInt,i)
	return(bitAnd(flags,testInt)==testInt)
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




#start a gfServer on an open port and return port
#nibDir: directory containing nib files to align against
#options: options to gfServer
#gfServer: path to gfServer program
#startPort: begin looking for open (no other R called gfServer using) ports on startPort and add 1 until finding one
#nibSuffix: nib files end in nibSuffix
#wait: number of 5 second intervals to wait for gfServer to finish starting before giving up
startBlat<-function(nibDir,options='',gfServer='gfServer',bit2=NULL,startPort=37900,nibSuffix='.nib',wait=120){
	port<-startPort
	counter<-1
	while(checkBlat(port)){
		port<-port+1
		if(counter>30)stop(simpleError('All 30 checked ports taken'))
		counter<-counter+1
	}
	if(!is.null(bit2))cmd<-sprintf('%s start localhost %d %s %s',gfServer,port,options,normalizePath(bit2))
	else cmd<-sprintf('%s start localhost %d %s %s/*%s',gfServer,port,options,normalizePath(nibDir),nibSuffix)
	message(cmd)
	system(cmd,wait=FALSE)
	counter<-1
	while(system(sprintf('%s status localhost %d >/dev/null',gfServer,port),ignore.stderr=TRUE)!=0){
		Sys.sleep(5)
		if(counter>wait){
			message("Blat took too long")
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
	
	cmd<-sprintf('%s localhost %d / %s %s %s',gfClient,port,faFile,gfClientOptions,outFile)
	message(cmd)
	errorCode<-system(cmd)
	if(errorCode!=0){
		message('Error running blat on port ',port, '. Trying again')
		errorCode<-system(cmd)
	}

	print(checkBlat(port))
	message('Code: ',errorCode,' on port ',port)
	message('Killing blat server on port ',port)
	killBlat(port)
}


#run blat from blat executable instead of gfserver/clien
#reads: named vector of sequences
#refs: named vector of references
#blat: path to blat program
#blatArgs: string of arguments for blat
#tmpDir: a directory to write tempfiles to
#outfile: file to write blat to 
#readFile: file to use instead of reads (will be deleted)
runBlatNoServer<-function(reads=NULL,refs,blatArgs='',outFile='out.blat',blat='blat',tmpDir=tempdir(),readFile=NULL,deleteFiles=!is.null(reads)){
	if(!file.exists(tmpDir))dir.create(tmpDir)
	if(is.null(reads)&is.null(readFile))stop(simpleError('Please provide reads or readFile'))
	if(is.null(readFile)){
		readFile<-sprintf('%s/read.fa',tmpDir)
		write.fa(names(reads),reads,readFile)
	}
	refFile<-sprintf('%s/refs.fa',tmpDir)
	write.fa(names(refs),refs,refFile)
	cmd<-sprintf('%s %s %s %s %s',blat,refFile,readFile,blatArgs,outFile)
	message(cmd)
	errorCode<-system(cmd)
	if(errorCode!=0)stop(simpleError('Problem running blat'))
	if(deleteFiles)file.remove(refFile,readFile,tmpDir)
	return(TRUE)
}	

#run blat from blat executable instead of gfserver/clien
#reads: named vector of sequences
#refs: named vector of references
#blat: path to blat program
#blatArgs: string of arguments for blat
#tmpDir: a directory to write tempfiles to
#nCore: number of cores to use
#outfile: file to write blat to (if ends in .gz then isGz defaults to true and writes to gzipped file)
multiRunBlatNoServer<-function(reads,refs,outFile,nCore=4,tmpDir=tempdir(),condense=TRUE,isGz=grepl('.gz$',outFile),sleepIncrement=1,...){
	prefix<-paste(sample(c(letters,LETTERS),20,TRUE),collapse='')
	library(parallel)
	if(!file.exists(tmpDir))dir.create(tmpDir)
	message('Preparing files')
	runFiles<-mapply(function(seqs,id){
		thisTmpDir<-sprintf('%s/work__%s__%d',tmpDir,prefix,id)
		dir.create(thisTmpDir)
		tmpFile<-sprintf('%s__%d.blat',sub('\\.blat(\\.gz)?$','',outFile),id)
		faFile<-sprintf('%s/reads.fa',thisTmpDir)
		write.fa(names(seqs),seqs,faFile)
		return(c(thisTmpDir,tmpFile,faFile))
	},split(reads,sort(rep(1:nCore,length.out=length(reads)))),1:nCore,SIMPLIFY=FALSE)
	message('Running blats')
	bigRun<-mclapply(runFiles,function(x){
		message('Starting ',x[2])
		runBlatNoServer(readFile=x[3],refs=refs,outFile=x[2],tmpDir=x[1],deleteFiles=TRUE,...)
		return(x[2])
	},mc.cores=nCore)

	if(any(!sapply(bigRun,file.exists)))stop(simpleError('Blat file missing'))
	if(isGz)outFile<-gzfile(outFile,open='w+')
	else outFile<-file(outFile,open='w+')
	counter<-0
	for(i in 1:length(bigRun)){
		message('Reading blat file ',i)
		blat<-readLines(bigRun[[i]])
		#take off header in later files
		if(i!=1){
			blat<-blat[-1:-5]
			append<-TRUE
		}else{
			append<-FALSE
		}
		writeLines(blat,sep="\n",con=outFile)
		file.remove(bigRun[[i]])
		counter<-counter+length(blat)
	}
	message('Wrote ',counter,' blat lines')
	close(outFile)
	return(outFile)
}

#make a 2bit file from one set of reads and blat another set against it
#reads: vector of query reads with names
#refs: vector of reference reads with names
#faToTwoBit: path to faToTwoBit program from blat
#tmpDir: directory to store work files
#...: additional arguments to run blat
blatReadsVsRefs<-function(reads,refs,outFile,faToTwoBit='faToTwoBit',tmpDir=sprintf('%s/%s',tempdir(),paste(sample(c(letters,LETTERS),20),collapse='')),...){
	dir.create(tmpDir)
	readFile<-sprintf('%s/read.fa',tmpDir)
	write.fa(names(reads),reads,readFile)
	refFile<-sprintf('%s/refs.fa',tmpDir)
	write.fa(names(refs),refs,refFile)
	twobitFile<-sprintf('%s/refs.2bit',tmpDir)
	system(sprintf('%s %s %s',faToTwoBit,refFile,twobitFile))
	runBlat(readFile,outFile=outFile,bit2=twobitFile,...)
	file.remove(twobitFile,refFile,readFile,tmpDir)
}

#run blat parallel (requires parallel package included in R 2.4.1)
#reads: vector of query reads with names
#refs: vector of query refs with names
#tmpDir: directory to store work files
#...:arguments for blatReadsVsRefs
multiBlatReadsVsRefs<-function(reads,refs,outFile,nCore=4,tmpDir=tempdir(),isGz=grepl('.gz$',outFile),condense=TRUE,gzipPath='gzip',...){
	prefix<-paste(sample(c(letters,LETTERS),20,TRUE),collapse='')
	library(parallel)
	if(!file.exists(tmpDir))dir.create(tmpDir)
	bigRun<-mclapply(mapply(function(x,y)list(x,y),split(reads,sort(rep(1:nCore,length.out=length(reads)))),1:nCore,SIMPLIFY=FALSE),function(x){
		tmpFile<-sprintf('%s__%d.blat',sub('\\.blat(\\.gz)?$','',outFile),x[[2]])
		thisTmpDir<-sprintf('%s/work__%s__%d',tmpDir,prefix,x[[2]])
		blatReadsVsRefs(x[[1]],refs,tmpFile,startPort=37900+x[[2]]*10,tmpDir=thisTmpDir,...) #likes to start on same port apparently
		if(isGz){
			system(sprintf('%s %s',gzipPath,tmpFile))
			tmpFile<-sprintf('%s.gz',tmpFile)
		}
		#print((tmpFile))
		return(tmpFile)
	},mc.cores=nCore)
	if(any(!sapply(bigRun,file.exists)))stop(simpleError('Blat file missing'))
	if(condense){
		if(isGz)outFile<-gzfile(outFile,open='w+')
		else outFile<-file(outFile,open='w+')
		counter<-0
		for(i in 1:length(bigRun)){
			blat<-readLines(bigRun[[i]])
			#take off header in later files
			if(i!=1){
				blat<-blat[-1:-5]
				append<-TRUE
			}else{
				append<-FALSE
			}
			writeLines(blat,sep="\n",con=outFile)
			file.remove(bigRun[[i]])
			counter<-counter+length(blat)
		}
		message('Wrote ',counter,' blat lines')
		close(outFile)
	}
	return(bigRun)
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

#seqs:sequences to be trimmed
#start: trim starts?
#end: trim ends?
#nonNStretch: number of nonNs required
#nStretch: delete until no stretch of Ns greater than this
trimNs<-function(seqs,start=TRUE,end=TRUE,nonNStretch=NULL,nStretch=10){
	regex<-sprintf('[ACTG]{%d,}',nonNStretch)
	if(start&!is.null(nonNStretch)){
		starts<-regexpr(sprintf('%s',regex),seqs)
		seqs[starts==-1]<-''
		seqs<-substring(seqs,starts)
	}
	if(end&!is.null(nonNStretch)){
		starts<-sapply(gregexpr(sprintf('%s',regex),seqs),tail,1)
		lengths<-sapply(gregexpr(sprintf('%s',regex),seqs),function(x)tail(attr(x,'match.length'),1))
		seqs[starts==-1]<-''
		seqs<-substring(seqs,1,starts+lengths-1)
	}
	if(!is.null(nStretch)){
		regex<-sprintf('N{%d,}',nStretch)		
		nStretches<-gregexpr(regex,seqs)
		seqs<-mapply(function(reg,seq){
			starts<-unique(c(reg,nchar(seq)+1))
			ends<-c(0,reg+attr(reg,'match.length')-1)
			diffs<-starts-ends-1
			select<-which.max(diffs)
			return(substring(seq,ends[select]+1,starts[select]-1))
		},nStretches,seqs)
	}
	return(seqs)
}

trimEnd<-function(seqs,revCompPrimer,trimmed=rep(FALSE,length(seqs)),minSubstring=8,location='end'){
	if(length(seqs)!=length(trimmed)) stop(simpleError('Already trimmed vector and seqs vector not same length'))
	if(!location %in% c('start','end'))stop(simpleError('Please specify location as start or end'))
	trimSeq<-seqs
	trimmedBak<-trimmed
	nCharPrimer<-nchar(revCompPrimer)
	for(i in nCharPrimer:minSubstring){
			  if(location=='end')thisRegex<-paste(substr(revCompPrimer,1,i),'$',sep='')
			  else thisRegex<-paste('^',substr(revCompPrimer,nCharPrimer-i+1,nCharPrimer),sep='')
			  selector<-!trimmed & grepl(thisRegex,seqs)
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
	if(adjustStart!=0)startsList<-lapply(startsList,function(x)as.numeric(x)+adjustStart)
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

#seq: sequence
#gaps: gap characters to remove
#return: string with no gap characters
removeGaps<-function(seq,gaps=c('*','-','.')){
	gsub(sprintf('[%s]+',paste(gaps,collapse='')),'',seq)
}

#seq: DNA sequence
#flowOrder: order of nucleotide flows
#outputLength: minimum output length
#return: vector of flows
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
	if(is.character(flow)&length(flow)==1)flow<-strsplit(flow,'\t')[[1]]
	flow<-as.numeric(flow)
	flow[flow<0]<-0
	chars<-rep(flowOrder,length.out=length(flow))
	output<-paste(rep(chars,round(as.numeric(flow))),collapse='')
	return(output)
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


#find and remove columns that are entirely gaps
removeGapColumns<-function(seqs,gapChars=c('-','.')){
	if(any(nchar(seqs)!=nchar(seqs[1])))stop(simpleError('All seqs not same length'))
	regex<-sprintf('[%s]',paste(gapChars,collapse=''))
	invRegex<-sprintf('[^%s]',paste(gapChars,collapse=''))
	#nonGap<-#unique(unlist(gregexpr(invRegex,seqs[1:min(1000,length(seqs))])))
	nonGap<-c()
	nSeq<-length(seqs)
	for(i in seq(1,length(seqs),500)){
		#if(i %% 1000==1)message(i)
		nonGap<-unique(c(nonGap,unlist(gregexpr(invRegex,seqs[i:min(i+499,nSeq)]))))
	}
	baseRanges<-index2range(nonGap)
	out<-apply(apply(baseRanges,1,function(x)substring(seqs,x[1],x[2])),1,paste,collapse='')
	if(any(nchar(out)!=nchar(out[1])))stop(simpleError('Output seqs not same length'))
	return(out)
}

#isGap<-apply(do.call(rbind,strsplit(scer$alignCut,'')),2,function(x)mean(x=='-'))>.5

#noGapSeq: sequence to expand
#seqLength: output sequence length
#isNotGap: indexes of non gap positions (should be >= nchar(noGapSeq))
addGaps<-function(noGapSeq,isNotGap,seqLength=max(isNotGap)){
	if(is.logical(isNotGap))isNotGap<-which(isNotGap)
	isNotGap<-isNotGap[isNotGap<=seqLength]
	if(length(isNotGap)<nchar(noGapSeq))noGapSeq<-substring(noGapSeq,1,length(isNotGap))
	outSeq<-paste(rep('-',seqLength),collapse='')
	for(i in 1:length(isNotGap))substring(outSeq,isNotGap[i],isNotGap[i])<-substring(noGapSeq,i,i)
	return(outSeq)
}

#forward:sequence for forward primer (ambigous bases ok)
#reverse:sequence for reverse primer (ambigous bases ok)
#seqs: gapped sequences to look for primer in
#refSeq: reference sequence for base counting
#groups: groupings for display next to sequence plots
#fName: forward primer name
#rName: reverse primer name
#outDir: directory to put output files in
#baseName: string to prepend on file names
plotPrimers<-function(forward,reverse,seqs,refSeq,groups,fName='Forward',rName='Reverse',outDir='.',baseName='',...){
	outDir<-sprintf('%s/%s',outDir,baseName)
	revPos<-findPrimer(seqs,reverse=reverse)
	forPos<-findPrimer(seqs,forward)
	cuts<-substring(seqs,forPos[1],revPos[2])
	primerCut1<-substring(seqs,forPos[1],forPos[2])
	primerCut2<-substring(seqs,revPos[1],revPos[2])
	notGap<-which(strsplit(refSeq,'')[[1]]!='-')

	gapFor<-addGaps(forward,which(apply(do.call(rbind,strsplit(primerCut1,'')),2,function(x)mean(x=='-'))<.5))
	gapRev<-addGaps(revComp(reverse),which(apply(do.call(rbind,strsplit(primerCut2,'')),2,function(x)mean(x=='-'))<.5))

	repNum<-ceiling(length(seqs)/20)
	gapRevRep<-expandAmbigous(gapRev)[[1]]
	gapRevRep<-rep(gapRevRep,length.out=max(repNum,length(gapRevRep)))
	gapForRep<-expandAmbigous(gapFor)[[1]]
	gapForRep<-rep(gapForRep,length.out=max(repNum,length(gapForRep)))

	gapCombo<-addGaps(sprintf('%s%s',forward,revComp(reverse)),notGap[(notGap>=forPos[1]&notGap<=forPos[2])|(notGap>=revPos[1]&notGap<=revPos[2])]-forPos[1]+1)
	gapComboRep<-gsub('-','.',expandAmbigous(gapCombo)[[1]])
	gapComboRep<-rep(gapComboRep,length.out=max(repNum,length(gapComboRep)))


	plotSeq(c(primerCut2,gapRevRep),sprintf('%s%s.png',outDir,rName),groups=c(groups,rep(sprintf('$%s',rName),length(gapRevRep))),refSeq=refSeq,xstart=gap2NoGap(refSeq,revPos[1]),...)
	plotSeq(c(primerCut1,gapForRep),sprintf('%s%s.png',outDir,fName),groups=c(groups,rep(sprintf('$%s',fName),length(gapForRep))),refSeq=refSeq,xstart=gap2NoGap(refSeq,forPos[1]),...)
	plotSeq(c(cuts,gapComboRep),sprintf('%s%s--%s.png',outDir,fName,rName),groups=c(groups,rep(sprintf('$%s--%s',fName,rName),length(gapComboRep))),refSeq=refSeq,xstart=gap2NoGap(refSeq,forPos[1]),...)
}

#forward:sequence for forward primer (ambigous bases ok)
#reverse:sequence for reverse primer (ambigous bases ok)
#seqs: gapped sequences to look for primer in
#padding: extra bases to include to each side
findPrimer<-function(seqs,forward=NULL,reverse=NULL,padding=0){
	if(is.null(forward)&&is.null(reverse))stop(simpleError('Please provide forward or reverse primer'))
	trims<-degap(seqs)

	if(!is.null(reverse)){
		revExpand<-expandAmbigous(revComp(reverse))[[1]]
		reverseStart<-sapply(trims,function(x)multiMismatch(revExpand,x)[2],USE.NAMES=FALSE)
		revEnd<-mostAbundant(mapply(noGap2Gap,seqs,reverseStart+nchar(reverse)-1+padding,USE.NAMES=FALSE))
		revStart<-mostAbundant(mapply(noGap2Gap,seqs,reverseStart-padding,USE.NAMES=FALSE))
		if(is.null(forward))return(as.numeric(c(revStart,revEnd)))
	}
	if(!is.null(forward)){
		forExpand<-expandAmbigous(forward)[[1]]
		forwardStart<-sapply(trims,function(x)multiMismatch(forExpand,x)[2],USE.NAMES=FALSE)
		forStart<-mostAbundant(mapply(noGap2Gap,seqs,forwardStart-padding,USE.NAMES=FALSE))
		forEnd<-mostAbundant(mapply(noGap2Gap,seqs,forwardStart+nchar(forward)-1+padding,USE.NAMES=FALSE))
		if(is.null(reverse))return(as.numeric(c(forStart,forEnd)))
	}

	return(as.numeric(c(forStart,revEnd)))
}

#condense any blocks of matches that aren't separated by condenseLimit bases
#start: vector (or string) of start locations
#block: vector (or string)  of match lengths
#condenseLimit: condense two neighboring matches together if separated by <=condenseLimit
#start2: a vector (or string) of start locations. do not condense unless both starts are <=condenseLimit
#synchronizeCondenseLimit: condense two neighboring matches together if both starts and starts2 have a gap here <=synchronizeCondenseLimit and gap2-gap1-1<condenseLimit
#invertCoords: Change from internal negative indexing to positive (inverting starts and ends)
condenseCoords<-function(start,block,condenseLimit=5,start2=NULL,synchronizeCondenseLimit=-Inf,invertCoords=TRUE){
	if(is.character(start))start<-as.numeric(strsplit(start,',')[[1]])
	if(is.character(block))block<-as.numeric(strsplit(block,',')[[1]])
	if(!is.null(start2)&&is.character(start2))start2<-as.numeric(strsplit(start2,',')[[1]])
	if(start[1]>start[length(start)])start<-0-start-block+1
  	ends<-start+block-1

	isDual<-!is.null(start2)
	if(isDual){
		if(start2[1]>start2[length(start2)])start2<-0-start2-block+1
		ends2<-start2+block-1
	}

	if(length(start)>1){
		gaps<-c(start[-1]-ends[-length(ends)]-1,Inf)
		if(isDual){
			gaps2<-c(start2[-1]-ends2[-length(ends2)]-1,Inf)
			condenseSelect<-which(((gaps<=condenseLimit&gaps2<=condenseLimit)|(gaps<=synchronizeCondenseLimit&gaps2<=synchronizeCondenseLimit&abs(gaps2-gaps)-1<condenseLimit)))
		}else{
			condenseSelect<-which((gaps<=condenseLimit))
		}
		if(any(condenseSelect)){
			start<-start[-(condenseSelect+1)]
			ends<-ends[-condenseSelect]
			if(isDual){
				start2<-start2[-(condenseSelect+1)]
				ends2<-ends2[-condenseSelect]
			}
		}
	}
	if(invertCoords){
		if(all(start<0)&&all(ends<0)){
			tmp<--start
			start<--ends
			ends<-tmp
		}
		if(isDual&&all(start2<0)&&all(ends2<0)){
			tmp<--start2
			start2<--ends2
			ends2<-tmp
		}
	}
	if(isDual)out<-data.frame('start'=start,'end'=ends,'start2'=start2,'end2'=ends2)
	else out<-data.frame('start'=start,'end'=ends)
	out<-out[order(out$start),]
	return(out)
}

#Convert index to interesting ranges joining small gaps 
#index: index to join up
#buffer: extra bit to include at start and end
#bufferMultiple: condense adjacent ranges within buffer * bufferMultiple + bufferAdd of each other
#bufferAdd: condense adjacent ranges within buffer * bufferMultiple + bufferAdd of each other
findInterestingRanges<-function(index,buffer=6,bufferMultiple=3,bufferAdd=0){
	if(is.logical(index))index<-which(index)
	ranges<-index2range(index)
	condensed<-condenseCoords(ranges$start,ranges$end-ranges$start+1,buffer*bufferMultiple+bufferAdd)
	condensed$start<-condensed$start-buffer
	condensed$end<-condensed$end+buffer
	return(condensed)
}

#convert starts and ends of matches to gaps between those matches
#x: a 2xX matrix/dataframe of start and end coords with columns $start and $end
findCoordGaps<-function(x){
	if(nrow(x)>1)return(data.frame('start'=x$end[-nrow(x)]+1,'end'=x$start[-1]-1))
	else return(data.frame('start'=-99999,'end'=-999999)[0,])
}

#condense several starts and ends to unique ranges
#x: a 2xX matrix/dataframe of start and end coords
#nullReturn: return NULL if NULL or 0 rows? else error
reduceCoords<-function(x,nullReturn=FALSE){
	if(nullReturn&is.null(x))return(NULL)
	allX<-do.call(rbind,x)
	if(nullReturn&nrow(allX)==0)return(NULL)
	if(nrow(allX)==1)return(allX)
	indices<-unique(unlist(apply(allX,1,function(x)x[1]:x[2])))
	return(index2range(indices))
}


#convenience function 
#x: coords to standardize to 0-1
#range: range to standardize to
standardizeCoords<-function(x,range)(x-range[1])/diff(range)

#add polygons connecting two horizontal lines showing matches between two sequences
#tStarts: comma separated string of start coords for first set
#qStarts: comma separated string of start coords for second set
#blockSizes: comma separated string of length of matches
#tRange: 2 element vector of range for first set
#qRange: 2 element vector of range for second set
#yPos: 2 element vector y positions to connect
#condenseLimit: connect matches not separated by more than condenseLimit
#...: extra arguments for polygon e.g. col
connectGenomes<-function(tStarts,qStarts,blockSizes,tRange,qRange,yPos=c(1,2),condenseLimit=5,...){
	tStart<-as.numeric(strsplit(tStarts,',')[[1]])
	qStart<-as.numeric(strsplit(qStarts,',')[[1]])
	blocks<-as.numeric(strsplit(blockSizes,',')[[1]])
	nBlocks<-length(tStart)
	if(nBlocks!=length(qStart)|nBlocks!=length(blocks))stop(simpleError('Start and block size numbers do not match up'))
	pieces<-condenseCoords(tStart,blocks,condenseLimit,start2=qStart)
	pieces$tStart<-standardizeCoords(pieces$start,tRange)
	pieces$tEnd<-standardizeCoords(pieces$end,tRange)
	pieces$qStart<-standardizeCoords(pieces$start2,qRange)
	pieces$qEnd<-standardizeCoords(pieces$end2,qRange)
	mapply(function(tS,qS,tE,qE,yPos)polygon(c(tS,qS,qE,tE),c(yPos[1:2],yPos[2:1]),...),pieces$tStart,pieces$qStart,pieces$tEnd,pieces$qEnd,MoreArgs=list(yPos))
	return(NULL)
}


#blat: a blat data.frame e.g. from readBlat
#tRange: 2 element vector of range for first set
#qRange: 2 element vector of range for second set
#names: 2 element vector of names to display on left axis
#main: title for plot
drawTwoCoords<-function(blat,tRange=c(1,blat$tSize[1]),qRange=c(1,blat$qSize[1]),names=blat[1,c('tName','qName')],main=''){
	plot(1,1,type='n',xlim=c(0,1),ylim=c(1,2),xaxt='n',yaxt='n',ylab='',xlab='',xaxs='i',yaxs='i',las=1,main=main)
	axis(2,1:2,names,las=1)
	prettyT<-pretty(tRange)
	axis(1,standardizeCoords(prettyT,tRange),prettyT)
	#segments(c(0,0),1:2,c(1,1),1:2)
	#prettyRead<-seq(readRange[1],readRange[2],100)
	prettyQ<-pretty(qRange)
	axis(3,standardizeCoords(prettyQ,qRange),prettyQ)
	#segments(prettyStandardRead,rep(1.98,length(prettyRead)),prettyStandardRead,rep(2.02,length(prettyRead)))
	apply(blat,1,function(x)connectGenomes(x['tStarts'],x['qStarts'],x['blockSizes'],tRange,qRange,col=ifelse(x['strand']=='+','#FF000044','#0000FF44'),border=NA,yPos=c(1,2)))
}

#look for primers in sequences
#seqs: dna sequences
#barcode: barcodes to look for
#primer: primers to look for 
#id: id to assign for barcode-primer matches (can be a dataframe with same # of rows as barcodes to return multiple columns)
#returns: vector of ids (or NA if no match) or data.frame of same columns as id (if id is a dataframe)
#example: fa<-cbind(fa,findPrimers(fa$seq,primers$barcode,primers$primerHybrid,primers[,c('primer','base','remove','tissue','sample','dir')]))
#example: fa$primer<-findPrimers(fa$seq,primers$barcode,primers$primerHybrid,primers$name) 
findPrimers<-function(seqs,barcode,primer=rep('',length(barcode)),id,vocal=FALSE){
	if(length(seqs)<1) stop(simpleError("No sequences given"))
	numPrimers<-length(primer)
	#is ids a data frame or a vector?
	#dfIds<-is.data.frame(id)
	if(!is.data.frame(id))id<-data.frame('id'=id,stringsAsFactors=FALSE)
	if(length(barcode)!=numPrimers|nrow(id)!=numPrimers) stop(simpleError("Length of barcode, primer and id not equal"))
	primerA454=toupper("gcctccctcgcgccatcag");
	primerB454=toupper("gccttgccagcccgctcag");
	seqs<-toupper(seqs)
	#Not sure this is actually going to find anything or not
	if(any(grep(paste('^',primerA454,sep=""),seqs))|any(grep(paste('^',primerB454,sep=""),seqs))) stop(simpleError("454 primers still attached to start of sequence"))
	regexprs<-paste('^',barcode,ambigous2regex(primer),sep="")
	print(regexprs)
	if(any(table(regexprs)>1)){
		print(table(regexprs)[table(regexprs)>1])
		browser()
		stop(simpleError('Found duplicate barcode-primers'))
	}
	fakeRow<-id[1,,drop=FALSE]
	fakeRow[1,]<-NA
	output<-fakeRow[rep(1,length(seqs)),,drop=FALSE]
	for (i in 1:numPrimers){
		#check for primer
		selector<-grep(regexprs[i],seqs)
		if(vocal)message('Primer ',id[i,1],' found ',length(selector),' sequences (',regexprs[i],')')
		if(any(selector)){
			output[selector,]<-id[rep(i,length(selector)),,drop=FALSE]
		}
	}
	#need to check for reaching the opposite primer?
	if(vocal)message('Found ',sum(!is.na(output[,1])),' sequences out of ',length(seqs),' possible')
	#reduce dimensions if a single column output
	if(ncol(output)==1)output<-output[,1]
	return(output)
}

#calculate the number of combinations possible from nGroups groups with groupsize members and n observations
#nGroups: number of equally sized groups
#groupSize: size of the equally sized groups
#n: number of observations pulled from the individuals
#returns: number of combinations containing at least one individual from each group
chooseAtLeastOneFromEach<-function(nGroups,groupSize,n){
	if(nGroups<=0|nGroups>n)return(0)
	if(nGroups>1){
		children<-sapply(1:(nGroups-1),function(x)chooseAtLeastOneFromEach(nGroups-x,groupSize,n))
		nChildren<-sapply(1:(nGroups-1),function(x)choose(nGroups,x))
	}else{
		nChildren<-0
		children<-0
	}
	answer<-choose(nGroups*groupSize,n)-sum(children*nChildren)
	return(answer)
}

#calculate the probability of observing observedX species from nGroups groups with groupsize members with n observations
#nGroups: number of equally sized groups
#groupSize: size of the equally sized groups
#n: number of observations pulled from the individuals
#observedX: number of species observed
#returns: probability of observing observedX species from nGroups groups with groupsize members with n observations
pRare<-function(nGroups,groupSize,n,observedX){
	choose(nGroups,observedX)*chooseAtLeastOneFromEach(observedX,groupSize,n)/choose(nGroups*groupSize,n)
}

#calculate the Wilson score interval http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
#nTrue: number of trues observed
#nFalse: number of falses observed
#alpha: error percentile
#returns: lower, upper bounds of the interval
wilsonInt<-function(nTrue,nFalse,alpha=.05){
	n<-nTrue+nFalse
	prop<-nTrue/n
	z<-pnorm(1-alpha/2)
	plusMinus<-z*sqrt(prop*(1-prop)/n+z^2/4/n^2)
	return((prop+1/2/n*z^2+c(-plusMinus,plusMinus))/(1+1/n*z^2))
}


#rainbow based on circle of lab color space (to equalize intensity I hope)
#n: number of colors desired
#start: start angle (in proportion of circle) in lab space
#end: end angle (in proportion of circle) in lab space
rainbow.lab<-function(n,start=1.5,end=-3,alpha=1,lightScale=.5,lightMultiple=1){
	#something is going crazy with R's implementation of lab
	#angles<-seq(start*2*pi,end*2*pi,length.out=n)
	#a<-sin(angles)*radius
	#b<-cos(angles)*radius
	#Lab<-cbind(luminance,a,b)
	#srgb<-convertColor(Lab,from="Lab",to="sRGB")
	#cols<-rgb(srgb[,1],srgb[,2],srgb[,3],alpha)
	#return(cols)
	l<-seq(1,.4,length.out=n)^lightScale*lightMultiple
	rgb<-cl2pix(seq(0,1,length.out=n),l,start=start,end=end,toColor=FALSE)
	return(rgb(rgb,alpha=alpha))
}



#http://davidad.net/colorviz/
#/* Convert from L*a*b* doubles to XYZ doubles
#* Formulas drawn from http://en.wikipedia.org/wiki/Lab_color_spaces
#*/
lab2xyz<-function(lab) {
	if(is.null(dim(lab)))matrix(lab,nrow=1)
	finv<-function(t) ifelse(t>6/29,t^3,3*(6/29)^2*(t-4/29))
	sl = (lab[,1]+0.16)/1.16
	ill = c(0.9643,1.00,0.8251) #D50
	out<-matrix(ill,nrow=nrow(lab),ncol=3,byrow=TRUE) * finv(cbind(sl,sl+lab[,2]/5,sl-lab[,3]/2))
	return(out)
}

#/* Convert from XYZ doubles to sRGB bytes
#* Formulas drawn from http://en.wikipedia.org/wiki/Srgb
#*/
xyz2rgb<-function(xyz) {
	rgb <- cbind(3.2406*xyz[,1] - 1.5372*xyz[,2] - 0.4986*xyz[,3], -0.9689*xyz[,1] + 1.8758*xyz[,2] + 0.0415*xyz[,3], 0.0557*xyz[,1] - 0.2040*xyz[,2] + 1.0570*xyz[,3])
	clip <- rgb>1|rgb<0
	if(any(clip)) {
		rgb[rgb>1]<-1
		rgb[rgb<0]<-0
	}
	#Uncomment the below to detect clipping by making clipped zones red.
	#if(clip) {rl=1.0;gl=bl=0.0;}
	correct<-function(cl) {
		a <- 0.055
		return (ifelse(cl<=0.0031308,12.92*cl,(1+a)*cl^(1/2.4)-a))
	}
	out<-correct(rgb)
	return(out)
}

#/* Convert from LAB doubles to sRGB bytes 
#* (just composing the above transforms)
#*/
lab2rgb<-function(lab){
  xyz<-lab2xyz(lab)
  xyz2rgb(xyz)
}

#/* Convert from a qualitative parameter c and a quantitative parameter l to a 24-bit pixel
 #* These formulas were invented by me to obtain maximum contrast without going out of gamut
 #* if the parameters are in the range 0-1
 #*/
cl2pix<-function(c, l,start=-3,end=4,toColor=TRUE) {
  L <- l #L of L*a*b*
  angle <- start+c*(end-start)
  r <- l*0.4+0.1 #~chroma
  a <- sin(angle+start)*r
  b <- cos(angle+start)*r
  out<-lab2rgb(cbind(L,a,b))
  if(toColor)out<-rgb(out)
  return(out)
}

# find coords for arrow plotting
#left: left coordinate of block
#right: right coordinate of block
#y: y position for middle of block
#arrowLength: arrow length in usr coords
#shaft: half of shaft thickness in usr coords
#point: half of arrow thickness in usr coords
#concat: concatenate multiple arrows into one data frame separated by NAs (ready for poly)?
arrow<-function(left,right,y,arrowLength=diff(par('usr')[1:2])*.05,shaft=.2,point=.4,concat=TRUE){
	if(any(left>right))stop(simpleError('Left border > right border of arrow'))
	arrowX<-right-arrowLength
	arrowX<-ifelse(arrowX<left,left+(right-left)/10,arrowX)
	coords<-mapply(function(left,right,y,arrowX){data.frame('x'=c(left,arrowX,arrowX,right,arrowX,arrowX,left),'y'=y+c(shaft,shaft,point,0,-point,-shaft,-shaft))},left,right,y,arrowX,SIMPLIFY=FALSE)
	if(concat)coords<-do.call(rbind,lapply(coords,function(x)return(rbind(x,c(NA,NA)))))
	return(coords)
}

#stack regions into smallest number of lines (using greedy algo)
#starts: vector of starts of regions (note can add buffer by substracting arbitrary spacer)
#ends: vector of ends of regions (note can add buffer by adding arbitrary spacer)
stackRegions<-function(starts,ends){
	nReg<-length(starts)
	if(nReg!=length(ends))stop(simpleError('Starts and ends must be same length'))
	startEnds<-data.frame('id'=1:nReg,start=starts,end=ends)
	startEnds<-startEnds[order(startEnds$start,startEnds$end),]
	lineNum<-rep(NA,nReg)
	linePos<-rep(min(starts)-1,nReg)
	for(i in 1:nReg){
		selectLine<-min(which(linePos<startEnds$start[i]))
		lineNum[i]<-selectLine
		linePos[selectLine]<-max(startEnds$end[i],linePos[selectLine])
	}
	return(lineNum[order(startEnds$id)])
}


#chr: vector of chromosomes
#start: vector of start coordinate 1-based
#end: vector of end coordinate 1-based
#strand: vector of strands
#chainFile: location of chain file for liftover
#liftoverBin: location of the liftOver executable
#return: 4 column data frame of coordinates
liftCoords<-function(chr,start,end,strand,chainFile,liftoverBin='liftOver',vocal=FALSE,removeNAs=TRUE){
	tmpFiles<-c(tempfile(),tempfile(),tempfile())
	y<-data.frame('chr'=chr,'start'=start,'end'=end,'strand'=strand,stringsAsFactors=FALSE)
	y$id<-1:nrow(y)
	writeLines(sprintf('%s\t%d\t%d\t%d\t%d\t%s',y$chr,y$start-1,y$end,y$id,1,strand),tmpFiles[1])
	cmd<-sprintf('%s %s %s %s %s',liftoverBin,tmpFiles[1],chainFile,tmpFiles[2],tmpFiles[3])
	if(vocal)message(cmd)
	returnCode<-system(cmd)
	if(vocal)message('Return: ',returnCode)
	lift<-read.table(tmpFiles[2],stringsAsFactors=FALSE)
	lift<-lift[,-5]
	colnames(lift)<-c('chr','start','end','id','strand')
	if(any(table(lift$id)>1))warning('Liftover created duplicates')
	y<-merge(y[,'id',drop=FALSE],lift,all.x=TRUE)
	y<-y[,colnames(y)!='id']
	y$start<-y$start+1
	selector<-is.na(y$chr)|is.na(y$start)|is.na(y$end)
	if(removeNAs&&any(selector)){
		message("Removing ",sum(selector)," of ",length(selector)," reads for failing to liftover")
		y<-y[!selector,]
	}
	return(y)
}

#seqs: sequences to get GC content of
#chars: characters to count (G,C by default)
#returns: proportion of GC in sequence
gcPercent<-function(seqs,chars=c('C','G')){
	regex<-sprintf('[%s]+',paste(chars,collapse=''))
	gcs<-nchar(gsub(regex,'',seqs,perl=TRUE))
	return(gcs/nchar(seqs))
}

#seqs: vector of strings, sequences to form a position weight matrix from (all same length)
#chars: allowed characters
#priors: additional counts to add to each column with (default is a lazy way to get a named 0 vector)
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

#seqs: vector of strings, sequences to score
#pwm: position weight matrix with uniqueBasesxseqLength dimensions
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

#line: line to convert to user coordinates
#axis: axis to do conversion on (1:4 same as axis, mtext command)
convertLineToUser<-function(line,axis=1){
	if(!(axis %in% 1:4))stop(simpleError('Undefined axis'))
	axisPair<-sort((c(axis-1,axis+1)%%4)+1)
	isHeight<-(axis%%2)==1
	isSecond<-axis>2
	thisMar<-par('mar')[axis]
	marWidth<-thisMar/sum(par('mar')[axisPair])*(par('fin')-par('pin'))[isHeight+1]
	widthPerLine<-marWidth/thisMar
	#find base line + add in if plot doesn't cover whole device e.g. par(mfrow=c(2,1))
	base<-ifelse(isSecond,par('fin')[isHeight+1]-widthPerLine*thisMar,widthPerLine*thisMar) + par('fig')[1+isHeight*2]*par('din')[isHeight+1]
	func<-if(isHeight)grconvertY else grconvertX
	out<-func(base+line*widthPerLine*ifelse(isSecond,1,-1),'inches','user')
	return(out)
}

#data.frame of sam flags
samFlags<-data.frame('short'=c('paired','properPair','unmapped','mateUnmapped','reverse','mateReverse','first','second','notPrimary','fail','dupe'),'desc'=c('read paired','read mapped in proper pair','read unmapped','mate unmapped','read reverse strand','mate reverse strand','first in pair','second in pair','not primary alignment','read fails platform/vendor quality checks','read is PCR or optical duplicate'),stringsAsFactors=FALSE)
samFlags$bit<-2^(0:(nrow(samFlags)-1))
rownames(samFlags)<-samFlags$short
#data.frame of amino acids
aminoAcids<-data.frame('codon'=c('UUU','UUC','UCU','UCC','UAU','UAC','UGU','UGC','UUA','UCA','UAA','UGA','UUG','UCG','UAG','UGG','CUU','CUC','CCU','CCC','CAU','CAC','CGU','CGC','CUA','CUG','CCA','CCG','CAA','CAG','CGA','CGG','AUU','AUC','ACU','ACC','AAU','AAC','AGU','AGC','AUA','ACA','AAA','AGA','AUG','ACG','AAG','AGG','GUU','GUC','GCU','GCC','GAU','GAC','GGU','GGC','GUA','GUG','GCA','GCG','GAA','GAG','GGA','GGG'),'abbr'=c('Phe','Phe','Ser','Ser','Tyr','Tyr','Cys','Cys','Leu','Ser','Ochre','Opal','Leu','Ser','Amber','Trp','Leu','Leu','Pro','Pro','His','His','Arg','Arg','Leu','Leu','Pro','Pro','Gln','Gln','Arg','Arg','Ile','Ile','Thr','Thr','Asn','Asn','Ser','Ser','Ile','Thr','Lys','Arg','Met','Thr','Lys','Arg','Val','Val','Ala','Ala','Asp','Asp','Gly','Gly','Val','Val','Ala','Ala','Glu','Glu','Gly','Gly'),'code'=c('F','F','S','S','Y','Y','C','C','L','S','X','X','L','S','X','W','L','L','P','P','H','H','R','R','L','L','P','P','Q','Q','R','R','I','I','T','T','N','N','S','S','I','T','K','R','M','T','K','R','V','V','A','A','D','D','G','G','V','V','A','A','E','E','G','G'),'name'=c('Phenylalanine','Phenylalanine','Serine','Serine','Tyrosine','Tyrosine','Cysteine','Cysteine','Leucine','Serine','Stop','Stop','Leucine','Serine','Stop','Tryptophan','Leucine','Leucine','Proline','Proline','Histidine','Histidine','Arginine','Arginine','Leucine','Leucine','Proline','Proline','Glutamine','Glutamine','Arginine','Arginine','Isoleucine','Isoleucine','Threonine','Threonine','Asparagine','Asparagine','Serine','Serine','Isoleucine','Threonine','Lysine','Arginine','Methionine','Threonine','Lysine','Arginine','Valine','Valine','Alanine','Alanine','Aspartic acid','Aspartic acid','Glycine','Glycine','Valine','Valine','Alanine','Alanine','Glutamic acid','Glutamic acid','Glycine','Glycine'),row.names=c('UUU','UUC','UCU','UCC','UAU','UAC','UGU','UGC','UUA','UCA','UAA','UGA','UUG','UCG','UAG','UGG','CUU','CUC','CCU','CCC','CAU','CAC','CGU','CGC','CUA','CUG','CCA','CCG','CAA','CAG','CGA','CGG','AUU','AUC','ACU','ACC','AAU','AAC','AGU','AGC','AUA','ACA','AAA','AGA','AUG','ACG','AAG','AGG','GUU','GUC','GCU','GCC','GAU','GAC','GGU','GGC','GUA','GUG','GCA','GCG','GAA','GAG','GGA','GGG'),stringsAsFactors=FALSE)



#http://life.nthu.edu.tw/~fmhsu/rasframe/SHAPELY.HTM
aminoColors<-data.frame('aa'=c("ASP","GLU","CYS","MET","LYS","ARG","SER","THR","PHE","TYR","ASN","GLN","GLY","LEU","VAL","ILE","ALA","TRP","HIS","PRO"),'col'=c("#E60A0A","#E60A0A","#E6E600","#E6E600","#145AFF","#145AFF","#FA9600","#FA9600","#3232AA","#3232AA","#00DCDC","#00DCDC","#EBEBEB","#0F820F","#0F820F","#0F820F","#C8C8C8","#B45AB4","#8282D2","#DC9682"),stringsAsFactors=FALSE)
tmp<-with(aminoAcids[!duplicated(aminoAcids[,c('code')]),],data.frame('letter'=code,row.names=toupper(abbr),stringsAsFactors=FALSE))
rownames(aminoColors)<-tmp[aminoColors$aa,'letter']
tmpAngles1<-cos(1+1:nrow(aminoColors)/nrow(aminoColors)*pi)
tmpAngles2<-sin(1:nrow(aminoColors)/nrow(aminoColors)*pi)
tmpAngles1<-tapply(tmpAngles1,aminoColors$col,c)
tmpAngles2<-tapply(tmpAngles2,aminoColors$col,c)
aminoColors$spreadCol<-ave(aminoColors$col,aminoColors$col,FUN=function(x){
	if(length(x)==1)return(x)
	spacer<-20
	angles1<-tmpAngles1[[x[1]]]
	angles2<-tmpAngles2[[x[1]]]
	spacing<-seq((length(x)-1)*-spacer,(length(x)-1)*spacer,length.out=length(x))
	lab<-convertColor(t(col2rgb(x[1])),from='sRGB',to='Lab',scale.in=255)[rep(1,length(x)),]
	lab[,'b']<-lab[,'b']-spacing*angles1
	lab[,'a.x']<-lab[,'a.x']+spacing*angles2
	rgbs<-convertColor(lab,from='Lab',to='sRGB')
	return(rgb(rgbs))
	#plot(1:nrow(aminoColors),col=aminoColors$col,lwd=30,xaxt='n');points(1:nrow(aminoColors),1:nrow(aminoColors)+1,col=aminoColors$spreadCol,lwd=20);axis(1,1:nrow(aminoColors),aminoColors$aa,las=2)
})
rm(tmpAngles1,tmpAngles2)
