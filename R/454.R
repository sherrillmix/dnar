#' Convert a DNA sequence into expected 454 flows
#'
#' @param seq DNA sequence
#' @param flowOrder order of nucleotide flows
#' @param outputLength minimum output length
#' @export
#' @return vector of flows
#' @examples
#' seq2flow('ACTAAAAAT')
seq2flow<-function(seq,flowOrder=c('T','A','C','G'),outputLength=NULL){
	seqSplit<-strsplit(seq,'')[[1]]
	if(any(!seqSplit %in% flowOrder))stop(simpleError('Unknown character found'))
	dif<-seqSplit!=c(seqSplit[-1],'DUMMY')
	strLengths<-diff(c(0,which(dif)))
	chars<-seqSplit[dif]
	nextChars<-c(chars[-1],NA)
	numFlows<-length(flowOrder)
	distMat<-matrix(NA,nrow=numFlows,ncol=numFlows,dimnames=list(flowOrder,flowOrder))
	distMat[,1]<-length(flowOrder):1
	if(numFlows>1)for(i in 2:numFlows)distMat[,i]<-c(distMat[numFlows,i-1],distMat[-numFlows,i-1])
	dists<-rep(NA,length(chars))
	for(i in flowOrder){
		dists[chars==i&is.na(nextChars)]<-NA
		dists[chars==i&!is.na(nextChars)]<-distMat[i,nextChars[chars==i&!is.na(nextChars)]]
	}
	dists<-c(which(flowOrder==chars[1]),dists[-length(dists)])
	if(is.null(outputLength))outputLength<-sum(dists)
	else outputLength<-max(outputLength,sum(dists))
	output<-rep(0,outputLength)
	output[cumsum(dists)]<-strLengths
	names(output)<-rep(flowOrder,length.out=length(output))
	return(output)
}

#' Find flow corresponding to a given base position
#'
#' @param flow vector of flow values (e.g. produced by seq2flow)
#' @param coords 1-based coordinates of desired positions
#' @export
#' @return flow number for each coord. Infinity if coordinate is outside sequence.
#' @examples
#' indexFlow(seq2flow('CCTTAA'),1:6)
#' indexFlow(seq2flow('GTGTG'),1:5)
indexFlow<-function(flow,coords){
	indices<-cumsum(flow)
	if(any(coords<1))stop(simpleError('Negative coords found'))
	output<-sapply(coords,function(x,y)return(min(which(y>=x),Inf)),indices)
	return(output)
}

#' Convert 454 flows to DNA sequences
#'
#' @param flow vector of flow values (e.g. produced by seq2flow)
#' @param flowOrder order of nucleotide flows
#' @export
#' @return DNA sequence corresponding to the flowgram
flow2seq<-function(flow,flowOrder=c('T','A','C','G')){
	if(is.character(flow)&length(flow)==1)flow<-strsplit(flow,'\t')[[1]]
	flow<-as.numeric(flow)
	flow[flow<0]<-0
	chars<-rep(flowOrder,length.out=length(flow))
	output<-paste(rep(chars,round(as.numeric(flow))),collapse='')
	return(output)
}

