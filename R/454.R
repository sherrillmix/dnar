#' Convert a DNA sequence into expected 454 flows
#'
#' @param seq DNA sequence
#' @param flowOrder order of nucleotide flows
#' @param outputLength minimum output length
#' @export
#' @return vector of flows
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

#' Find flow corresponding to a given base position
#'
#' @param flow vector of flow values (e.g. produced by seq2flow)
#' @param coords bp coordinates of desired positions
#' @export
#' @return flow number for each coord
indexFlow<-function(flow,coords){
	indices<-cumsum(flow)
	output<-sapply(coords,function(x,y)return(min(which(y>=x))),indices)
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


