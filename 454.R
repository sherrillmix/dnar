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


