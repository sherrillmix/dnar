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

#' Break 454 file(s) into separate files for each sample
#'
#' @param names names of reads e.g. >DRDR12A125
#' @param samples sample ID for each above read
#' @param barcodes list of vectors of barcodes or barcode/primers with each entry indexed by sample
#' @param outDir directory to put seperate sffs
#' @param sffDir directory to look for sffs
#' @param baseName prepend to sample names 
#' @param sffBinDir location of sfffile binary
#' @param vocal if TRUE output status messages
#' @export
#' @return NULL
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
	return(NULL)
}
