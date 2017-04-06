#' Start a blat server
#'
#' @param nibDir directory containing nib files to align against
#' @param options options to gfServer
#' @param gfServer path to gfServer program
#' @param bit2 a bit2 file alternative to nib directory
#' @param startPort begin looking for open (no other R called gfServer using) ports on startPort and add 1 until finding one
#' @param nibSuffix nib files end in nibSuffix
#' @param wait number of 5 second intervals to wait for gfServer to finish starting before giving up
#' @export
#' @return port of blat server
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

#' Check if blat is running on port
#' 
#' @param port port to check if blat is running on
#' @export
#' @return TRUE if blat running on port otherwise FALSE
checkBlat<-function(port){
	return(system(sprintf('pgrep -f "gfServer *start *localhost *%d">/dev/null',port))==0)
}

#' Start a blat server, run blat on file and kill blat server
#' @param faFile path to fasta file for alignment
#' @param outFile file to write blat results to
#' @param gfClient path to gfClient program
#' @param gfClientOptions arguments for gfClient
#' @param ... arguments for startBlat
#' @export
#' @return NULL
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
	return(NULL)
}


#' Run blat from blat executable instead of gfserver/client
#'
#' @param reads named vector of sequences
#' @param refs named vector of references
#' @param blat path to blat program
#' @param blatArgs string of arguments for blat
#' @param tmpDir a directory to write tempfiles to
#' @param outFile file to write blat to
#' @param readFile file to use instead of reads (will be deleted)
#' @export
#' @return NULL
runBlatNoServer<-function(reads=NULL,refs,blatArgs='',outFile='out.blat',blat='blat',tmpDir=tempdir(),readFile=NULL){
	if(!file.exists(tmpDir))dir.create(tmpDir,recursive=TRUE)
	if(is.null(reads)&is.null(readFile))stop(simpleError('Please provide reads or readFile'))
	tmpReadFile<-sprintf('%s/read.fa',tmpDir)
	if(is.null(readFile)){
		write.fa(names(reads),reads,tmpReadFile)
	}else{
		file.symlink(readFile,tmpReadFile)
	}
	refFile<-sprintf('%s/refs.fa',tmpDir)
	write.fa(names(refs),refs,refFile)
	cmd<-sprintf('%s %s %s %s %s',blat,refFile,tmpReadFile,blatArgs,outFile)
	message(cmd)
	errorCode<-system(cmd)
	if(errorCode!=0)stop(simpleError('Problem running blat'))
	file.remove(refFile,tmpReadFile,tmpDir)
	return(NULL)
}

#' Run blat from parallel blat executables
#'
#' @param reads named vector of sequences
#' @param refs named vector of references
#' @param nCore number of cores to use
#' @param tmpDir a directory to write tempfiles to
#' @param outFile file to write blat to (if ends in .gz then isGz defaults to true and writes to gzipped file)
#' @param isGz if TRUE write to gz file else normal file
#' @param ... additional arguments to runBlatNoServer
#' @export
#' @return NULL
multiRunBlatNoServer<-function(reads,refs,outFile,nCore=4,tmpDir=tempdir(),isGz=grepl('.gz$',outFile),...){
	prefix<-paste(sample(c(letters,LETTERS),20,TRUE),collapse='')
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
	bigRun<-parallel::mclapply(runFiles,function(x){
		message('Starting ',x[2])
		runBlatNoServer(readFile=x[3],refs=refs,outFile=x[2],tmpDir=x[1],...)
		return(x[2])
	},mc.cores=nCore)
    if(any(isError(bigRun))){
        print(bigRun)
        stop(simpleError("Error in running blat"))
    }
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
		}
		writeLines(blat,sep="\n",con=outFile)
		file.remove(bigRun[[i]])
		counter<-counter+length(blat)
	}
	message('Wrote ',counter,' blat lines')
	close(outFile)
	return(NULL)
}

#' Make a 2bit file from one set of reads and blat another set against it
#'
#' @param reads vector of query reads with names
#' @param refs vector of reference reads with names
#' @param outFile file to write blat to (if ends in .gz then isGz defaults to true and writes to gzipped file). If NULL then output to temp file, read in and return.
#' @param faToTwoBit path to faToTwoBit program from blat
#' @param tmpDir directory to store working files
#' @param ... additional arguments to run blat
#' @export
#' @return NULL if outFile defined or an alignment data.frame if outFile is NULL
blatReadsVsRefs<-function(reads,refs,outFile=NULL,faToTwoBit='faToTwoBit',tmpDir=tempfile('dir'),...){
	dir.create(tmpDir)
  if(is.null(outFile)){
    isTmp<-TRUE
    outFile<-tempfile()
  }else{
    isTmp<-FALSE
  }
  if(is.null(names(reads)))names(reads)<-1:length(reads)
	readFile<-sprintf('%s/read.fa',tmpDir)
	write.fa(names(reads),reads,readFile)
	refFile<-sprintf('%s/refs.fa',tmpDir)
  if(is.null(names(refs)))names(refs)<-1:length(refs)
	write.fa(names(refs),refs,refFile)
	twobitFile<-sprintf('%s/refs.2bit',tmpDir)
	system(sprintf('%s %s %s',faToTwoBit,refFile,twobitFile))
	runBlat(readFile,outFile=outFile,bit2=twobitFile,...)
	file.remove(twobitFile,refFile,readFile,tmpDir)
  if(!isTmp)return(NULL)
  out<-read.blat(outFile)
  file.remove(outFile)
  return(out)
}

#' Run blatReadsVsRefs in parallel
#' @param reads vector of query reads with names
#' @param refs vector of query refs with names
#' @param outFile file to write blat to (if ends in .gz then isGz defaults to true and writes to gzipped file)
#' @param nCore number of cores to use
#' @param tmpDir directory to store work files
#' @param isGz if TRUE write to gz file else normal file
#' @param ... additional arguments for blatReadsVsRefs
#' @export
#' @return NULL
multiBlatReadsVsRefs<-function(reads,refs,outFile,nCore=4,tmpDir=tempdir(),isGz=grepl('.gz$',outFile),...){
	prefix<-paste(sample(c(letters,LETTERS),20,TRUE),collapse='')
	if(!file.exists(tmpDir))dir.create(tmpDir)
	bigRun<-parallel::mclapply(mapply(function(x,y)list(x,y),split(reads,sort(rep(1:nCore,length.out=length(reads)))),1:nCore,SIMPLIFY=FALSE),function(x){
		tmpFile<-sprintf('%s__%d.blat',sub('\\.blat(\\.gz)?$','',outFile),x[[2]])
		thisTmpDir<-sprintf('%s/work__%s__%d',tmpDir,prefix,x[[2]])
		blatReadsVsRefs(x[[1]],refs,tmpFile,startPort=37900+x[[2]]*10,tmpDir=thisTmpDir,...) #likes to start on same port apparently
		return(tmpFile)
	},mc.cores=nCore)
	if(any(!sapply(bigRun,file.exists)))stop(simpleError('Blat file missing'))
	#condense files
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
	return(NULL)
}

#' Kill blat running on a given port
#'
#' @param port pkill blat running on port port
#' @export
#' @return NULL
killBlat<-function(port){
	if(!checkBlat(port)){
		stopError('Blat does not appear to be running on port ',port)
	}
	code<-system(sprintf('pkill -f "gfServer *start *localhost *%d">/dev/null',port))
	if(checkBlat(port)){
		stopError('Could not kill blat on port ',port)
	}
	return(NULL)
}

#' Run liftover to convert coordinates from one genome version to another
#'
#' @param chr vector of chromosomes
#' @param start vector of start coordinate 1-based
#' @param end vector of end coordinate 1-based
#' @param strand vector of strands
#' @param chainFile location of chain file for liftover
#' @param liftoverBin location of the liftOver executable
#' @param vocal if TRUE then message extra information
#' @param removeNAs if TRUE then remove positions that did could not be lifted over
#' @export
#' @return four column data frame of chr, start, end and strand of coordinates
liftCoords<-function(chr,start,end,strand,chainFile,liftoverBin='liftOver',vocal=FALSE,removeNAs=TRUE){
	tmpFiles<-c(tempfile(),tempfile(),tempfile())
	y<-data.frame('chr'=chr,'start'=start,'end'=end,'strand'=strand,stringsAsFactors=FALSE)
	y$id<-1:nrow(y)
	writeLines(sprintf('%s\t%d\t%d\t%d\t%d\t%s',y$chr,y$start-1,y$end,y$id,1,strand),tmpFiles[1])
	cmd<-sprintf('%s %s %s %s %s',liftoverBin,tmpFiles[1],chainFile,tmpFiles[2],tmpFiles[3])
	if(vocal)message(cmd)
	returnCode<-system(cmd)
	if(vocal)message('Return: ',returnCode)
	lift<-utils::read.table(tmpFiles[2],stringsAsFactors=FALSE)
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


