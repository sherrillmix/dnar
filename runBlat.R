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
runBlatNoServer<-function(reads=NULL,refs,blatArgs='',outFile='out.blat',blat='blat',tmpDir=tempfile(),readFile=NULL,deleteFiles=!is.null(reads),gzFile=grepl('\\.gz$',outFile)){
	if(!file.exists(tmpDir))dir.create(tmpDir,recursive=TRUE)
	if(is.null(reads)&is.null(readFile))stop(simpleError('Please provide reads or readFile'))
	if(is.null(readFile)){
		readFile<-sprintf('%s/read.fa',tmpDir)
		write.fa(names(reads),reads,readFile)
	}
	refFile<-sprintf('%s/refs.fa',tmpDir)
	write.fa(names(refs),refs,refFile)
	if(gzFile)outFile<-sub('\\.gz$','',outFile)
	cmd<-sprintf('%s %s %s %s %s',blat,refFile,readFile,blatArgs,outFile)
	message(cmd)
	errorCode<-system(cmd)
	if(errorCode!=0)stop(simpleError('Problem running blat'))
	if(deleteFiles)file.remove(refFile,readFile,tmpDir)
	if(gzFile){
		system(sprintf('gzip %s',outFile))
		if(!file.exists(sprintf('%s.gz',outFile)))stop(simpleError('Problem gzipping file'))
	}
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
multiRunBlatNoServer<-function(reads,refs,outFile,nCore=4,tmpDir=tempdir(),condense=TRUE,isGz=grepl('.gz$',outFile),sleepIncrement=1,runFilter=NULL,...){
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
    if(any(isError(bigRun))){
        print(bigRun)
        stop(simpleError("Error in running blat"))
    }
	if(any(!sapply(bigRun,file.exists)))stop(simpleError('Blat file missing'))
	if(!is.null(runFilter)){
		bigRun<-mclapply(bigRun,function(x,runFilter){
			newFile<-sprintf('%s__2',x)
			message('Filtering ',x)
			system(sprintf('cat %s|%s>%s',x,runFilter,newFile))
			file.remove(x)
			return(newFile)
		},runFilter,mc.cores=nCore)
        if(any(isError(bigRun))){
            print(bigRun)
            stop(simpleError("Error in filtering blat"))
        }
	}
	if(isGz)outFile<-gzfile(outFile,open='w+')
	else outFile<-file(outFile,open='w+')
	counter<-0
	for(i in 1:length(bigRun)){
		message('Reading blat file ',i)
		blat<-readLines(bigRun[[i]])
		#take off header in later files
		if(i!=1&is.null(runFilter)){
			blat<-blat[-1:-5]
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
