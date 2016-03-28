#' convenience function to resize R console window
#'
#' @export
adjustWindow<-function()options(width=as.integer(Sys.getenv('COLUMNS')))

#' Convenience function to list objects by size
#'
#' @param env the environment to list
#' @export
#' @return object sizes in decreasing order
object.sizes<-function(env=.GlobalEnv)sort(sapply(ls(envir=env),function(x)object.size(get(x))),decreasing=TRUE)

#' Generate an error
#'
#' Convenience function for stop(simpleError()
#' @param ... Strings to be concatenated into error message
#' @export
#' @return Generates stopping error before returning
stopError<-function(...){
	stop(simpleError(paste(...,sep='')))
}

#' Get the conservative edge of a confidence interval
#'
#' @param boundaries upper and lower values
#' @param base value, e.g. 0 or 1
#' @export
#' @return single numberic. if the boundaries overlap the base then return base otherwise return value closest to base
conservativeBoundary<-function(boundaries,base=0){
	boundaries<-sort(boundaries)
	return(ifelse(all(boundaries>base),boundaries[1],ifelse(all(boundaries<base),boundaries[2],base)))
}


#' Convenience function for lagging with NA
#'
#' @param x a vector to be lagged
#' @param lag an integer for how much to lag x. lag greater than zero moves values left. lag less than zero moves values right
#' @param fill the value to fill at the start or end  of the output lagged vector
#' @export
#' @return a vector the same length as x 
lagNA<-function(x,lag=1,fill=NA){
	out<-x	
	if(lag>0)out<-c(out[-(1:lag)],rep(fill,lag))
	if(lag<0)out<-c(rep(fill,abs(lag)),out[-(length(x)-(1:lag)+1)])
	return(out)
}

#' Calculating a moving average/max/min statistic over a vector
#'
#' @param vec to average over
#' @param statFunc function to use
#' @param spacer how many neighbors to the left and right to average
#' @export
#' @return a vector the same length as vec
movingStat<-function(vec,statFunc=max,spacer=2){
	n<-length(vec)
	sapply(1:n,function(x)statFunc(vec[max(1,x-spacer):min(n,x+spacer)]))
}

#' Convenience function to check all args are same length
#'
#' @param ... arguments to check the lengths of 
#' @export
#' @return logical indicating whether all input are the same length
allSameLength<-function(...){
	args<-list(...)
	ns<-sapply(args,length)
	return(all(ns==ns[1]))
}


#' Convenience function to check whether any  args have NAs
#' 
#' @param ... arguments to check if contain NA
#' @export
#' @return logical indicating whether any argument contains NA
anyNa<-function(...){
	args<-list(...)
	hasNa<-sapply(args,function(x)any(is.na(x)))
	return(any(hasNa))
}

#' Cache an operation to disk or load if cache file present
#'
#' @param cacheFile File location to save data to 
#' @param operation A function taking ... arguments and returning object to be stored
#' @param ... Arguments for operation function
#' @param OVERWRITE If FALSE throw an error if hash of ... changes from cached values. If TRUE redo operation and overwrite cache without asking
#' @param VOCAL If TRUE report on status of caching
#' @param EXCLUDE Vector of names of arguments to exclude from md5 digest comparison (for very large arguments)
#' @export
#' @return Output from operation function with ... arguments
cacheOperation<-function(cacheFile,operation,...,OVERWRITE=FALSE,VOCAL=TRUE,EXCLUDE=NULL){
	#avoid evaluation EXCLUDEd args until necessary since they're probably big
	unevalArgs<-match.call(expand.dots=FALSE)$'...'
	varSelector<-if(is.null(names(unevalArgs)))rep(TRUE,length(unevalArgs)) else !names(unevalArgs) %in% EXCLUDE
	allArgs<-lapply(unevalArgs[varSelector],eval)
	#try to prevent function scope from changing things
	md5<-digest::digest(lapply(allArgs,function(x)if(is.function(x))deparse(x)else x))
	if(file.exists(cacheFile)){
		tmp<-new.env()
		load(cacheFile,envir=tmp)
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

#' Check predictions of glm by cross validation
#' @param model a glm to cross validate
#' @param K number of pieces to split data into
#' @param nCores number of cores to use
#' @param subsets predefined subset of same length as data specifying groupings
#' @param vocal echo progress indicator
#' @export
#' @return a dataframe giving predictions for each data point with columns prediction and subsetId
cv.glm.par<-function(model,K=nrow(thisData),nCores=1,subsets=NULL,vocal=TRUE){
	modelCall<-model$call
	thisData<-eval(modelCall$data)
	n<-nrow(thisData)
	if(is.null(subsets))subsets<-split(1:n,sample(rep(1:K,length.out=n)))
	preds<-parallel::mclapply(subsets,function(outGroup){
		if(vocal)cat('.')
		subsetData<-thisData[-outGroup,,drop=FALSE]
		predData<-thisData[outGroup,,drop=FALSE]
		thisModel<-modelCall
		thisModel$data<-subsetData
		return(predict(eval(thisModel),predData))
	},mc.cores=nCores)
	pred<-unlist(preds)[order(unlist(subsets))]
	subsetId<-rep(1:K,sapply(subsets,length))[order(unlist(subsets))]
	return(data.frame(pred,subsetId))
}

#' A function to break an mclapply into parts, save to disk and run independently
#' 
#' @param x vector/list to apply over
#' @param mc.cores number of cores to use
#' @param applyFunc function to apply
#' @param extraCode character vector of setup code (each command self-contained within a cell or concatenated with ;)
#' @param ... arguments for mclapply
#' @param nSplits number of splits to make (if > mc.cores then R will restart more frequently)
#' @export
#' @return the concatenated outputs from applyFunc 
cleanMclapply<-function(x,mc.cores,applyFunc,...,extraCode='',nSplits=mc.cores){
	if(nSplits<mc.cores)nSplits<-mc.cores
	splits<-unique(round(seq(0,length(x),length.out=nSplits+1)))
	if(length(splits)<nSplits+1)nSplits<-length(splits)-1 #not enough items to fill so set lower
	dotVars<-match.call(expand.dots=FALSE)$'...'
	extraArgs<-lapply(dotVars,eval)
	files<-c()
	outFiles<-c()
	scriptFiles<-c()
	logFiles<-c()
	for(ii in 1:nSplits){
		message("Writing core ",ii," data")
		thisInRdat<-tempfile()
		thisRScript<-tempfile()
		thisOutRdat<-tempfile()
		thisLog<-tempfile()
		THISDATA__<-x[(splits[ii]+1):splits[ii+1]]
		SAVEDATA__<-c('THISDATA__'=list(THISDATA__),'applyFunc'=applyFunc,extraArgs)
		save(SAVEDATA__,file=thisInRdat)
		script<-sprintf('load("%s")\nwith(SAVEDATA__,{%s})\nout<-with(SAVEDATA__,lapply(THISDATA__,applyFunc%s%s))\nsave(out,file="%s");',thisInRdat,paste(extraCode,collapse=';'),ifelse(length(extraArgs)>0,',',''),paste(names(extraArgs),names(extraArgs),sep='=',collapse=','),thisOutRdat)
		writeLines(script,thisRScript)
		outFiles<-c(outFiles,thisOutRdat)		
		scriptFiles<-c(scriptFiles,thisRScript)		
		logFiles<-c(logFiles,thisLog)		
	}
	message("Running")
	message("Logs: ",paste(logFiles,collapse=', '))
	exitCode<-parallel::mclapply(mapply(c,scriptFiles,logFiles,SIMPLIFY=FALSE),function(x){out<-system(sprintf("R CMD BATCH --no-save --no-restore %s %s",x[1],x[2]));cat('.');return(out)},mc.cores=mc.cores)
	cat('\n')
	if(any(exitCode!=0)){message('Problem running multi R code');browser()}
	message("Loading split outputs")
	out<-do.call(c,lapply(outFiles,function(outFile){load(outFile);return(out)}))
	return(out)
}


#' Find the  order of one vector in another
#'
#' @param query values to be sorted in target order
#' @param target order for query to be sorted into
#' @param strict if true error if query not in target, if false append unknown queries to end of target
#' @param orderFunc function to use for ordering
#' @param ... arguments for order
#' @export
#' @return a vector the same length as query specifying the order of query elements in target
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


#' Loop through list, get unique names and make sure every element has those names
#'
#' @param x list to loop through
#' @param namesList names to make sure every element has (and delete extras)
#' @param fill value to insert in missing elements
#' @export
#' @return list with missing elements filled in and elements not in namesList removed
fillList<-function(x,namesList=unique(unlist(lapply(x,names))),fill=NA){
	output<-lapply(x,function(x){
		x<-x[names(x) %in% namesList]
		x[namesList[!namesList %in% names(x)]]<-fill
		return(x)
	})	
	return(output)
}


#' Convert from axis line to user coordinates
#'
#' @param line line to convert to user coordinates
#' @param axis axis to do conversion on (1:4 same as axis, mtext command)
#' @export
#' @return position in user coordinates
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

#' Convert from user coordinates to axis line
#'
#' @param usr user coordinate to convert to line
#' @param axis axis to do conversion on (1:4 same as axis, mtext command)
#' @export
#' @return axis line position
convertUserToLine<-function(usr,axis=1){
	if(!(axis %in% 1:4))stop(simpleError('Undefined axis'))
	axisPair<-sort((c(axis-1,axis+1)%%4)+1)
	isHeight<-(axis%%2)==1
	isSecond<-axis>2
	thisMar<-par('mar')[axis]
	marWidth<-thisMar/sum(par('mar')[axisPair])*(par('fin')-par('pin'))[isHeight+1]
	widthPerLine<-marWidth/thisMar
	#find base line + add in if plot doesn't cover whole device e.g. par(mfrow=c(2,1))
	plotWidth<-par('fin')[isHeight+1]
	func<-if(isHeight)grconvertY else grconvertX
	usrInches<-func(usr,'user','inches')
	base<-ifelse(isSecond,par('fin')[isHeight+1]-widthPerLine*thisMar,widthPerLine*thisMar) + par('fig')[1+isHeight*2]*par('din')[isHeight+1]
	out<-(usrInches-base)/widthPerLine*ifelse(isSecond,1,-1)
	return(out)
}

#' Find coords for arrow plotting
#'
#' @param left left coordinate of block
#' @param right right coordinate of block
#' @param y y position for middle of block
#' @param arrowLength arrow length in usr coords
#' @param shaft half of shaft thickness in usr coords
#' @param point half of arrow thickness in usr coords
#' @param concat logical indicating whether to concatenate multiple arrows into one data frame separated by NAs (e.g. ready for poly)
#' @export
#' @return a data.frame with columns x and y specifying coordinates for arrows
arrow<-function(left,right,y,arrowLength=diff(par('usr')[1:2])*.05,shaft=.2,point=.4,concat=TRUE){
	if(any(left>right))stop(simpleError('Left border > right border of arrow'))
	arrowX<-right-arrowLength
	arrowX<-ifelse(arrowX<left,left+(right-left)/10,arrowX)
	coords<-mapply(function(left,right,y,arrowX){data.frame('x'=c(left,arrowX,arrowX,right,arrowX,arrowX,left),'y'=y+c(shaft,shaft,point,0,-point,-shaft,-shaft))},left,right,y,arrowX,SIMPLIFY=FALSE)
	if(concat)coords<-do.call(rbind,lapply(coords,function(x)return(rbind(x,c(NA,NA)))))
	return(coords)
}

#' Stack regions into smallest number of lines (using greedy algorithm)
#'
#' @param starts vector of starts of regions (note can add buffer by substracting arbitrary spacer)
#' @param ends vector of ends of regions (note can add buffer by adding arbitrary spacer)
#' @export
#' @return vector of row ids
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


#' Calculate the Wilson score interval http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
#'
#' @param nTrue number of trues observed
#' @param nFalse number of falses observed
#' @param alpha error percentile
#' @export
#' @return two element vector giving the lower and upper bounds of the interval 
wilsonInt<-function(nTrue,nFalse,alpha=.05){
	n<-nTrue+nFalse
	prop<-nTrue/n
	z<-pnorm(1-alpha/2)
	plusMinus<-z*sqrt(prop*(1-prop)/n+z^2/4/n^2)
	return((prop+1/2/n*z^2+c(-plusMinus,plusMinus))/(1+1/n*z^2))
}



