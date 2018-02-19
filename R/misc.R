#' Convenience function to check for errors in a list of results e.g. from mclapply
#'
#' @param x vector or list to check for errors
#' @export
#' @return logical vector of length(x) specifying if each element of x was an error
#' @seealso \code{\link[parallel]{mclapply}}
#' @examples
#' isError(c(list(1:3),list(simpleError('Error!'))))
isError<-function(x){
    unlist(lapply(x,function(y)inherits(y,'simpleError')|inherits(y,'try-error')))
}

#' Convenience function for selecting multiple elements from a matrix by rows and columns position
#' 
#' @param rows row coordinate of matrix
#' @param cols col coorinate of matrix
#' @param mat Matrix of interest
#' @param returnIndex If TRUE return the one dimensional index for the items otherwise return the selected elements of the matrix
#' @export
#' @return A vector of indices if returnIndex is TRUE or a vector of the selected matrix elements
#' @examples
#' indexMatrix(c(1,2,3),c(1,2,1),matrix(1:9,nrow=3,byrow=TRUE))
indexMatrix<-function(rows,cols,mat,returnIndex=FALSE){
  mat<-as.matrix(mat)
  if(is.character(rows)){tmp<-1:nrow(mat);names(tmp)<-rownames(mat);rows<-tmp[rows]}
  if(is.character(cols)){tmp<-1:ncol(mat);names(tmp)<-colnames(mat);cols<-tmp[cols]}
  if(any(is.na(rows))||any(is.na(cols)))stop(simpleError('NAs in indices'))
  if(length(rows)!=length(cols))stop(simpleError("rows and cols different lengths"))
  if(min(rows)<1|max(rows)>nrow(mat))stop(simpleError("rows outside matrix bounds"))
  if(min(cols)<1|max(cols)>ncol(mat))stop(simpleError("cols outside matrix bounds"))
  index<-(cols-1)*nrow(mat)+rows
  if(returnIndex)return(as.vector(index)) #as.vector to remove names
  else return(as.vector(mat[index]))
}

#' Convenience function for selecting multiple elements from a matrix by x,y position
#' convenience function to resize R console window
#' @param width width to adjust window to
#'
#' @export
adjustWindow<-function(width=as.integer(Sys.getenv('COLUMNS')))if(!is.na(width))options(width=width)

#' Convenience function to list objects by size
#'
#' @param env the environment to list
#' @export
#' @return object sizes in decreasing order
#' object.sizes()
object.sizes<-function(env=.GlobalEnv)sort(sapply(ls(envir=env),function(x)utils::object.size(get(x,env))),decreasing=TRUE)

#' Generate an error
#'
#' Convenience function for stop(simpleError())
#' @param ... Strings to be concatenated into error message
#' @export
#' @return Generates stopping error before returning
#' tryCatch(stopError('Error: we had ',sum(1:10),' problems'),error=function(e)print(e))
stopError<-function(...){
  stop(simpleError(paste(...,sep='')))
}

#' Get the conservative edge of a confidence interval
#'
#' @param boundaries upper and lower values
#' @param base value, e.g. 0 or 1
#' @export
#' @return single numeric. if the boundaries overlap the base then return base otherwise return value closest to base
#' @examples
#' conservativeBoundary(c(1,100),10)
#' conservativeBoundary(c(20,100),10)
#' conservativeBoundary(c(-100,-20),0)
#' mat<-matrix(1:4,nrow=2)
#' fish<-fisher.test(mat)
#' conservativeBoundary(fish$conf.int,1)
#' mat2<-matrix(c(10,1,1,10),nrow=2)
#' fish2<-fisher.test(mat2)
#' conservativeBoundary(fish2$conf.int,1)
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
#' @examples
#' lagNA(1:10,3)
#' lagNA(1:10,-3)
#' lagNA(1:10,-3,-99)
lagNA<-function(x,lag=1,fill=NA){
  out<-x  
  if(abs(lag)>=length(x))return(rep(fill,length(x)))
  if(lag>0)out<-c(out[-(1:lag)],rep(fill,lag))
  if(lag<0)out<-c(rep(fill,abs(lag)),out[-(length(x)+(-1:lag)+1)])
  return(out)
}

#' Calculating a moving average/max/min statistic over a vector
#'
#' @param vec to average over
#' @param statFunc function to use
#' @param spacer how many neighbors to the left and right to average
#' @export
#' @return a vector the same length as vec
#' @examples
#' movingStat(1:10,mean)
#' movingStat(1:10,max)
#' movingStat(1:20,min)
#' movingStat(rnorm(100),mean,5)
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
#' @param recursive A logical indicating whether to recursively look inside arguments for NAs
#' @export
#' @return logical indicating whether any argument contains NA
anyArgsNA<-function(...,recursive=FALSE){
  args<-list(...)
  hasNa<-anyNA(args,recursive=recursive)
  return(hasNa)
}

#' Cache an operation to disk or load if cache file present
#'
#' Cache an operation to disk to save computation time in the future. A hash of the function and arguments is generated to avoid erroneous caching but any global variables or arguments excluded by EXCLUDE are not checked and are vulnerable erroneous output if they are changed.
#' 
#' @param cacheFile File location to save data to 
#' @param operation A function taking ... arguments and returning object to be stored
#' @param ... Arguments for operation function
#' @param OVERWRITE If FALSE throw an error if hash of ... changes from cached values. If TRUE redo operation and overwrite cache without asking
#' @param VOCAL If TRUE report on status of caching
#' @param EXCLUDE Vector of names of arguments to exclude from md5 digest comparison (for very large arguments)
#' @export
#' @return Output from operation function with ... arguments
#' @examples
#' cache<-tempfile()
#' cacheOperation(cache,mean,1:10)
#' cacheOperation(cache,mean,1:10)
#' cacheOperation(cache,mean,x=1:20,OVERWRITE=TRUE)
#' cacheOperation(cache,mean,x=1:20,OVERWRITE=TRUE,EXCLUDE='x')
#' #Note that EXCLUDEd arguments are not checked and can generate false output
#' y<-2
#' cacheOperation(cache,function(x)sum(x^y),1:10,OVERWRITE=TRUE)
#' #Note that global variables can generate false output
#' y<-3
#' cacheOperation(cache,function(x)sum(x^y),1:10,OVERWRITE=TRUE)
#' cacheOperation(cache,mean,x=1:100,OVERWRITE=TRUE,EXCLUDE='x')
cacheOperation<-function(cacheFile,operation,...,OVERWRITE=FALSE,VOCAL=TRUE,EXCLUDE=NULL){
  parent<-parent.frame()
  #avoid evaluation EXCLUDEd args until necessary since they're probably big
  unevalArgs<-match.call(expand.dots=FALSE)$'...'
  varSelector<-if(is.null(names(unevalArgs)))rep(TRUE,length(unevalArgs)) else !names(unevalArgs) %in% EXCLUDE
  allArgs<-lapply(unevalArgs[varSelector],eval,envir=parent)
  #try to prevent function scope from changing things
  md5<-digest::digest(lapply(c(operation,allArgs),function(x)if(is.function(x))deparse(x)else x))
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
  out<-do.call(operation,allArgs,envir=parent)
  save(md5,out,file=cacheFile)
  return(out)
}

#' Check predictions of glm by cross validation
#'
#' @param model a glm to cross validate
#' @param data data for the model
#' @param K number of pieces to split data into
#' @param nCores number of cores to use
#' @param subsets predefined subset of same length as data specifying groupings
#' @param vocal echo progress indicator
#' @export
#' @return a dataframe giving predictions for each data point with columns prediction and subsetId
#' @examples
#' z<-data.frame('x'=1:11,y=c(rnorm(10,1:10),100))
#' cv.glm.par(glm(y~x,data=z))
cv.glm.par<-function(model,data=eval(modelCall$data),K=nrow(data),nCores=1,subsets=NULL,vocal=FALSE){
  modelCall<-model$call
  n<-nrow(data)
  if(is.null(subsets))subsets<-split(1:n,sample(rep(1:K,length.out=n)))
  preds<-parallel::mclapply(subsets,function(outGroup){
    if(vocal)cat('.')
    subsetData<-data[-outGroup,,drop=FALSE]
    predData<-data[outGroup,,drop=FALSE]
    thisModel<-modelCall
    thisModel$data<-subsetData
    return(stats::predict(eval(thisModel),predData))
  },mc.cores=nCores)
  pred<-unlist(preds)[order(unlist(subsets))]
  subsetId<-rep(1:K,sapply(subsets,length))[order(unlist(subsets))]
  out<-data.frame(pred,subsetId)
  rownames(out)<-NULL
  return(out)
}

#' A function to break an mclapply into parts, save to disk and run independently
#' 
#' @param x vector/list to apply over
#' @param applyFunc function to apply
#' @param ... arguments for mclapply
#' @param mc.cores number of cores to use
#' @param extraCode character vector of setup code (each command self-contained within a cell or concatenated with ;)
#' @param nSplits number of splits to make (if > mc.cores then R will restart more frequently and load less data at once)
#' @param VOCAL if TRUE then output status information
#' @param envir environment to look for variables in
#' @export
#' @return the concatenated outputs from applyFunc 
#' cleanMclapply(1:10,2,sqrt)
#' cleanMclapply(1:10,2,function(x,y)x^y,y=3)
cleanMclapply<-function(x,applyFunc,...,mc.cores=parallel::detectCores(),extraCode='',nSplits=mc.cores,VOCAL=TRUE,envir=.GlobalEnv){
  #otherwise global variables can get pulled along with function environment
  environment(applyFunc)<-envir
  if(nSplits<mc.cores)nSplits<-mc.cores
  splits<-unique(round(seq(0,length(x),length.out=nSplits+1)))
  if(length(splits)<nSplits+1)nSplits<-length(splits)-1 #not enough items to fill so set lower
  dotVars<-match.call(expand.dots=FALSE)$'...'
  extraArgs<-lapply(dotVars,eval,envir=envir)
  if(is.null(names(extraArgs)))names(extraArgs)<-rep('',length(extraArgs))
  names(extraArgs)[names(extraArgs)=='']<-sprintf('POSITIONAL_%d__',1:sum(names(extraArgs)==''))
  files<-c()
  outFiles<-c()
  scriptFiles<-c()
  logFiles<-c()
  allInRdat<-tempfile()
  if(VOCAL)message('Writing extra arguments')
  EXTRAARGS__<-extraArgs
  save(EXTRAARGS__,file=allInRdat)
  for(ii in 1:nSplits){
    if(VOCAL)message("Writing core ",ii," data")
    thisInRdat<-tempfile()
    thisRScript<-tempfile()
    thisOutRdat<-tempfile()
    thisLog<-tempfile()
    THISDATA__<-x[(splits[ii]+1):splits[ii+1]]
    SAVEDATA__<-c('THISDATA__'=list(THISDATA__),'applyFunc'=applyFunc)
    save(SAVEDATA__,file=thisInRdat)
    script<-sprintf(
      'load("%s")\nload("%s")\n%s\nout<-with(EXTRAARGS__,with(SAVEDATA__,lapply(THISDATA__,applyFunc%s%s)))\nsave(out,file="%s")',
      allInRdat,
      thisInRdat,
      paste(extraCode,collapse=';'),
      ifelse(length(extraArgs)>0,',',''),
      paste(ifelse(grepl('POSITIONAL_[0-9]+__',names(extraArgs)),
        names(extraArgs),
        paste(names(extraArgs),names(extraArgs),sep='=')
      ),collapse=','),
      thisOutRdat
    )
    writeLines(script,thisRScript)
    outFiles<-c(outFiles,thisOutRdat)    
    scriptFiles<-c(scriptFiles,thisRScript)    
    logFiles<-c(logFiles,thisLog)    
  }
  if(VOCAL)message("Running")
  if(VOCAL)message("Logs: ",paste(logFiles,collapse=', '))
  exitCode<-parallel::mclapply(mapply(c,scriptFiles,logFiles,SIMPLIFY=FALSE),function(x){out<-system(sprintf("R CMD BATCH --no-save --no-restore %s %s",x[1],x[2]));if(VOCAL)cat('.');return(out)},mc.cores=mc.cores)
  if(VOCAL)cat('\n')
  if(any(exitCode!=0)){
    if(VOCAL)message(paste(utils::tail(readLines(logFiles[exitCode!=0][1]),30),collapse='\n'))
    stop(simpleError('Problem running multi R code'))
  }
  if(VOCAL)message("Loading split outputs")
  out<-do.call(c,lapply(outFiles,function(outFile){if(VOCAL)cat('.');load(outFile);return(out)}))
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
#' @examples
#' orderIn(c('z','c','a','t'),letters)
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
#' @examples
#' fillList(list(c('a'=1,'b'=2),c('c'=3)))
#' fillList(list(c('a'=1,'b'=2),c('c'=3)),namesList='a')
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
#' @examples
#' plot(1:10)
#' points(rep(5,5),convertLineToUser(0:4,1),xpd=NA)
#' points(convertLineToUser(0:4,2),rep(5,5),xpd=NA)
convertLineToUser<-function(line,axis=1){
  if(!(axis %in% 1:4))stop(simpleError('Undefined axis'))
  axisPair<-sort((c(axis-1,axis+1)%%4)+1)
  isHeight<-(axis%%2)==1
  isSecond<-axis>2
  thisMar<-graphics::par('mar')[axis]
  marWidth<-thisMar/sum(graphics::par('mar')[axisPair])*(graphics::par('fin')-graphics::par('pin'))[isHeight+1]
  widthPerLine<-marWidth/thisMar
  #find base line + add in if plot doesn't cover whole device e.g. graphics::par(mfrow=c(2,1))
  base<-ifelse(isSecond,graphics::par('fin')[isHeight+1]-widthPerLine*thisMar,widthPerLine*thisMar) + graphics::par('fig')[1+isHeight*2]*graphics::par('din')[isHeight+1]
  func<-if(isHeight)graphics::grconvertY else graphics::grconvertX
  out<-func(base+line*widthPerLine*ifelse(isSecond,1,-1),'inches','user')
  return(out)
}

#' Convert from user coordinates to axis line
#'
#' @param usr user coordinate to convert to line
#' @param axis axis to do conversion on (1:4 same as axis, mtext command)
#' @export
#' @return axis line position
#' @examples
#' plot(1:10)
#' mtext('Test1',1,convertUserToLine(0,1))
#' points(5.5,0,xpd=NA)
#' mtext('Test2',2,convertUserToLine(0,2))
#' points(0,5.5,xpd=NA)
convertUserToLine<-function(usr,axis=1){
  if(!(axis %in% 1:4))stop(simpleError('Undefined axis'))
  axisPair<-sort((c(axis-1,axis+1)%%4)+1)
  isHeight<-(axis%%2)==1
  isSecond<-axis>2
  thisMar<-graphics::par('mar')[axis]
  marWidth<-thisMar/sum(graphics::par('mar')[axisPair])*(graphics::par('fin')-graphics::par('pin'))[isHeight+1]
  widthPerLine<-marWidth/thisMar
  #find base line + add in if plot doesn't cover whole device e.g. graphics::par(mfrow=c(2,1))
  plotWidth<-graphics::par('fin')[isHeight+1]
  func<-if(isHeight)graphics::grconvertY else graphics::grconvertX
  usrInches<-func(usr,'user','inches')
  base<-ifelse(isSecond,graphics::par('fin')[isHeight+1]-widthPerLine*thisMar,widthPerLine*thisMar) + graphics::par('fig')[1+isHeight*2]*graphics::par('din')[isHeight+1]
  out<-(usrInches-base)/widthPerLine*ifelse(isSecond,1,-1)
  return(out)
}

#' Label axis with slanted labels
#' 
#' @param side an integer specifying which side to draw axis on. 1=bottom, 2=left, 3=top and 4=right.
#' @param at locations for ticks and labels
#' @param labels character vector of tick labels
#' @param srt a numeric specifying string rotation in degrees
#' @param location an integer specifying the line to start the text
#' @param adj one or two values in [0, 1] which specify the x (and optionally y) adjustment of the labels.  On most devices values outside that interval will also work.
#' @param axisArgs a list of additional arguments for axis
#' @param ... additional arguments to text
#' @return NULL
#' @export
#' @examples
#' par(mar=c(8,8,8,8))
#' plot(1:10,xlab='',ylab='',xaxt='n',yaxt='n')
#' labels<-c('A label','Another label','A longer longer label','A really really\nlong label')
#' slantAxis(1,seq(2,8,2),labels)
#' slantAxis(2,seq(2,8,2),labels)
#' slantAxis(3,seq(2,8,2),labels)
#' slantAxis(4,seq(2,8,2),labels,srt=-30,cex=.8,axisArgs=list(col.ticks='red'),lwd=2)
slantAxis<-function(side,at,labels=at,srt=ifelse(side %in% c(1,4),-45,45),location=1.2,adj=ifelse(side==2,1,0),axisArgs=list(),...){
  
  do.call(graphics::axis,c(list(side,at,label=FALSE),axisArgs))
  if(side %in% c(1,3)){
    graphics::text(at, convertLineToUser(location,side), srt = srt, adj = adj, labels = labels, xpd = TRUE,...)
  }else{
    graphics::text(convertLineToUser(location,side),at, srt = srt, adj = adj, labels = labels, xpd = TRUE,...)
  }
  return(invisible(NULL))
}

#' Find coords for plotting horizontal arrows
#'
#' @param tail x coordinate for tail of the arrow
#' @param head x coordinate for the head of the arrow
#' @param y y position for middle of block
#' @param arrowLength arrow length in usr coords
#' @param shaft half of shaft thickness in usr coords
#' @param point half of arrow thickness in usr coords
#' @param concat logical indicating whether to concatenate multiple arrows into one data frame separated by NAs (e.g. ready for poly)
#' @export
#' @return a data.frame with columns x and y specifying coordinates for arrows
#' @examples
#' plot(1:10)
#' polygon(arrow(1:5,3:7*1.5,2:6))
#' polygon(arrow(3:6,1:4/1.5,7:10))
#' polygon(arrow(7:10,7:10+0:3*.1,7:10))
#' polygon(arrow(6:9,6:9-0:3*.1,7:10))
arrow<-function(tail,head,y,arrowLength=diff(graphics::par('usr')[1:2])*.05,shaft=.2,point=.4,concat=TRUE){
  isLeft<-tail<head
  arrowX<-head+arrowLength*ifelse(isLeft,-1,1)
  arrowX<-ifelse((arrowX<tail&isLeft)|(arrowX>tail&!isLeft),tail+(head-tail)/10,arrowX)
  coords<-mapply(function(tail,head,y,arrowX){data.frame('x'=c(tail,arrowX,arrowX,head,arrowX,arrowX,tail),'y'=y+c(shaft,shaft,point,0,-point,-shaft,-shaft))},tail,head,y,arrowX,SIMPLIFY=FALSE)
  if(concat)coords<-do.call(rbind,lapply(coords,function(x)return(rbind(x,c(NA,NA)))))
  return(coords)
}

#' Stack regions into smallest number of lines (using greedy algorithm)
#'
#' @param starts vector of starts of regions (note can add buffer by substracting arbitrary spacer)
#' @param ends vector of ends of regions (note can add buffer by adding arbitrary spacer)
#' @export
#' @return vector of row ids
#' @examples
#' stackRegions(1:10,1:10+2)
#' stackRegions(1:10,1:10+5)
stackRegions<-function(starts,ends){
  if(any(ends<starts))stop(simpleError('Ends less than starts'))
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


#' Calculate the Wilson score interval 
#'
#' @param nTrue single number of trues observed
#' @param nFalse single number of falses observed
#' @param alpha error percentile
#' @references \url{http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval}
#' @export
#' @return two element vector giving the lower and upper bounds of the interval 
#' @examples
#' wilsonInt(10,100)
#' wilsonInt(0,100)
#' wilsonInt(100,1)
wilsonInt<-function(nTrue,nFalse,alpha=.05){
  nTrue<-nTrue[1]
  nFalse<-nFalse[1]
  if(nTrue<0|nFalse<0)stop(simpleError('Counts less than zero'))
  n<-nTrue+nFalse
  prop<-nTrue/n
  z<-stats::qnorm(1-alpha/2)
  plusMinus<-z*sqrt(prop*(1-prop)/n+z^2/4/n^2)
  return((prop+1/2/n*z^2+c(-plusMinus,plusMinus))/(1+1/n*z^2))
}


#' Convenience function for picking first most abundant of a set of values
#'
#' @param values A vector of items
#' @export
#' @return First most abundant item as a string
#' @examples
#' mostAbundant(c(1:10,1))
#' mostAbundant(c('d',rep(letters,10)))
mostAbundant<-function(values){
  tmp<-table(values)
  return(names(tmp)[which.max(tmp)])
}

#' Convenience function for converting 1-d table to named vector
#'
#' @param tab A one-dimensional table e.g. the output from \code{table}
#' @export
#' @return a named vector with elements corresponding to the counts in table
#' @examples
#' table2vector(1:10)
table2vector<-function(tab){
  if(length(dim(tab))>1)stop(simpleError('table2vector only works with 1 dimensional tables'))
  return(structure(as.vector(tab),.Names=names(tab)))
}

#' Convenience function to escape characters for insertion into regex brackets
#'
#' Note to escape a '\code{]}' using this you'll have to use \code{perl=TRUE} in the downstream function
#'
#' @param regexChars a vector of single characters to escape if in escapeChars
#' @param escapeChars a vector of which characters to escape
#' @export
#' @return a vector of escaped characters
#' @examples
#' escapeRegexBracketChars(c(']','\\','-','A','B'))
escapeRegexBracketChars<-function(regexChars,escapeChars=c('^','-','[',']','\\')){
  return(sprintf('%s%s',ifelse(regexChars %in% escapeChars,'\\',''),regexChars))
}

#' Convenience function to fill values down a column
#'
#' Takes a vector of values and fills any empty strings or NAs with the last non-empty value.
#'
#' @param x a vector to be filled
#' @param emptyStrings a vector of strings to be considered empty. Defaults to empty string and NA.
#' @param errorIfFirstEmpty if TRUE then throw an error if the first entry is empty
#' @export
#' @return a vector of the same length as x with empty entries filled
#' @examples
#' fillDown(c(1:5,NA,NA,6,NA,7))
#' fillDown(c('a','c','d','Z',NA),'Z')
fillDown<-function(x,emptyStrings=c(NA,''),errorIfFirstEmpty=TRUE){
  #depending on %in% to catch NAs if necessary
  isEmpty<-x %in% emptyStrings
  if(isEmpty[1]&errorIfFirstEmpty)stop(simpleError('First value empty'))
  #if first is empty and we don't want errors then have to just fill down from it anyway
  isEmpty[1]<-FALSE
  ids<-1:length(x)
  ids[isEmpty]<-0
  ids<-cummax(ids)
  return(x[ids])
}


#' Alternative with() function with explicit naming
#'
#' An alternative to the \code{with} function in base R where \code{data} is represented by a temporary variable within a new environment rather than evaluating \code{expr} directly within \code{data}. This can help explicitly show which variables are contained within \code{data} from those outside.
#'
#' @param ... data to be used in expression. Argument names are used as temporary variable names. So \code{xx=data.frame('a'=1:10)} will result in a data.frame called \code{xx} with column \code{a} in the evaluating environment.
#' @param expr expression to evaluate. if expr is not defined by name then the last unnamed argument is used.
#' @export
#' @return the value of the evaluated \code{expr}
#' @examples
#' d<-3:12
#' longNameDataFrame<-data.frame('a'=1:10,'b'=2:11)
#' with(longNameDataFrame,a+b+d)
#' withAs(xx=longNameDataFrame,xx$a+xx$b+d)
#' anotherDf<-data.frame('c'=2:12)
#' withAs(xx=longNameDataFrame,zz=anotherDf,xx$a+xx$b+zz$c)
withAs<-function(...,expr=NULL){
  parent<-parent.frame()
  env<-new.env(parent=parent)
  dotVars<-match.call(expand.dots=FALSE)$'...'
  if(missing(expr)){
    expr<-dotVars[[length(dotVars)]]
    dotVars<-dotVars[-length(dotVars)]
  }
  if(is.null(names(dotVars))||any(names(dotVars)==''))stop(simpleError('Unassigned variables passed to withAs'))
  #make sure to get parent above and not within the mapply
  mapply(function(as,val)assign(as,eval(val,envir=parent),env),names(dotVars),dotVars)
  return(eval(substitute(expr), env))
}


#' Add a pretty log axis to a plot
#'
#' Adds a log axis to the current plot with \code{log='x'} or \code{log='y'} with ticks and labels specified with powers of 10
#'
#' @param side an integer specifying which side of the plot the axis is to be drawn on.  The axis is placed as follows: 1=below, 2=left, 3=above and 4=right. 
#' @param exponent If true label sides using exponents of 10. Otherwise just use normal numbers
#' @param addExtra If true add additional labels at 5 or 2 if there are few labels 
#' @param minorTcl Length of the minor unlabelled ticks. NA to suppress
#' @param axisMin The minimum to extend ticks to
#' @param offset An integer giving the amount to offset labels from their true position. For example, offset=1 if the plotted data was generated by adding to 1 to avoid 0s in the original data.
#' @param col.ticks colors for the tick marks
#' @param ... Additional arguments to \code{axis}
#' @export
#' @return An invisible list containing the minor and major axis label positions
#' @examples
#' plot(1:1000,log='xy',xaxt='n',yaxt='n')
#' logAxis(las=1)
#' logAxis(1,exponent=FALSE)
logAxis<-function(side=2,exponent=TRUE,addExtra=!exponent,minorTcl=-.2,axisMin=-Inf,offset=0,col.ticks='black',...){
  if(side %in% c(2,4)) parX<-sort(graphics::par('usr')[3:4])
  else parX<-sort(graphics::par('usr')[1:2])
  minX<-max(10^parX[1],axisMin)
  maxX<-10^parX[2]
  if(log10(maxX)-log10(minX)>400)stop(simpleError('Huge range in logged axis'))
  allTicks<-unlist(lapply(floor(log10(minX)):ceiling(log10(maxX)),function(x)1:9*10^x))
  allTicks<-allTicks[allTicks<=maxX & allTicks>=minX]
  graphics::axis(side,allTicks+offset,rep('',length(allTicks)),tcl=minorTcl,col.ticks=col.ticks)
  if(ceiling(log10(minX))<=floor(log10(maxX)))prettyY<-seq(ceiling(log10(minX)),floor(log10(maxX)),1)
  else prettyY<-c()
  graphics::axis(side,10^prettyY+offset,rep('',length(prettyY)),tcl=minorTcl*2,col.ticks=col.ticks)
  if(length(prettyY)>7)prettyY<-pretty(prettyY)
  if(length(prettyY)==0)prettyY<-c(ceiling(log10(minX)),floor(log10(maxX)))
  if(addExtra){
    origPretty<-prettyY
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(5),origPretty-log10(10/5)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(2),origPretty-log10(10/2)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(3),origPretty-log10(10/3)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(7),origPretty-log10(10/7)))
  }
  if(exponent){
    if(any(prettyY%%1!=0))labs<-sapply(prettyY,function(x)as.expression(bquote(.(10^(x%%1))%*%10^.(floor(x)))))
    else labs<-ifelse(prettyY==0,1,sapply(prettyY,function(x)as.expression(bquote(10^.(floor(x))))))
  }
  else labs<-10^prettyY
  graphics::axis(side,10^prettyY+offset,labs,col.ticks=col.ticks,...)
  return(invisible(list('minor'=allTicks,'major'=10^prettyY)))
}

#' Add an inset color scale to a plot
#'
#' Adds an inset color scale to a plot without disrupting future plotting. Especially useful with \code{heatmap} or \code{image} plots.
#'
#' @param breaks a set of finite numeric breakpoints for the colours: must have one more breakpoint than colour and be in increasing order.
#' @param col a list of colors such as that generated by 'rainbow', 'heat.colors', 'topo.colors', 'terrain.colors' or similar functions.
#' @param insetPos a four element numeric vector giving the position of the bottom, left, top and right edges of the scale in nfc coordinates (0 = bottom/left of figure, 1 = top/right of figure)
#' @param main a character string giving label for color scale
#' @param at positions for labels on the axis
#' @param labels labels for the axis
#' @param offset add a small offset to rectangles to suppress pdf viewers showing a gap at exact border between rectangles
#' @export
#' @seealso \code{\link[stats]{heatmap}}, \code{\link[graphics]{image}}
#' @return NULL
#' @examples
#' dists<-as.matrix(dist(sort(runif(40))))
#' breaks<-seq(0,1,.01)
#' cols<-rev(heat.colors(length(breaks)-1))
#' heatmap(dists,col=cols)
#' insetScale(breaks,cols,main='Distance')
insetScale<-function(breaks,col,insetPos=c(.025,.015,.04,.25),main='',offset=1e-3,at=NULL,labels=NULL){
  if(length(breaks)!=length(col)+1)stop('Number of breaks must be one more than colors')
  insetPos<-c(graphics::grconvertY(insetPos[1],'nfc','user'),graphics::grconvertX(insetPos[2],'nfc','user'),graphics::grconvertY(insetPos[3],'nfc','user'),graphics::grconvertX(insetPos[4],'nfc','user'))
  breakPos<-((breaks)-(min(breaks)))/max((breaks)-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  #add a bit of offset to avoid pdf viewers displaying breaks between exact rectangle border meeting
  offsetPos<-breakPos[-1]+c(rep(offset*diff(range(breakPos)),length(breakPos)-2),0)
  graphics::rect(breakPos[-length(breakPos)],insetPos[1],offsetPos,insetPos[3],col=col,xpd=NA,border=NA)
  graphics::rect(insetPos[2],insetPos[1],insetPos[4],insetPos[3],xpd=NA)
  if(is.null(at)){
    at<-pretty(breaks)
    at<-at[at<=max(breaks)&at>=min(breaks)]
  }
  if(is.null(labels))labels<-at
  convertPos<-(at-(min(breaks)))/((max(breaks))-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  graphics::segments(convertPos,insetPos[1],convertPos,insetPos[1]-diff(insetPos[c(1,3)])*.1,xpd=NA)
  graphics::text(convertPos,insetPos[1]-diff(insetPos[c(1,3)])*.175,labels,xpd=NA,adj=c(.5,1),cex=.85)
  graphics::text(mean(insetPos[c(2,4)]),insetPos[3]+diff(insetPos[c(1,3)])*.45,main,xpd=NA,adj=c(.5,0))
  invisible(NULL)
}

