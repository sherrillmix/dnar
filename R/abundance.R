#' Calculate Shannon diversity of a vector
#'
#' @param x Vector of counts or proportions
#' @param base Base of logarithm
#' @param standardize If TRUE divide x by sum x
#' @export
#' @return Shannon diversity
#' @examples
#' shannon(1:10)
#' shannon(1:10,base=2)
shannon<-function(x,base=exp(1),standardize=TRUE){
   if(any(x<0))stop(simpleError('Please give positive values of x'))
   x<-x[x>0]
   if(standardize)prop<-x/sum(x)
   else prop<-x
   return(-sum(prop*log(prop,base)))
}


#' Calculate Rao diversity of a vector
#'
#' @param x Vector of counts or proportions of N elements
#' @param dist Distance matrix of NxN elements
#' @export
#' @return Rao diversity
#' rao(1:3,as.matrix(dist(1:3)))
rao<-function(x,dist){
  if(any(dim(dist)!=length(x)))stop(simpleError('Distance and proportion vectors do not match'))
  zeros<-x==0
  x<-x[!zeros]
  dist<-dist[!zeros,!zeros]
  props<-x/sum(x)
  propTable<-outer(props,props)
  return(sum(propTable*dist))
}

#' Calculate Jensen Shannon divergence between two probability distributions
#'
#' @param x Vector of counts or proportions of N elements
#' @param y Second vector of counts or proportions of N elements
#' @param base Base of logarithm
#' @export
#' @return Jensen-Shannon divergence
#' @examples
#' jensenShannon(1:4,4:1)
jensenShannon<-function(x,y,base=2){
  if(length(x)!=length(y))stop(simpleError('Length of x and y differ'))
  propX<-x/ifelse(sum(x)==0,1,sum(x))
  propY<-y/ifelse(sum(y)==0,1,sum(y))
  out<-shannon((propX+propY)/2,base,FALSE)-shannon(propX,base,FALSE)/2-shannon(propY,base,FALSE)/2
  #if(method=='shannon')out<-shannon((propX+propY)/2,base,FALSE)-shannon(propX,base,FALSE)/2-shannon(propY,base,FALSE)/2
  #if(method=='kullback')out<-.5*kullback(propX,(propY+propX)/2,base,FALSE)+.5*kullback(propY,(propY+propX)/2,base,FALSE)
  #if(round(shannonCalc,5)!=round(kullbackCalc,5))stop(simpleError('Problem calculating shannon jensen'))
  return(out)
}

#' Calculate Kullback-Leibler divergence between two probability distributions
#'
#' Undefined if zero counts/proportions
#'
#' @param x Vector of counts or proportions of N elements
#' @param y Second vector of counts or proportions of N elements
#' @param base Base of logarithm
#' @param standardize If TRUE divide x and y by sum(x) and sum(y)
#' @export
#' @return Kullback-Leibler divergence
#' @examples
#' kullback(1:3,3:1)
kullback<-function(x,y,base=2,standardize=TRUE){
  selector<-y<=0&x<=0
  x<-x[!selector];y<-y[!selector]
  if(standardize){
    propX<-x/sum(x)
    propY<-y/sum(y)
  }else{
    propX<-x
    propY<-y
  }
  sum(propX*log(propX/propY,base))
}

#' Calculate chao diversity index
#'
#' @param counts A vector of counts with one entry per "species"
#' @export
#' @return Chao index
#' @examples
#' chao(c(1,1,1,2,3,5,100))
chao<-function(counts){
  counts<-counts[counts>0]
  return(length(counts)+sum(counts==1)*(sum(counts==1)-1)/2/(sum(counts==2)+1))
}

#' Calculate rarefaction using boostrapping
#'
#' @param counts Counts of species
#' @param samples Vector of numbers of draws to calculate rarefaction at
#' @param reps How many random samples to take at each step
#' @param chaoAdjust If TRUE calculate chao-predicted species number on each random draw
#' @param replace If TRUE sample with replacements. If FALSE sample without replacement.
#' @param minObs Minimum number of counts to be counted as present e.g. minobs=2 discards singletons (ignored if chaoAdjust=TRUE)
#' @param statFunc Function to apply to counts e.g. mean or median
#' @param ... Additional arguments to statFunc
#' @export
#' @return Dataframe of calculated quantiles with rownames of the number of samples drawn
#' @seealso \code{\link{rareEquation}}, \code{\link{chao}}
#' @examples
#' rarefy(1:20,reps=100)
#' rarefy(1:20,reps=100,chaoAdjust=TRUE)
rarefy<-function(counts,samples=unique(round(sum(counts)*seq(.1,1,.1))),reps=1000,chaoAdjust=FALSE,replace=FALSE,minObs=1,statFunc=stats::quantile,...){
  species<-rep(1:length(counts),counts)
  output<-lapply(samples,function(sample,reps,species){
    numSpecies<-sapply(1:reps,function(rep,species,sample,chaoAdjust){
      thisSpecies<-sample(species,sample,replace=replace)
      if(chaoAdjust){
        return(chao(table(thisSpecies)))
      } else{
        speciesTable<-table(thisSpecies)
        return(length(speciesTable[speciesTable>=minObs]))
      }
    },species,sample,chaoAdjust)
    estimate<-statFunc(numSpecies,...)
    return(estimate)
  },reps,species)
  output<-as.data.frame(do.call(rbind,output))
  rownames(output)<-samples
  return(output)
}

#' Calculate rarefaction using formula
#'
#' @param counts Counts of species
#' @param samples Vector of numbers of draws to calculate rarefaction at
#' @param minObs Minimum number of counts to be counted as present e.g. minObs=2 discards singletons
#' @export
#' @return Vector with predicted rarefied counts for entries and the sampled number of reads for names
#' @seealso \code{\link{rarefy}}
#' @examples
#' rareEquation(1:20)
rareEquation<-function(counts,samples=unique(round(sum(counts)*seq(.1,1,.1))),minObs=1){
  counts<-counts[counts>0]
  if(length(samples)>1){
    out<-sapply(samples,function(xx)rareEquation(counts,xx,minObs=minObs))
  }else{
    expectMisses<-do.call(cbind,lapply(0:(minObs-1),function(ii)exp(lchoose(sum(counts)-counts,samples-ii)+lchoose(counts,ii)-lchoose(sum(counts),samples))))
    expectMiss<-apply(expectMisses,1,sum)
    out<-sum(1-expectMiss)
  }
  names(out)<-samples
  return(out)
}

#' Calculate the probability of observing observedX species from nGroups groups with groupsize members with n observations
#'
#' @param nGroups number of equally sized groups
#' @param n number of observations pulled from the individuals
#' @param groupSize size of the equally sized groups
#' @param observedX number of species observed
#' @return probability of observing observedX species from nGroups groups with groupsize members with n observations
#' pRare(4,10,10,10)
pRare<-function(observedX,n,nGroups,groupSize){
  choose(nGroups,observedX)*chooseAtLeastOneFromEach(n,observedX,groupSize)/choose(nGroups*groupSize,n)
}

#' Calculate the number of combinations possible from nGroups groups with groupsize members and n observations
#'
#' @param n number of observations pulled from the individuals
#' @param nGroups number of equally sized groups
#' @param groupSize size of the equally sized groups
#' @export
#' @seealso \code{\link{pRare}}
#' @return number of combinations containing at least one individual from each group
#' @examples
#' chooseAtLeastOneFromEach(10,10,2)
chooseAtLeastOneFromEach<-function(n,nGroups,groupSize){
  if(nGroups<=0|groupSize<=0)return(0)
  if(nGroups>n)return(0)
  if(nGroups>1){
    children<-sapply(1:(nGroups-1),function(x)chooseAtLeastOneFromEach(n,nGroups-x,groupSize))
    nChildren<-sapply(1:(nGroups-1),function(x)choose(nGroups,x))
  }else{
    nChildren<-0
    children<-0
  }
  answer<-choose(nGroups*groupSize,n)-sum(children*nChildren)
  return(answer)
}



#' Calculate the shared branch length from a matrix of taxonomic labels
#'
#' @param xx a matrix with a row for each read and a column for each taxonomic assignment
#' @param yy a matrix with a row for each read and a column for each taxonomic assignment
#' @param weighted logical indicating whether to sum the differences in reads for each branch or if FALSE to look only at presence absence
#' @param checkUpstream a logical indicating whether to make sure taxa with the same name are counted separetely if their upstream taxa differ e.g. two species Bos taurus and Repipta taurus. If you are sure this is not a problem, then computation can be reduced by skipping this step.
#' @return the proportion of the shared branches 
#' @export
#' @examples
#' x<-matrix(c('a','b','a','c'),ncol=2,byrow=TRUE)
#' y<-matrix(c('a','b','a','d'),ncol=2,byrow=TRUE)
#' unifracMatrix(x,x)
#' unifracMatrix(x,y)
unifracMatrix<-function(xx,yy,weighted=TRUE,checkUpstream=TRUE){
  n<-ncol(xx)
  xx[is.na(xx)]<-'__NAFILLER__'
  yy[is.na(yy)]<-'__NAFILLER__'
  if(ncol(yy)!=n)stop('Number of taxonomic ranks not the same')
  if(checkUpstream){
    pastedX<-t(apply(xx,1,function(taxas)Reduce(function(hh,tt)paste(hh,tt,sep='_|_'),taxas,accumulate=TRUE)))
    pastedY<-t(apply(yy,1,function(taxas)Reduce(function(hh,tt)paste(hh,tt,sep='_|_'),taxas,accumulate=TRUE)))
  }else{
    pastedX<-xx
    pastedY<-yy
  }
  dists<-do.call(rbind,lapply(n:1,function(ii){
    xTaxa<-table(pastedX[,ii])
    xTaxa<-xTaxa/sum(xTaxa)
    yTaxa<-table(pastedY[,ii])
    yTaxa<-yTaxa/sum(yTaxa)
    #slightly slower
    #merged<-merge(data.frame('name'=names(xTaxa),'x'=as.numeric(xTaxa)),data.frame('name'=names(yTaxa),'y'=as.numeric(yTaxa)),all=TRUE)
    allNames<-unique(c(names(yTaxa),names(xTaxa)))
    merged<-cbind('x'=xTaxa[allNames],'y'=yTaxa[allNames])
    merged[is.na(merged)]<-0
    rownames(merged)<-allNames
    if(weighted){
      dists<-apply(merged[,c('x','y'),drop=FALSE],1,diff)
      out<-c(sum(abs(dists)),2)
    }else{
      isDiff<-sum(apply(merged[,c('x','y'),drop=FALSE]>0,1,function(zz)zz[1]!=zz[2]))
      out<-c(isDiff,nrow(merged))
    }
    return(out)
  }))
  return(sum(dists[,1])/sum(dists[,2]))
}

