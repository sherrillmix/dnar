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
jensenShannon<-function(x,y,base=2){
	propX<-x/ifelse(sum(x)==0,1,sum(x))
	propY<-y/ifelse(sum(y)==0,1,sum(y))
	shannonCalc<-shannon((propX+propY)/2,base,FALSE)-shannon(propX,base,FALSE)/2-shannon(propY,base,FALSE)/2
	kullbackCalc<-.5*kullback(propX,(propY+propX)/2,base,FALSE)+.5*kullback(propY,(propY+propX)/2,base,FALSE)
	if(round(shannonCalc,5)!=round(kullbackCalc,5))stop(simpleError('Problem calculating shannon jensen'))
	return(kullbackCalc)
}

#' Calculate Kullback-Leibler divergence between two probability distributions
#'
#' @param x Vector of counts or proportions of N elements
#' @param y Second vector of counts or proportions of N elements 
#' @param base Base of logarithm
#' @param standardize If TRUE divide x and y by sum(x) and sum(y)
#' @export
#' @return Kullback-Leibler divergence
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

#' Calculate chao diversity index
#'
#' @param counts A vector of counts with one entry per "species"
#' @export
#' @return Chao index
chao<-function(counts){
	counts<-counts[counts>0]
	return(length(counts)+sum(counts==1)*(sum(counts==1)-1)/2/(sum(counts==2)+1))
}

#' Calculate rarefaction using boostrapping
#' 
#' @param species Ids of species
#' @param counts Corresponding counts of species
#' @param samples Vector of numbers of draws to calculate rarefaction at
#' @param reps How many random samples to take at each step
#' @param quants Quantiles to return
#' @param chaoAdjust If TRUE calculate chao-predicted species number on each random draw
#' @param debug If TRUE display debugging information
#' @param minCount Remove and species with less than minCount counts
#' @param replaceSpecies If TRUE sample with replacements. If FALSE sample without replacement.
#' @export
#' @return Dataframe of calculated quantiles with rownames of the number of samples drawn
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
		estimate<-stats::quantile(numSpecies,quants)
		return(estimate)
	},reps,species,debug)
	output<-do.call(rbind,output)
	return(list(output,samples))
}

#' Calculate rarefaction using formula
#'
#' @param sample Vector of numbers of individuals per species
#' @param step Size of sampling steps for rarefaction
#' @param maxN The maximum number of counts to rarefy to
#' @export
#' @return Matrix with two columns; number of draws and rarefaction value
quickRare<-function(sample,step=10,maxN=sum(sample)){
	sampleSize<-20;
	steps<-unique(c(seq(step,maxN,step),maxN))
	output<-sapply(steps,function(x)rareEquation(sample,x))
	return(data.frame('rare'=output,'sampleN'=steps))
}

#' Calculate rarefaction of single sample using equation
#'
#' @param speciesCounts Vector of counts for each "species" e.g. c(10,100,5)
#' @param sampleSize Single value of number of draws from sample 
#' @export
#' @return Rarefaction value
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



#' Calculate the probability of observing observedX species from nGroups groups with groupsize members with n observations
#'
#' @param nGroups number of equally sized groups
#' @param groupSize size of the equally sized groups
#' @param n number of observations pulled from the individuals
#' @param observedX number of species observed
#' @return probability of observing observedX species from nGroups groups with groupsize members with n observations
pRare<-function(nGroups,groupSize,n,observedX){
	choose(nGroups,observedX)*chooseAtLeastOneFromEach(observedX,groupSize,n)/choose(nGroups*groupSize,n)
}

#' Calculate the number of combinations possible from nGroups groups with groupsize members and n observations
#'
#' @param nGroups number of equally sized groups
#' @param groupSize size of the equally sized groups
#' @param n number of observations pulled from the individuals
#' @export
#' @return number of combinations containing at least one individual from each group
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


