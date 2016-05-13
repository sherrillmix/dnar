#' Rainbow based on circle of lab color space
#'
#' @param n number of colors desired
#' @param start start angle (in proportion of circle) in lab space
#' @param end end angle (in proportion of circle) in lab space
#' @param alpha transparency (0 for transparent, 1 for opaque)
#' @param lightScale exponent to scale lightness by 
#' @param lightMultiple multiple to scale lightness by
#' @export 
#' @return vector of color strings
#' @examples 
#' cols<-rainbow.lab(100)
#' plot(1:100,cex=3,col=cols,lwd=5)
#' cols<-rainbow.lab(100,lightScale=1)
#' plot(1:100,col=cols,lwd=5)
rainbow.lab<-function(n,start=1.5,end=-3,alpha=1,lightScale=0,lightMultiple=.7){
	#something is going crazy with R's implementation of lab
	#angles<-seq(start*2*pi,end*2*pi,length.out=n)
	#a<-sin(angles)*radius
	#b<-cos(angles)*radius
	#Lab<-cbind(luminance,a,b)
	#srgb<-convertColor(Lab,from="Lab",to="sRGB")
	#cols<-rgb(srgb[,1],srgb[,2],srgb[,3],alpha)
	#return(cols)
	l<-seq(1,.4,length.out=n)^lightScale*lightMultiple
	rgb<-cl2pix(seq(0,1,length.out=n),l,start=start,end=end,toColor=FALSE)
	return(grDevices::rgb(rgb,alpha=alpha))
}


#http://davidad.net/colorviz/
#/* Convert from L*a*b* doubles to XYZ doubles
#* Formulas drawn from http://en.wikipedia.org/wiki/Lab_color_spaces
#*/
lab2xyz<-function(lab) {
	if(is.null(dim(lab)))matrix(lab,nrow=1)
	finv<-function(t) ifelse(t>6/29,t^3,3*(6/29)^2*(t-4/29))
	sl = (lab[,1]+0.16)/1.16
	ill = c(0.9643,1.00,0.8251) #D50
	out<-matrix(ill,nrow=nrow(lab),ncol=3,byrow=TRUE) * finv(cbind(sl,sl+lab[,2]/5,sl-lab[,3]/2))
	return(out)
}

#http://davidad.net/colorviz/
#/* Convert from XYZ doubles to sRGB bytes
#* Formulas drawn from http://en.wikipedia.org/wiki/Srgb
#*/
xyz2rgb<-function(xyz) {
	rgb <- cbind(3.2406*xyz[,1] - 1.5372*xyz[,2] - 0.4986*xyz[,3], -0.9689*xyz[,1] + 1.8758*xyz[,2] + 0.0415*xyz[,3], 0.0557*xyz[,1] - 0.2040*xyz[,2] + 1.0570*xyz[,3])
	clip <- rgb>1|rgb<0
	if(any(clip)) {
		rgb[rgb>1]<-1
		rgb[rgb<0]<-0
	}
	#Uncomment the below to detect clipping by making clipped zones red.
	#if(clip) {rl=1.0;gl=bl=0.0;}
	correct<-function(cl) {
		a <- 0.055
		return (ifelse(cl<=0.0031308,12.92*cl,(1+a)*cl^(1/2.4)-a))
	}
	out<-correct(rgb)
	return(out)
}

#http://davidad.net/colorviz/
#/* Convert from LAB doubles to sRGB bytes 
#* (just composing the above transforms)
#*/
lab2rgb<-function(lab){
  xyz<-lab2xyz(lab)
  xyz2rgb(xyz)
}

#http://davidad.net/colorviz/
#/* Convert from a qualitative parameter c and a quantitative parameter l to a 24-bit pixel
 #* These formulas were invented by davidad to obtain maximum contrast without going out of gamut
 #* if the parameters are in the range 0-1
 #*/
cl2pix<-function(c, l,start=-3,end=4,toColor=TRUE) {
  L <- l #L of L*a*b*
  angle <- start+c*(end-start)
  r <- l*0.4+0.1 #~chroma
  a <- sin(angle+start)*r
  b <- cos(angle+start)*r
  out<-lab2rgb(cbind(L,a,b))
  if(toColor)out<-grDevices::rgb(out)
  return(out)
}


