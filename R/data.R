
#data.frame of sam flags
samFlags<-data.frame('short'=c('paired','properPair','unmapped','mateUnmapped','reverse','mateReverse','first','second','notPrimary','fail','dupe'),'desc'=c('read paired','read mapped in proper pair','read unmapped','mate unmapped','read reverse strand','mate reverse strand','first in pair','second in pair','not primary alignment','read fails platform/vendor quality checks','read is PCR or optical duplicate'),stringsAsFactors=FALSE)
samFlags$bit<-2^(0:(nrow(samFlags)-1))
rownames(samFlags)<-samFlags$short
#data.frame of amino acids
aminoAcids<-data.frame('codon'=c('UUU','UUC','UCU','UCC','UAU','UAC','UGU','UGC','UUA','UCA','UAA','UGA','UUG','UCG','UAG','UGG','CUU','CUC','CCU','CCC','CAU','CAC','CGU','CGC','CUA','CUG','CCA','CCG','CAA','CAG','CGA','CGG','AUU','AUC','ACU','ACC','AAU','AAC','AGU','AGC','AUA','ACA','AAA','AGA','AUG','ACG','AAG','AGG','GUU','GUC','GCU','GCC','GAU','GAC','GGU','GGC','GUA','GUG','GCA','GCG','GAA','GAG','GGA','GGG'),'abbr'=c('Phe','Phe','Ser','Ser','Tyr','Tyr','Cys','Cys','Leu','Ser','Ochre','Opal','Leu','Ser','Amber','Trp','Leu','Leu','Pro','Pro','His','His','Arg','Arg','Leu','Leu','Pro','Pro','Gln','Gln','Arg','Arg','Ile','Ile','Thr','Thr','Asn','Asn','Ser','Ser','Ile','Thr','Lys','Arg','Met','Thr','Lys','Arg','Val','Val','Ala','Ala','Asp','Asp','Gly','Gly','Val','Val','Ala','Ala','Glu','Glu','Gly','Gly'),'code'=c('F','F','S','S','Y','Y','C','C','L','S','X','X','L','S','X','W','L','L','P','P','H','H','R','R','L','L','P','P','Q','Q','R','R','I','I','T','T','N','N','S','S','I','T','K','R','M','T','K','R','V','V','A','A','D','D','G','G','V','V','A','A','E','E','G','G'),'name'=c('Phenylalanine','Phenylalanine','Serine','Serine','Tyrosine','Tyrosine','Cysteine','Cysteine','Leucine','Serine','Stop','Stop','Leucine','Serine','Stop','Tryptophan','Leucine','Leucine','Proline','Proline','Histidine','Histidine','Arginine','Arginine','Leucine','Leucine','Proline','Proline','Glutamine','Glutamine','Arginine','Arginine','Isoleucine','Isoleucine','Threonine','Threonine','Asparagine','Asparagine','Serine','Serine','Isoleucine','Threonine','Lysine','Arginine','Methionine','Threonine','Lysine','Arginine','Valine','Valine','Alanine','Alanine','Aspartic acid','Aspartic acid','Glycine','Glycine','Valine','Valine','Alanine','Alanine','Glutamic acid','Glutamic acid','Glycine','Glycine'),row.names=c('UUU','UUC','UCU','UCC','UAU','UAC','UGU','UGC','UUA','UCA','UAA','UGA','UUG','UCG','UAG','UGG','CUU','CUC','CCU','CCC','CAU','CAC','CGU','CGC','CUA','CUG','CCA','CCG','CAA','CAG','CGA','CGG','AUU','AUC','ACU','ACC','AAU','AAC','AGU','AGC','AUA','ACA','AAA','AGA','AUG','ACG','AAG','AGG','GUU','GUC','GCU','GCC','GAU','GAC','GGU','GGC','GUA','GUG','GCA','GCG','GAA','GAG','GGA','GGG'),stringsAsFactors=FALSE)



#http://life.nthu.edu.tw/~fmhsu/rasframe/SHAPELY.HTM
aminoColors<-data.frame('aa'=c("ASP","GLU","CYS","MET","LYS","ARG","SER","THR","PHE","TYR","ASN","GLN","GLY","LEU","VAL","ILE","ALA","TRP","HIS","PRO"),'col'=c("#E60A0A","#E60A0A","#E6E600","#E6E600","#145AFF","#145AFF","#FA9600","#FA9600","#3232AA","#3232AA","#00DCDC","#00DCDC","#EBEBEB","#0F820F","#0F820F","#0F820F","#C8C8C8","#B45AB4","#8282D2","#DC9682"),stringsAsFactors=FALSE)
tmp<-with(aminoAcids[!duplicated(aminoAcids[,c('code')]),],data.frame('letter'=code,row.names=toupper(abbr),stringsAsFactors=FALSE))
rownames(aminoColors)<-tmp[aminoColors$aa,'letter']
tmpAngles1<-cos(1+1:nrow(aminoColors)/nrow(aminoColors)*pi)
tmpAngles2<-sin(1:nrow(aminoColors)/nrow(aminoColors)*pi)
tmpAngles1<-tapply(tmpAngles1,aminoColors$col,c)
tmpAngles2<-tapply(tmpAngles2,aminoColors$col,c)
aminoColors$spreadCol<-ave(aminoColors$col,aminoColors$col,FUN=function(x){
	if(length(x)==1)return(x)
	spacer<-20
	angles1<-tmpAngles1[[x[1]]]
	angles2<-tmpAngles2[[x[1]]]
	spacing<-seq((length(x)-1)*-spacer,(length(x)-1)*spacer,length.out=length(x))
	lab<-convertColor(t(col2rgb(x[1])),from='sRGB',to='Lab',scale.in=255)[rep(1,length(x)),]
	lab[,'b']<-lab[,'b']-spacing*angles1
	lab[,'a.x']<-lab[,'a.x']+spacing*angles2
	rgbs<-convertColor(lab,from='Lab',to='sRGB')
	return(rgb(rgbs))
	#plot(1:nrow(aminoColors),col=aminoColors$col,lwd=30,xaxt='n');points(1:nrow(aminoColors),1:nrow(aminoColors)+1,col=aminoColors$spreadCol,lwd=20);axis(1,1:nrow(aminoColors),aminoColors$aa,las=2)
})
rm(tmpAngles1,tmpAngles2)
