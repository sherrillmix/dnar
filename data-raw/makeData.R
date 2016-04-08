compressSave<-function(...,file){
	save(...,file=file)
	tools::resaveRdaFiles(file)
}

#' Base codes used to indicate ambiguous bases
ambiguousBaseCodes<-c(
	'R'='AG',
	'Y'='CT',
	'M'='AC',
	'K'='GT',
	'S'='CG',
	'W'='AT',
	'H'='ACT',
	'B'='CGT',
	'V'='ACG',
	'D'='AGT',
	'N'='ACGT'
)
reverseAmbiguous<-structure(names(ambiguousBaseCodes),.Names=ambiguousBaseCodes)
compressSave(ambiguousBaseCodes,file='data/ambiguousBaseCodes.RData')
compressSave(reverseAmbiguous,file='data/reverseAmbiguous.RData')




#data.frame of sam flags
samFlags<-data.frame('short'=c('paired','properPair','unmapped','mateUnmapped','reverse','mateReverse','first','second','notPrimary','fail','dupe'),'desc'=c('read paired','read mapped in proper pair','read unmapped','mate unmapped','read reverse strand','mate reverse strand','first in pair','second in pair','not primary alignment','read fails platform/vendor quality checks','read is PCR or optical duplicate'),stringsAsFactors=FALSE)
samFlags$bit<-2^(0:(nrow(samFlags)-1))
rownames(samFlags)<-samFlags$short
compressSave(samFlags,file='data/samFlags.RData')


#data.frame of amino acids
aminoAcids<-data.frame('codon'=c('UUU','UUC','UCU','UCC','UAU','UAC','UGU','UGC','UUA','UCA','UAA','UGA','UUG','UCG','UAG','UGG','CUU','CUC','CCU','CCC','CAU','CAC','CGU','CGC','CUA','CUG','CCA','CCG','CAA','CAG','CGA','CGG','AUU','AUC','ACU','ACC','AAU','AAC','AGU','AGC','AUA','ACA','AAA','AGA','AUG','ACG','AAG','AGG','GUU','GUC','GCU','GCC','GAU','GAC','GGU','GGC','GUA','GUG','GCA','GCG','GAA','GAG','GGA','GGG'),'abbr'=c('Phe','Phe','Ser','Ser','Tyr','Tyr','Cys','Cys','Leu','Ser','Xxx','Xxx','Leu','Ser','Xxx','Trp','Leu','Leu','Pro','Pro','His','His','Arg','Arg','Leu','Leu','Pro','Pro','Gln','Gln','Arg','Arg','Ile','Ile','Thr','Thr','Asn','Asn','Ser','Ser','Ile','Thr','Lys','Arg','Met','Thr','Lys','Arg','Val','Val','Ala','Ala','Asp','Asp','Gly','Gly','Val','Val','Ala','Ala','Glu','Glu','Gly','Gly'),'code'=c('F','F','S','S','Y','Y','C','C','L','S','X','X','L','S','X','W','L','L','P','P','H','H','R','R','L','L','P','P','Q','Q','R','R','I','I','T','T','N','N','S','S','I','T','K','R','M','T','K','R','V','V','A','A','D','D','G','G','V','V','A','A','E','E','G','G'),'name'=c('Phenylalanine','Phenylalanine','Serine','Serine','Tyrosine','Tyrosine','Cysteine','Cysteine','Leucine','Serine','Stop (ochre)','Stop (opal)','Leucine','Serine','Stop (amber)','Tryptophan','Leucine','Leucine','Proline','Proline','Histidine','Histidine','Arginine','Arginine','Leucine','Leucine','Proline','Proline','Glutamine','Glutamine','Arginine','Arginine','Isoleucine','Isoleucine','Threonine','Threonine','Asparagine','Asparagine','Serine','Serine','Isoleucine','Threonine','Lysine','Arginine','Methionine','Threonine','Lysine','Arginine','Valine','Valine','Alanine','Alanine','Aspartic acid','Aspartic acid','Glycine','Glycine','Valine','Valine','Alanine','Alanine','Glutamic acid','Glutamic acid','Glycine','Glycine'),row.names=c('UUU','UUC','UCU','UCC','UAU','UAC','UGU','UGC','UUA','UCA','UAA','UGA','UUG','UCG','UAG','UGG','CUU','CUC','CCU','CCC','CAU','CAC','CGU','CGC','CUA','CUG','CCA','CCG','CAA','CAG','CGA','CGG','AUU','AUC','ACU','ACC','AAU','AAC','AGU','AGC','AUA','ACA','AAA','AGA','AUG','ACG','AAG','AGG','GUU','GUC','GCU','GCC','GAU','GAC','GGU','GGC','GUA','GUG','GCA','GCG','GAA','GAG','GGA','GGG'),stringsAsFactors=FALSE)
compressSave(aminoAcids,file='data/aminoAcids.RData')




