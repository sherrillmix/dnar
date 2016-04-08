#' Ambiguous base codes
#'
#' @docType data
#' @format A vector with names corresponding to ambiguous base codes and items giving concatenated strings of corresponding bases
#' @references \url{https://en.wikipedia.org/wiki/Nucleic_acid_notation}
#' @source system.file("data-raw", "makeData.R", package = "dnar")
"ambiguousBaseCodes"

#' Reverse lookup for ambiguous base codes
#'
#' @docType data
#' @format A vector with names giving sortedconcatenated strings of bases for ambiguous base codes with items giving the corresponding ambiguous base code
#' @references \url{https://en.wikipedia.org/wiki/Nucleic_acid_notation}
#' @source system.file("data-raw", "makeData.R", package = "dnar")
"reverseAmbiguous"

#' SAM/BAM flags
#'
#' @docType data
#' @format A data frame with 64 rows one for each codon and columns:
#' \describe{
#'   \item{short}{short name for SAM flags}
#'   \item{desc}{longer description of SAM flags}
#'   \item{bit}{bit corresponding to the SAM flag}
#' }
#' @references \url{https://samtools.github.io/hts-specs/SAMv1.pdf}
#' @source system.file("data-raw", "makeData.R", package = "dnar")
"samFlags"

#' Amino acid for RNA codons
#'
#' @docType data
#' @format A data frame with 64 rows, one for each RNA codon, with 4 columns:
#' \describe{
#'   \item{codon}{3 RNA bases coding for an amino acid (DNA Ts will be coded as RNA Us)}
#'   \item{abbr}{3 letter abbreviation for amino acid (stop codons are coded as Xxx)}
#'   \item{code}{single letter code for amino acid. Stop codons are coded as X}
#'   \item{name}{name of the amino acids. Stop codons given their color in parenthesis}
#' }
#' @references \url{https://en.wikipedia.org/wiki/Genetic_code}
#' @source system.file("data-raw", "makeData.R", package = "dnar")
"aminoAcids"

