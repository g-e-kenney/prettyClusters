#' Trim .fa files based on a list of sequence IDs
#'
#' Simple function for trimming .fa files (based on simplified headers), using a list of sequence IDs.
#' @param seqInput Filename for the FASTA-formatted file containing protein sequences.  Text ("seqInput.fa")
#' @param trimInput Filename for the list of desired sequences, as a single column with "gene_oid" as a header.  Text ("trimInput.txt")
#' @param seqName Filename identifier to add do your file name, along with the date, trim marker, etc.  Text ("genE_subset")
#' @return Output files and data frame with trimmed sequences.
#' @export
#' @importFrom utils read.csv write.csv
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @examples
#' \dontrun{
#' trimFastaOut <- trimFasta(seqInput = "seqInput.fa", 
#'                           trimInput="trimInput.txt", 
#'                           seqName = "genE")
#' }
#' 
trimFasta <- function(seqInput = seqInput, trimInput=trimInput, seqName = seqName) {
  fileDate <- format(Sys.Date(),format="%Y%m%d")
  outputSeqsFile <- paste(fileDate,"_trimmed_",seqName,".fa",sep="")
  proteinSeqs <- seqinr::read.fasta(file=seqInput, seqtype="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
  trimList <- read.csv(file=trimInput, header=TRUE, sep=",", stringsAsFactors=FALSE)
  seqNames <- names(proteinSeqs)
  seqTable <- data.frame()
  reList <- list()
  for (i in 1:length(seqNames)) {
    seqTable[i,1] <- seqNames[i]
    seqTable[i,2] <- proteinSeqs[[i]]
  }
  colnames(seqTable) <- c("gene_oid", "sequence")
  trimmedSeqs <- seqTable %>% dplyr::filter(.data$gene_oid %in% trimList$gene_oid)
  for (i in 1:length(trimmedSeqs[,1])) {
    reList[[trimmedSeqs[i,1]]] <- trimmedSeqs[i,2] 
  }
  seqinr::write.fasta(reList, names=names(reList),file.out=outputSeqsFile)
  print(".fa file trimmed according to your list of desired sequences.")
  return(trimmedSeqs)
}
