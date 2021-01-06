#' Function for trimming IMG gene metadata and neighbor metadata to match EFI-EST-generated representative nodes
#'
#' After SSN generation via EFI-EST, this function ties together the EFI-EST gene IDs with gene IDs in existing metadata files and produces trimmed metadata and sequence files for the genes and their neighbors using the list of representative nodes at the desired cutoff.
#' @param imgGenes Filename for the (edited) IMG metadata file for genes of interest.  Text ("geneFile.txt")
#' @param imgNeighbors Filename for the (edited) IMG metadata file for neighboring genes.  Text ("neighborFile.txt")
#' @param imgGeneSeqs Filename for the FASTA amino acid sequence file for genes of interest.  Text ("geneFile.fa")
#' @param imgNeighborSeqs Filename for the FASTA amino acid sequence file for neighboring genes.  Text ("neighborFile.fa")
#' @param geneName Name of gene of interest for filenames. Text ("genE")
#' @param efiFullMetadata Filename for the Cytoscape node table for the full EFI-EST network. Text ("fullnetwork.txt")
#' @param efiFinalMetadata Filename for the Cytoscape node table for the EFI-EST network at the final repnode %ID cutoff. Text ("finalnetwork.txt")
#' @return List with the trimmed gene and neighbor metadata
#' @export
#' @examples
#' repnodeTrimOut <- repnodeTrim(imgGenes = "geneFile.txt", imgNeighbors="neighborFile.txt", imgGeneSeqs = "geneSeqs.fa", imgNeighborSeqs = "neighborSeqs.fa", geneName = "genE", efiFullMetadata = "efiFullMetadata.csv", efiFinalMetadata = "efiFinalMetadata.csv")
#'
repnodeTrim <- function(imgGenes = imgGenes, imgNeighbors=imgNeighbors, imgGeneSeqs = imgGeneSeqs, imgNeighborSeqs = imgNeighborSeqs, geneName = geneName, efiFullMetadata = efiFullMetadata, efiFinalMetadata = efiFinalMetadata) {
    # leaving trimming of the gene seq as optional - 
  if(exists(x="imgGenes") == FALSE | exists(x="imgNeighbors") == FALSE | exists(x="geneName") == FALSE | exists(x="imgGeneSeqs") == FALSE | exists(x="imgNeighborSeqs") == FALSE | exists(x="efiFullMetadata") == FALSE | exists(x="efiFinalMetadata") == FALSE) {
    print("Missing a required term (genes and neighbors metadata files, or name of gene of interest")
    return(0)
  }
  ## filenames
  fileDate <- format(Sys.Date(),format="%Y%m%d")  
  fileName <- paste(fileDate,"_repnodeTrim_",geneName,sep="")
  fullGenesFile <- paste(fileName,"_imgGenes.txt",sep="")
  fullNeighborsFile <- paste(fileName,"_imgNeighbors.txt",sep="")
  repGenesFile <- paste(fileName,"_repnodeGenes.txt",sep="")
  repNeighborsFile <- paste(fileName,"_repnodeNeighbors.txt",sep="")
  repGeneSeqsFile <- paste(fileName,"_repnodeGeneSeqs.fa",sep="")
  repNeighborSeqsFile <- paste(fileName,"_repnodeNeighborSeqs.fa",sep="")
  ## importing data
  if(typeof(imgGenes) == "character") {
    imgGenesTemp <-read.csv(imgGenes, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  } else {
    imgGenesTemp <- imgGenes
  }
  if(typeof(imgNeighbors) == "character") {
    imgNeighborsTemp <- read.csv(imgNeighbors, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  } else {
    imgNeighborsTemp <- imgNeighbors
  }
  imgNeighborSeqs <- seqinr::read.fasta(file=imgNeighborSeqs, seqtyp="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
  imgGeneSeqs <- seqinr::read.fasta(file=imgGeneSeqs, seqtyp="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
  efiFullNodes <- read.csv(efiFullMetadata, header=TRUE, stringsAsFactors=FALSE)
  efiFinalNodes <- read.csv(efiFinalMetadata, header=TRUE, stringsAsFactors=FALSE)
  ## starting to reformat
  efiFullNodes <- efiFullNodes %>% dplyr::rename(gene_oid=Description)
  efiFullNodes <- efiFullNodes %>% dplyr::rename(efi_oid=name)
  efiFinalNodes <- efiFinalNodes %>% dplyr::rename(gene_oid=Description)
  efiFinalNodes <- efiFinalNodes %>% dplyr::rename(efi_oid=name)
  efiFullNodes$gene_oid <- as.character(efiFullNodes$gene_oid)
  imgGenesTemp$gene_oid <- as.character(imgGenesTemp$gene_oid)
  imgNeighborsTemp$gene_oid <- as.character(imgNeighborsTemp$gene_oid)
  imgNeighborsTemp$source_gene_oid <- as.character(imgNeighborsTemp$source_gene_oid)
  ## adding EFI-generated faux uniprot IDs to the metadata
  efiGenesData <- imgGenesTemp %>% dplyr::left_join(efiFullNodes[,c("gene_oid","efi_oid")], by="gene_oid")
  ## what repnode does a gene have, and is it a repnode?
  repnodeList <- efiFinalNodes$shared.name
  efiGenesData <- efiGenesData %>% dplyr::mutate(isRepnode = ifelse(efi_oid %in% repnodeList,TRUE,FALSE))
  for (i in 1:length(efiGenesData$gene_oid)) {
    if (is.na(efiGenesData$efi_oid[i])) {
      next
    } else {
      efiGenesData$repnodeIs[i] <- efiFinalNodes$shared.name[grep(efiGenesData$gene_oid[i],efiFinalNodes$gene_oid)]
      next
    }
  }
  repGenesTrimmed <- efiGenesData %>% dplyr::filter(isRepnode == TRUE)
  efiFullNodes <- efiFullNodes %>% dplyr::rename(source_gene_oid=gene_oid)
  efiNeighborsData <- imgNeighborsTemp %>% dplyr::left_join(efiFullNodes[,c("source_gene_oid","efi_oid")], by="source_gene_oid")
  efiNeighborsData <- efiNeighborsData %>% dplyr::mutate(isRepnode = ifelse(efi_oid %in% repnodeList,TRUE,FALSE))
  repNeighborsTrimmed <- efiNeighborsData %>% dplyr::filter(isRepnode == TRUE)
  efiFullNodes <- efiFullNodes %>% dplyr::rename(gene_oid=source_gene_oid)
  ## exporting the datatables as full and repnode-only datases
  write.table(efiGenesData, file=fullGenesFile, row.names=F, col.names = T, sep="\t")
  write.table(efiNeighborsData, file=fullNeighborsFile, row.names=F, col.names = T, sep="\t")
  write.table(repGenesTrimmed, file=repGenesFile, row.names=F, col.names = T, sep="\t")
  write.table(repNeighborsTrimmed, file=repNeighborsFile, row.names=F, col.names = T, sep="\t")        
  ## trimming the sequence files for neighboring genes
  nSeqNames <- names(imgNeighborSeqs)
  nSeqTable <- data.frame()
  nReList <- list()
  for (i in 1:length(nSeqNames)) {
    nSeqTable[i,1] <- nSeqNames[i]
    nSeqTable[i,2] <- imgNeighborSeqs[[i]]
  }
  colnames(nSeqTable) <- c("gene_oid", "sequence")
  nRepSeqs <- nSeqTable %>% dplyr::filter(gene_oid %in% repNeighborsTrimmed$gene_oid)
  for (i in 1:length(nRepSeqs[,1])) {
    nReList[[nRepSeqs[i,1]]] <- nRepSeqs[i,2] 
  }
  repNeighborSeqs <- nReList
  write.fasta(repNeighborSeqs, names=names(repNeighborSeqs),file.out=repNeighborSeqsFile)
  ## trimming the sequence files for the repnodes themselves
  gSeqNnames <- names(imgGeneSeqs)
  gSeqTable <- data.frame()
  gReList <- list()
  for (i in 1:length(gSeqNnames)) {
    gSeqTable[i,1] <- gSeqNnames[i]
    gSeqTable[i,2] <- imgGeneSeqs[[i]]
  }
  colnames(gSeqTable) <- c("gene_oid", "sequence")
  gRepSeqs <- gSeqTable %>% dplyr::filter(gene_oid %in% repGenesTrimmed$gene_oid)
  for (i in 1:length(gRepSeqs[,1])) {
    gReList[[gRepSeqs[i,1]]] <- gRepSeqs[i,2] 
  }  
  repGeneSeqs <- gReList
  write.fasta(repGeneSeqs, names=names(repGeneSeqs),file.out=repGeneSeqsFile)
  print("Repnode metadata files generated.")           
  return(list(repGenesTrimmed=repGenesTrimmed,repNeighborsTrimmed=repNeighborsTrimmed))
}