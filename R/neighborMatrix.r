#' Subfunction to matrix-ify data about protein presence in a given genomic neighborhood
#'
#' This function uses binary present/absent data for various protein families in the proximity of a set of genes of interest and matrixifies it.
#' @param imgGenesTrimmed Data frame from analyzeNeighbors containing metadata for genes of interest
#' @param neighborBinary Data frame from analyzeNeighbors/neighborHere binary data regarding presence of protein families in close proximity
#' @param familyList Data from analyzeNeighbors/neighborCatalog listing common nearby protein families
#' @param geneName Name of gene of interest as string
#' @return Matrix object
#' @export
#' @importFrom utils read.csv write.csv write.table read.table
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @examples
#' \dontrun{
#' neighborMatrixOut <- neighborMatrix(imgGenesTrimmed=imgGenesTrimmed, 
#'                                     neighborBinary = neighborBinary, 
#'                                     familyList = familyList, 
#'                                     geneName = "genE")
#' }
#'
neighborMatrix <- function(imgGenesTrimmed = imgGenesTrimmed,
                           neighborBinary = neighborBinary,
                           familyList = familyList,
                           geneName = geneName) {
  fileDate <- format(Sys.Date(),format="%Y%m%d")
  fileName <- paste(fileDate,"_neighborMatrix_",geneName,sep="")
  finalCSV <- paste(fileName,".csv",sep="")
  allFalse <- list()
  familyNum <- length(familyList)
  for (i in 1:familyNum) {
    allFalse <- list(c(allFalse, FALSE))
  }
  matrixColNames <- c("gene_oid", familyList)
  matrixColNum <- length(matrixColNames)
  ## getting the list of ids for the gene of interest
  centralGenes <- unique(imgGenesTrimmed$gene_oid)
  centralGeneNum <- length(centralGenes)
  neighborMatrixData <- data.frame(matrix(NA, nrow=centralGeneNum, ncol=matrixColNum))
  colnames(neighborMatrixData) <- matrixColNames
  ## going through the gene ids and connecting the binary there/not there values for families of interests from the neighboring genes to the gene of interest
  colTally <- matrixColNum -1
  for (i in 1:centralGeneNum) {
    neighborMatrixData$gene_oid[i] <- centralGenes[i]
    ## getting the mini-set of neighbors of that source gene id
    nearestBinary <- neighborBinary %>% dplyr::filter(.data$source_gene_oid == centralGenes[i])
    for (j in 1:colTally) {
      if(any(grepl(1,nearestBinary[[matrixColNames[j+1]]])) == TRUE) {
        neighborMatrixData[[matrixColNames[j+1]]][i] <- 1
      } else {    
        neighborMatrixData[[matrixColNames[j+1]]][i] <- 0
      }
    }
  }
  write.csv(neighborMatrixData, finalCSV)
  print("Matrix of binary family values per gene of interest created.")
  return(neighborMatrixData)
}

