#' Subfunction to quantify protein family presence in the proximity of genes of interest
#'
#' This subfunction identifies whether members are a given protein family are present in the genomic neighborhood of a gene of interest
#' @param imgNeighborsTrimmed Data frame from analyzeNeighbors containing metadata for neighbors of genes of interest
#' @param familyList Data from from analyzeNeighbors/neighborCatalog containing a list of protein families for further analysis
#' @param geneName Name of gene of interest as string
#' @param coreGeneName Name of gene of interest as string, without any suffix.
#' @return Table of binary data indicating presence of absence of a protein family in the proximity of a gene of interest
#' @export
#' @importFrom utils read.csv write.csv write.table read.table
#' @importFrom data.table :=
#' @importFrom magrittr %>%
#' @examples
#' neighborHereOut <- neighborHere(imgNeighborsTrimmed=imgNeighborsTrimmed, familyList = familyList, geneName = geneName, coreGeneName = coreGeneName) 
#'
neighborHere <- function(imgNeighborsTrimmed = imgNeighborsTrimmed, familyList = familyList, geneName = geneName, coreGeneName = coreGeneName) {
  fileDate <- format(Sys.Date(),format="%Y%m%d")
  fileName <- paste(fileDate,"_neighborHere_",geneName,sep="")
  finalCSV <- paste(fileName,"_neighborBinary.csv",sep="")
                                        # gene_oid and source_gene_oid start neighborBinary
imgNeighborsTrimmed <- imgNeighborsTrimmed %>% dplyr::distinct()
  nbTemp <- list(gene_oid = as.character(imgNeighborsTrimmed$gene_oid), source_gene_oid = as.character(imgNeighborsTrimmed$source_gene_oid))
  neighborBinary <- data.frame(nbTemp, stringsAsFactors = FALSE)
  famNum <- length(familyList)
  geneNum <- length(unique(imgNeighborsTrimmed$gene_oid))
  ## loop through the families
  for (i in 1:famNum){
    tempFam <- vector()
    tempFamName <-familyList[i]
    tempFamName <- as.character(tempFamName)
    ## loop through the neighbor-genes
    for (j in 1:geneNum){
      pfam <- grepl(familyList[i], imgNeighborsTrimmed$Pfam[j])
      tigrfam <- grepl(familyList[i], imgNeighborsTrimmed$Tigrfam[j])
      iprfam <- grepl(familyList[i], imgNeighborsTrimmed$InterPro[j])
      hypofam <- grepl(familyList[i], imgNeighborsTrimmed$Hypofam[j])
      ## need to alter how this step works to include img terms
      ## since they have no prefix, they could have a false hit in a pfam/tigrfam
      ##imgfam <- grepl(familyList[i], imgNeighborsTrimmed$IMG.Term[j])
      if (any(c(pfam,tigrfam,iprfam,hypofam))==TRUE) {
        tempFam[j] <- 1
      } else {
        tempFam[j] <- 0
      }
    }
    ## add the data for that family as a new column in neighborBinary
    neighborBinary <- neighborBinary %>% dplyr::mutate(!! tempFamName := tempFam)
  }
                                        # export the table
  write.csv(neighborBinary, file=finalCSV)
  print("Binary family values assigned to neighboring genes.")
  return(neighborBinary)
}
