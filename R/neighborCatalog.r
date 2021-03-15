#' Subfunction to catalog protein families in IMG metadata file
#'
#' This function identifies and quantifies protein families found in genomic neighborhoods (as defined by IMG metadata)
#' @param imgGenesTrimmed Data frame from analyzeNeighbors containing metadata for genes of interest
#' @param imgNeighborsTrimmed Data frame from analyzeNeighbors containing metadata for neighbors of genes of interest
#' @param geneName Name of gene of interest as string
#' @param coreGeneName Name of gene of interest as string, without any suffix.
#' @param neighborThreshold Abundance of protein family required to include it in further analysis (0-1)
#' @return List of relevant protein families
#' @export
#' @examples
#' neighborCatalogOut <- neighborCatalog(imgGenesTrimmed=imgGenesTrimmed, imgNeighborsTrimmed = imgNeighborsTrimmed, neighborNumber=neighborNumber, includeGene=includeGene, geneName=geneName, coreGeneName = coreGeneName)
#'
neighborCatalog <- function(imgGenesTrimmed = imgGenesTrimmed, imgNeighborsTrimmed = imgNeighborsTrimmed, geneName = geneName, neighborThreshold = neighborThreshold, coreGeneName = coreGeneName) { 
  fileDate <- format(Sys.Date(),format="%Y%m%d")
  fileName <- paste(fileDate,"_neighborCatalog_",geneName,sep="")
  finaltxtname <- paste(fileName,"_family-abundance.txt",sep="")
    ## identifying extant families
  numGenes <- length(unique(imgNeighborsTrimmed$source_gene_oid))
  pfamFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$Pfam, "pfam")==TRUE)
  pfamFams <- pfamFams$Pfam
  allFams <- pfamFams
  tigrFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$Tigrfam, "TIGR")==TRUE)
  tigrFams <- tigrFams$Tigrfam
  allFams <- append(allFams, tigrFams)
  iprFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$InterPro, "IPR")==TRUE)
  iprFams <- iprFams$InterPro
  allFams <- append(allFams, iprFams)
  hypoFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$Hypofam, "hypo")==TRUE)
  hypoFams <- hypoFams$Hypofam
  allFams <- append(allFams, hypoFams)
  ## to implement: IMG fams
  ## using the stupid XML traces since the IDs are number-only
  ##imgFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$IMG.Term, "br")==TRUE)
  ##imgFams <- imgFams$IMG.Term
  ##allFams <- append(allFams, imgFams)
  uniquePfam <- unique(pfamFams)
  uniqueTigrfam <- unique(tigrFams)
  uniqueIprfam <- unique(iprFams)
  uniqueHypofam <- unique(hypoFams)
  ##uniqueIMGfam <- unique(imgFams)
  uniqueFams <- list()
  ## finding unique pfams and dealing with the punctuation
  for (i in 1:length(uniquePfam)) {
    tempSplit <- unlist(strsplit(uniquePfam[i],"[[:punct:][:space:]]+"))
    tempSplit <- tempSplit[stringr::str_detect(tempSplit, "pfam")]
    uniqueFams <- append(uniqueFams, tempSplit)
  }
    ## finding unique tigrfams and dealing with the punctuation
  for (i in 1:length(uniqueTigrfam)) {
    tempSplit <- unlist(strsplit(uniqueTigrfam[i],"[[:punct:][:space:]]+"))
    tempSplit <- tempSplit[stringr::str_detect(tempSplit, "TIGR") ]
    uniqueFams <- append(uniqueFams, tempSplit)
  }
  ## finding unique InterPro fams and dealing with the punctuation
  for (i in 1:length(uniqueIprfam)) {
    tempSplit <- unlist(strsplit(uniqueIprfam[i],"[[:punct:][:space:]]+"))
    tempSplit <- tempSplit[stringr::str_detect(tempSplit, "IPR") ]
    uniqueFams <- append(uniqueFams, tempSplit)
  }
  ## finding unique hypofams, if any
  for (i in 1:length(uniqueHypofam)) {
    ## no weirdo punctuation here!
    tempSplit <- unlist(strsplit(uniqueHypofam[i],"[[:space:]]+"))
    tempSplit <- tempSplit[stringr::str_detect(tempSplit, "hypo") ]
    uniqueFams <- append(uniqueFams, tempSplit)
  }
  ## finding unique IMGfams, if any
  ##for (i in 1:length(uniqueIMGfam)) {
    ## haven't seen multi-annotations for this yet, but....
    ##tempSplit <- unlist(strsplit(uniqueIMGfam[i],"[[:punct:][:space:]]+"))
    ##tempSplit <- tempSplit[stringr::str_detect(tempSplit, "[[:digit:]][[:digit:]][[:digit:]][[:digit:]]") ]
    ##uniqueFams <- append(uniqueFams, tempSplit)
  ##}
  uniqueFams <- unlist(unique(uniqueFams))
  uniqueFamNum <- length(uniqueFams)
  tempFamList <- list(fams=unlist(uniqueFams), abund=sprintf("",1:uniqueFamNum))
  uDataFams <- data.frame(tempFamList, stringsAsFactors=FALSE)
  trimmedFams <- list()
    ## calculating abundance
    ## and adding families that pass to the list
  for (j in 1:length(uDataFams[,1]))  {
    tempAbund <- grep(uDataFams[j,1], allFams)
    uDataFams[j,2] <- length(tempAbund)
    if (length(tempAbund)/numGenes >= neighborThreshold)  {
      if(exists(x="trimmedFams")) {
        trimmedFams <- append(trimmedFams, uDataFams[j,1])
        next
      } else   {
        trimmedFams <-  uDataFams[j,1]
        next
      }
    } else {
      next
    }
  }
  trimmedFams <- unlist(trimmedFams)
  write.table(uniqueFams, file=finaltxtname, sep="\t", quote=FALSE)
  familyList <- trimmedFams
  print("Protein families in neighboring genes catalogued.")
  return(familyList)
}

