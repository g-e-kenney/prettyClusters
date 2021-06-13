#' Subfunction to catalog protein families in IMG metadata file
#'
#' This function identifies and quantifies protein families found in genomic neighborhoods (as defined by IMG metadata)
#' @param imgGenesTrimmed Data frame from analyzeNeighbors containing metadata for genes of interest. Character string, required.
#' @param imgNeighborsTrimmed Data frame from analyzeNeighbors containing metadata for neighbors of genes of interest. Character string, required.
#' @param geneName Name of gene of interest as string. Character string, required.
#' @param neighborThreshold Abundance of protein family required to include it in further analysis. Number (0-1)
#' @param useInterPro Should InterPro families be used for analysis? Boolean, defaults to FALSE.
#' @param useHypo Should hypothetical protein families from prepNeighbors be used for analysis? Boolean, defaults to TRUE.
#' @return List of relevant protein families
#' @export
#' @importFrom magrittr %>% 
#' @importFrom utils read.csv write.csv write.table read.table data
#' @examples
#' \dontrun{
#' neighborCatalogOut <- neighborCatalog(imgGenesTrimmed=imgGenesTrimmed, 
#'                                       imgNeighborsTrimmed = imgNeighborsTrimmed, 
#'                                       neighborNumber=neighborNumber, 
#'                                       includeGene=includeGene, 
#'                                       geneName=geneName)
#' }
#' 
neighborCatalog <- function(imgGenesTrimmed = imgGenesTrimmed,
                             imgNeighborsTrimmed = imgNeighborsTrimmed,
                             geneName = geneName,
                             neighborThreshold = neighborThreshold,
                             useInterPro = useInterPro,
                             useHypo = useHypo) { 
  fileDate <- format(Sys.Date(),format="%Y%m%d")
  fileName <- paste(fileDate,"_neighborCatalog_",geneName,sep="")
  fileSummary <- paste(fileName,"_familySummary.txt",sep="")
  fileList <- paste(fileName,"_familyList.txt",sep="")
    ## identifying extant families
  numGenes <- length(unique(imgNeighborsTrimmed$source_gene_oid))
  pfamFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$Pfam, "pfam")==TRUE)
  pfamFams <- pfamFams$Pfam
  uniquePfam <- unique(pfamFams)
  allFams <- pfamFams
  tigrFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$Tigrfam, "TIGR")==TRUE)
  tigrFams <- tigrFams$Tigrfam
  uniqueTigrfam <- unique(tigrFams)
  allFams <- append(allFams, tigrFams)
  if (useInterPro == TRUE) {
      iprFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$InterPro, "IPR")==TRUE)
      iprFams <- iprFams$InterPro
      uniqueIprfam <- unique(iprFams)
      allFams <- append(allFams, iprFams)
  }
  if (useHypo == TRUE) {
      hypoFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$Hypofam, "hypo")==TRUE)
      hypoFams <- hypoFams$Hypofam
      uniqueHypofam <- unique(hypoFams)
      allFams <- append(allFams, hypoFams)
  }
  ## to implement: IMG fams
  ## using the stupid XML traces since the IDs are number-only
  ##imgFams <- imgNeighborsTrimmed %>% dplyr::filter(stringr::str_detect(imgNeighborsTrimmed$IMG.Term, "br")==TRUE)
  ##imgFams <- imgFams$IMG.Term
  ##allFams <- append(allFams, imgFams)
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
  if (useInterPro == TRUE) {
      ## finding unique InterPro fams and dealing with the punctuation
      for (i in 1:length(uniqueIprfam)) {
          tempSplit <- unlist(strsplit(uniqueIprfam[i],"[[:punct:][:space:]]+"))
          tempSplit <- tempSplit[stringr::str_detect(tempSplit, "IPR") ]
          uniqueFams <- append(uniqueFams, tempSplit)
      }
  }
  if (useHypo == TRUE) {
      ## finding unique hypofams, if any
      for (i in 1:length(uniqueHypofam)) {
          ## no weirdo punctuation here!
          tempSplit <- unlist(strsplit(uniqueHypofam[i],"[[:space:]]+"))
          tempSplit <- tempSplit[stringr::str_detect(tempSplit, "hypo") ]
          uniqueFams <- append(uniqueFams, tempSplit)
      }
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
  tempFamList <- list(fams=unlist(uniqueFams), abund=vector(mode = "character", length = uniqueFamNum))
  uDataFams <- data.frame(tempFamList, stringsAsFactors=FALSE)
  trimmedFams <- list()
    ## calculating abundance
    ## and adding families that pass to the list
  for (j in 1:length(uDataFams[,1]))  {
    tempAbund <- grep(uDataFams[j,1], allFams)
    uDataFams[j,2] <- length(tempAbund)
  }
  if (any(is.na(uDataFams))) {
    ## we don't want to include 
    uDataFams <- uDataFams[-which(is.na(uDataFams)),]
  }
  uDataFams$pAbund <- as.numeric(uDataFams$abund) / numGenes
  commonFams <- uDataFams %>% dplyr::filter(.data$pAbund >= neighborThreshold)
  colnames(commonFams) <- c("imgNames", "abund", "pAbund")
  data(famInfoList)
  commonFams <- dplyr::left_join(commonFams, famInfoList)
  commonFams[which(!(commonFams$imgNames %in% famInfoList$imgNames)),4:6] <- ""
  familyList <- commonFams$imgNames
  write.table(commonFams, file=fileSummary, row.names=FALSE, sep="\t", quote=FALSE)
  write.table(uniqueFams, file=fileList, row.names=FALSE, sep="\t", quote=FALSE)
  print("Protein families in neighboring genes catalogued.")
  return(familyList)
}

