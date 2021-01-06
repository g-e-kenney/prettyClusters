#' Umbrella function that identifies sets of related genome neighborhoods in IMG-derived data
#'
#' This function encompasses several subfunctions that catalog and quantify protein families in the neighborhood of genes of interest and use that data to identify genomic neighborhoods that are themselves similar.  
#' @param imgGenesTrimmed What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param imgNeighborsTrimmed What is the file with the metadata for neighbors of your gene of interest? Filename as string ("filename.txt")
#' @param efiRepnodes Are values here post-EFI-EST? T/F value, defaults to FALSE.
#' @param neighborThreshold In what percentage of neighborhoods must a protein family show up to be of interest?  Number, defaults to 0.025
#' @param geneName What is the name of your gene? Gene name as string ("genE")
#' @param autoClust Should clusters be automatically identified? T/F, defaults to TRUE.
#' @param clustMethod What method should be used to identify these clusters? String ("tidygraph" or "pvclust")
#' @param alphaVal What alpha value cutoff should be used for pvclust?  Number, defaults to 0.95
#' @param bootStrap How many bootstrap rounds for pvclust? Integer, defaults to 10
#' @return Updated metadata and misc. figures and files en route
#' @export
#' @examples
#' analyzeNeighborsOutput <- analyzeNeighbors <- function(imgGenes="repnodeGenes.txt", imgNeighbors = "repnodeNeighbors.txt", geneName = "genE") 
#'
analyzeNeighbors <- function(imgGenes = imgGenes, imgNeighbors = imgNeighbors, efiRepnodes = FALSE, neighborThreshold = 0.025, geneName = geneName, autoClust = TRUE, clustMethod = "tidygraph", alphaVal = 0.95, bootStrap= 10)  {

  if(exists(x="imgGenes") == FALSE | exists(x="imgNeighbors") == FALSE | exists(x="geneName") == FALSE) {
    print("Missing a required term (gene and neighbor metadata files or name of gene of interest)")
    return(0)
  }
  fileDate <- format(Sys.Date(),format="%Y%m%d")
  if(typeof(imgGenes) == "character") {
    imgGenesTrimmed <-read.csv(imgGenes, header=TRUE, sep="\t", stringsAsFactors=FALSE )
  } else {
    imgGenesTrimmed <- imgGenes
  }
  if(typeof(imgNeighbors) == "character")  {
    imgNeighborsTrimmed <- read.csv(imgNeighbors, header=TRUE, sep="\t", stringsAsFactors=FALSE )
  } else {
    imgNeighborsTrimmed <- imgNeighbors
  }
                                        # a purely cosmetic option to track whether input data consists of EFI repnodes or not
                                        # will help for keeping track of ten zillion output files
  if (efiRepnodes == TRUE) {
    coreGeneName <- geneName
    geneName <- paste(geneName, "_repnodes",sep="")
  } else {
    geneName <- geneName
    coreGeneName <- geneName
  }
                                        # process imgNeighborsTrimmed here: if genes-of-interest are not present in it, add in imgGeneMetadata. trim extra columns, adjust formatting if needed.
  neighborCatalogOut <- neighborCatalog(imgNeighborsTrimmed = imgNeighborsTrimmed, imgGenesTrimmed=imgGenesTrimmed, geneName = geneName, neighborThreshold = neighborThreshold, coreGeneName = coreGeneName)
                                        # output: familyList (list of neighboring families), familyAbundance, familyAbundance.txt)
  neighborHereOut <- neighborHere(imgNeighborsTrimmed = imgNeighborsTrimmed, familyList = neighborCatalogOut, geneName = geneName, coreGeneName = coreGeneName)
                                        # input: neighbor metadata table, family abundance file, gene of interest, fileDate)
                                        # output: neighborBinary, neighborBinary (.txt)
  neighborMatrixOut <- neighborMatrix(imgGenesTrimmed = imgGenesTrimmed, neighborBinary = neighborHereOut, familyList = neighborCatalogOut, geneName = geneName, coreGeneName = coreGeneName)
                                        # input: imgGenesTrimmed, the neighbors binary file, the gene of interest, the file date
                                        # output: neighborMatrixData
  neighborClustersOut <- neighborClusters(imgGenesTrimmed = imgGenesTrimmed, imgNeighborsTrimmed = imgNeighborsTrimmed, geneName = geneName, neighborMatrixData = neighborMatrixOut, autoClust = autoClust, clustMethod = clustMethod, alphaVal = alphaVal, bootStrap = bootStrap, coreGeneName = coreGeneName)
                                        # input: matrix of genes/families, metadata file for genes of interest, gene of interest, auto-generated date for file names
                                        # output: neighborMatrixClustered, neighborMatrixClustered (.txt), clusterHeatmap (.pdf, png), imgGenesClustered (imgGenes + clusters)  & .txt, efiNodeDataClustered.txt (efiNodeData + clusters) & .txt
                                        # put here - trimming imgNeighborsTrimmed to the right input format
  print ("Genome neighborhood clusters & associated diagrams generated.")
  return(neighborClustersOut)
}    

    
