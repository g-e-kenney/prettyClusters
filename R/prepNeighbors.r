#' Function for preliminary analysis of neighboring genes from IMG
#'
#' This function uses IMG metadatasets for genes of interest and their genomic neighbors to do some QC and identify subfamilies of hypothetical proteins.
#' @param imgGenesTemp What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param imgNeighborsTemp What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param geneSeqs What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param neighborSeqs What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param neighborsContext File from neighborGenerate tying genes of interest to their neighbors. Filename as string ("filename.txt")
#' @param geneName What is the name of your gene? String ("genE")
#' @param efiRepnodes Does the input dataset consist of EFI repnodes?  T/F, defaults to FALSE.
#' @param neighborNumber How many neighbors do you want to look at on each side of the gene? Integer.
#' @param neighborThreshold What percent of gene clusters should a protein family occur in to be of interest? Number, defaults to 0.05
#' @param hypoAnalysis Should hypothetical proteins be clustered and analyzeD? T/F value, defaults to TRUE.
#' @param clustMethod What clustering method should be used? String value ("tidygraph", "pvclust"), defaults to "tidygraph".
#' @param sysTerm What terminal are you using? String value ("wsl", "nix"), no default.
#' @param numThreads How many threads should processes use?  Number depends on your processor, defaults to 1 to be safe.
#' @param alphaVal What alpha value cutoff should be used for pvclust? Number from 0-1, defaults to 0.95
#' @param bootStrap How many bootstrap rounds does pvclust get? Integer, defaults to 10.
#' @param pidCutoff Below what percent ID should edges be deleted? Number from 1-100, defaults to 35.
#' @param trimShortClusters Should gene clusters with fewer than the minimum gene number be visualized? T/F value.
#' @return List with trimmed metadata sets for both genes of interest and neighboring genes (additional files generated en route)
#' @export
#' @examples 
#' prepNeighborsOutput <- prepNeighbors(imgGenes="geneFile.txt", imgNeighbors="neighborsFile.txt", geneSeqs="geneSeqs.fa", neighborSeqs="neighborSeqs.fa", neighborsContext = "context.txt", geneName="genE",  neighborNumber=10, sysTerm="wsl")
#'
prepNeighbors <- function(imgGenes = imgGenes, imgNeighbors = imgNeighbors, geneSeqs = geneSeqs, neighborSeqs = neighborSeqs, neighborsContext = neighborsContext, geneName = geneName, neighborNumber = neighborNumber, sysTerm = sysTerm, efiRepnodes = FALSE, neighborThreshold = 0.025, hypoAnalysis = TRUE, clustMethod = "tidygraph", numThreads = 1, alphaVal = 0.95, bootStrap = 10, pidCutoff = 35, trimShortClusters = TRUE)  {
                                        # starting stuff
  if(exists(x="imgGenes") == FALSE | exists(x="imgNeighbors") == FALSE | exists(x="geneSeqs") == FALSE | exists(x="neighborSeqs") == FALSE | exists(x="neighborNumber") == FALSE | exists(x="geneName") == FALSE | exists(x="sysTerm") == FALSE) {
    print("Missing a required term (gene and neighbor metadata and sequence files, context file from generateNeighbors, name of gene of interest, number of neighbors, and system terminal)")
    return(0)
  } 
  fileDate <- format(Sys.Date(),format="%Y%m%d")
                                                # a purely cosmetic option to track whether input data consists of EFI repnodes or not
    ## i feel like we don't wanna risk bad sequences becomming repnodes
    ## thus running this pre-EFI-EST
  if (efiRepnodes == TRUE) {
    geneName <- paste(geneName, "_EFI",sep="")
  } else {
    geneName <- geneName
  }
                                        #        continue with existing data or upload new data
  imgCols <- list("gene_oid", "Locus.Tag", "Gene.Product.Name", "Genome.ID", "Genome.Name", "Gene.Symbol", "GenBank.Accession", "Chromosome", "Start.Coord", "End.Coord", "Strand", "DNA.Sequence.Length..bp.", "Amino.Acid.Sequence.Length..aa.", "Transmembrane.Helices", "Signal.Peptides", "Scaffold.ID", "Scaffold.External.Accession", "Scaffold.Name", "Scaffold.GC..", "COG", "Pfam", "Tigrfam", "SMART.ID", "SUPERFam.ID", "CATH.FunFam.ID", "Enzyme", "KO", "IMG.Term")
# decided it's easier just to explicitly re-import at this stage
  if (typeof(imgGenes) == "character" && typeof(neighborsContext) == "character") {
## by default, 
    imgGenesData <-as.data.frame(read.csv(imgGenes, header=TRUE, sep="\t" )[,1:45])
    imgGenesData <- imgGenesData %>% dplyr::mutate_all(~ replace_na(.x, ""))
    imgGenesData <- imgGenesData[names(imgGenesData) %in% imgCols]
    imgNeighborsContext <-as.data.frame(read.csv(neighborsContext, header=TRUE, sep="\t" , stringsAsFactors=FALSE))
  } else {
    print("Please load an imgGenesData textfile (can be a raw IMG file) and an imgNeighborsContext file.")
    return(0)
  }
  if(typeof(imgNeighbors) == "character")  {
    imgNeighborsFull <- as.data.frame(read.csv(imgNeighbors, header=TRUE, sep="\t", stringsAsFactors=FALSE)[,1:45])
    imgNeighborsData <- imgNeighborsFull %>% dplyr::mutate_all(~ replace_na(.x, ""))
    imgNeighborsData <- imgNeighborsFull[names(imgNeighborsFull) %in% imgCols]
  } else {
    print("Please load an imgNeighborsTemp textfile (can be a raw IMG metadata textfile).")
    return(0)
  }
# make sure the neighbors metadata is un-annoying, and loading more sequences
  ndVector <- imgNeighborsData$gene_oid
  ndVector <- as.vector(ndVector)
  rowidx <- order(ndVector,decreasing=FALSE,na.last=TRUE)
  imgNeighborsData <- imgNeighborsData[rowidx,,drop=TRUE]
  imgNeighborsContext$gene_oid <- as.numeric(imgNeighborsContext$gene_oid)
  imgNeighborsContext$source_gene_oid <- as.numeric(imgNeighborsContext$source_gene_oid)
  imgNeighborsContext$source_scaffold_id <- as.numeric(imgNeighborsContext$source_scaffold_id)
  ncVector <-imgNeighborsContext$gene_oid
  ncVector <- as.vector(ncVector)
  rowidx2 <- order(ncVector,decreasing=FALSE,na.last=TRUE)
  imgNeighborsContext <- imgNeighborsContext[rowidx2,,drop=TRUE]
  imgNeighborsContext$gene_oid <- as.list(imgNeighborsContext$gene_oid)
  imgNeighborSeqs <- seqinr::read.fasta(file=neighborSeqs, seqtyp="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
  imgGeneSeqs <- seqinr::read.fasta(file=geneSeqs, seqtyp="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
  ## if analyzing hypothetical proteins, we follow both subroutines; otherwise, just the basic trimming function
  if (hypoAnalysis == TRUE) {
    imgNeighborsData <- neighborHypothetical(imgGenesData = imgGenesData, imgNeighborsData = imgNeighborsData, imgNeighborSeqs=imgNeighborSeqs, geneName = geneName, alphaVal = alphaVal, bootStrap = bootStrap, sysTerm = sysTerm, numThreads = numThreads, clustMethod = clustMethod, pidCutoff = pidCutoff)
	  	  ## output: imgNeighborsData (updated), hypoClusters.txt, hypoClusterX.fa (per cluster), hypoClusterXaln.fa (per cluster), hypoSettings.txt
    neighborTrimOutput <- neighborTrim(imgNeighborsData = imgNeighborsData, imgGenesData = imgGenesData, imgNeighborsContext = imgNeighborsContext, trimShortClusters = trimShortClusters, neighborNumber = neighborNumber, geneName = geneName, imgNeighborSeqs=imgNeighborSeqs, imgGeneSeqs=imgGeneSeqs)
    ## output: imgNeighborsTrimmed (includes source_gene_oid column), imgGenesTrimmed
  } else {
    neighborTrimOutput <- neighborTrim(imgNeighborsData = imgNeighborsData, imgGenesData = imgGenesData, imgNeighborsContext = imgNeighborsContext, trimShortClusters = trimShortClusters, neighborNumber = neighborNumber, geneName = geneName, imgNeighborSeqs=imgNeighborSeqs, imgGeneSeqs=imgGeneSeqs)
      ## output: imgNeighborsTrimmed (includes source_gene_oid column), imgGenesTrimmed
  }
  imgGenesTrimmed <- neighborTrimOutput$imgGenesTrimmed
  imgNeighborsTrimmed <- neighborTrimOutput$imgNeighborsTrimmed
  return(list(imgNeighborsTrimmed=imgNeighborsTrimmed,imgGenesTrimmed=imgGenesTrimmed))
}