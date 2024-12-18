#' Function for preliminary analysis of neighboring genes from IMG
#'
#' This function uses IMG metadatasets for genes of interest and their genomic neighbors to do some QC and identify subfamilies of hypothetical proteins.
#' @param imgGenes What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param imgNeighbors What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param geneSeqs What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param neighborSeqs What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param neighborsContext File from neighborGenerate tying genes of interest to their neighbors. Filename as string ("filename.txt")
#' @param geneName What is the name of your gene? String ("genE")
#' @param efiRepnodes Does the input dataset consist of EFI repnodes?  T/F, defaults to FALSE.
#' @param neighborNumber How many neighbors do you want to look at on each side of the gene? Integer.
#' @param trimShortClusters Should gene clusters with fewer than the minimum neighbor number be removed? T/F value, defaults to FALSE.
#' @param hypoAnalysis Should hypothetical proteins be clustered and analyzed (requires BLAST)? T/F value, defaults to FALSE.
#' @param sysTerm If running hypoAnalysis, what terminal are you using? String value ("wsl", "nix"), defaults to "nix".
#' @param numThreads How many threads should processes use?  Number depends on your processor, defaults to 1 to be safe.
#' @param neighborThreshold What percent of gene clusters should a protein family occur in to be of interest? Number, defaults to 0.05
#' @param clustMethod What clustering method should be used? String value ("tidygraph", "pvclust"), defaults to "tidygraph".
#' @param alphaVal What alpha value cutoff should be used for pvclust? Number from 0-1, defaults to 0.95
#' @param bootStrap How many bootstrap rounds does pvclust get? Integer, defaults to 10.
#' @param pidCutoff Below what percent ID should edges be deleted? Number from 1-100, defaults to 35.
#' @param pepScreen Should subgroups of peptides be identified (annotated or not)?  T/F, defaults to FALSE.
#' @param pepMax Maximum size (in aa) for peptides in pepScreen. Number, defaults to 150.
#' @param matchLength BLAST matches need to above this fraction of of the length of the smaller seq. Number from 0-1, defaults to .65.
#' @param alnClust Should MAFFT alignments be made of members of a hypothetical protein cluster? T/F, defaults to FALSE.
#' @param hmmClust Should HMM models be made for a given hypothetical protein cluster? T/F, defaults to FALSE.
#' @return List with trimmed metadata sets for both genes of interest and neighboring genes (additional files generated en route)
#' @export
#' @importFrom utils read.csv write.csv write.table read.table
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @examples 
#' \dontrun{
#' prepNeighborsOutput <- prepNeighbors(imgGenes="geneFile.txt", 
#'                                      imgNeighbors="neighborsFile.txt", 
#'                                      geneSeqs="geneSeqs.fa", 
#'                                      neighborSeqs="neighborSeqs.fa", 
#'                                      neighborsContext = "context.txt", 
#'                                      geneName="genE",  
#'                                      neighborNumber=10, 
#'                                      sysTerm="wsl")
#' }
#'
prepNeighbors <- function(imgGenes = imgGenes,
                          imgNeighbors = imgNeighbors,
                          geneSeqs = geneSeqs,
                          neighborSeqs = neighborSeqs,
                          neighborsContext = neighborsContext,
                          geneName = geneName,
                          neighborNumber = neighborNumber,
                          efiRepnodes = FALSE,
                          trimShortClusters = FALSE,
                          hypoAnalysis = FALSE,
                          sysTerm = "nix",
                          clustMethod = "tidygraph",
                          numThreads = 1,
                          neighborThreshold = 0.025,
                          alphaVal = 0.95,
                          bootStrap = 10,
                          pidCutoff = 35,
                          pepScreen = FALSE,
                          pepMax = 150,
                          matchLength = 0.65,
                          alnClust = FALSE,
                          hmmClust = FALSE)  {
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
    ## continue with existing data or upload new data
    ## NOTE:  this does not require InterPro info but the next bits will import it if it is there.
    imgCols <- list("gene_oid",
                    "Locus.Tag",
                    "Gene.Product.Name",
                    "Genome.ID",
                    "Genome.Name",
                    "Gene.Symbol",
                    "GenBank.Accession",
                    "Chromosome",
                    "Start.Coord",
                    "End.Coord",
                    "Strand",
                    "DNA.Sequence.Length..bp.",
                    "Amino.Acid.Sequence.Length..aa.",
                    "Transmembrane.Helices",
                    "Signal.Peptides",
                    "Scaffold.ID",
                    "Scaffold.External.Accession",
                    "Scaffold.Name",
                    "Scaffold.GC..",
                    "COG",
                    "Pfam",
                    "Tigrfam",
                    "SMART.ID",
                    "SUPERFam.ID",
                    "CATH.FunFam.ID",
                    "Enzyme",
                    "KO",
                    "IMG.Term")
    ## given that gene_oids and so on are long numbers
    ## we're just gonna kill scientific notation
    options(scipen = 999)
    ## decided it's easier just to explicitly re-import at this stage
    if (typeof(imgGenes) == "character" && typeof(neighborsContext) == "character") {
        ## by default, 
        imgGenesFull <- as.data.frame(read.csv(imgGenes, header=TRUE, sep="\t", stringsAsFactors=FALSE))
        imgGenesData <- imgGenesFull[names(imgGenesFull) %in% imgCols]
        imgGenesData <- imgGenesData %>% dplyr::mutate_all(~ tidyr::replace_na(as.character(.x), "")) 
        ## if interpro exists add it
        if (any(grepl("InterPro", colnames(imgGenesFull)))) {
            imgGenesData$InterPro <- imgGenesFull$InterPro
        } else {
            imgGenesData$InterPro <- ""
        }  
        imgNeighborsContext <-as.data.frame(read.csv(neighborsContext, header=TRUE, sep="\t", stringsAsFactors=FALSE))
    } else {
        print("Please load an imgGenesData textfile (can be a raw IMG file) and an imgNeighborsContext file.")
        return(0)
    }
    if(typeof(imgNeighbors) == "character")  {
        imgNeighborsFull <- as.data.frame(read.csv(imgNeighbors, header=TRUE, sep="\t", stringsAsFactors=FALSE))
        imgNeighborsData <- imgNeighborsFull[names(imgNeighborsFull) %in% imgCols]
        imgNeighborsData <- imgNeighborsData %>% dplyr::mutate_all(~ tidyr::replace_na(as.character(.x), ""))
        ## adding InterPro and Hypofam columns if we don't have them
        ## if we don't run hypothetical protein or peptide detection the latter will just be blank
        ## if we don't run incorpIprScan and/or the IMG doesn't start more consistently adding InterPro values the former will be blank
        if (any(grepl("InterPro", colnames(imgNeighborsFull)))) {
            imgNeighborsData$InterPro <- imgNeighborsFull$InterPro
        } else {
            imgNeighborsData$InterPro <- ""
        }
        if (any(grepl("Hypofam", colnames(imgNeighborsFull)))) {
            imgNeighborsData$Hypofam <- imgNeighborsFull$Hypofam
        } else {
            imgNeighborsData$Hypofam <- ""
        }
    } else {
        print("Please load an imgNeighborsTemp textfile (can be a raw IMG metadata textfile).")
        return(0)
    }
    ## let's take a minute and get rid of quotation marks in the metadata that can break things later
    ## Thank you annotations for streptomycin 3''-kinase aka streptomycin 3"-kinase
    ## Why must you make everything terrible
    imgGenesData <- imgGenesData %>% dplyr::mutate_all(stringr::str_replace_all, "\"", "\'\'")
    imgNeighborsData <- imgNeighborsData  %>% dplyr::mutate_all(stringr::str_replace_all, "\"", "\'\'")
    ## make sure the metadata is un-annoying, i.e. sorted
    imgNeighborsData$gene_oid <- as.numeric(imgNeighborsData$gene_oid)
    ndVector <- as.vector(imgNeighborsData$gene_oid)
    rowidx <- order(ndVector,decreasing=FALSE,na.last=TRUE)
    imgNeighborsData <- imgNeighborsData[rowidx,,drop=TRUE]
    ## sorting the context data
    imgNeighborsContext$gene_oid <- as.numeric(imgNeighborsContext$gene_oid)
    imgNeighborsContext$source_gene_oid <- as.numeric(imgNeighborsContext$source_gene_oid)
    imgNeighborsContext$source_scaffold_id <- as.numeric(imgNeighborsContext$source_scaffold_id)
    ncVector <- as.vector(imgNeighborsContext$gene_oid)
    rowidx2 <- order(ncVector,decreasing=FALSE,na.last=TRUE)
    imgNeighborsContext <- imgNeighborsContext[rowidx2,,drop=TRUE]
    imgNeighborsContext$gene_oid <- as.list(imgNeighborsContext$gene_oid)
    ## sorting the gene data
    imgGenesData$gene_oid <- as.numeric(imgGenesData$gene_oid)
    gdVector <- as.vector(imgGenesData$gene_oid)
    rowidx3 <- order(gdVector,decreasing=FALSE,na.last=TRUE)
    imgGenesData <- imgGenesData[rowidx3,,drop=TRUE]  
    ## import amino acid sequences
    imgNeighborSeqs <- seqinr::read.fasta(file=neighborSeqs, seqtyp="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
    imgGeneSeqs <- seqinr::read.fasta(file=geneSeqs, seqtyp="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
    ## if analyzing hypothetical proteins, we follow both subroutines; otherwise, just the basic trimming function
    if (hypoAnalysis == TRUE) {
        imgNeighborsData <- neighborHypothetical(imgGenesData = imgGenesData,
                                                 imgNeighborsData = imgNeighborsData,
                                                 imgNeighborSeqs=imgNeighborSeqs,
                                                 geneName = geneName,
                                                 alphaVal = alphaVal,
                                                 bootStrap = bootStrap,
                                                 sysTerm = sysTerm,
                                                 numThreads = numThreads,
                                                 clustMethod = clustMethod,
                                                 pidCutoff = pidCutoff,
                                                 screenPep = FALSE,
                                                 alnClust = alnClust,
                                                 hmmClust = hmmClust,
                                                 matchLength = matchLength)
    }
    ## output: imgNeighborsData (updated), hypoClusters.txt, hypoClusterX.fa (per cluster), hypoClusterXaln.fa (per cluster), hypoSettings.txt
    if (pepScreen == TRUE)  {
        imgNeighborsData <- neighborHypothetical(imgGenesData = imgGenesData,
                                                 imgNeighborsData = imgNeighborsData,
                                                 imgNeighborSeqs=imgNeighborSeqs,
                                                 geneName = geneName,
                                                 alphaVal = alphaVal,
                                                 bootStrap = bootStrap,
                                                 sysTerm = sysTerm,
                                                 numThreads = numThreads,
                                                 clustMethod = clustMethod,
                                                 pidCutoff = pidCutoff,
                                                 screenPep = TRUE,
                                                 pepMax = pepMax,
                                                 alnClust = alnClust,
                                                 hmmClust = hmmClust,
                                                 matchLength = matchLength)
    }
    ## output: imgNeighborsData (updated), pepClusters.txt, pepClusterX.fa (per cluster), pepClusterXaln.fa (per cluster), pepSettings.txt
    ## note: this will always trim scaffold mismatches
    ## if trimShortClusters is false, it will stop here.
    ## if true, it will, well, trim short clusters too.
    neighborTrimOutput <- neighborTrim(imgNeighborsData = imgNeighborsData,
                                            imgGenesData = imgGenesData,
                                            imgNeighborsContext = imgNeighborsContext,
                                            trimShortClusters = trimShortClusters,
                                            neighborNumber = neighborNumber,
                                            geneName = geneName,
                                            imgNeighborSeqs=imgNeighborSeqs,
                                            imgGeneSeqs=imgGeneSeqs)
    ## output: imgNeighborsTrimmed (includes source_gene_oid column), imgGenesTrimmed
    imgGenesTrimmed <- neighborTrimOutput$imgGenesTrimmed
    imgNeighborsTrimmed <- neighborTrimOutput$imgNeighborsTrimmed
    return(list(imgNeighborsTrimmed=imgNeighborsTrimmed,imgGenesTrimmed=imgGenesTrimmed))
}
