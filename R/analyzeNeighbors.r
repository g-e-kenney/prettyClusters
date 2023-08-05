#' Umbrella function that identifies sets of related genome neighborhoods in IMG-derived data
#'
#' This function encompasses several subfunctions that catalog and quantify protein families in the neighborhood of genes of interest and use that data to identify genomic neighborhoods that are themselves similar.  
#' @param imgGenes What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param imgNeighbors What is the file with the metadata for neighbors of your gene of interest? Filename as string ("filename.txt")
#' @param efiRepnodes Are values here post-EFI-EST? Boolean, defaults to FALSE.
#' @param neighborThreshold In what percentage of neighborhoods must a protein family show up to be of interest?  Number, defaults to 0.025
#' @param geneName What is the name of your gene? Gene name as string ("genE")
#' @param useInterPro Should InterPro families be used in neighborhood analyses? Boolean, defaults to FALSE
#' @param useHypo Should hypothetical protein families from prepNeighbors be used in analyses?  Boolean, defaults to FALSE
#' @param autoClust Should clusters be automatically identified? Boolean, defaults to FALSE.
#' @param clustMethod What method should be used to identify these clusters? String ("tidygraph" or "pvclust")
#' @param alphaVal What alpha value cutoff should be used for pvclust?  Number, defaults to 0.95
#' @param bootStrap How many bootstrap rounds for pvclust? Integer, defaults to 10
#' @param tgCutoff What sort of edge similarity should be kept for tidygraph? Number (0-1), defaults to 0.65.
#'
#' @return Updated metadata and misc. figures and files en route
#' @export
#'
#' @importFrom magrittr %>% 
#' @importFrom utils read.csv write.csv write.table read.table
#'
#' @examples
#' \dontrun{
#' analyzeNeighborsOutput <- analyzeNeighbors(imgGenes="repnodeGenes.txt", 
#'                                            imgNeighbors = "repnodeNeighbors.txt", 
#'                                            geneName = "genE", 
#'                                            tgCutoff = 0.6, 
#'                                            efiRepnodes = FALSE) 
#' }
#'
analyzeNeighbors <- function(imgGenes = imgGenes,
                             imgNeighbors = imgNeighbors,
                             efiRepnodes = FALSE,
                             neighborThreshold = 0.025,
                             geneName = geneName,
                             useInterPro = FALSE,
                             useHypo = FALSE,
                             autoClust = FALSE,
                             clustMethod = "tidygraph",
                             alphaVal = 0.95,
                             bootStrap= 10,
                             tgCutoff = 0.6)  {
    if(exists(x="imgGenes") == FALSE | exists(x="imgNeighbors") == FALSE | exists(x="geneName") == FALSE) {
        print("Missing a required term (gene and neighbor metadata files or name of gene of interest)")
        return(0)
    }
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    imgCols <- list("gene_oid",
                    "source_gene_oid",
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
    ## now let's import data
    if(typeof(imgGenes) == "character") {
        imgGenesTemp <-read.csv(imgGenes, header=TRUE, sep="\t", stringsAsFactors=FALSE )
    } else {
        imgGenesTemp <- imgGenes
    }
    if(typeof(imgNeighbors) == "character")  {
        imgNeighborsTemp <- read.csv(imgNeighbors, header=TRUE, sep="\t", stringsAsFactors=FALSE )
    } else {
        imgNeighborsTemp <- imgNeighbors
    }
    ## some cleanup for the input files, JUST IN CASE
    ## even though they should be fine if coming out of prepNeighbors or repnodeTrim
    ## nevertheless, adds the two columns that may not be in standard IMG metadata files
    imgNeighborsTemp <- imgNeighborsTemp %>% dplyr::mutate_all(~ tidyr::replace_na(as.character(.x), ""))
    imgNeighborsTrimmed <- imgNeighborsTemp[names(imgNeighborsTemp) %in% imgCols]
    ## add in source_gene_oid
    imgNeighborsTrimmed$source_gene_oid <- imgNeighborsTemp$source_gene_oid
    imgNeighborsTrimmed$source_scaffold_id <- imgNeighborsTemp$source_scaffold_id
    ## if interpro exists add it
    if (any(grepl("InterPro", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$InterPro <- imgNeighborsTemp$InterPro
    } else {
        imgNeighborsTrimmed$InterPro <- ""
    }
    if (any(grepl("Hypofam", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$Hypofam <- imgNeighborsTemp$Hypofam
    } else {
        imgNeighborsTrimmed$Hypofam <- ""
    }     
    imgGenesTemp <- imgGenesTemp %>% dplyr::mutate_all(~ tidyr::replace_na(as.character(.x), ""))
    imgGenesTrimmed <- imgGenesTemp[names(imgGenesTemp) %in% imgCols]
    ## if interpro exists add it
    if (any(grepl("InterPro", colnames(imgGenesTemp)))) {
        imgGenesTrimmed$InterPro <- imgGenesTemp$InterPro
    } else {
        imgGenesTrimmed$InterPro <- ""
    }
    if (any(grepl("Hypofam", colnames(imgGenesTemp)))) {
        imgGenesTrimmed$Hypofam <- imgGenesTemp$Hypofam
    } else {
        imgGenesTrimmed$Hypofam <- ""
    } 
    ## a purely cosmetic option to track whether input data consists of EFI repnodes or not
    ## will help for keeping track of ten zillion output files
    ## and while we're at it we'll keep 
    if (efiRepnodes == TRUE) {
        coreGeneName <- geneName
        geneName <- paste(geneName, "_repnodes",sep="")
        if (any(grepl("efi_oid", colnames(imgGenesTemp)))) {
            imgGenesTrimmed$efi_oid <- imgGenesTemp$efi_oid
        }
        if (any(grepl("efi_oid", colnames(imgNeighborsTemp)))) {
            imgNeighborsTrimmed$efi_oid <- imgNeighborsTemp$efi_oid
        }
    } else {
        geneName <- geneName
        coreGeneName <- geneName
        if (any(grepl("efi_oid", colnames(imgGenesTemp)))) {
            imgGenesTrimmed$efi_oid <- imgGenesTemp$efi_oid
        }
        if (any(grepl("isRepnode", colnames(imgGenesTemp)))) {
            imgGenesTrimmed$isRepnode <- imgGenesTemp$isRepnode
        }
        if (any(grepl("repnodeIs", colnames(imgGenesTemp)))) {
            imgGenesTrimmed$repnodeIs <- imgGenesTemp$repnodeIs
        }    
        if (any(grepl("efi_oid", colnames(imgNeighborsTemp)))) {
            imgNeighborsTrimmed$efi_oid <- imgNeighborsTemp$efi_oid
        }
        if (any(grepl("isRepnode", colnames(imgNeighborsTemp)))) {
            imgNeighborsTrimmed$isRepnode <- imgNeighborsTemp$isRepnode
        }
    }
                                        # process imgNeighborsTrimmed here: if genes-of-interest are not present in it, add in imgGeneMetadata. trim extra columns, adjust formatting if needed.
    neighborCatalogOut <- neighborCatalog(imgNeighborsTrimmed = imgNeighborsTrimmed,
                                          imgGenesTrimmed=imgGenesTrimmed,
                                          geneName = geneName,
                                          neighborThreshold = neighborThreshold,
                                          useInterPro = useInterPro,
                                          useHypo = useHypo)
                                        # output: familyList (list of neighboring families), familyAbundance, familyAbundance.txt)
    neighborHereOut <- neighborHere(imgNeighborsTrimmed = imgNeighborsTrimmed,
                                    familyList = neighborCatalogOut,
                                    geneName = geneName)
                                        # input: neighbor metadata table, family abundance file, gene of interest, fileDate)
                                        # output: neighborBinary, neighborBinary (.txt)
    neighborMatrixOut <- neighborMatrix(imgGenesTrimmed = imgGenesTrimmed,
                                        neighborBinary = neighborHereOut,
                                        familyList = neighborCatalogOut,
                                        geneName = geneName)
                                        # input: imgGenesTrimmed, the neighbors binary file, the gene of interest, the file date
                                        # output: neighborMatrixData
    neighborClustersOut <- neighborClusters(imgGenesTrimmed = imgGenesTrimmed,
                                            imgNeighborsTrimmed = imgNeighborsTrimmed,
                                            geneName = geneName,
                                            neighborMatrixData = neighborMatrixOut,
                                            autoClust = autoClust,
                                            clustMethod = clustMethod,
                                            alphaVal = alphaVal,
                                            bootStrap = bootStrap,
                                            coreGeneName = coreGeneName,
                                            tgCutoff=tgCutoff)
                                        # input: matrix of genes/families, metadata file for genes of interest, gene of interest, auto-generated date for file names
                                        # output: neighborMatrixClustered, neighborMatrixClustered (.txt), clusterHeatmap (.pdf, png), imgGenesClustered (imgGenes + clusters)  & .txt, efiNodeDataClustered.txt (efiNodeData + clusters) & .txt
                                        # put here - trimming imgNeighborsTrimmed to the right input format
    print ("Genome neighborhood clusters & associated diagrams generated.")
    return(neighborClustersOut)
}    


