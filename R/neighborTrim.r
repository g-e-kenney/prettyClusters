#' Subfunction to clean up IMG metadata files
#'
#' Given a desired genome neighborhood and lists of genes of interest and their neighbors, this function removes any genes not on the same scaffold as their ostensible neighbors and trims truncated clusters, removing them from further analysis.
#' @param imgGenesData Data frame from hypoProteins or neighborPrep with metadata for your gene of interest? 
#' @param imgNeighborsData Data frame from hypoProteins or neighborPrep with metadata for neighbors of your gene of interest.
#' @param imgGeneSeqs Data frame from hypoProteins or neighborPrep with sequence information for your gene of interest.
#' @param imgNeighborSeqs Data frame from hypoProteins or neighborPrep with sequence information for neighbors of your gene of interest.
#' @param imgNeighborsContext Data frame from hypoProteins or neighborPrep tying genes of interest to their neighbors.
#' @param geneName Character string with your gene of interest's name.
#' @param neighborNumber How many neighbors do you want to look at on each side of the gene? Integer greater than 1.
#' @param trimShortClusters Should gene clusters with fewer than the minimum gene number be visualized? Boolean, defaults to TRUE in parent..
#' @return List with trimmed metadata sets for both genes of interest and neighboring genes (additional files generated en route)
#' @export
#' @importFrom utils read.csv write.csv write.table read.table
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @examples 
#' \dontrun{
#' neighborTrimOutput <- neighborTrim(imgNeighborsData = imgNeighborsData, 
#'                                    imgGenesData = imgGenesData, 
#'                                    imgNeighborsContext = imgNeighborsContext, 
#'                                    trimShortClusters = TRUE, 
#'                                    neighborNumber = 10, 
#'                                    geneName = "genE", 
#'                                    imgNeighborSeqs=imgNeighborSeqs, 
#'                                    imgGeneSeqs=imgGeneSeqs)
#' }
#'
neighborTrim <- function(imgNeighborsData = imgNeighborsData,
                         imgGenesData = imgGenesData,
                         imgNeighborsContext = imgNeighborsContext,
                         trimShortClusters = trimShortClusters,
                         neighborNumber = neighborNumber,
                         geneName = geneName,
                         imgNeighborSeqs=imgNeighborSeqs,
                         imgGeneSeqs=imgGeneSeqs) {
                                        # let's make things
    trimmedGenes <- list()
    notInMetadata <- list(gene_oid = 0)
    wrongScaffold <- list(gene_oid = 0)
    smallNeighborhood <- list(gene_oid = 0)
    nMatch <- list()
    nSeqTable <- data.frame()
    gSeqTable <- data.frame()
    nReList <- list()    
    gReList <- list()   
    ## filenaming nonsense
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    fileName <- paste(fileDate,"_","neighborTrim_",geneName,sep="")
    fileNameiNT <- paste(fileName,"_neighbors.txt",sep="")
    fileNameiNCT <- paste(fileName,"_neighborsContext.txt",sep="")
    fileNameiGT <- paste(fileName,"_genes.txt",sep="")
    fileNameiGST <- paste(fileName,"_geneSeqs.fa",sep="")
    fileNameiNST <- paste(fileName,"_neighborSeqs.fa",sep="")
    ## and moving on
    print("Performing quality control on neighbors.")
    uniqueNeighbors <- unique(imgNeighborsData$gene_oid)
    numNeighbors <- length(imgNeighborsData$gene_oid)
    ## deals with the occasional disappearance of gene_oids from IMG metadata
    ## yeah i dunno why this happens
    ## also handles duplicate neighbor entries (e.g. from when you have two GoIs (or two fragments of one GoI) in one cluster)
    imgNeighborsContext <- data.table::as.data.table(imgNeighborsContext)
    if(length(imgNeighborsData$gene_oid) != length(imgNeighborsContext$gene_oid)) {
        print("You have different numbers of genes in your current neighbor metadata and your original neighbor list; adjusting the latter to match the former.") 
    }
    imgNeighborsTrimmed <- dplyr::filter(imgNeighborsData, .data$gene_oid %in% uniqueNeighbors)
    imgNeighborsContextTrimmed <- dplyr::filter(imgNeighborsContext, .data$gene_oid %in% uniqueNeighbors)
    imgNeighborSeqsTrimmed <- imgNeighborSeqs[which(names(imgNeighborSeqs) %in% uniqueNeighbors)]
    imgGeneSeqsTrimmed <- imgGeneSeqs
    imgGenesTrimmed <- imgGenesData
                                        # fasta to data frame because the seqinR format is a wonky list
    nSeqNames <- names(imgNeighborSeqsTrimmed)
    for (i in 1:length(nSeqNames)) {
        nSeqTable[i,1] <- nSeqNames[i]
        nSeqTable[i,2] <- imgNeighborSeqsTrimmed[[i]]
    }
    colnames(nSeqTable) <- c("gene_oid", "sequence")
    gSeqNnames <- names(imgGeneSeqsTrimmed)
    for (i in 1:length(gSeqNnames)) {
        gSeqTable[i,1] <- gSeqNnames[i]
        gSeqTable[i,2] <- imgGeneSeqsTrimmed[[i]]
    }
    colnames(gSeqTable) <- c("gene_oid", "sequence")
                                        # check if neighbor scaffold equals source scaffold; if not, removes those genes
                                        # because we don't know they are actually neighbors
    ## note that this removes the need for scaffold checks later on (e.g. in prettyClusterDiagrams)
    uniqueNeighborsTrimmed <- unique(imgNeighborsTrimmed$gene_oid)
    numNeighborsTrimmed <- length(imgNeighborsTrimmed$gene_oid)
                                        #imgNeighborsContextTrimmed$gene_oid <- as.character(imgNeighborsContextTrimmed$gene_oid)
    multiHit <- list(gene_oid=0)
                                        # just in case trimming screwed up the order or something? shouldn't but...
                                        # this is pretty much deprecated but 
    ##for (i in 1:numNeighborsTrimmed) {
    ##    geneIdxCon <- grep(imgNeighborsTrimmed$gene_oid[i], imgNeighborsContextTrimmed$gene_oid)
    ##  geneIdxCon <- as.numeric(unlist(geneIdxCon))
    ##if (length(geneIdxCon)>=2) {append(multiHit, geneIdxCon)}
    ##if (as.character(imgNeighborsTrimmed$Scaffold.ID[i]) != as.character(imgNeighborsContextTrimmed$source_scaffold_id[geneIdxCon])) {
    ##      wrongScaffold <- append(wrongScaffold, as.character(imgNeighborsTrimmed$gene_oid[i]))
    ##  } else {
    ##  wrongScaffold <- wrongScaffold
    ##    }
    ##}
    imgNeighborsContextTrimmed$gene_oid <- as.character(imgNeighborsContextTrimmed$gene_oid)    
    imgNeighborsTrimmed$gene_oid <- as.character(imgNeighborsTrimmed$gene_oid)
    imgNeighborsTrimmed <- imgNeighborsTrimmed %>% dplyr::left_join(imgNeighborsContextTrimmed, by="gene_oid")
    wrongScaffold <- imgNeighborsTrimmed %>% dplyr::filter(.data$Scaffold.ID != .data$source_scaffold_id)
    wrongScaffold <- wrongScaffold$gene_oid
    imgNeighborsTrimmed <-imgNeighborsTrimmed %>% dplyr::filter(!.data$gene_oid %in% wrongScaffold)
    imgNeighborsContextTrimmed <- imgNeighborsContextTrimmed %>% dplyr::filter(!.data$gene_oid %in% wrongScaffold)
    nSeqTable <- nSeqTable %>% dplyr::filter(!.data$gene_oid %in% wrongScaffold)
    print("Scaffold mismatches trimmed.")
                                        # trimming gene clusters that are truncated
    ## it's possible to not remove the short gene clusters, and appropriate in some circumstances
    ## but keeping them can result in under-weighting of things a little further out in larger gene clusters, when doing cluster analysis
    ## counts the number of hits for source genes;
    ## if under 2*neighborNumber + 1 (i.e. your window on each side of your GoI), add source_gene_oid to remove list
    if (trimShortClusters == TRUE) {
        smallNeighborhood <- list()
        smallNeighbors <- list()
        ## let's get only unique GoI listings
        sourceIDs <- unique(imgNeighborsContextTrimmed$source_gene_oid)
        ## let's define our minimum neighborhood size (and include our GoI in it)
        minNeighbors <- neighborNumber * 2 + 1
        for (i in 1:length(sourceIDs)) {
            ## this originally had the as.character as unSeq - why?
            ## smallMatch <- imgNeighborsContextTrimmed[which(imgNeighborsContextTrimmed$source_gene_oid == as.character(sourceIDs[i])),]
            ## because we can get double the gene count if our GoI is broken in two, let's 
            smallIndex <- grep(as.character(sourceIDs[i]), imgNeighborsContextTrimmed$source_gene_oid)
            ## if we have fewer than our minNeighbor number, add it to the list of too-small neighborhoods...
            if (length(smallIndex) != minNeighbors) {
                for (j in 1:length(smallIndex)) {
                    smallNeighborhood <- append(smallNeighborhood, imgNeighborsContextTrimmed$source_gene_oid[smallIndex[j]])
                    smallNeighbors <- append(smallNeighbors, imgNeighborsContextTrimmed$gene_oid[smallIndex[j]])
                } 
            } else  {
                next
            }    
        }
        ## now we remove everything on the shortlist from the datafiles
        imgNeighborsTrimmed <- imgNeighborsTrimmed %>%  dplyr::filter(!.data$gene_oid %in% smallNeighbors)
        imgNeighborsContextTrimmed <-imgNeighborsContextTrimmed %>% dplyr::filter(!.data$gene_oid %in% smallNeighbors)
        imgGenesTrimmed <-imgGenesTrimmed %>% dplyr::filter(!.data$gene_oid %in% smallNeighborhood)
        gSeqTable <- gSeqTable %>% dplyr::filter(!.data$gene_oid %in% smallNeighborhood)
        nSeqTable <- nSeqTable %>% dplyr::filter(!.data$gene_oid %in% smallNeighbors)
        print("Truncated neighborhoods trimmed.")
    }
    ## exporting metadata tables
    ## note: note keeping quote=FALSE on the metadata, since that can have weird punctuation
    ## and since we handled quotation marks on input
    ## the punctuation can stay marked as strings
    write.table(imgNeighborsTrimmed, file=fileNameiNT, row.names=FALSE, col.names = TRUE, sep="\t")
    ##    imgNeighborsContextTrimmed <- as.character(imgNeighborsContextTrimmed)
    write.table(imgNeighborsContextTrimmed, file=fileNameiNCT, row.names=FALSE, col.names = TRUE, sep="\t", quote=FALSE)
    write.table(imgGenesTrimmed, file=fileNameiGT, row.names=FALSE, col.names = TRUE, sep="\t")
    ## fasta un-wrangling & export
    for (i in 1:length(nSeqTable[,1])) {
        nReList[[nSeqTable[i,1]]] <- nSeqTable[i,2] 
    }
    imgNeighborSeqsTrimmed <- nReList
    for (i in 1:length(gSeqTable[,1])) {
        gReList[[gSeqTable[i,1]]] <- gSeqTable[i,2] 
    }
    imgGeneSeqsTrimmed <- gReList
    seqinr::write.fasta(imgGeneSeqsTrimmed, names=names(imgGeneSeqsTrimmed), file.out=fileNameiGST)
    seqinr::write.fasta(imgNeighborSeqsTrimmed, names=names(imgNeighborSeqsTrimmed), file.out=fileNameiNST)
    neighborTrimOutput <- list(imgNeighborsTrimmed = imgNeighborsTrimmed, imgGenesTrimmed = imgGenesTrimmed)
    print("Genome neighborhoods trimmed.")
    return(neighborTrimOutput)
}
