#' Subfunction to identify subgroups of hypothetical proteins in larger IMG datasets
#'
#' Given a desired genome neighborhood and lists of genes of interest and their neighbors, along with their protein sequences, this program uses all-by-all-blast & a clustering method of choice to identify subsets of hypothetical proteins that may represent meaningfully associated proteins.
#' @param iprScanSource Data frame from neighborPrep with metadata for your genes of interest. Character string, required.
#' @param imgNeighborsSource Data frame from neighborPrep with metadata for neighbors of your genes of interest. Character string, required.
#' @param geneName Character string with your gene of interest's name.  Required.
#' @param addPfam Should Pfam assignments be added to the metadata table? Boolean, defaults to TRUE.
#' @param addTigrfam Should TIGRfam assignments be added to the metadata table? Boolean, defaults to TRUE.
#' @param addIPRfam Should InterPro family assignments be added to the metadata table? Boolean, defaults to TRUE.
#' @return Updated metadata for neighboring genes (additional files generated en route)
#' @export
#' @importFrom magrittr %>% 
#' @importFrom utils read.csv write.csv write.table read.table
#' @examples 
#' \dontrun{
#' incorpIprScanOut <- incorpIprScan(iprScanSource = "iprScan.txt", 
#'                                   imgNeighborsSource = "fauxNeighborsData.txt", 
#'                                   geneName = "genE", 
#'                                   addPfam = TRUE, 
#'                                   addTigrfam = TRUE)  
#' }
#' 
incorpIprScan <- function(iprScanSource = iprScanSource, imgNeighborsSource = imgNeighborsSource, geneName = geneName, addPfam = TRUE, addTigrfam = TRUE, addIPRfam = TRUE) { 
    ## import the IMG neighbor data
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
    ## go away scientific notation
    options(scipen = 999)
    imgNeighborsDataTemp <- as.data.frame(read.csv(imgNeighborsSource, header=TRUE, sep="\t" ))
    imgNeighborsDataTemp <- imgNeighborsDataTemp %>% dplyr::mutate_all(~ tidyr::replace_na(as.character(.x), ""))
    imgNeighborsData <- imgNeighborsDataTemp[names(imgNeighborsDataTemp) %in% imgCols]
    if (any(grepl("InterPro", colnames(imgNeighborsDataTemp)))) {
        imgNeighborsData$InterPro <- imgNeighborsDataTemp$InterPro
    }   else {
        imgNeighborsData$Interpro <- ""
    }
    ## import the InterProScan annotations
    ## note that this assumes a basic InterProScan run (including InterPro family output but no GO terms or whatever
    ## if those don't show up at the very end, this could be an issue for some people
    iprScan <- as.data.frame(read.csv(iprScanSource, header=FALSE, sep="\t" ))[,1:13]
    colnamesIPR <-c("gene_oid", "md5", "aa", "sig.analysis", "sig.accession", "sig.description", "start", "stop", "score", "status", "date", "interpro.accession", "interpro.description")
                                        #    iprScan <- iprScan %>% dplyr::select(colnamesIPR)
                                        #    iprScan <- iprScan %>% dplyr::mutate_all(~ tidyr::replace_na(as.character(.x), ""))
    colnames(iprScan) <- colnamesIPR
    ## if we want to add pfam annotation
    if (addPfam == TRUE) {
        ## Find the Pfam, Tigrfam, and InterPro hits
        pfamLocs <- which(iprScan$sig.analysis == "Pfam")
        iprPfams <- iprScan[pfamLocs,]
        havePfams <- which(imgNeighborsData$gene_oid %in% iprPfams$gene_oid)
        ## add Pfam data to the neighbor data table
        ## also IMG-format the pfams
        for (i in 1:length(havePfams)) {
            pfamIdx <- grep(imgNeighborsData$gene_oid[havePfams[i]],iprPfams$gene_oid)
            imgNeighborsData$Pfam[havePfams[i]] <-  stringr::str_c(unlist(iprPfams$sig.accession[pfamIdx]), collapse=" ")
            imgNeighborsData$Pfam[havePfams[i]] <- gsub("PF","pfam",imgNeighborsData$Pfam[havePfams[i]])
        }
    }
    
    ## if we want to add tigrfam annotation
    if (addTigrfam == TRUE) {
        tigrLocs <- which(iprScan$sig.analysis == "TIGRFAM")
        iprTigrfams <-  iprScan[tigrLocs,]
        haveTigrfams <- which(imgNeighborsData$gene_oid %in% iprTigrfams$gene_oid)
        ## add Tigrfam data to the neighbor data table
        for (i in 1:length(haveTigrfams)) {
            tigrfamIdx <- grep(imgNeighborsData$gene_oid[haveTigrfams[i]],iprTigrfams$gene_oid)
            imgNeighborsData$Tigrfam[haveTigrfams[i]] <-  stringr::str_c(unlist(iprTigrfams$sig.accession[tigrfamIdx]), collapse=" ")
        }
    }
    ## if we want to add InterPro annotation
    if (addIPRfam == TRUE) {
        iprLocs <- grep("IPR",iprScan$interpro.accession)
        iprIPRs <- iprScan[iprLocs,]
        haveIPRs <- which(imgNeighborsData$gene_oid %in% iprIPRs$gene_oid)
        ## add InterPro data to the neighbor data table
        for (i in 1:length(haveIPRs)) {
            iprIdx <- grep(imgNeighborsData$gene_oid[haveIPRs[i]],iprIPRs$gene_oid)
            imgNeighborsData$InterPro[haveIPRs[i]] <-  stringr::str_c(unlist(iprIPRs$interpro.accession[iprIdx]), collapse=" ")
        }
    }
    ## save
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    fileName <- paste(fileDate,"_incorpIprScan_",geneName,"_neighborData.txt",sep="")
    write.table(imgNeighborsData, file=fileName, row.names=FALSE, col.names = TRUE, sep="\t")
    print("InterProScan annotations added!")
    return(imgNeighborsData)
}

