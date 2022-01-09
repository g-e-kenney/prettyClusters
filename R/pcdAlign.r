#' Subfunction to automatically align BGCs for visualization.
#'
#' Aligns BGCs for visualization, centering around a gene (family) of interest.
#' @param finalGeneSets Partially processed data frame from IMG or analyzeNeighbors with metadata for your gene (family) of interest.  Has already been through the initial prettyClusterDiagrams input process.
#' @param standAlone Specifies whether this dataset has been through the rest of the prettyClusters pipeline
#' @param uniqueBGCs List of unique BGCs in the dataset
#' @param coreGeneName Name of gene (family) of interest - mostly important if this is being run in standAlone mode
#' @return dataframe containing updated metadata with relative positioning information for the final gene cluster diagrams
#' @export
#' @importFrom utils read.csv write.csv write.table read.table
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @examples 
#' \dontrun{
#' pcdAlign <- pcdAlign(finalGeneSets = finalGeneSets, 
#'                                    standAlone = standAlone, 
#'                                    uniqueBGCs = uniqueBGCs,
#'									  coreGeneName, coreGeneName)
#' }
#'
pcdAlign <- function(finalGeneSets = finalGeneSets, standAlone = standAlone, uniqueBGCs = uniqueBGCs, coreGeneName = coreGeneName) {
	if (standAlone == FALSE) {
            uniqueGOIs <- unique(finalGeneSets$source_gene_oid)
            for (i in 1:length(uniqueGOIs)) {
                aligningGenes <- data.frame()
                aligningGenes <- dplyr::filter(finalGeneSets, .data$source_gene_oid == uniqueGOIs[i])
                core <- which(aligningGenes$gene_oid == uniqueGOIs[i])
                ## neither of these failure cases should be likely since we are the source_gene_oid
                ## so even pseudogenes shouldn't screw this up
                ## case:  gene is missing???
                if (length(core)==0) {next}
                ## case: the gene appears twice as a GOI (shouldn't happen, though it could occur twice in a neighborhood
                ## we can combine this with case: one gene copy
                ## because for the purposes of fixing directions, we can deal easily as long as both copies point the same way
                ## and just choose the first if not
                if (length(core)>=1) {
                    if (all(aligningGenes$direction[core]==1 || aligningGenes$direction[core[1]] == 1)) {
                        if (exists(x="dirGeneSets") == FALSE) {
                            dirGeneSets <- data.frame(aligningGenes, stringsAsFactors = FALSE)
                        } else {
                            dirGeneSets <- rbind(dirGeneSets, aligningGenes)
                        }                    
                    } else if (all(aligningGenes$direction[core]==-1 || aligningGenes$direction[core[1]] == -1)) {
                        last <- length(aligningGenes$end)
                        endgene <- aligningGenes$end[last]
                        newend <- abs(endgene - aligningGenes$start)
                        aligningGenes$start <- abs(endgene - aligningGenes$end)
                        aligningGenes$end <- newend
                        aligningGenes$direction <- aligningGenes$direction * -1
                        if (exists(x="dirGeneSets") == FALSE) {
                            dirGeneSets <- data.frame(aligningGenes, stringsAsFactors = FALSE)
                        } else {
                            dirGeneSets <- rbind(dirGeneSets, aligningGenes)
                        }
                    }
                } 
                rm(aligningGenes)
            }
            processed <- dirGeneSets
            print("Gene clusters aligned around the genes of interest.")
        } else if (standAlone == TRUE) {
        ## i.e. we cannot rely on having source_gene_oid
            for (i in 1:length(uniqueBGCs)) {
                aligningGenes <- data.frame()
                aligningGenes <- dplyr::filter(finalGeneSets, .data$bgc == uniqueBGCs[i])
                core <- which(aligningGenes$gene == coreGeneName)
                ## see above section for notes - note that we have to fall back into sorting things by BGC here
                ## which makes the "more than one core gene" thing more likely, annoyingly
                if (length(core)==0) {next}
                if (length(core)>=1) {
                    if (all(aligningGenes$direction[core]==1 || aligningGenes$direction[core[1]] == 1)) {
                        if (exists(x="dirGeneSets") == FALSE) {
                            dirGeneSets <- data.frame(aligningGenes, stringsAsFactors = FALSE)
                        } else {
                            dirGeneSets <- rbind(dirGeneSets, aligningGenes)
                        }                    
                    } else if (all(aligningGenes$direction[core]==-1 || aligningGenes$direction[core[1]] == -1)) {
                        last <- length(aligningGenes$end)
                        endgene <- aligningGenes$end[last]
                        newend <- abs(endgene - aligningGenes$start)
                        aligningGenes$start <- abs(endgene - aligningGenes$end)
                        aligningGenes$end <- newend
                        aligningGenes$direction <- aligningGenes$direction * -1
                        if (exists(x="dirGeneSets") == FALSE) {
                            dirGeneSets <- data.frame(aligningGenes, stringsAsFactors = FALSE)
                        } else {
                            dirGeneSets <- rbind(dirGeneSets, aligningGenes)
                        }
                    }
                } 
                rm(aligningGenes)
            }
            processed <- dirGeneSets
            print("Gene clusters aligned around the family of genes of interest.")                     
        } 
		return(processed)
	} 