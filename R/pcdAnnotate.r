#' Subfunction to flag genes for highlighting in prettyClusterDiagrams
#'
#' Given IMG-formatted metadata for genomic neighborhoods and a user-provided guide file that defines the families that will be highlighted in the BGC diagrams, this function identifies which genes are hits for these user-provided annotation values.
#' @param imgGenes Data frame from IMG or analyzeNeighbors with metadata for your gene (family) of interest.
#' @param geneSets Partially processed data frame from IMG or analyzeNeighbors with metadata for your gene (family) of interest.  Has already been through the initial prettyClusterDiagrams input process.
#' @param coreGeneList List of genes of interest (your core gene family)
#' @param coreGeneName Name assigned to genes of interest.
#' @param annotationGuide User-provided table with annotation guidelines.
#' @return Data frame with custom assignments for prettyClusterDiagrams based on protein family criteria
#' @export
#' @importFrom utils read.csv write.csv write.table read.table
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @examples 
#' \dontrun{
#' pcdAnnotate <- pcdAnnotate(imgGenes = imgGenes, 
#'                                    geneSets = geneSets, 
#'                                    coreGeneList = coreGeneList, 
#'                                    coreGeneName = coreGeneName, 
#'                                    annotationGuide = annotationGuide)
#' }
#'
 pcdAnnotate <- function(imgGenes = imgGenes, 
						geneSets = geneSets, 
						coreGeneList = coreGeneList, 
						coreGeneName = coreGeneName, 
						annotationGuide = annotationGuide) {
		coreGeneList <- imgGenes$gene_oid
        for (i in 1:length(geneSets$gene)) {
            tempGeneID <- geneSets$gene_oid[i]
            if (any(grepl(tempGeneID, coreGeneList)) == TRUE) {
                geneSets$gene[i] <- as.character(coreGeneName)
            } else {
                next
            }
        }
        geneSets$gene <- as.character(geneSets$gene)
        geneSets$scaffold <- as.character(geneSets$scaffold)
        ## slightly a kludge -  genes where a complete match is required write over previous annotations
        reorderGenes1 <- annotationGuide %>% dplyr::filter(!.data$Requirement == "and")
        reqGenes <- annotationGuide %>% dplyr::filter(.data$Requirement == "and")
        namedGenes <- rbind(reorderGenes1, reqGenes)
        ## bride of the slight kludge - fused genes _also_ write over previous annotations
        reorderGenes2 <- namedGenes %>% dplyr::filter(!.data$Fusion == "yes")
        fusedGenes <- namedGenes %>% dplyr::filter(.data$Fusion == "yes")
        namedGenes <- rbind(reorderGenes2, fusedGenes)
        namedGenes <- namedGenes %>% dplyr::rename(IMGfam = .data$IMG.Term)
        ## "" or NA becomes "none"
        ## then none will become a list element
        namedGenes$Pfam <- namedGenes$Pfam %>% tidyr::replace_na("none")
        namedGenes$Tigrfam <- namedGenes$Tigrfam %>% tidyr::replace_na("none")
        namedGenes$InterPro <- namedGenes$InterPro %>% tidyr::replace_na("none")
        namedGenes$Hypofam <- namedGenes$Hypofam %>% tidyr::replace_na("none")
        namedGenes$IMGfam <- namedGenes$IMGfam %>% tidyr::replace_na("none")
        ## son of the bride of the kludge
        if ("IMGfam" %in% colnames(namedGenes) == FALSE) {namedGenes$IMGfam <- "none"}
        if ("Hypofam" %in% colnames(namedGenes) == FALSE) {namedGenes$Hypofam <- "none"}
        ## the return of the son of the bride of the kludge
        namedGenes$Pfam[namedGenes$Pfam==""]<-"none"
        namedGenes$Tigrfam[namedGenes$Tigrfam==""]<-"none"
        namedGenes$InterPro[namedGenes$InterPro==""]<-"none"      
        namedGenes$IMGfam[namedGenes$IMGfam==""]<-"none"
        namedGenes$Hypofam[namedGenes$Hypofam==""]<-"none"
        annotList <- namedGenes$geneSymbol
        for (i in 1:length(geneSets$gene)) {
            foundMe <- "no"
            foundIn <- "nowhere"
            tempRealData <- geneSets[i,] %>% dplyr::select(.data$Pfam, .data$Tigrfam, .data$InterPro, .data$IMGfam, .data$Hypofam)
            ## entirely hypothetical proteins are easy to flag - no annotations.
            if (all(tempRealData=="none") == TRUE) {
                geneSets$gene[i] <- "hyp"
                next
            } else if (geneSets$gene[i] != "" && any(grepl(geneSets$gene[i],namedGenes$geneSymbol)==TRUE)) {
                ## if it's already annotated according to our key, let's move on
                next
            } else {    
                for (j in 1:length(annotList)) {
                    tempVector1 <- vector()
                    tempVector2 <- vector()
                    tempList1 <- list()
                    tempList2 <- list()
                    tempList3 <- list()
                    ## first list has IMGfam and Hypofam families - both are specific and singular, so any non-none match is legit
                    tempList1 <- strsplit(namedGenes$IMGfam[j], " ")
                    ## adding the dollar sign here makes the regex work, ugh what a kludge
                    ## this is only an issue for Hypofams due to the way the family names/numbers are generated
                    ## might be possible to fix that back in neighborHypothetical?
                    tempList1 <- append(tempList1, paste(strsplit(namedGenes$Hypofam[j], " "),"$",sep=""))
                    ## second list is Pfam and TIGRfam and InterPro - here you can have multiple hits, and we may or may not consider the absence of a family meaningful
                    tempList2 <- strsplit(namedGenes$Pfam[j], " ")
                    tempList2 <- append(tempList2, strsplit(namedGenes$Tigrfam[j], " "))
                    tempList2 <- append(tempList2, strsplit(namedGenes$InterPro[j], " "))
                    ## we do not want "none" to be a possible match when looking at annotation rules for which "any" is an option (thus tempList3), or for commonly absent annotations (IMG.Term, Hypofam)
                    ## but we do want to be able to look for it in some situations (e.g. members of a superfamily but not a subfamily), thus tempList2 retains it
                    tempList1 <- unlist(tempList1[tempList1!="none$"&&tempList1!="none"])
                    tempList2 <- unlist(tempList2) 
                    tempList3 <- unlist(tempList2[tempList2!="none"])
                    ## matches between data and criteria
                    tempVector1 <- sapply(lapply(tempList1,grepl,tempRealData[,4:5]),any)
                    tempVector2 <- sapply(lapply(tempList2,grepl,tempRealData[,1:3]),any)
                    tempVector3 <- sapply(lapply(tempList3,grepl,tempRealData[,1:3]),any)   
                    if (length(tempList1) != 0 && any(tempVector1)) {
                        ## there are annotations and they match the most specific terms (IMG Term or HypoFam)
                        geneSets$gene[i] <- namedGenes$geneSymbol[j]
                        foundMe <- "yes"
                        ## this is a strong match and generally should not be overwritten
                        foundIn <- "exact"
                                  #      print("exit A")
                        next
                    } else if (namedGenes$Requirement[j]=="all" && all(tempVector2) && length(tempList1) != 0 && any(tempVector1) != FALSE) {
                        ## all annotations had to match the requirements and all did (including requirements for the absence of a specific annotation AND HypoFam and/or IMG.Term)
                        geneSets$gene[i] <- namedGenes$geneSymbol[j]
                        foundMe <- "yes"
                        ## this is a strong match and generally should not be ovewritten
                        foundIn <- "exact"
                                 #      print("exit B")
                        next
                    } else if (namedGenes$Requirement[j]=="all" && all(tempVector2) && length(tempList1) == 0) {
                        ## all annotations had to match the requirements and all did (including requirements for the absence of a specific annotation)
                        ## but we didn't have any IMG.Term/HypoFam options, so we're working off of tempList2 (TIGRfam, Pfam, InterPro)
                        geneSets$gene[i] <- namedGenes$geneSymbol[j]
                        foundMe <- "yes"
                        ## this is still a strong match and generally should not be ovewritten
                        foundIn <- "exact"
                                #       print("exit C")
                        next
                    } else if (namedGenes$Requirement[j]=="all" && all(tempVector3)  && namedGenes$RequireNone[j] == "no") {
                        ## all annotations had to match the requirement, and all listed annotations were found, but there may have been extra annotations
                        ## generally in a non-required column (e.g. no InterPro requirements were provided, but we got InterPro annotations)
                        ## BUT we don't actively care about the _absence_ of annotations, so this is still probably a good match
                        ## prevent a less specific match from overwriting a better one, if there is one
                        if (namedGenes$Fusion[j] == "yes") {
                            ## this is a strong match, and because this is a fusion gene, we can overwrite previous matches for one of the components
                            geneSets$gene[i] <- namedGenes$geneSymbol[j]
                            foundMe <- "yes"
                            foundIn <- "exact"
                               #          print("exit D1")
                            next
                        } else if (foundMe != "yes" || foundIn != "exact") {
                            ## this is a strong match, and we don't think we already have one
                            geneSets$gene[i] <- namedGenes$geneSymbol[j]
                            foundMe <- "yes"
                            foundIn <- "exact"
                                        # print("exit D2")
                            next
                        }  
                    } else if (namedGenes$Requirement[j]=="all" && all(tempVector3)   && namedGenes$RequireNone[j] == "yes") {
                        ## all annotations had to match the requirement, and all that weren't "none" did, but we actually cared about "none"
                        ## so the extra annotations may make this a suboptimal match
                        ## prevent a less specific match from overwriting a better one, if there is one
                        if (namedGenes$Fusion[j] == "yes") {
                            ## this is pretty decent match, and because this is a fusion gene, we can overwrite previous matches for one of the components
                            ## but we still want to flag it for the searcher because we did care about "none"
                            geneSets$gene[i] <- paste(namedGenes$geneSymbol[j],"*",sep="")
                            foundMe <- "maybe"
                                         # print("exit E1")
                            next
                        } else if (foundMe != "yes" || foundIn != "exact") {
                            ## we don't have a better answer yet and this is a decent match  
                            ## let's flag this for the end-user so that they can investigate if necessary
                            geneSets$gene[i] <- paste(namedGenes$geneSymbol[j]," \u2020",sep="")
                            foundMe <- "maybe"
                                        # print("exit E2")
                            next
                        }  
                    } else if (namedGenes$Requirement[j]=="all" && any(tempVector3)) {
                        ## all annotations had to match the requirement but only some did, and no fusion proteins were expected
                        ## prevent a less specific match from overwriting a better one, if there is one
                        if (foundMe != "yes" || foundIn != "exact") {
                            ## note that here, we don't want fusion proteins to overwrite: "any" match in tempVector3 could mean matching just 1 component               
                            if (namedGenes$Fusion[j] == "no") {
                                geneSets$gene[i] <- paste(namedGenes$geneSymbol[j]," \u2020",sep="")
                                ## this is a weaker match and it's possible we could overwrite it
                                foundMe <- "maybe"
                                       #  print("exit F")
                                next
                            }
                        }
                    } else if (namedGenes$Requirement[j] =="any" && length(tempList3) != 0 && any(tempVector3)) {
                        ## there were annotations that matched (and any match was allowed)
                        if (foundIn == "exact") {
                            ## an exact match was already found, this is probably not better, move on
                                      # print("exit G1")
                            next
                        } else {
                            geneSets$gene[i] <- namedGenes$geneSymbol[j]
                            ## probably we found it? (we weren't being too picky)
                            foundMe <- "maybe"
                                      #  print("exit G2")
                            next
                        }
                    } else {
                        if (j != length(annotList) && foundMe == "no") {
                            ## no matches were found on this round but it isn't the last, so let's call it other for now
                            ## note that this is by default overwriteable  
                            geneSets$gene[i] <- "other"
                            foundMe <- "maybe"
                                      # print("exit H1")
                            next
                        } else if (j == length(annotList) && foundMe == "no") {
                            ## no matches were found and this is the last chance, so let's finalize on other
                            geneSets$gene[i] <- "other"
                            foundMe <- "yes"
                                      #  print("exit H2")
                            next
                        }
                    } 
                } 
            }
        }
        print("Genes automatically annotated based on the annotationGuideFile file.")
		return(geneSets)
    }