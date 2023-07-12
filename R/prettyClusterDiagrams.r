#' Umbrella function generates pretty gene cluster diagrams
#'
#' This function generates gene cluster diagrams incorporating annotation and clustering info from several sources.
#' @param imgGenesFile What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param imgNeighborsFile What is the file with the metadata for neighbors of your gene of interest? Filename as string ("filename.txt")
#' @param annotationGuideFile User-defined annotation fine, assigning gene names via protein families.  Filename as string ("filename.txt")
#' @param geneName What is the name of your gene? Gene name as string ("genE")
#' @param efiRepnodes Are values here post-EFI-EST? T/F value, defaults to FALSE
#' @param neighborNumber How many genes are expected on either side of the gene of interest? Integer (e.g. 10)
#' @param annotateGenes Should the genes be auto-annotated using the geneFormat file? T/F, defaults to TRUE.
#' @param standAlone Should we assume the data has not been through the full pipeline? T/F, defaults to FALSE.
#' @param markClusters Should cluster labels be included in the final figure? T/F, defaults to FALSE.
#' @param autoColor If TRUE, generates palettes automatically. If FALSE, requires hex-formatted colors per gene in the last column of the geneFormat data. Defaults to TRUE.
#' @param colorType What palette family should be used? String, defaults to "viridis"
#' @param paletteInput What specific palette within that family should be used? String, defaults to "plasma"
#' @param showScaffold Should the scaffold be included in name labels? T/F, defaults to FALSE
#' @param alignToCore Should the core gene be centered in the diagrams? T/F, defaults to TRUE
#' @param labelGenes Should all genes be labeled in the diagrams? T/F, defaults to FALSE
#' @param subclusterDiagrams Should additional figures be made for each group of gene clusters? T/F, defaults to FALSE. 
#' @param everyScale Should every single gene cluster have a scale added? T/F, defaults to FALSE.
#' @param makeScale Should a 1 kb-delineated scale be included in the figures? T/F, defaults to TRUE.
#' @return Gene diagrams as .png and .pdf files.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @importFrom utils capture.output combn read.csv write.csv write.table
#' @importFrom stats dist end hclust na.omit start setNames
#' @importFrom grDevices cairo_pdf cairo_ps col2rgb colorRampPalette dev.off hsv pdf rgb2hsv
#' @examples
#' \dontrun{
#' prettyClusterOutput <- prettyClusterDiagrams(imgGenesFile = "genesFile.txt", 
#'                                              imgNeighborsFile = "neighborsFile.txt", 
#'                                              annotationGuideFile = "formatFile.txt", 
#'                                              geneName = "genE", 
#'                                              neighborNumber = 10)
#' }
#'  
prettyClusterDiagrams <- function(imgGenesFile = imgGenesFile,
                                  imgNeighborsFile = imgNeighborsFile,
                                  annotationGuideFile = annotationGuideFile,
                                  geneName = geneName,
                                  neighborNumber = neighborNumber,
                                  efiRepnodes = FALSE,
                                  annotateGenes = TRUE,
                                  standAlone = FALSE,
                                  markClusters = FALSE,
                                  autoColor = TRUE,
                                  colorType = "viridis",
                                  paletteInput = "plasma",
                                  showScaffold = FALSE,
                                  alignToCore = TRUE,
                                  labelGenes = FALSE,
                                  subclusterDiagrams = FALSE,
                                  everyScale = FALSE,
                                  makeScale = TRUE) { 
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    if (efiRepnodes == TRUE) {
        coreGeneName <- geneName
        geneName <- paste(geneName, "_repnodes",sep="")
    } else {
        coreGeneName <- geneName
    }
    fileName <- paste(fileDate,"_prettyClusters_",geneName,sep="")
    finalpngname <- paste(fileName,"_with-axes.png",sep="")
    finalpdfname <- paste(fileName,"_with-axes.pdf",sep="")
    finalpdfnameX <- paste(fileName,"_no-axes.pdf",sep="")
    finalpngnameX <- paste(fileName,"_no-axes.png",sep="")
    annotationFileName <- paste(fileName, "_annotation.txt",sep="")
    imgCols <- list("gene_oid","Locus.Tag", "Start.Coord","End.Coord", "Strand", "Gene.Symbol", "Scaffold.Name", "Scaffold.ID", "Genome.Name","Genome.ID","Pfam", "Tigrfam", "InterPro")
    ggCols <- list("gene_oid", "locus_tag", "start", "end", "strand","gene","scaffold", "scaffoldID", "genome","genomeID","Pfam", "Tigrfam", "InterPro", "clustNum", "clustOrd","source_gene_oid")
                                        # geneset and annotation processing
    annotationGuide <- read.csv(annotationGuideFile, header=TRUE,  sep = ",", stringsAsFactors = FALSE)
    annotationGuide$Pfam[which(is.na(annotationGuide$Pfam))] <- ""
    annotationGuide$Tigrfam[which(is.na(annotationGuide$Tigrfam))] <- ""
    annotationGuide$InterPro[which(is.na(annotationGuide$InterPro))] <- ""  
    annotationGuide$Hypofam[which(is.na(annotationGuide$Hypofam))] <- ""
    annotationGuide$IMG.Term[which(is.na(annotationGuide$IMG.Term))] <- ""
    annotationGuide$Color[which(is.na(annotationGuide$Color))] <- ""
    if(typeof(imgGenesFile) == "character") {
        imgGenes <-read.csv(imgGenesFile, header=TRUE, sep="\t", stringsAsFactors=FALSE )
    } else {
        imgGenes <- imgGenesFile
    }
    if(typeof(imgNeighborsFile) == "character")  {
        imgNeighbors <- read.csv(imgNeighborsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE )
    } else {
        imgNeighbors <- imgNeighborsFile
    }
    print("Opening data files.")
    ## img metadata conversion - filters out things other than the core info because why do we need it
    inputColNames <- colnames(imgNeighbors)
    inputColNames <- stringr::str_replace_all(inputColNames, c(" " = "." , "," = "" ))
    colnames(imgNeighbors) <- inputColNames
    geneSets <- imgNeighbors[names(imgNeighbors) %in% imgCols]
    geneSets <- geneSets %>% dplyr::rename(locus_tag = .data$Locus.Tag)
    geneSets <- geneSets %>% dplyr::rename(scaffold = .data$Scaffold.Name)
    geneSets <- geneSets %>% dplyr::rename(scaffoldID = .data$Scaffold.ID)
    geneSets <- geneSets %>% dplyr::rename(start = .data$Start.Coord)
    geneSets <- geneSets %>% dplyr::rename(end = .data$End.Coord)
    geneSets <- geneSets %>% dplyr::rename(strand = .data$Strand)
    geneSets <- geneSets %>% dplyr::rename(gene = .data$Gene.Symbol)
    geneSets <- geneSets %>% dplyr::rename(genome = .data$Genome.Name)
    geneSets <- geneSets %>% dplyr::rename(genomeID = .data$Genome.ID)
    geneSets <- data.table::as.data.table(geneSets, stringsAsFactors=FALSE)
    ## this is easier than dealing with different possible IMG-like metadata inputs
    if ("InterPro" %in% colnames(imgNeighbors)) {
        geneSets$InterPro <- as.character(imgNeighbors$InterPro)
    } else {
        geneSets$InterPro <- "none"
    }
    if ("IMG.Term" %in% colnames(imgNeighbors)) {
        geneSets$IMGfam <- as.character(imgNeighbors$IMG.Term)
    } else {
        geneSets$IMGfam <- "none"
    }
    if ("Hypofam" %in% colnames(imgNeighbors)) {
        geneSets$Hypofam <- as.character(imgNeighbors$Hypofam)
    } else {
        geneSets$Hypofam <- "none"
    }
    if ("clustOrd" %in% colnames(imgNeighbors)) {
        geneSets$clustOrd <- as.character(imgNeighbors$clustOrd)
    } else {
        geneSets$clustOrd <- "none"
    }
    if ("clustNum" %in% colnames(imgNeighbors)) {
        geneSets$clustNum <- as.character(imgNeighbors$clustNum)
    } else {
        geneSets$clustNum <- "none"
    }
    if ("source_gene_oid" %in% colnames(imgNeighbors)) {
        geneSets$source_gene_oid <- as.character(imgNeighbors$source_gene_oid)
    } else {
        geneSets$source_gene_oid <- "none"
    }
    geneSets <- geneSets[order(geneSets$gene_oid),]
    ## if this is being run as a standalone figure generator, there possibly aren't any source gene IDs
    if (standAlone == TRUE) {
        ## if labels are manually added, let's respect them
        coreGeneList <- imgGenes$gene_oid
        coreGeneList <- sort(coreGeneList)
        for (m in 1:length(coreGeneList)) {
                ## we need to figure out what GOIs have overlapping neighborhoods - if any
            multiCore <- coreGeneList %>% dplyr::between(coreGeneList[m]-neighborNumber, coreGeneList[m]+neighborNumber) %>% which()
            realMultis <- coreGeneList[multiCore]
            if (length(multiCore) == 1) {
                ## if we have only one GOI locally we can handle this easily
                neighborSubset <- geneSets %>% dplyr::filter(.data$gene_oid >= coreGeneList[m] - neighborNumber & .data$gene_oid <= coreGeneList[m] + neighborNumber)
                neighborList <- neighborSubset$gene_oid
                neighborSubset$source_gene_oid[neighborSubset$gene_oid %in% neighborList] <- coreGeneList[m]
                if (exists(x="fixGeneSets") == FALSE) {
                    fixGeneSets <- data.frame(neighborSubset, stringsAsFactors = FALSE)
                } else {
                    fixGeneSets <- rbind(fixGeneSets, neighborSubset)
                } 
            } else {
                ## if we have 3+ GOIs locally this is gonna be hideous and algorithmic fixes might break in a truncated gene cluster
                ## let's bail and make the user fix it
                if (length(realMultis) > 2)  {
                    goiError <- paste("There are 3 or more GOIs sharing a neighborhood with ",
                                      coreGeneList[m],
                                      ". Please address this manually in the source_gene_oid column of the neighborMetadata."
                                     ,sep="")
                    print(goiError)
                    next
                }                             
                ## if we have 2 GOIs in one area and this is the first time we've hit this area, let's handle it
                if (length(multiCore) > 1 && coreGeneList[m] == coreGeneList[multiCore[1]]) {
                    neighborSubset <- geneSets %>% dplyr::filter(.data$gene_oid >= coreGeneList[multiCore[1]] - neighborNumber & .data$gene_oid <= coreGeneList[multiCore[length(multiCore)]] + neighborNumber)
                    neighborList <- neighborSubset$gene_oid
                    neighborSubset$index <- seq_len(nrow(neighborSubset)) %% 2
                    neighborSubset$index[1] <- 0
                    neighborSubset$index[length(neighborList)] <- 1
                    neighborSubset$source_gene_oid[neighborSubset$index == 0] <- coreGeneList[multiCore[1]]
                    neighborSubset$source_gene_oid[neighborSubset$index == 1] <- coreGeneList[multiCore[2]]
                    neighborSubset <- neighborSubset[,-"index"]
                    if (exists(x="fixGeneSets") == FALSE) {
                        fixGeneSets <- data.frame(neighborSubset, stringsAsFactors = FALSE)
                    } else {
                        fixGeneSets <- rbind(fixGeneSets, neighborSubset)
                    } 
                } else if (length(multiCore) > 1 && coreGeneList[m] != coreGeneList[multiCore[1]]) {
                    ## in this case we should have already handled the area
                    next
                }
            }
        }
        geneSets <- fixGeneSets
        print("Source gene IDs added to neighbor metadata.")
    }
    ## dealing with blank cells & misc metadata processing
    geneSets$Pfam[geneSets$Pfam==""]<-"none"
    geneSets$Tigrfam[geneSets$Tigrfam==""]<-"none"
    geneSets$InterPro[geneSets$InterPro==""]<-"none"
    geneSets$IMGfam[geneSets$IMGfam==""]<-"none"
    geneSets$Hypofam[geneSets$Hypofam==""]<-"none"
    geneSets$direction <- ifelse(geneSets$strand == "+", 1, -1)
    ## dealing with NAs
    geneSets$Pfam <- geneSets$Pfam %>% tidyr::replace_na("none")
    geneSets$Tigrfam <- geneSets$Tigrfam %>% tidyr::replace_na("none")
    geneSets$InterPro <- geneSets$InterPro %>% tidyr::replace_na("none")
    geneSets$Hypofam <- geneSets$Hypofam %>% tidyr::replace_na("none")
    geneSets$IMGfam <- geneSets$IMGfam %>% tidyr::replace_na("none")
    ## misc 
    rowidx <- order(geneSets$gene_oid)
    geneSets <- geneSets[rowidx,,drop=FALSE]
    ## in case of stupid factor stuff
    geneSets$gene <- as.character(geneSets$gene)
    geneSets$genome <- as.character(geneSets$genome)
    geneSets$scaffold <- as.character(geneSets$scaffold)
    geneSets$strand <- as.character(geneSets$strand)
                                        # annotation - if chosen - begins here
    ## Adding annotations acording to the highlightGenes file IF annotateGenes is true
    ## By default requires Pfam and Tigrfam, can accept Hypofam and IMG Terms
    ## Deals with (known sets of) fused genes and multiple options per category
    ## With un-specified fused genes, you'll just end up with whichever annotation category comes later, ugh
    ## Also lets you specify matching some or all terms (the latter for subgroups)
    if (annotateGenes == TRUE) {
        ## this looks slightly different in the standalone version
        ## because we can't assume there's a source_gene_oid column
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
            } else if (geneSets$gene[i] != "" && length(which(namedGenes$geneSymbol == geneSets$gene[i]))>0) {
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
                    ## we may eventually want to change this for hypofam since there could be different subgroups
                    ## in the interim just use the most specific
                    tempList1 <- strsplit(namedGenes$IMGfam[j], " ")
                    ## adding the dollar sign here makes the regex work, ugh what a kludge
                    ## this is only an issue for Hypofams due to the way the family names/numbers are generated
                    ## i.e. we don't want hypofam_2 to dredge up hypofam_20
                    ## might be possible to fix that back in neighborHypothetical?
                    if (grepl(" ",namedGenes$Hypofam[j])==FALSE) {
                        ## this is for a case where we have only one hypofam term
                        tempList1 <- append(tempList1, paste(strsplit(namedGenes$Hypofam[j], " "),"$",sep=""))
                    } else {
                        ## in this case, we listed multiple hypofam terms
                        tempList1 <- append(tempList1, paste(unlist(strsplit(namedGenes$Hypofam[j], " ")),"$",sep=""))
                    }
                    ## second list is Pfam and TIGRfam and InterPro - here you can have multiple hits, and we may or may not consider the absence of a family meaningful
                    tempList2 <- strsplit(namedGenes$Pfam[j], " ")
                    tempList2 <- append(tempList2, strsplit(namedGenes$Tigrfam[j], " "))
                    tempList2 <- append(tempList2, strsplit(namedGenes$InterPro[j], " "))
                    ## we do not want "none" to be a possible match when looking at annotation rules for which "any" is an option (thus tempList3), or for commonly absent annotations (IMG.Term, Hypofam)
                    ## but we do want to be able to look for it in some situations (e.g. members of a superfamily but not a subfamily), thus tempList2 retains it
                    tempList1 <- unlist(tempList1[tempList1 != "none$" & tempList1 != "none"])
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
                                     #   print("exit A")
                        next
                    } else if (namedGenes$Requirement[j]=="all" && all(tempVector2) && length(tempList1) != 0 && any(tempVector1) != FALSE) {
                        ## all annotations had to match the requirements and all did (including requirements for the absence of a specific annotation AND HypoFam and/or IMG.Term)
                        geneSets$gene[i] <- namedGenes$geneSymbol[j]
                        foundMe <- "yes"
                        ## this is a strong match and generally should not be ovewritten
                        foundIn <- "exact"
                                     #  print("exit B")
                        next
                    } else if (namedGenes$Requirement[j]=="all" && all(tempVector2) && length(tempList1) == 0) {
                        ## all annotations had to match the requirements and all did (including requirements for the absence of a specific annotation)
                        ## but we didn't have any IMG.Term/HypoFam options, so we're working off of tempList2 (TIGRfam, Pfam, InterPro)
                        geneSets$gene[i] <- namedGenes$geneSymbol[j]
                        foundMe <- "yes"
                        ## this is still a strong match and generally should not be ovewritten
                        foundIn <- "exact"
                                    #   print("exit C")
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
                                     #   print("exit D1")
                            next
                        } else if (foundMe != "yes" || foundIn != "exact") {
                            ## this is a strong match, and we don't think we already have one
                            geneSets$gene[i] <- namedGenes$geneSymbol[j]
                            foundMe <- "yes"
                            foundIn <- "exact"
                                     #    print("exit D2")
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
                                     #     print("exit E1")
                            next
                        } else if (foundMe != "yes" || foundIn != "exact") {
                            ## we don't have a better answer yet and this is a decent match  
                            ## let's flag this for the end-user so that they can investigate if necessary
                            geneSets$gene[i] <- paste(namedGenes$geneSymbol[j]," \u2020",sep="")
                            foundMe <- "maybe"
                                     #    print("exit E2")
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
                                     #    print("exit F")
                                next
                            }
                        }
                    } else if (namedGenes$Requirement[j] =="any" && length(tempList3) != 0 && any(tempVector3)) {
                        ## there were annotations that matched (and any match was allowed)
                        if (foundIn == "exact") {
                            ## an exact match was already found, this is probably not better, move on
                                     #  print("exit G1")
                            next
                        } else {
                            geneSets$gene[i] <- namedGenes$geneSymbol[j]
                            ## probably we found it? (we weren't being too picky)
                            foundMe <- "maybe"
                                     #   print("exit G2")
                            next
                        }
                    } else {
                        if (j != length(annotList) && foundMe == "no") {
                            ## no matches were found on this round but it isn't the last, so let's call it other for now
                            ## note that this is by default overwriteable  
                            geneSets$gene[i] <- "other"
                            foundMe <- "maybe"
                                    #   print("exit H1")
                            next
                        } else if (j == length(annotList) && foundMe == "no") {
                            ## no matches were found and this is the last chance, so let's finalize on other
                            geneSets$gene[i] <- "other"
                            foundMe <- "yes"
                                     #   print("exit H2")
                            next
                        }
                    } 
                } 
            }
        }
        print("Genes automatically annotated based on the annotationGuideFile file.")
    }
    ## this should be true independent of auto-annotation
    geneSets$gene[which(geneSets$gene_oid %in% imgGenes$gene_oid)] <- coreGeneName
    ## if this is being run as a standalone figure generator, we can't assume scaffolds have been trimmed to avoid wrong-scaffold neighbors
    if (standAlone == TRUE) {
        initScaff <- unique(geneSets$scaffoldID)
        for (m in 1:length(initScaff)) {
            scaffTest <- geneSets %>% dplyr::filter(.data$scaffoldID == initScaff[m])
            if (any(grepl(coreGeneName,scaffTest$gene)) == TRUE) {
                next
            } else {
                geneSets <- geneSets %>% dplyr::filter(!.data$scaffoldID == initScaff[m])
            }
        }
        print("Genes on incorrect scaffolds trimmed.")
    }
    ## occasionally img bizarrely lacks a genome ID and/or name???
    ## i dunno man?  still beats dealing with NCBI/WGS though
    ## anyway, this uses the scaffold name and ID to generate faux IMG names/IDs for the genome
    for (i in 1:length(geneSets$genomeID)) {
        if (is.na(geneSets$genome[i]) || geneSets$genome[i] == "") {
            tempGName <- geneSets$scaffold[i]
            tempGName <- sub(" [:].*$","",tempGName)
            geneSets$genome[i] <- tempGName
            ## not going to overlap with existing numbering schemes
            geneSets$genomeID[i] <- geneSets$scaffoldID[i]+10000000000
        } else {
            next
        }
        print("Absent genome names addressed.")
    }
    ## similarly, sometimes genome names are duplicated
    ## usually when people "helpfully" reannotate and resubmit genomes or whatever
    ## this deals with duplicated genome names ane makes them v1, v2, etc.
    maxNeighbor <- neighborNumber*2+1
    goiSeek <- grep(coreGeneName,geneSets$gene)
    if (any(duplicated(geneSets$genome[goiSeek]))) {
        dupes <- which(duplicated(geneSets$genome[goiSeek]) == TRUE)
        for (i in 1:length(dupes)) {
            ## because genome names can have horrible characters for regex purposes...      
            grepName <- geneSets$genome[goiSeek][dupes[i]]
            grepName <- gsub("\\\\", "\\\\\\", grepName)   
            grepName <- gsub("\\(", "\\\\(", grepName)
            grepName <- gsub("\\)", "\\\\)", grepName)
            grepName <- gsub("\\.", "\\\\.", grepName)
            grepName <- gsub("\\|", "\\\\|", grepName)
            grepName <- gsub("\\+", "\\\\+", grepName)           
            grepName <- gsub("\\?", "\\\\?", grepName)   
            grepName <- gsub("\\^", "\\\\^", grepName)  
            grepName <- gsub("\\$", "\\\\$", grepName)  
            grepName <- gsub("\\[", "\\\\[", grepName)
            grepName <- gsub("\\]", "\\\\]", grepName)
            grepName <- gsub("\\{", "\\\\{", grepName)    
            grepName <- gsub("\\}", "\\\\}", grepName)   
            grepName <- gsub("\\-", "\\\\-", grepName)        
            ## I don't think =, ;, :, or / should cause issues but here's a note in case they do
            dupeList <- grep(grepName,geneSets$genome)
            dupeGenomes <- unique(geneSets$genomeID[dupeList])
            if (length(dupeGenomes) == 1)  {
                ## i.e. there are multiple bgcs within a single genome, which will be dealt with later
                next
            } else {
                ## there are multiple genomes, so they get version numbers
                for (k in 1:length(dupeGenomes)) {
                    dupeGenomeLoc <- grep(dupeGenomes[k], geneSets$genomeID)
                    newGName <- geneSets$genome[dupeGenomeLoc[1]]
                    newGName <- paste(newGName, " v",k,sep="")
                    geneSets$genome[dupeGenomeLoc] <- newGName
                }
            }
        }
        print("Resquenced genome name duplication addressed.")
    }
    ## dealing with species that have multiple BGCs with your gene type of interest
    ## initially it looks for places where the max and min coords are < 50kb apart
    ## if they are within 50 kb, they're considered contiguous gene clusters, with multiple copies and/or pseudogenes as the likely explanation
    ## (why 50 kb? might be huge NRPS megagenes or something)
    ## if they are further, they get more analysis
    ## unless a big jump comes between 2 genes (25k!), a bgc end/beginning is not registered
    ## names are adjusted as necessary to add multiple BGCs and thus not break gggenes
    contigs <- unique(geneSets$scaffoldID)
    species <- unique(geneSets$genomeID)
    multiBGC <- data.frame()
    multicoord <- list()
    for(k in 1:length(species)) {
        tempSpecies <- geneSets %>% dplyr::filter(.data$genomeID == species[k])
        specContigs <- unique(tempSpecies$scaffoldID)
        for (j in 1:length(specContigs))  {
            tempContig <- tempSpecies %>% dplyr::filter(.data$scaffoldID == specContigs[j])
            ## using a rough "within 50kb" guideline to catch multiple gene clusters on the same scaffold
            if (max(tempContig$start) <= min(tempContig$start) + 50000 && k == 1 && j == 1) {
                tempContig$bgc <- paste(tempContig$genome[1], " (bgc ",j,")",sep="")
                multiBGC <- tempContig
                next
            } else if  (max(tempContig$start) <= min(tempContig$start) + 50000 && ((j >= 2) || (k >= 2))) {
                tempContig$bgc <- paste(tempContig$genome[1], " (bgc ",j,")",sep="")
                multiBGC <- rbind(multiBGC, tempContig)
                next
            } else {
                ## ok so we failed our 50 kb, now what?
                for (m in 1:length(tempContig$start)) {
                    if (m == 1) {
                        tempCoord <- tempContig$start[m]
                    } else if (m != length(tempContig$start)) {
                        tempDiff <- tempContig$start[m] - tempContig$start[m-1]
                        ## this should catch huge-ass gene clusters
                        if (abs(tempDiff) > 25000) {
                            tempCoord <- append(tempCoord, tempContig$start[m-1])
                            tempCoord <- append(tempCoord, tempContig$start[m])
                        } else  {
                            next
                        }
                    } else if (m == length(tempContig$start)) {
                        tempCoord <- append(tempCoord, tempContig$start[m])
                    }
                }
                ## if we have at least 4 coordinates - 2x cluster start and 2x cluster end - on one scaffold...
                if (length(tempCoord) >= 4) {
                    numBGC <- length(tempCoord)/2
                    ## this'll number BGCs that show up on the same scaffold (e.g. Genus species BGC 1.1)
                    for (g in 1:length(tempContig$start)) {
                        if (g == 1) {
                            whichBGC <- 1
                            tempContig$bgc[g] <- paste(tempContig$genome[g], " (bgc ",j,".",whichBGC,")",sep="")
                        } else {
                            if (dplyr::between(tempContig$start[g], min(tempCoord[c(whichBGC*2-1,whichBGC*2)]), max(tempCoord[c(whichBGC*2-1,whichBGC*2)])) == TRUE) {
                                whichBGC <- whichBGC
                                tempContig$bgc[g] <- paste(tempContig$genome[g], " (bgc ",j,".",whichBGC,")",sep="")
                            } else  {
                                whichBGC <- whichBGC + 1
                                tempContig$bgc[g] <- paste(tempContig$genome[g], " (bgc ",j,".",whichBGC,")",sep="")
                            }
                        }
                    }
                    if (j == 1 && k == 1) {
                        multiBGC <- tempContig
                    } else {
                        multiBGC <- rbind(multiBGC, tempContig)
                    }
                } else {
                    if (j == 1 && k == 1) {
                        tempContig$bgc <- paste(tempContig$genome[1], " (bgc ",j,")",sep="")
                        multiBGC <- tempContig
                    } else {
                        tempContig$bgc <- paste(tempContig$genome[1], " (bgc ",j,")",sep="")
                        multiBGC <- rbind(multiBGC, tempContig) 
                    } 
                } 
            }  
        }
    }    
    uniBGCs <- unique(multiBGC$bgc)    
    print("Gene cluster naming fixed in species with multiple clusters.")
    ## this normalizes positioning. 
    ## we do not fundamentally care that much about absolute genome position as a default here
    ## nucleotide position number is arbitrary to begin with and meaningless when comparing across wildly different species
    ## so let's make every cluster start at the beginning
    ## some of this may end up in flux as i try to nuke the extra x-axes.
    maxval <- 0
    for (i in 1:length(uniBGCs)) {
        tempPositioning <- dplyr::filter(multiBGC, .data$bgc == uniBGCs[i])
        ## JUST IN CASE
        tempPositioning <- tempPositioning[order(tempPositioning$start),]
        zeroval <- tempPositioning$start[1]
        tempPositioning$start <- tempPositioning$start - zeroval
        tempPositioning$end <- tempPositioning$end - zeroval
        genenum <- length(tempPositioning$end)
        finalnum <- tempPositioning$end[genenum]
        if (finalnum >= maxval) {
            maxval <- finalnum
        }
        if (i == 1) {
            finalGeneSets <- data.frame(tempPositioning, stringsAsFactors = FALSE)
        } else {
            finalGeneSets <- rbind(finalGeneSets, tempPositioning)
        }
        rm(tempPositioning)
    }
    finalGeneSets$start <- as.numeric(finalGeneSets$start)
    finalGeneSets$end <- as.numeric(finalGeneSets$end)
    finalGeneSets$direction <- as.numeric(finalGeneSets$direction)
    print("Genes repositioned on a uniform relative scale.")
    ## finally, quantification of gene types
    ## adds percent of BGCs that have a specific gene type as a caption
    ## skips hypothetical genes, MGEs, and non-specified genes
    ## also adds a title because why not
    ## note that unless both ORF identification and annotation are entirely correct, the present-in-X%-of-BGCs count may be an undercount
    ## to-do: cluster-based abundance values (i.e. in cluster one, genA is there 25% of the time, but it's there 95% of the time in cluster two)
    if (annotateGenes == TRUE)  {
        finalUniqueBGCs <- unique(finalGeneSets$bgc)
        numNamedGenes <- length(namedGenes$gene)
        wrapGenes <- ceiling(numNamedGenes/2)
        namedGenes$num <- 0
        namedGenes$percent <- 0
        for (i in 1:length(finalUniqueBGCs))  {
            forGrep <- finalUniqueBGCs[i]
            ## oof look at that regex backslash pileup
            forGrep <- gsub("\\\\", "\\\\\\", forGrep)   
            forGrep <- gsub("\\(", "\\\\(", forGrep)
            forGrep <- gsub("\\)", "\\\\)", forGrep)
            forGrep <- gsub("\\.", "\\\\.", forGrep)
            forGrep <- gsub("\\|", "\\\\|", forGrep)
            forGrep <- gsub("\\+", "\\\\+", forGrep)           
            forGrep <- gsub("\\?", "\\\\?", forGrep)   
            forGrep <- gsub("\\^", "\\\\^", forGrep)  
            forGrep <- gsub("\\$", "\\\\$", forGrep)  
            forGrep <- gsub("\\[", "\\\\[", forGrep)
            forGrep <- gsub("\\]", "\\\\]", forGrep)
            forGrep <- gsub("\\{", "\\\\{", forGrep)    
            forGrep <- gsub("\\}", "\\\\}", forGrep)   
            forGrep <- gsub("\\-", "\\\\-", forGrep) 
            bgcLoc <- grep(forGrep,finalGeneSets$bgc) 
            for (j in 1:length(namedGenes$geneSymbol))  {
                hitNum <- length(grep(namedGenes$geneSymbol[j], finalGeneSets$gene[bgcLoc]))
                if (hitNum >= 1) {
                    namedGenes$num[j] <- namedGenes$num[j] + 1
                } else {
                    namedGenes$num[j] <- namedGenes$num[j] + 0
                }
            }
        }
        captionText <- "**gene abundance:** "
        for (k in 1:numNamedGenes) {
            if ( namedGenes$geneSymbol[k] == "hyp" || namedGenes$geneSymbol[k] == "mge") {
                next
            } else  {
                if (k != numNamedGenes) {
                    namedGenes$percent[k] <- round(namedGenes$num[k] / length(finalUniqueBGCs), digits=2)
                    printPercent <- namedGenes$percent[k]*100
                    tempGeneSymbol <- namedGenes$geneSymbol[k]
                    captionText <- paste(captionText,"*",tempGeneSymbol,"* (",printPercent,"%), ",sep="")
                } else  {
                    namedGenes$percent[k] <- round(namedGenes$num[k] / length(finalUniqueBGCs), digits=2)
                    printPercent <- namedGenes$percent[k]*100
                    tempGeneSymbol <- namedGenes$geneSymbol[k]
                    captionText <- paste(captionText,"*",tempGeneSymbol,"* (",printPercent,"%).",sep="")
                }
            }
        }
        figureTitle <- paste("<b>genomic neighborhood diagrams for *",coreGeneName,"*</b> (\u00b1", neighborNumber," genes displayed)",sep="")
        print("Gene abundance calculated based on annotationGuideFile-derived annotations.")
    } else {
        figureTitle <- paste("<b>genomic neighborhood diagrams for *",coreGeneName,"*</b> (\u00b1", neighborNumber," genes displayed)",sep="")
        captionText <- "No automated annotation and gene quantification for these diagrams."
    }
    uniqueBGCs <- unique(finalGeneSets$bgc)
    ## aligning gene clusters to your gene of interest
    ## without alignToCore selected, expect all clusters to be left-aligned, and in many situations, this is not preferred
    ## also without this they won't necessarily all point the same direction
    if (alignToCore == TRUE) {
        if (standAlone == FALSE) {
            uniqueGOIs <- unique(finalGeneSets$source_gene_oid)
            for (i in 1:length(uniqueGOIs)) {
                aligningGenes <- data.frame()
                aligningGenes <- dplyr::filter(finalGeneSets, .data$source_gene_oid == uniqueGOIs[i])
                core <- which(aligningGenes$gene_oid == uniqueGOIs[i])
                ## neither of these failure cases should be likely since we have the source_gene_oid
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
            rm(dirGeneSets)
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
            rm(dirGeneSets)
        } 
    } else {
        ## no aligning to core
        processed <- finalGeneSets
    }
    ## palette generation - needed or not?
    if (autoColor == FALSE) {
        ## add bit to sort annotationGuide by geneSymbol
        ## then export finalColors and allGeneTypes
        geneTypes <- unique(processed$gene)
        finalColors <- annotationGuide$Color
        allGeneTypes <- annotationGuide$geneSymbol
        ## if we added "hyp" in annotation and it wasn't an autogenerated color
        if(length(which(geneTypes=="hyp"))>0 && length(which(allGeneTypes=="hyp"))<1) {
            allGeneTypes <- append(allGeneTypes, "hyp")
            finalColors <- append(finalColors, "#FFFFFF")
        }
        if(length(which(geneTypes=="mge"))>0 && length(which(allGeneTypes=="mge"))<1) {
            allGeneTypes <- append(allGeneTypes, "mge")
            finalColors <- append(finalColors, "#858585")
        }
        if(length(which(geneTypes=="other"))>0 && length(which(allGeneTypes=="other"))<1) {
            allGeneTypes <- append(allGeneTypes, "other")
            finalColors <- append(finalColors, "#DEDEDE")
        }  
##        if(length(which(geneTypes=="tal"))>0 && length(which(allGeneTypes=="tal"))<1) {
##            allGeneTypes <- append(allGeneTypes, "tal")
##            finalColors <- append(finalColors, "#000000")
##        }
        alphaGeneTypes <- sort(allGeneTypes)
        for (k in 1:length(alphaGeneTypes)) {
            geneIdx <- which(allGeneTypes==alphaGeneTypes[k])
            if (k == 1) {
                alphaColors <- finalColors[geneIdx]
            } else {
                alphaColors <- append(alphaColors, finalColors[geneIdx])
            }
        }    
        allGeneTypes <- alphaGeneTypes
        finalColors <- alphaColors
        print("Palettes generated.")
    } else {
        ## this will color genes according the the genes you specified
        ## in its default form it has distinct colors for hypothetical proteins, unspecified proteins, and mobile genetic elements
        ## it will also add a 2nd light/dark color shade if there are > 10 named genes
        ## standardizing on a minimal palette set so that all get colorRamped similarly
        ## fish!
        if (colorType == "fishualize") {
            if(requireNamespace("fishualize", quietly=TRUE)) {
                tempPalette <- fishualize::fish(5,option=paletteInput)
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "ghibli") {
            if(requireNamespace("ghibli", quietly=TRUE)) {
                tempPalette <- ghibli::ghibli_palette(paletteInput)[1:7]
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "lisa") { 
            if(requireNamespace("lisa", quietly=TRUE)) {
                tempPalette <- lisa::lisa_palette(paletteInput, 5)[1:5]
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "nord") {
            if(requireNamespace("nord", quietly=TRUE)) {
                tempPalette <- nord::nord(paletteInput,5)      
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "rtist") {
            if(requireNamespace("rtist", quietly=TRUE)) {
                tempPalette <- rtist::rtist_palette(paletteInput,5)[1:5]      
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "scico") {
            if(requireNamespace("scico", quietly=TRUE)) {
                tempPalette <- scico::scico(5,palette=paletteInput)      
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "wesanderson") {
            if(requireNamespace("wesanderson", quietly=TRUE)) {
                tempPalette <- wesanderson::wes_palette(paletteInput,5)      
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "viridis") {tempPalette <- viridis::viridis(5,option=paletteInput)}
        ## now tailoring palettes to our gene set
        geneTypes <- unique(processed$gene)
        colorNum <- length(geneTypes)
        geneTypes <- geneTypes[stringr::str_order(geneTypes)]
        allGeneTypes <- geneTypes
        notMe <- which(geneTypes =="hyp")
        notMeEither <- which(geneTypes == "mge")
        norMe <- which(geneTypes == "other")
        if (any(notMe>0, notMeEither>0, norMe>0)) {
            geneTypes <- geneTypes[-c(notMe,notMeEither,norMe)]
        }
        numGeneTypes <- length(geneTypes)
        halfGeneTypes <- numGeneTypes
        if(numGeneTypes >= 10) {
            halfGeneTypes <- ceiling(numGeneTypes/2)
        }
        if (halfGeneTypes == numGeneTypes) {
            ## keep a reasonably bright version if possible
            geneColors <- colorRampPalette(tempPalette)(numGeneTypes)
        } else if (halfGeneTypes != numGeneTypes) {
            if (colorType == "ghibli") {
                ## it has built-in medium, light, and dark options
                paletteDark <- gsub("Medium","Dark",paletteInput)
                paletteDark <- ghibli::ghibli_palette(paletteDark)
                paletteLight <- gsub("Medium","Light",paletteInput)
                paletteLight <- ghibli::ghibli_palette(paletteLight)
                darkColors <- colorRampPalette(paletteDark)(halfGeneTypes)
                lightColors <- colorRampPalette(paletteLight)(halfGeneTypes)
                geneColors <- c(darkColors, lightColors)
            } else {
                ## for the rest we roll our own, i guess
                ## these are sorta arbitrary functions
                geneColors <- colorRampPalette(tempPalette)(halfGeneTypes)
                darkColors <- rgb2hsv(col2rgb(geneColors))
                darkColors["s",] <- abs(0.7*darkColors["s",])
                darkColors["v",] <- abs(darkColors["v",] - 0.4*(darkColors["v",]))
                darkColors <- apply(darkColors,2,function(x) hsv(x[1], x[2], x[3]))
                lightColors <- rgb2hsv(col2rgb(geneColors))
                lightColors["s",] <- abs(0.4*lightColors["s",])
                lightColors["v",] <- abs(lightColors["v",] + 0.7*(1-lightColors["v",]))
                lightColors <- apply(lightColors,2,function(x) hsv(x[1], x[2], x[3]))
                geneColors <- c(darkColors, lightColors)
            }
        }       
        ## this specifically codes (unidentified) hypothetical proteins, mobile genetic elements, and known proteins not on the annotation list
        ## and reinserts them into the palette, so now it doesn't sacrifice pretty colors for them either
        ## as a default: hyp. = white, MGE = dark grey, other = middle gray
        ## update, i have done this in a dumb and ludicrously overcomplicated way and should trim it someday
        if (length(notMe)>=1) {
            ## if we have hypothetical proteins
            if (notMe == 1) {
                tempColors <- "#FFFFFF"
            } else  {
                tempColors <- geneColors[1:(notMe-1)]
                tempColors <- append(tempColors, "#FFFFFF")
            }
            ## ok if we have mge & other as well
            if ((length(notMeEither)>=1) && (length(norMe)>=1)) {
                ## print("case A")
                ## if the mge is right after hyp vs. if it's further down
                if (notMe == (notMeEither-1))  {
                    tempColors <- append(tempColors, "#858585")
                } else {
                    temp2Colors <- geneColors[(notMe):(notMeEither-2)]
                    tempColors <- append(tempColors, temp2Colors)
                    tempColors <- append(tempColors, "#858585")
                }
                ## if "other" is right after MGE, vs. if its further down
                if (notMeEither == (norMe-1)) {
                    tempColors <- append(tempColors, "#DEDEDE")
                } else {
                    temp3Colors <- geneColors[(notMeEither-1):(norMe-3)]
                    tempColors <- append(tempColors, temp3Colors)            
                    tempColors <- append(tempColors, "#DEDEDE")
                }
                if (geneColors[(norMe-3)] != geneColors[length(geneColors)] && length(geneColors > 1)) {
                    temp4Colors <-geneColors[(norMe-2):length(geneColors)]
                    tempColors <- append(tempColors, temp4Colors)
                } else if (length(geneColors) == 1 && notMe == 1 && length(geneColors)>0) {
                    ## if we hit this point, we already assigned hyp, mge, and other, and we didn't assign a color, let's add it
                    if (length(tempColors) < 4) {
                        tempColors <- append(tempColors, geneColors)  
                    }  
                }
            } else if (length(notMeEither)>=1) {
                ## print("case B")
                ## if there is an MGE (and a hyp) but not an other
                ## if it's right after hyp vs. if it's further down
                if (notMe == (notMeEither-1)) {
                    tempColors <- append(tempColors, "#858585")
                } else {
                    temp5Colors <- geneColors[(notMe):(notMeEither-2)]
                    tempColors <- append(tempColors, temp5Colors)            
                    tempColors <- append(tempColors, "#858585")
                }
                if (geneColors[(notMeEither-2)] != geneColors[length(geneColors)] && length(geneColors) > 1) {
                    temp6Colors <-geneColors[(notMeEither-1):length(geneColors)]
                    tempColors <- append(tempColors, temp6Colors)
                } else if (length(geneColors) == 1 && notMe == 1 && length(geneColors)>0) {
                    ## if we hit this point, we assigned hyp and mge and not the one color
                    if (length(tempColors) < 3) {
                        tempColors <- append(tempColors, geneColors)  
                    }    
                }             
            } else if (length(norMe)>=1) {
                ## print("case C")
                ## if there is an MGE (and a hyp) but not an other
                ## if it's right after hyp vs. if it's further down
                if (notMe == (norMe-1)) {
                    tempColors <- append(tempColors, "#DEDEDE")
                } else {
                    temp7Colors <- geneColors[(notMe):(norMe-2)]
                    tempColors <- append(tempColors, temp7Colors)            
                    tempColors <- append(tempColors, "#DEDEDE")
                }
                if (geneColors[(norMe-2)] != geneColors[length(geneColors)] && length(geneColors) > 1) {
                    temp8Colors <-geneColors[(norMe-1):length(geneColors)]
                    tempColors <- append(tempColors, temp8Colors)
                } else if (length(geneColors) == 1 && notMe == 1 && length(geneColors)>0) {
                    ## if we hit this point, we assigned hyp and other and not the one color
                    if (length(tempColors) < 3) {
                        tempColors <- append(tempColors, geneColors)  
                    }  
                }                 
            } else {
                ## print("case D")
                        ## we've only got hypotheticals
                if (notMe == 1) {
                    ## we can only have assigned hyp so far
                    if (length(geneColors)>0) {
                        tempColors <- append(tempColors, geneColors)
                    }
                } else {
                ## we have already assigned a color + hyp
                    if (length(geneColors)>1) {
                        temp9Colors <-geneColors[(notMe):length(geneColors)]
                        tempColors <- append(tempColors, temp9Colors)
                    }
                }    
            }
        } else {
            ## if we have entered a brave new world where there are no hypothetical proteins
            ## if we have both mge + other
            if ((length(notMeEither)>=1 && length(norMe)>=1)) {
                ## print("case E")
                ## if MGE is first
                if (notMeEither == 1) {
                    tempColors <- "#858585"
                } else  {
                    tempColors <- geneColors[1:(notMeEither-1)]
                    tempColors <- append(tempColors, "#858585")
                }
                ## if MGE is immediately followed by other
                if (notMeEither == (norMe-1)) {
                    tempColors <- append(tempColors, "#DEDEDE")
                } else {
                    temp2Colors <- geneColors[(notMeEither):(norMe-2)]
                    tempColors <- append(tempColors, temp2Colors)            
                    tempColors <- append(tempColors, "#DEDEDE")
                }
                if (geneColors[(norMe-2)] != geneColors[length(geneColors)] && length(geneColors)>1) {
                    temp3Colors <- geneColors[(norMe-1):length(geneColors)]
                    tempColors <- append(tempColors, temp3Colors)
                } else if (length(geneColors) == 1 && notMeEither == 1) {
                    ## if we hit this point, we assigned mge and other and not the one color
                    if (length(tempColors) < 3) {
                        tempColors <- append(tempColors, geneColors)  
                    }    
                }   
            } else if (length(notMeEither)>=1) {
                ## print("case F")
                ## if there is an MGE but no other
                if (notMeEither == 1) {
                    tempColors <- "#858585"
                    tempColors <- append(tempColors, geneColors)
                } else  {
                    tempColors <- geneColors[1:(notMeEither-1)]
                    tempColors <- append(tempColors, "#858585")
                    if (geneColors[(notMeEither-1)] != geneColors[length(geneColors)] && length(geneColors)>1) {
                        temp3Colors <- geneColors[(notMeEither):length(geneColors)]
                        tempColors <- append(tempColors, temp3Colors)
                    } else if (length(geneColors) == 1 && notMeEither == 1) {
                    ## if we hit this point, we assigned mge and other and not the one color
                        if (length(tempColors) < 2) {
                            tempColors <- append(tempColors, geneColors)  
                        }    
                    }   
                }              
            } else if (length(norMe)>=1) {
                ## print("case G")
                ## if there is not an MGE or hyp, but there is still "other"
                if (norMe == 1) {
                    tempColors <- "#DEDEDE"
                    tempColors <- append(tempColors, geneColors)
                } else  {
                    tempColors <- geneColors[1:(norMe-1)]
                    tempColors <- append(tempColors, "#DEDEDE")
                    if (geneColors[(norMe-1)] != geneColors[length(geneColors)] && length(geneColors)>1) {
                        temp3Colors <- geneColors[(norMe):length(geneColors)]
                        tempColors <- append(tempColors, temp3Colors)
                    } else if (length(geneColors) == 1 && norMe == 1) {
                    ## if we hit this point, we assigned mge and other and not the one color
                        if (length(tempColors) < 2) {
                            tempColors <- append(tempColors, geneColors)  
                        }    
                    }   
                }              
            } else {
                ## print("case H")
                ## if literally everything is annotated somehow
                tempColors <- geneColors
            }
        }
        finalColors <- tempColors
        print("Palettes generated.")
    }
    scale_fill_genes <- function(...) {
            ggplot2:::manual_scale(
                          'fill', 
                          values = setNames(finalColors, allGeneTypes), 
                          ...)
    }  
    ## if markclusters is chosen, this will add cluster numbers to the sequence IDs
    ## note that this adds the sequence number (ordered by clustering) and the cluster number
    ## these are generated in analyzeNeighbors (order is, at least; cluster number can be generated manually too)
    ## so it's also gonna sort the sequences by cluster order, making it easy to match with the heatmap
    if (markClusters == TRUE && standAlone == FALSE) {
        clustNumTemp <- unique(processed$clustNum)
        if (length(which(is.na(clustNumTemp)))>0) {
            clustNumTemp <- sort(clustNumTemp[-which(is.na(clustNumTemp))])
        } else {
            clustNumTemp <- sort(clustNumTemp)
        }    
        clustActualColors <- viridis::viridis(length(clustNumTemp))
        clustColors <- as.data.frame(list("clustNum"=clustNumTemp, "clustColors"=clustActualColors))
        processed$clustNum[which(is.na(processed$clustNum))]<-"none"
        mcGenes <- unique(processed$source_gene_oid)
        for (i in 1:length(mcGenes)) {
            tempSpeciesLoc <- grep(mcGenes[i], processed$source_gene_oid)
            tempLabel <- processed$bgc[tempSpeciesLoc[1]]
            if (showScaffold == TRUE) {
                tempLabel <- paste(tempLabel, ", scaffold ", processed$scaffoldID[tempSpeciesLoc[1]],sep="")
            }    
            tempSeqID <- as.numeric(processed$clustOrd[tempSpeciesLoc[1]]) + 100000
            tempSeqID <- as.character(tempSeqID)
            tempClustID <- processed$clustNum[tempSpeciesLoc[1]]
            tempColorID <- clustColors$clustColor[which(clustColors$clustNum %in% tempClustID)]
            if (any(grepl("none",tempClustID),is.na(tempClustID)) == TRUE)  {
                extName <- paste("<span style = 'font-size:9pt; color:#9F9F9F'>sequence ",
                                 tempSeqID,
                                 "</span><br>",
                                 "<span style = 'font-size:9pt; color:#000000'>",
                                 tempLabel,
                                 "</span>",
                                 sep="") 
            } else {
                ## this will match coloring with the coloring from the heatmap in analyzeNeighbors
                extName <- paste("<span style = 'font-size:9pt; color:#9F9F9F'>sequence ",
                                 tempSeqID,
                                 ", </span><span style = 'font-size:9pt; color:",
                                 tempColorID,
                                 "'>**cluster ",
                                 numbers2words(as.numeric(tempClustID)),
                                 "**</span><br>",
                                 "<span style = 'font-size:9pt; color:#000000'>",
                                 tempLabel,
                                 "</span>",
                                 sep="")
            }
            processed$bgc[tempSpeciesLoc] <- extName 
        }
        ## while we're at it, let's make sure we're sorted for export
        ## in this case, we're going by clustOrd (i.e. display order) and then gene_oid
        ## and uh let's make sure that clustOrd is an integer here
        processed$clustOrd <- as.numeric(processed$clustOrd)
        processed <- processed[order(processed$clustOrd, processed$gene_oid),]
                                        #    clustOrdNums <- unique(processed$clustOrd)
                                        #    clustOrdNums <- clustOrdNums[order(clustOrdNums)]  
                                        #    processed$fact <- factor(processed$clustOrd, levels=clustOrdNums))
        print("Cluster number labels added to diagram, and diagrams ordered by cluster order.")
    } else {
        ## if markClusters is not chosen, there's no fancy coloring, and we'll stick with the default sorting by cluster name - which is gonna default to the genome name
        ## but we might still want to add the scaffold
        ## what unique bgc regions are there?  
        mcGenes <- unique(processed$source_gene_oid)
        for (i in 1:length(mcGenes)) {
            tempSpeciesLoc <- grep(mcGenes[i], processed$source_gene_oid)
            tempLabel <- processed$bgc[tempSpeciesLoc[1]]
            if (showScaffold == TRUE) {
                tempLabel <- paste(tempLabel, ", scaffold ", processed$scaffoldID[tempSpeciesLoc[1]],sep="")
            }    
            ## this'll give the same font size etc. as the fancier cluster-colored version
            extName <- paste("<span style = 'font-size:9pt; color:#000000'>", tempLabel, "</span>",sep="")
            processed$bgc[tempSpeciesLoc] <- extName 
        }
        ## while we're at it, let's make sure we're sorted for export
        ## in this case, we're going by BGC (i.e. display order) and then gene_oid
        processed <- processed[order(processed$bgc, processed$gene_oid),]
                                        #      genomeNames <- unique(processed$genome)
                                        #      genomeNames <- genomeNames[order(genomeNames)]
                                        #      processed$fact <- factor(processed$genome, levels=genomeNames))
        print("Diagrams ordered by species name.")
    }
    ## let's save our data
    write.table(processed, annotationFileName, row.names=FALSE,sep="\t", quote=FALSE)  
    ## this function defaults to making a single scale (with 1 kb tickmarks)
    ## by dint of making a fake BGC
    ## there has to be a better way to do this, but messing with facets is screwing up gene alignment
    if (makeScale == TRUE) {
        num1k <- floor(max(as.double(processed$end))/1000)
        fakeSource <- as.character(9000000000 + num1k/2)
        fakeOrd <- as.character(max(as.numeric(processed$clustOrd))+1)
        fakeClust <- "0"
        ## assuming that a given dataset has a gene type has caused problems before
        ## this is overkill but ensures that we don't do that
        ## but starts with the color schemes least likely to cause surprising random problems
        ## this shouldn't even be visible since these are faux-tick marks, not legit genes, but.
        if (any(grep("other",unique(processed$gene))==TRUE)) {
            fakeType <- "other"
        } else if (any(grep("hyp",unique(processed$gene))==TRUE)) {  
            fakeType <- "hyp"
        } else if (any(grep("mge",unique(processed$gene))==TRUE)) {  
            fakeType <- "mge"
        } else {
            fakeType <- unique(processed$gene)[1]
        }
        if (markClusters == TRUE) {
            for (m in 1:num1k) {
                fakeNum <- as.double(9000000000 + m - 1)
                fakeStart <- as.double(0 + (m-1)*1000)
                fakeStop <- as.double(1 + (m-1)*1000)
                tempRow <- c(fakeNum, "","","scale (kb)",fakeType,fakeStart, fakeStop,"+","","","","","","","",fakeOrd, fakeClust, fakeSource, "1", "scale (kb)")
                processed <- rbind(processed, tempRow)
            }
        } else {
             for (m in 1:num1k) {
                fakeNum <- as.double(9000000000 + m - 1)
                fakeStart <- as.double(0 + (m-1)*1000)
                fakeStop <- as.double(1 + (m-1)*1000)
                tempRow <- c(fakeNum, "","","scale (kb)",fakeType,fakeStart, fakeStop,"+","","","","","","","",fakeSource, "1", "scale (kb)")
                processed <- rbind(processed, tempRow)
            }           
        }    
        processed$gene_oid <- as.double(processed$gene_oid)
        processed$start <- as.double(processed$start)
        processed$end <- as.double(processed$end)
        processed$direction <- as.double(processed$direction)
    }
    ## actually making pretty gene clusters
    ## labelGenes establishes whether genes are, well, labeled (can get busy on big sets of clusters)
    ## alignToCore establishes whether your gene of interest centers clusters (otherwise you get left align)
    ## X-axes are currently all or nothing - superfluous when things are scaled, but you do want one for publication.
    ## so a version with and without the axes for all clusters is also generated as a .pdf so that you can have a properly scaled vector scale
    ## anyway starting with the barebones diagrams
    if (alignToCore == TRUE) {
        ## note that where there are 2+ GOIs, this just centers them (it's looking for coreGeneName)
        dummies <- gggenes::make_alignment_dummies(processed, ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$bgc, id = .data$gene), on = coreGeneName)
    }
    if (labelGenes == TRUE ) {
        if (alignToCore == TRUE) {
            clusterDiagram <- ggplot2::ggplot(processed, ggplot2::aes(xmin = .data$start,
                                                                      xmax = .data$end,
                                                                      y = .data$bgc,
                                                                      fill = .data$gene,
                                                                      forward = .data$direction,
                                                                      label="gene_oid")) +
                gggenes::geom_gene_arrow(arrowhead_height=ggplot2::unit(3,"mm"), arrowhead_width=ggplot2::unit(1,"mm")) +
                ggplot2::geom_blank(data=dummies, ggplot2::aes(forward = 1))
        } else  {
            clusterDiagram <- ggplot2::ggplot(processed, ggplot2::aes(xmin = .data$start,
                                                                      xmax = .data$end,
                                                                      y = .data$bgc,
                                                                      fill = .data$gene,
                                                                      forward = .data$direction,
                                                                      label="gene_oid")) +
                gggenes::geom_gene_arrow(arrowhead_height=ggplot2::unit(3,"mm"), arrowhead_width=ggplot2::unit(1,"mm"))
        }
    }
    if (labelGenes == FALSE ) {
        if (alignToCore == TRUE) {      
            clusterDiagram <- ggplot2::ggplot(processed, ggplot2::aes(xmin = .data$start,
                                                                      xmax = .data$end,
                                                                      y = .data$bgc,
                                                                      fill = .data$gene,
                                                                      forward = .data$direction)) +
                gggenes::geom_gene_arrow(arrowhead_height=ggplot2::unit(3,"mm"), arrowhead_width=ggplot2::unit(1,"mm")) +
                ggplot2::geom_blank(data=dummies, ggplot2::aes(forward = 1))
        } else  {
            clusterDiagram <- ggplot2::ggplot(processed, ggplot2::aes(xmin = .data$start,
                                                                      xmax = .data$end,
                                                                      y = .data$bgc,
                                                                      fill = .data$gene,
                                                                      forward = .data$direction)) +
                gggenes::geom_gene_arrow(arrowhead_height=ggplot2::unit(3,"mm"), arrowhead_width=ggplot2::unit(1,"mm"))       
        }
    }
    ## adding details - colors, captions, etc.
    clusterDiagram <- clusterDiagram +
        ggplot2::facet_wrap(~ .data$bgc, ncol = 1, scales = "free") +
        ggplot2::scale_fill_manual(values = finalColors) +
        ggplot2::labs(title=figureTitle, caption=captionText, y="**genomic neighborhoods**", x=NULL) +
        gggenes::theme_genes()
    ## defining the diagram with axes
    xPlusDiagram <- clusterDiagram +
        ggplot2::theme(plot.title = ggtext::element_textbox_simple(size=14,
                                                                   padding=ggplot2::margin(10,10,10,10),
                                                                   margin=ggplot2::margin(5,0,5,0),
                                                                   r=ggplot2::unit(8,"pt"),
                                                                   fill="#EEEEEE",
                                                                   halign=0.5,
                                                                   valign=0.5,
                                                                   width = grid::unit(0.75, "npc")),
                       plot.title.position = "plot",
                       axis.title.y = ggtext::element_markdown(size=12),
                       axis.text.y = ggtext::element_markdown(size=9),
                       plot.caption = ggtext::element_textbox_simple(lineheight=1.2,
                                                                     size=11,
                                                                     padding=ggplot2::margin(10,10,10,10),
                                                                     margin=ggplot2::margin(5,0,5,0),
                                                                     r=ggplot2::unit(8,"pt"),
                                                                     fill="#D9D9D9",
                                                                     halign=0.5,
                                                                     valign=0.5,
                                                                     width = grid::unit(0.75, "npc")),
                       plot.caption.position = "plot")
    ## if the ratio of BGCs to gene types is lower than a cutoff, we probably need to make sure the legend has columns for legibility
    if (length(unique(processed$bgc))/length(finalColors) >= 0.5) {
        xPlusDiagram <- xPlusDiagram + ggplot2::guides(fill=ggplot2::guide_legend(ncol=1, bycol=TRUE)) 
    } else {
        xPlusDiagram <- xPlusDiagram + ggplot2::guides(fill=ggplot2::guide_legend(ncol=2, bycol=TRUE)) 
    }
    ## get rid of the extra axes (but means no scale on-figure, for now)
    xMinusDiagram <- xPlusDiagram + ggplot2::theme(axis.line.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),  axis.ticks.x=ggplot2::element_blank())
    ## this gets rid of that ggplot2 grid error by scaling sizing to the number of gene clusters being visualized
    ## also adapts width if you have chosen the longest names
    print("Cluster diagrams generated.")
    ## some rough guidelines to come up with vaguely reasonable dimensions
    modheight <- length(uniqueBGCs)/2
    if (exists(x="showScaffold") && showScaffold==TRUE) {
        modwidth <- modheight/2
    } else if (length(uniqueBGCs)/2 <=50) {
        modwidth <- modheight  
    } else {  
        modwidth <- modheight/3
    }
    ## handling really small genesets that will otherwise cause grid.call errors
    if (length(uniqueBGCs) < 5) {
        modwidth <- 15
        if (modheight < 5) { modheight <- 15 }
    }
        
    ## exports a .pdf and a .png version
    if (everyScale == TRUE) {
        ggplot2::ggsave(file=finalpdfname, plot=xPlusDiagram, device="pdf", height=modheight, width=modwidth, limitsize=FALSE)
        ggplot2::ggsave(file=finalpngname, plot=xPlusDiagram, device="png", height=modheight, width=modwidth, limitsize=FALSE)
    }
    ggplot2::ggsave(file=finalpdfnameX, plot=xMinusDiagram, device="pdf", height=modheight, width=modwidth, limitsize=FALSE)
    ggplot2::ggsave(file=finalpngnameX, plot=xMinusDiagram, device="png", height=modheight, width=modwidth, limitsize=FALSE)
    print("Diagrams for complete set of gene clusters generated.")
    ## this way you have a copy of the annotated dataset and can skip auto-annotation if desired 
    if (subclusterDiagrams == TRUE)   {
        scaleBar <- processed %>% dplyr::filter(.data$bgc == "scale (kb)")
        forClusters <- processed %>% dplyr::filter(.data$bgc != "scale (kb)")
        clustList <- unique(forClusters$clustNum)
        ## taking it one cluster at a time
        for (i in 1:length(clustList))  {
            clustFileNamePDF <- paste(fileName, "_cluster_",clustList[i],".pdf",sep="")
            processedCluster <- forClusters %>% dplyr::filter(.data$clustNum == clustList[i])         
            ## let's update the title and annotations to apply only to this cluster
            subFigureTitle <- paste("<b>genomic neighborhood diagrams for *",coreGeneName,"*, cluster ",clustList[i],"</b> (\u00b1", neighborNumber," genes displayed)",sep="")
            if (annotateGenes == TRUE)  {
                uniqueClusterBGCs <- unique(processedCluster$bgc)
                subNamedGenes <- namedGenes
                numNamedGenes <- length(subNamedGenes$gene)
                wrapGenes <- ceiling(numNamedGenes/2)
                subNamedGenes$num <- 0
                subNamedGenes$percent <- 0
                for (j in 1:length(uniqueClusterBGCs))  {
                    bgcLoc <- which(processedCluster$bgc == uniqueClusterBGCs[j])
                    for (m in 1:length(subNamedGenes$geneSymbol))  {
                        hitNum <- length(grep(subNamedGenes$geneSymbol[m], processedCluster$gene[bgcLoc]))
                        if (hitNum >= 1) {
                            subNamedGenes$num[m] <- subNamedGenes$num[m] + 1
                        } else {
                            subNamedGenes$num[m] <- subNamedGenes$num[m] + 0
                        }
                    }
                }
                subCaptionText <- paste("**gene abundance in cluster ",clustList[i],":** ",sep="")
                for (k in 1:numNamedGenes) {
                    if ( subNamedGenes$geneSymbol[k] == "hyp" || subNamedGenes$geneSymbol[k] == "mge") {
                        next
                    } else  {
                        if (k != numNamedGenes) {
                            subNamedGenes$percent[k] <- round(subNamedGenes$num[k] / length(uniqueClusterBGCs), digits=2)
                            printPercent <- subNamedGenes$percent[k]*100
                            tempGeneSymbol <- subNamedGenes$geneSymbol[k]
                            subCaptionText <- paste(subCaptionText,"*",tempGeneSymbol,"* (",printPercent,"%), ",sep="")
                        } else  {
                            subNamedGenes$percent[k] <- round(subNamedGenes$num[k] / length(uniqueClusterBGCs), digits=2)
                            printPercent <- subNamedGenes$percent[k]*100
                            tempGeneSymbol <- subNamedGenes$geneSymbol[k]
                            subCaptionText <- paste(subCaptionText,"*",tempGeneSymbol,"* (",printPercent,"%).",sep="")
                        }
                    }
                }    
            } else {
                subCaptionText <- "No automated annotation and gene quantification for these diagrams."
            }
            ## appending the scale bar?
            processedCluster <- rbind(processedCluster,scaleBar)
            ## let's trim the palette for what actually ends up being visualized here 
            ## not trimming has weird consequences
            ## (this is incidentally why we had to add the scalebar back now)
            someGeneTypes <- allGeneTypes[which(allGeneTypes %in% processedCluster$gene)]
            someColors <- finalColors[which(allGeneTypes %in% processedCluster$gene)]
            scale_fill_genes2 <- function(...) {
                ggplot2:::manual_scale(
                          'fill', 
                          values = setNames(someColors, someGeneTypes), 
                          ...)
            }
            if (alignToCore == TRUE) {      
                subDummies <- gggenes::make_alignment_dummies(processedCluster,
                                                              ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$bgc, id = .data$gene), on = coreGeneName)
                subclusterDiagram <- ggplot2::ggplot(processedCluster,
                                                     ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$bgc, fill = .data$gene, forward = .data$direction)) +
                    gggenes::geom_gene_arrow(arrowhead_height=ggplot2::unit(3,"mm"), arrowhead_width=ggplot2::unit(1,"mm")) +
                    ggplot2::geom_blank(data=subDummies, ggplot2::aes(forward = 1))
            } else  {
                subclusterDiagram <- ggplot2::ggplot(processedCluster,
                                                     ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$bgc, fill = .data$gene, forward = .data$direction)) +
                    gggenes::geom_gene_arrow(arrowhead_height=ggplot2::unit(3,"mm"), arrowhead_width=ggplot2::unit(1,"mm"))
            }
            subclusterDiagram <- subclusterDiagram +
                ggplot2::facet_wrap(~ .data$bgc, ncol = 1, scales = "free") +
                scale_fill_genes2() +
                ggplot2::labs(title=subFigureTitle, caption=subCaptionText, y="**genomic neighborhoods**", x=NULL) +
                gggenes::theme_genes() +
                ggplot2::theme(axis.line.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),  axis.ticks.x=ggplot2::element_blank())
            subModHeight <- (length(unique(processedCluster$bgc)) + 2)/2
            if (length(unique(processedCluster$bgc)) > 10) {
                subclusterDiagram <- subclusterDiagram + ggplot2::theme(plot.title = ggtext::element_textbox_simple(size=14,
                                                                           padding=ggplot2::margin(10,10,10,10),
                                                                           margin=ggplot2::margin(5,0,5,0),
                                                                           r=ggplot2::unit(8,"pt"),
                                                                           fill="#EEEEEE",
                                                                           halign =0.5,
                                                                           valign =0.5,
                                                                           width = grid::unit(0.75, "npc")),
                               plot.title.position = "plot",
                               axis.title.y = ggplot2::element_blank(),
                               axis.text.y = ggtext::element_markdown(size=9),
                               plot.caption = ggtext::element_textbox_simple(lineheight=1.2,
                                                                             size=11,
                                                                             padding=ggplot2::margin(10,10,10,10),
                                                                             margin=ggplot2::margin(5,0,5,0),
                                                                             r=ggplot2::unit(8,"pt"),
                                                                             fill="#D9D9D9",
                                                                             halign=0.5,
                                                                             valign=0.5,
                                                                             width=grid::unit(0.75, "npc")),
                               plot.caption.position = "plot") +
                               ggplot2::guides(fill=ggplot2::guide_legend(ncol=1, bycol=TRUE)) 
                } else {
                    subclusterDiagram <- subclusterDiagram + ggplot2::theme(plot.title = ggtext::element_textbox_simple(size=14,
                                                                           padding=ggplot2::margin(10,10,10,10),
                                                                           margin=ggplot2::margin(5,0,5,0),
                                                                           r=ggplot2::unit(8,"pt"),
                                                                           fill="#EEEEEE",
                                                                           halign=0.5,
                                                                           valign=0.5,
                                                                           width=grid::unit(0.75, "npc")),
                               plot.title.position = "plot",
                               axis.title.y = ggtext::element_markdown(size=12),
                               axis.text.y = ggtext::element_markdown(size=9),
                               plot.caption = ggtext::element_textbox_simple(lineheight=1.2,
                                                                             size=11,
                                                                             padding=ggplot2::margin(10,10,10,10),
                                                                             margin=ggplot2::margin(5,0,5,0),
                                                                             r=ggplot2::unit(8,"pt"),
                                                                             fill="#D9D9D9",
                                                                             halign=0.5,
                                                                             valign=0.5,
                                                                             width=grid::unit(0.75, "npc")),
                               plot.caption.position = "plot") +
                               ggplot2::guides(fill=ggplot2::guide_legend(ncol=2, bycol=TRUE)) 
                }
            ggplot2::ggsave(file=clustFileNamePDF, plot=subclusterDiagram, device="pdf", height=subModHeight, width=modwidth, limitsize=FALSE)
        }
        print("Diagrams for individual subgroups of gene clusters generated.")
    }   
    print("Pretty gene clusters illustrated.")
    ## they'd better be pretty after all of this
    return(processed)
}
