#' Umbrella function generates pretty gene cluster diagrams
#'
#' This function generates gene cluster diagrams incorporating annotation and clustering info from several sources.
#' @param imgGenes What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")
#' @param imgNeighbors What is the file with the metadata for neighbors of your gene of interest? Filename as string ("filename.txt")
#' @param geneFormat User-defined annotation fine, assigning gene names via protein families.  Filename as string ("filename.txt")
#' @param geneName What is the name of your gene? Gene name as string ("genE")
#' @param efiRepnodes Are values here post-EFI-EST? T/F value, defaults to FALSE
#' @param neighborNumber How many genes are expected on either side of the gene of interest? Integer (e.g. 10)
#' @param annotateGenes Should the genes be auto-annotated using the geneFormat file? T/F, defaults to TRUE.
#' @param standAlone Should we assume the data has not been through the full pipeline? T/F, defaults to FALSE.
#' @param markClusters Should cluster labels be included in the final figure? T/F, defaults to FALSE.
#' @param autoColor If TRUE, generates palettes automatically. If FALSE, requires hex-formatted colors per gene in the last column of the geneFormat data. Defaults to TRUE.
#' @param colorType What palette family should be used? String, defaults to "nord"
#' @param paletteInput What specific palette within that family should be used? String, defaults to "aurora"
#' @param showScaffold Should the scaffold be included in name labels? T/F, defaults to FALSE
#' @param alignToCore Should the core gene be centered in the diagrams? T/F, defaults to TRUE
#' @param labelGenes Should all genes be labeled in the diagrams? T/F, defaults to FALSE
#' @return Gene diagrams as .png and .pdf files.
#' @export
#' @examples
#' prettyClusterOutput <- prettyClusterDiagrams(imgGenes = "genesFile.txt", imgNeighbors = "neighborsFile.txt", geneFormat = "formatFile.txt", geneName = "genE", neighborNumber = 10)
#' 
prettyClusterDiagrams <- function(imgGenes = imgGenes, imgNeighbors = imgNeighbors, geneFormat = geneFormat, geneName = geneName, efiRepnodes = FALSE, neighborNumber = neighborNumber, annotateGenes = TRUE, standAlone = FALSE, markClusters = FALSE, autoColor = TRUE, colorType = "nord", paletteInput = "aurora", showScaffold = FALSE, alignToCore=TRUE, labelGenes = FALSE) { 
  fileDate <- format(Sys.Date(),format="%Y%m%d")
  fileName <- paste(fileDate,"_prettyClusters_",geneName,sep="")
  finalpngname <- paste(fileName,"_with-axes.png",sep="")
  finalpdfname <- paste(fileName,"_with-axes.pdf",sep="")
  finalpdfnameX <- paste(fileName,"_no-axes.pdf",sep="")
  finalpngnameX <- paste(fileName,"_no-axes.png",sep="")
  annotationFileName <- paste(fileName, "_annotation.txt",sep="")
  imgCols <- list("gene_oid","Start.Coord","End.Coord", "Strand", "Gene.Symbol", "Scaffold.Name", "Scaffold.ID", "Genome.Name","Genome.ID","Pfam", "Tigrfam", "clustNum", "clustOrd","source_gene_oid")
  ggCols <- list("gene_oid", "start", "end", "strand","gene","scaffold", "scaffoldID", "genome","genomeID","Pfam", "Tigrfam", "clustNum", "clustOrd","source_gene_oid")
                                       # geneset and annotation processing
  geneFormat <- read.table(geneFormat, header=TRUE,  sep = "\t", stringsAsFactors = FALSE)
  geneFormat$Pfam[which(is.na(geneFormat$Pfam))] <- ""
  geneFormat$Tigrfam[which(is.na(geneFormat$Tigrfam))] <- ""
  geneFormat$Hypofam[which(is.na(geneFormat$Hypofam))] <- ""
  geneFormat$IMG.Term[which(is.na(geneFormat$IMG.Term))] <- ""
  geneFormat$Color[which(is.na(geneFormat$Color))] <- ""
  if(typeof(imgGenes) == "character") {
    imgGenes <-read.csv(imgGenes, header=TRUE, sep="\t", stringsAsFactors=FALSE )
  } else {
    imgGenes <- imgGenes
  }
  colnames(imgGenes)[1] <- "gene_oid"
  if(typeof(imgNeighbors) == "character")  {
    imgNeighbors <- read.csv(imgNeighbors, header=TRUE, sep="\t", stringsAsFactors=FALSE )
  } else {
    imgNeighbors <- imgNeighbors
  }
  colnames(imgNeighbors)[1] <- "gene_oid"
  if (efiRepnodes == TRUE) {
    coreGeneName <- geneName
    geneName <- paste(geneName, "_repnodes",sep="")
  } else {
    geneName <- geneName
    coreGeneName <- geneName
  }
  print("Opening data files.")
  ## img metadata conversion - filters out things other than the core info because why do we need it
  inputColNames <- colnames(imgNeighbors)
  inputColNames <- stringr::str_replace_all(inputColNames, c(" " = "." , "," = "" ))
  colnames(imgNeighbors) <- inputColNames
  geneSets <- imgNeighbors[names(imgNeighbors) %in% imgCols]
  geneSets <- geneSets %>% dplyr::rename(scaffold = Scaffold.Name)
  geneSets <- geneSets %>% dplyr::rename(scaffoldID = Scaffold.ID)
  geneSets <- geneSets %>% dplyr::rename(start = Start.Coord)
  geneSets <- geneSets %>% dplyr::rename(end = End.Coord)
  geneSets <- geneSets %>% dplyr::rename(strand = Strand)
  geneSets <- geneSets %>% dplyr::rename(gene = Gene.Symbol)
  geneSets <- geneSets %>% dplyr::rename(genome = Genome.Name)
  geneSets <- geneSets %>% dplyr::rename(genomeID = Genome.ID)
  geneSets <- data.table::as.data.table(geneSets, stringsAsFactors=FALSE)
  ## this is easier than dealing with different possible IMG-like metadata inputs
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
  ## dealing with blank cells & misc metadata processing
  geneSets$Pfam[geneSets$Pfam==""]<-"none"
  geneSets$Tigrfam[geneSets$Tigrfam==""]<-"none"
  geneSets$IMGfam[geneSets$IMGfam==""]<-"none"
  geneSets$Hypofam[geneSets$Hypofam==""]<-"none"
  geneSets$direction <- ifelse(geneSets$strand == "+", 1, -1)
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
    reorderGenes1 <- geneFormat %>% dplyr::filter(!Requirement == "and")
    reqGenes <- geneFormat %>% dplyr::filter(Requirement == "and")
    namedGenes <- rbind(reorderGenes1, reqGenes)
    ## bride of the slight kludge - fused genes _also_ write over previous annotations
    reorderGenes2 <- namedGenes %>% dplyr::filter(!Fusion == "yes")
    fusedGenes <- namedGenes %>% dplyr::filter(Fusion == "yes")
    namedGenes <- rbind(reorderGenes2, fusedGenes)
    namedGenes <- namedGenes %>% dplyr::rename(IMGfam = IMG.Term)
    ## "" or NA becomes "none"
    ## then none will become a list element
    namedGenes$Pfam <- namedGenes$Pfam %>% tidyr::replace_na("none")
    namedGenes$Tigrfam <- namedGenes$Tigrfam %>% tidyr::replace_na("none")
    namedGenes$Hypofam <- namedGenes$Hypofam %>% tidyr::replace_na("none")
    namedGenes$IMGfam <- namedGenes$IMGfam %>% tidyr::replace_na("none")
    ## son of the bride of the kludge
    if ("IMGfam" %in% colnames(namedGenes) == FALSE) {namedGenes$IMGfam <- "none"}
    if ("Hypofam" %in% colnames(namedGenes) == FALSE) {namedGenes$Hypofam <- "none"}
    ## the return of the son of the bride of the kludge
    namedGenes$Pfam[namedGenes$Pfam==""]<-"none"
    namedGenes$Tigrfam[namedGenes$Tigrfam==""]<-"none"
    namedGenes$IMGfam[namedGenes$IMGfam==""]<-"none"
    namedGenes$Hypofam[namedGenes$Hypofam==""]<-"none"
    annotList <- namedGenes$geneSymbol
    for (i in 1:length(geneSets$gene)) {
      foundMe <- "no"
      foundIn <- "nowhere"
      tempRealData <- geneSets[i,] %>% dplyr::select(Pfam, Tigrfam, IMGfam, Hypofam)
      ## entirely hypothetical proteins are easy to flag - no annotations
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
          ## first list is IMGfam and Hypofam - both are specific and singular, so any non-none match is legit
          tempList1 <- strsplit(namedGenes$IMGfam[j], " ")
          tempList1 <- append(tempList1, strsplit(namedGenes$Hypofam[j], " "))
          ## second list is Pfam and TIGRfam - here you can have multiple hits, and "none" can be meaningful (e.g. for superfamilies)
          tempList2 <- strsplit(namedGenes$Pfam[j], " ")
          tempList2 <- append(tempList2, strsplit(namedGenes$Tigrfam[j], " "))
                ## we want "none" to be an element but _only_ for "and", it should not be a match for "any"
          tempList1 <- unlist(tempList1[tempList1!="none"])
          tempList2 <- unlist(tempList2)
          tempList3 <- unlist(tempList2[tempList2!="none"])
          tempVector1 <- sapply(lapply(tempList1,grepl,tempRealData[,3:4]),any)
          tempVector2 <- sapply(lapply(tempList2,grepl,tempRealData[,1:2]),any)
          tempVector3 <- sapply(lapply(tempList3,grepl,tempRealData[,1:2]),any)   
          if (length(tempList1) != 0 && any(tempVector1)) {
            ## there are annotations and they match the most specific terms (IMG Term or HypoFam)
            geneSets$gene[i] <- namedGenes$geneSymbol[j]
            foundMe <- "yes"
            ## this is a strong match and generally should not be overwritten
            foundIn <- "exact"
            #print("exit A")
            next
          } else if (namedGenes$Requirement[j]=="all" && all(tempVector2) && length(tempList1) != 0 && any(tempVector1) != FALSE) {
            ## all annotations had to match the requirements and all did
            geneSets$gene[i] <- namedGenes$geneSymbol[j]
            foundMe <- "yes"
            ## this is a strong match and generally should not be ovewritten
            foundIn <- "exact"
            #print("exit B")
            next
          } else if (namedGenes$Requirement[j]=="all" && any(tempVector3)  && namedGenes$Fusion[j] != "yes") {
            ## all annotations had to match the requirement but only some did, and no fusion proteins were expected
            ## prevent a less specific match from overwriting a better one, if there is one
            if (foundMe != "yes" || foundIn != "exact") {
              geneSets$gene[i] <- paste(namedGenes$geneSymbol[j],"*",sep="")
              ## this is a weaker match and it's possible we could overwrite it
              foundMe <- "maybe"
              # print("exit C")
              next
            }
          } else if (namedGenes$Requirement[j] =="any" && length(tempList3) != 0 && any(tempVector3)) {
            ## there were annotations that matched (and any match was allowed)
            if (foundIn == "exact") {
              ## an exact match was already found, move on
              #print("exit D1")
              next
            } else {
              geneSets$gene[i] <- namedGenes$geneSymbol[j]
              ## probably we found it
              foundMe <- "yes"
              #print("exit D2")
              next
            }
          } else {
            if (j != length(annotList) && foundMe == "no") {
              ## no matches were found on this round but it isn't the last, so let's call it other for now
              geneSets$gene[i] <- "other"
              foundMe <- "maybe"
              #print("exit E1")
              next
            } else if (j == length(annotList) && foundMe == "no") {
              ## no matches were found and this is the last chance, so let's finalize on other
              geneSets$gene[i] <- "other"
              foundMe <- "yes"
              #print("exit E2")
              next
            }
          } 
        } 
      }
    }
    print("Genes automatically annotated based on the geneFormat file.")
  }
  ## this should be true independent of auto-annotation
  geneSets$gene[which(geneSets$gene_oid %in% imgGenes$gene_oid)] <- coreGeneName
  ## if this is being run as a standalone figure generator, we can't assume scaffolds have been trimmed to avoid wrong-scaffold neighbors
  if (standAlone == TRUE) {
    initScaff <- unique(geneSets$scaffoldID)
    for (m in 1:length(initScaff)) {
      scaffTest <- geneSets %>% dplyr::filter(scaffoldID == initScaff[m])
      if (any(grepl(coreGeneName,scaffTest$gene)) == TRUE) {
        next
      } else {
        geneSets <- geneSets %>% dplyr::filter(!scaffoldID == initScaff[m])
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
      dupeList <- grep(geneSets$genome[goiSeek][dupes[i]],geneSets$genome)
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
    ## if they are within 50 kb, they're considered contiguous (might be huge NRPS megagenes or something)
    ## if they are further, they get more analysis
    ## unless a big jump comes between 2 genes (25k!), a bgc end/beginning is not registered
    ## names are adjusted as necessary to add multiple BGCs and thus not break gggenes
  contigs <- unique(geneSets$scaffoldID)
  species <- unique(geneSets$genomeID)
  multiBGC <- data.frame()
  multicoord <- list()
  for(k in 1:length(species)) {
    tempSpecies <- geneSets %>% dplyr::filter(genomeID == species[k])
    specContigs <- unique(tempSpecies$scaffoldID)
    for (j in 1:length(specContigs))  {
      tempContig <- tempSpecies %>% dplyr::filter(scaffoldID == specContigs[j])
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
              if (between(tempContig$start[g], min(tempCoord[c(whichBGC*2-1,whichBGC*2)]), max(tempCoord[c(whichBGC*2-1,whichBGC*2)])) == TRUE) {
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
  bgcScaffolds <- unique(multiBGC$scaffoldID)    
  print("Gene cluster naming fixed in species with multiple clusters.")
  ## this normalizes positioning. 
  ## we do not fundamentally care that much about absolute genome position as a default here
  ## nucleotide position number is arbitrary to begin with and meaningless when comparing across wildly different species
  ## so let's make every cluster start at the beginning
  ## some of this may end up in flux as i try to nuke the extra x-axes.
  maxval <- 0
  for (i in 1:length(bgcScaffolds)) {
    tempPositioning <- dplyr::filter(multiBGC, scaffoldID == bgcScaffolds[i])
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
      ## oof look at that regex
      forGrep <-sub("\\(","\\\\\\(",finalUniqueBGCs[i])
      forGrep <-sub("\\)","\\\\\\)",forGrep)
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
    figureTitle <- paste("<b>genomic neighborhood diagrams for *",coreGeneName,"*</b>",sep="")
    print("Gene abundance calculated based on geneFormat-derived annotations.")
  } else {
    figureTitle <- paste("<b>genomic neighborhood diagrams for *",coreGeneName,"*</b>",sep="")
    captionText <- "no automated annotation and gene quantification for these diagrams."
  }
  uniqueBGCs <- unique(finalGeneSets$bgc)
  ## aligning gene clusters to your gene of interest
  ## without this, expect all clusters to be left-aligned, and in many situations, this is not preferred 
  if (alignToCore == TRUE) {
    # dirGeneSets <- NULL
    for (i in 1:length(uniqueBGCs)) {
      aligningGenes <- data.frame()
      aligningGenes <- dplyr::filter(finalGeneSets, bgc == uniqueBGCs[i])
      core <- which(aligningGenes$gene == coreGeneName)
      if (length(core)==0) {next}
      if (aligningGenes$direction[core] == 1) {
        if (exists(x="dirGeneSets") == FALSE) {
          dirGeneSets <- data.frame(aligningGenes, stringsAsFactors = FALSE)
        } else {
          dirGeneSets <- rbind(dirGeneSets, aligningGenes)
        }
      } else {
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
      rm(aligningGenes)
    }
    processed <- dirGeneSets
    print("Gene clusters aligned around the family of genes of interest.")
  } else {
    processed <- finalGeneSets
  }
  ## palette generation - needed or not?
  if (autoColor == FALSE) {
    finalColors <- geneFormat$Color
  } else {
  ## this will color genes according the the genes you specified
  ## in its default form it has distinct colors for hypothetical proteins, unspecified proteins, and mobile genetic elements
  ## it will also add a 2nd light/dark color shade if there are > 10 named genes
  ## standardizing on a minimal palette set so that all get colorRamped similarly
  ## fish!
    if (colorType == "fishualize") {tempPalette <- fishualize::fish(5,option=paletteInput)}
    if (colorType == "ghibli") {tempPalette <- ghibli::ghibli_palette(paletteInput)[1:7]}
    if (colorType == "lisa") { tempPalette <- lisa::lisa_palette(paletteInput, 5)[1:5]} 
    if (colorType == "nord") {tempPalette <- nord::nord(paletteInput,5)}
    if (colorType == "rtist") {tempPalette <- rtist::rtist_palette(paletteInput,5)[1:5]}
    if (colorType == "scico") {tempPalette <- scico::scico(5,palette=paletteInput)}
    if (colorType == "viridis") {tempPalette <- viridis::viridis(5,option=paletteInput)}
    if (colorType == "wesanderson") {tempPalette <- wesanderson::wes_palette(paletteInput,5)}
  ## now tailoring palettes to our gene set
    geneTypes <- unique(processed$gene)
    colorNum <- length(geneTypes)
    geneTypes <- geneTypes[stringr::str_order(geneTypes)]
    notMe <- which(geneTypes =="hyp")
    notMeEither <- which(geneTypes == "mge")
    norMe <- which(geneTypes == "other")
    geneTypes <- geneTypes[-c(notMe,notMeEither,norMe)]
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
    tempColors <- geneColors[1:(notMe-1)]
    tempColors <- append(tempColors, "#FFFFFF")
    if (notMe == (notMeEither-1))  {
      tempColors <- append(tempColors, "#888888")
    } else {
      temp2Colors <- geneColors[(notMe):(notMeEither-2)]
      tempColors <- append(tempColors, temp2Colors)
      tempColors <- append(tempColors, "#888888")
    }
    if (notMeEither == (norMe-1)) {
      tempColors <- append(tempColors, "#DEDEDE")
    } else {    
      temp3Colors <- geneColors[(notMeEither-1):(norMe-3)]
      tempColors <- append(tempColors, temp3Colors)
      tempColors <- append(tempColors, "#EEEEEE")
    }
    if (geneColors[(norMe-3)] != geneColors[length(geneColors)]) {
      temp4Colors <-geneColors[(norMe-2):length(geneColors)]
      tempColors <- append(tempColors, temp4Colors)
    }
    finalColors <- tempColors
    print("Palettes generated.")
  }
  ## if markclusters is chosen, this will add cluster numbers to the sequence IDs
  ## note that this adds the sequence number (ordered by clustering) and the cluster number
  ## these are generated in analyzeNeighbors (order is, at least; cluster number can be generated manually too)
  ## so it's also gonna sort the sequences by cluster order, making it easy to match with the heatmap
  if (markClusters == TRUE) {
    clustActualColors <- viridis::viridis(length(unique(processed$clustNum)))
    clustNumTemp <- unique(processed$clustNum)
    clustNumTemp[which(is.na(clustNumTemp))] <- 0
    clustNumTemp <- sort(clustNumTemp)
    clustColors <- as.data.frame(list("clustNum"=clustNumTemp, "clustColors"=clustActualColors))
    processed$clustNum[processed$clustNum==""]<-"(none)"
    mcGenes <- unique(processed$source_gene_oid)
    for (i in 1:length(mcGenes)) {
      tempSpeciesLoc <- grep(mcGenes[i], processed$source_gene_oid)
      tempLabel <- processed$bgc[tempSpeciesLoc[1]]
      tempSeqID <- as.numeric(processed$clustOrd[tempSpeciesLoc[1]]) + 100000
      tempSeqID <- as.character(tempSeqID)
      tempClustID <- processed$clustNum[tempSpeciesLoc[1]]
      tempColorID <- clustColors$clustColor[grep(tempClustID,clustColors$clustNum)]
      if (grepl("none",tempClustID) == TRUE)  {
        extName <- paste("<span style = 'font-size:9pt; color:#9F9F9F'>sequence ",tempSeqID,"</span><br>","<span style = 'font-size:9pt; color:#000000'>",tempLabel, "</span>",sep="")
      } else {
        ## this will match coloring with the coloring from the heatmap in analyzeNeighbors
        extName <- paste("<span style = 'font-size:9pt; color:",tempColorID,"'>sequence ",tempSeqID,", **cluster ",numbers2words(as.numeric(tempClustID)),"**</span><br>","<span style = 'font-size:9pt; color:#000000'>",tempLabel,"</span>",sep="")
      }
      for (j in 1:length(tempSpeciesLoc)) {
         processed$bgc[tempSpeciesLoc[j]] <- extName
      }
    }
    print("Cluster number labels added to diagram.")
  }
    ## actually making pretty gene clusters
    ## labelGenes establishes whether genes are, well, labeled (can get busy on big sets of clusters)
    ## alignToCore establishes whether your gene of interest centers clusters (otherwise you get left align)
    ## X-axes are currently all or nothing - superfluous when things are scaled, but you do want one for publication.
    ## so a version with and without the axes for all clusters is also generated as a .pdf so that you can have a properly scaled vector scale
  if (labelGenes == TRUE) {
    if (alignToCore == TRUE) {
    dummies <- gggenes::make_alignment_dummies(processed, aes(xmin = start, xmax = end, y = bgc, id = gene), on = coreGeneName)
    clusterDiagram <- ggplot2::ggplot(processed, aes(xmin = start, xmax = end, y = bgc, fill = gene, forward = direction, label="gene_oid")) +
      gggenes::geom_gene_arrow(arrowhead_height=unit(3,"mm"), arrowhead_width=unit(1,"mm")) +
      ggplot2::geom_blank(data=dummies, aes(forward = 1)) +
      ggplot2::facet_wrap(~ bgc, ncol = 1, scales = "free") +
      ggplot2::scale_fill_manual(values = finalColors) +
      ggplot2::labs(title=figureTitle, caption=captionText, y="**genomic neighborhoods**", x=NULL) +
      gggenes::theme_genes()
    xPlusDiagram <- clusterDiagram +
      ggplot2::theme(plot.title = ggtext::element_textbox_simple(size=14,  padding=margin(10,10,10,10), margin=margin(5,0,5,0), r=grid::unit(8,"pt"), fill="#EEEEEE", halign=0.5, valign=0.5), plot.title.position = "plot", axis.title.y = ggtext::element_markdown(size=12), axis.text.y = ggtext::element_markdown(size=9), plot.caption = ggtext::element_textbox_simple(lineheight=1.2, size=11, padding=margin(10,10,10,10), margin=margin(5,0,5,0), r=grid::unit(8,"pt"), fill="#D9D9D9", halign=0.5, valign=0.5), plot.caption.position = "plot")
    } else  {
      clusterDiagram <- ggplot2::ggplot(processed, aes(xmin = start, xmax = end, y = bgc, fill = gene, forward = direction, label="gene_oid")) +
        gggenes::geom_gene_arrow(arrowhead_height=unit(3,"mm"), arrowhead_width=unit(1,"mm")) +
        ggplot2::facet_wrap(~ bgc, ncol = 1, scales = "free") +
        ggplot2::scale_fill_manual(values = finalColors) +
        ggplot2::labs(title=figureTitle, caption=captionText, y="**genomic neighborhoods**", x=NULL) +
        gggenes::theme_genes()
      xPlusDiagram <- clusterDiagram +
        ggplot2::theme(plot.title = ggtext::element_textbox_simple(size=14,  padding=margin(10,10,10,10), margin=margin(5,0,5,0), r=grid::unit(8,"pt"), fill="#EEEEEE", halign=0.5, valign=0.5), plot.title.position = "plot", axis.title.y = ggtext::element_markdown(size=12), axis.text.y = ggtext::element_markdown(size=9), plot.caption = ggtext::element_textbox_simple(lineheight=1.2, size=11, padding=margin(10,10,10,10), margin=margin(5,0,5,0), r=grid::unit(8,"pt"), fill="#D9D9D9", halign=0.5, valign=0.5), plot.caption.position = "plot")
    }
  } else  {
    if (alignToCore == TRUE) {
      dummies <- gggenes::make_alignment_dummies(processed, aes(xmin = start, xmax = end, y = bgc, id = gene), on = coreGeneName)
      clusterDiagram <- ggplot2::ggplot(processed, aes(xmin = start, xmax = end, y = bgc, fill = gene, forward = direction)) +
        gggenes::geom_gene_arrow(arrowhead_height=unit(3,"mm"), arrowhead_width=unit(1,"mm")) +
        ggplot2::geom_blank(data=dummies, aes(forward = 1)) +
        ggplot2::facet_wrap(~ bgc, ncol = 1, scales = "free") +
        ggplot2::scale_fill_manual(values = finalColors) +
        ggplot2::labs(title=figureTitle, caption=captionText, y="**genomic neighborhoods**", x=NULL) +
        gggenes::theme_genes()
      xPlusDiagram <- clusterDiagram +
        ggplot2::theme(plot.title = ggtext::element_textbox_simple(size=14,  padding=margin(10,10,10,10), margin=margin(5,0,5,0), r=grid::unit(8,"pt"), fill="#EEEEEE", halign=0.5, valign=0.5), plot.title.position = "plot", axis.title.y = ggtext::element_markdown(size=12), axis.text.y = ggtext::element_markdown(size=9), plot.caption = ggtext::element_textbox_simple(lineheight=1.2, size=11, padding=margin(10,10,10,10), margin=margin(5,0,5,0), r=grid::unit(8,"pt"), fill="#D9D9D9", halign=0.5, valign=0.5), plot.caption.position = "plot")
    } else  {
      clusterDiagram <- ggplot2::ggplot(processed, aes(xmin = start, xmax = end, y = bgc, fill = gene, forward = direction)) +
        gggenes::geom_gene_arrow(arrowhead_height=unit(3,"mm"), arrowhead_width=unit(1,"mm")) +
        ggplot2::facet_wrap(~ bgc, ncol = 1, scales = "free") +
        ggplot2::scale_fill_manual(values = finalColors) +
        ggplot2::labs(title=figureTitle, caption=captionText, y="**genomic neighborhoods**", x=NULL) +
        gggenes::theme_genes()
      xPlusDiagram <- clusterDiagram +
        ggplot2::theme(plot.title = ggtext::element_textbox_simple(size=14,  padding=margin(10,10,10,10), margin=margin(5,0,5,0), r=grid::unit(8,"pt"), fill="#EEEEEE", halign=0.5, valign=0.5), plot.title.position = "plot", axis.title.y = ggtext::element_markdown(size=12), axis.text.y = ggtext::element_markdown(size=9), plot.caption = ggtext::element_textbox_simple(lineheight=1.2, size=11, padding=margin(10,10,10,10), margin=margin(5,0,5,0), r=grid::unit(8,"pt"), fill="#D9D9D9", halign=0.5, valign=0.5), plot.caption.position = "plot")
    }
  }
    ## get rid of the extra axes (but means no scale on-figure, for now)
  xMinusDiagram <- xPlusDiagram + ggplot2::theme(axis.line.x=element_blank(), axis.text.x=element_blank(),  axis.ticks.x=element_blank())
    ## this gets rid of that ggplot2 grid error by scaling sizing to the number of gene clusters being visualized
    ## also adapts width if you have chosen the longest names
  print("Cluster diagrams generated.")
  modheight <- length(uniqueBGCs)/2
  if (exists(x="showScaffolds") && showScaffolds==TRUE) {
    modwidth <- modheight/2
  } else {
    modwidth <- modheight/3
  }
  ## exports a .pdf and a .png version
  ggplot2::ggsave(file=finalpdfname, plot=xPlusDiagram, device="pdf", height=modheight, width=modwidth, limitsize=FALSE)
  ggplot2::ggsave(file=finalpngname, plot=xMinusDiagram, device="png", height=modheight, width=modwidth, limitsize=FALSE)
  ggplot2::ggsave(file=finalpdfnameX, plot=xMinusDiagram, device="pdf", height=modheight, width=modwidth, limitsize=FALSE)
  ggplot2::ggsave(file=finalpngnameX, plot=xMinusDiagram, device="png", height=modheight, width=modwidth, limitsize=FALSE)
  ## this way you have a copy of the annotated dataset and can skip auto-annotation if desired
  write.table(processed, annotationFileName, row.names=FALSE,sep="\t", quote=FALSE)  
  print("Very pretty gene clusters exported.")
    ## they'd better be after all of this
  return(processed)
}
