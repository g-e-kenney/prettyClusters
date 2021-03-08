#' Making a fake IMG metadata file from a GenBank file
#'
#' This function takes GenBank files as input and outputs fake IMG-style metadata files for genes of interest and for their neighbors, along with sequences etc.  Replaces generateNeighbors for GenBank input.
#' @param dataFolder location of folder containing GenBank files to be processed. Character string.
#' @param neighborNum number of neighbors surrounding the gene of interest. Integer.
#' @param goiListInput filename of file containing list of locus tags for gene of interest.  Character string.
#' @param geneName name of family of genes of interest (for filenaming purposes).  Character string.
#' @param removeDupes remove duplicate gene entries if they show up.  Boolean.
#' @return
#' @export
#' @examples
#' gbToIMGOutput <- gbToIMG(dataFolder="/data/here", neighborNum=5, goiListInput = "goiList.txt", geneName="genE", removeDupes=TRUE)
#'
gbToIMG <- function(dataFolder=dataFolder, neighborNum = 10, goiListInput = goiListInput, geneName="genE", removeDupes=TRUE, scaffoldGenBase =3000000000, genomeGenBase=4000000000)  {
    ## gettign the list of .gb files in the folder
    fileList <- list.files(path=dataFolder, pattern="*\\.gb$", full.names=TRUE, recursive=FALSE)
    scaffNum <- length(fileList)
    goiList <- read.csv(goiListInput, header=TRUE, sep = "\t",stringsAsFactors=FALSE)
    colnames(goiList) <- c("locus_tag")
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    fileName <- paste(fileDate,"_gb2img_",geneName,sep="")
    fauxNeighborDataFile <- paste(fileName, "_neighborData.txt",sep="")
    fauxGeneDataFile <- paste(fileName, "_geneData.txt",sep="")
    fauxNeighborSeqsFile <- paste(fileName, "_neighborSeqs.fa",sep="")
    fauxGeneSeqsFile <- paste(fileName, "_geneSeqs.fa",sep="")
    fauxContextFile <- paste(fileName, "_neighborContext.txt",sep="")
    if (exists(x="scaffoldGenBase")==FALSE) {
        scaffoldGenBase <- 30000000000
    }
    if (exists(x="genomeGenBase")==FALSE) {
        genomeGenBase <- 40000000000
    }   
    for (i in 1:scaffNum) {
        ## these ID number ranges aren't really used by IMG IDs at this point
        scaffoldGenID <- scaffoldGenBase + i*10000
        genomeGenID <- genomeGenBase  + i*10000
        ## file input
        gbTemp <- genbankr::readGenBank(fileList[i], ret.seq = FALSE)
        ## skipping files that have no gene annotations (using genes, exons, and cds as a proxy for "gene annotations").
        ## will probably fail at someone's stupid non-standard file where all the info is in mRNA type entries or something
        if (length(genbankr::genes(gbTemp)) == 0 && length(genbankr::cds(gbTemp))==0 && length(genbankr::exons(gbTemp))==0) {
            tempWarning <- paste("There appear to be no annotations for ", fileList[i],". Moving on.", sep="")
            print(tempWarning)
            next
        }
        ## getting temporary gene info, privileging CDS (likely to have translations) over exons and genes
        ## but going with whatever has the most columns of metadata
        if (length(GenomicRanges::mcols(genbankr::cds(gbTemp))) > length(GenomicRanges::mcols(genbankr::genes(gbTemp))))  {
            ## antiSmash GenBank files are Special
            if (length(GenomicRanges::mcols(genbankr::genes(gbTemp)))==0 && any(grepl("aSProdPred",colnames(GenomicRanges::mcols(genbankr::cds(gbTemp)))))==TRUE)  {
                tempWarning <- paste("This is probably an antiSmash-annotated GenBank file and it will probably have issues: ", fileList[i],". Moving on.", sep="")
                print(tempWarning)
            } else  {
                genesTemp <- genbankr::cds(gbTemp)
            }
        } else if (length(GenomicRanges::mcols(genbankr::exons(gbTemp))) > length(GenomicRanges::mcols(genbankr::genes(gbTemp)))) {
            genesTemp <- genbankr::exons(gbTemp)
        } else {
            genesTemp <- genbankr::genes(gbTemp)
        }
        ## locate GOI and subset the genomic neighborhood around it
        ## note that this assumes we're working with a locus tag to identify the genes
        ## but will back up to gene_id, protein_id, and transcript_id just in case
        if (length(which(colnames(GenomicRanges::mcols(genesTemp))=="locus_tag"))!=0) {
            goiHits <- which(genesTemp$locus_tag %in% goiList$locus_tag)
            if (length(goiHits) >= 1) {
                ## this should allow for multiple genes of interest on one contig
                for (j in 1:length(goiHits))  {
                    if (goiHits[j] - neighborNum < 1) {
                        startIdx <- 1
                    } else {
                        startIdx <- goiHits[j] - neighborNum
                    }
                    if (goiHits[j] + neighborNum > length(genesTemp$locus_tag)) {
                        endIdx <- length(genesTemp$locus_tag)
                    } else {
                        endIdx <- goiHits[j] + neighborNum
                    }
                    if (j == 1) {
                        neighborSubset <- genesTemp[startIdx:endIdx,]
                    } else {
                        neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
                    }
                }
            }
        }
        if (exists(x="neighborSubset")==FALSE && length(which(colnames(GenomicRanges::mcols(genesTemp))=="gene_id"))!=0)  {
            goiHits <- which(genesTemp$gene_id %in% goiList$locus_tag)
            if (length(goiHits) >= 1) {
                for (j in 1:length(goiHits))  {
                    if (goiHits[j] - neighborNum < 1) {
                        startIdx <- 1
                    } else {
                        startIdx <- goiHits[j] - neighborNum
                    }
                    if (goiHits[j] + neighborNum > length(genesTemp$gene_id)) {
                        endIdx <- length(genesTemp$gene_id)
                    } else {
                        endIdx <- goiHits[j] + neighborNum
                    }
                    if (j == 1) {
                        neighborSubset <- genesTemp[startIdx:endIdx,]
                    } else {
                        neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
                    }
                }
            }
        }
        if (exists(x="neighborSubset")==FALSE && length(which(colnames(GenomicRanges::mcols(genesTemp))=="transcript_id"))!=0)  {
            goiHits <- which(genesTemp$transcript_id %in% goiList$locus_tag)
            if (length(goiHits) >= 1) {
                for (j in 1:length(goiHits))  {
                    if (goiHits[j] - neighborNum < 1) {
                        startIdx <- 1
                    } else {
                        startIdx <- goiHits[j] - neighborNum
                    }
                    if (goiHits[j] + neighborNum > length(genesTemp$transcript_id)) {
                        endIdx <- length(genesTemp$transcript_id)
                    } else {
                        endIdx <- goiHits[j] + neighborNum
                    }
                    if (j == 1) {
                        neighborSubset <- genesTemp[startIdx:endIdx,]
                    } else {
                        neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
                    }
                }
            }
        }
        if (exists(x="neighborSubset")==FALSE && length(which(colnames(GenomicRanges::mcols(genesTemp))=="protein_id"))!=0)  {
             goiHits <- which(genesTemp$protein_id %in% goiList$locus_tag)
            if (length(goiHits) >= 1) {
                for (j in 1:length(goiHits))  {
                    if (goiHits[j] - neighborNum < 1) {
                        startIdx <- 1
                    } else {
                        startIdx <- goiHits[j] - neighborNum
                    }
                    if (goiHits[j] + neighborNum > length(genesTemp$protein_id)) {
                        endIdx <- length(genesTemp$protein_id)
                    } else {
                        endIdx <- goiHits[j] + neighborNum
                    }
                    if (j == 1) {
                        neighborSubset <- genesTemp[startIdx:endIdx,]
                    } else {
                        neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
                    }
                }
            }
        }
        if (exists(x="neighborSubset")==FALSE) {
            tempWarning <- paste("No gene of interest found in locus_tag, gene_id, transcript_id, or protein_id for: ", fileList[i],". Moving on.", sep="")
            print(tempWarning)
            next
        }
        ## starting the fake IMG metadata table
        ## begin with the locus tag or gene ID
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="locus_tag"))!=0) {
            fauxIMG <- neighborSubset$locus_tag
        } else if  (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="gene_id"))!=0)  {
            fauxIMG <- neighborSubset$gene_id
        } else if  (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="transcript_id"))!=0)  {
            fauxIMG <- neighborSubset$transcript_id
        } else if  (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="protein_id"))!=0)  {
            fauxIMG <- neighborSubset$protein_id
        } else {
            tempWarning <- paste("No locus_tag, gene_id, transcript_id, or protein_id in gene info for: ", fileList[i],". Moving on.", sep="")
            print(tempWarning)
            next
        }
        fauxIMG <- as.data.frame(fauxIMG)
        colnames(fauxIMG) <- "Locus.Tag"
        ## making the fake gene_oid
        ## numbering is based on the neighborhood subset so it's not going to be consistent if you redo this with different size neighborhoods
        fauxIMG$gene_oid <-  as.numeric(rownames(fauxIMG)) + scaffoldGenID
        ## while we're here, let's quickly make the neighborsContext file
        contextTable <-  fauxIMG$gene_oid
        contextTable <- as.data.frame(contextTable)
        colnames(contextTable) <- c("gene_oid")
        contextTable$source_scaffold_id <- scaffoldGenID
        locGOI <- which(fauxIMG$Locus.Tag %in% goiList$locus_tag)
        contextTable$source_gene_oid <- ""
        if (length(locGOI) == 1)  {
            contextTable$source_gene_oid <- fauxIMG$gene_oid[locGOI]
        } else {
            ## in case we have multiple genes of interest, le sigh
            for (j in 1:length(locGOI)) {
                minTemp <- fauxIMG$gene_oid[locGOI[j]] - neighborNum
                maxTemp <- fauxIMG$gene_oid[locGOI[j]] + neighborNum
                contextTable <- contextTable  %>% dplyr::mutate(source_gene_oid = ifelse (gene_oid >= minTemp & gene_oid <= maxTemp,fauxIMG$gene_oid[locGOI[j]], source_gene_oid))
            }
        }
        ## Gene Product
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="product"))!=0)  {
            fauxIMG$Gene.Product.Name <- neighborSubset$product
        } else {
            fauxIMG$Gene.Product.Name <- ""
        }
        ## Genome ID
        fauxIMG$Genome.ID <- genomeGenID
        ## Genome Name - likely to actually be the accession
        fauxIMG$Genome.Name <-  as.character(unlist(GenomeInfoDb::genome(neighborSubset)))
        ## Gene Symbol - going through less-useless failure modes before returning empty values.  Starts with gene, then tries protein_id and gene_id.
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="gene"))!=0)  {
            fauxIMG$Gene.Symbol <- neighborSubset$gene
        } else if  (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="protein_id"))!=0)  {
            fauxIMG$Gene.Symbol <- neighborSubset$protein_id
        } else if  (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="gene_id"))!=0)  {
            fauxIMG$Gene.Symbol <- neighborSubset$gene_id
        } else {
            fauxIMG$Gene.Symbol <- ""
        }
        ## GenBank Accession - using protein IDs here because they fit more poorly elsewhere
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="protein_id"))!=0)  {
            fauxIMG$GenBank.Accession <- neighborSubset$protein_id
        } else {
            fauxIMG$GenBank.Accession <- ""
        }
        ## Chromosome - not likely to be useful and generally empty in IMG anyway
        fauxIMG$Chromosome <- ""
        ## Start Coord - should exist as long as genes do
        fauxIMG$Start.Coord <- GenomicRanges::start(neighborSubset)
        ## End Coord - should exist as long as genes do
        fauxIMG$End.Coord <- GenomicRanges::end(neighborSubset)
        ## Strand - should exist as long as genes do
        fauxIMG$Strand <- unlist(as.data.frame(GenomicRanges::strand(neighborSubset)))
        ## DNA Sequence Length (bp) - should exist as long as genes do
        fauxIMG$DNA.Sequence.Length..bp. <- GenomicRanges::end(neighborSubset) - GenomicRanges::start(neighborSubset) + 1
        ## Amino Acid Sequence Length (aa) - should exist as long as genes do
        fauxIMG$Amino.Acid.Sequence.Length..aa. <- (fauxIMG$DNA.Sequence.Length..bp. - 3)/3
        ## Locus Type
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="type"))!=0)  {
            fauxIMG$Locus.Type <- neighborSubset$type
        } else {
            fauxIMG$Locus.Type <- ""
        }
        ## Is Pseudogene
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="pseudo"))!=0)  {
            fauxIMG$Is.Pseudogene <- neighborSubset$pseudo
        } else {
            fauxIMG$Is.Pseudogene <- ""
        }
        ## Is Obsolete, Is Partial Gene, Add Date, Is Public - unlikely to auto-gen in any meaningful way
        fauxIMG$Is.Obsolete <- "No"
        fauxIMG$Is.Partial.Gene <- ""
        fauxIMG$Add.Date <- ""
        fauxIMG$Is.Public <- ""
        ## Transmembrane Helices - just in case, probably via TMHMM
        fauxIMG$Transmembrane.Helices <- ""
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="inference"))!=0 && any(grepl("TMHMM", neighborSubset$inference)))  {
            fauxIMG$Transmembrane.Helices <- ""
            sigLoc <- grep("TMHMM", neighborSubset$inference)
            fauxIMG$Transmembrane.Helices[sigLoc] <- "yes"
        } else if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="note"))!=0 && any(grepl("TMHMM", neighborSubset$note)))  {
            fauxIMG$Transmembrane.Helices <- ""
            sigLoc <- grep("TMHMM", neighborSubset$note)
            fauxIMG$Transmembrane.Helices[sigLoc] <- "yes"
        } else if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="db_xref"))!=0 && any(grepl("TMHMM", neighborSubset$db_xref)))  {
            fauxIMG$Transmembrane.Helices <- ""
            sigLoc <- grep("TMHMM", neighborSubset$db_xref)
            fauxIMG$Transmembrane.Helices[sigLoc] <- "yes"
        } else {
            fauxIMG$Transmembrane.Helices <- ""
        }             
        ## Signal Peptides - just in case, probably via SignalP if they exist
        fauxIMG$Signal.Peptides <- ""
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="inference"))!=0 && any(grepl("SignalP", neighborSubset$inference))) {
            fauxIMG$Signal.Peptides <- ""
            sigLoc <- grep("SignalP", neighborSubset$inference)
            fauxIMG$Signal.Peptides[sigLoc] <- "yes"
        } else if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="note"))!=0 && any(grepl("SignalP", neighborSubset$note)))  {
            fauxIMG$Signal.Peptides <- ""
            sigLoc <- grep("SignalP", neighborSubset$note)
            fauxIMG$Signal.Peptides[sigLoc] <- "yes"
        } else if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="db_xref"))!=0 && any(grepl("SignalP", neighborSubset$db_xref)))  {
            fauxIMG$Signal.Peptides <- ""
            sigLoc <- grep("SignalP", neighborSubset$db_xref)
            fauxIMG$Signal.Peptides[sigLoc] <- "yes"
        } else {
            fauxIMG$Signal.Peptides <- ""
        }
        ## Scaffold ID
        fauxIMG$Scaffold.ID <- scaffoldGenID
        ## Scaffold External Accession - using the GenBank accession. Should always exist.
        fauxIMG$Scaffold.External.Accession <- as.character(GenomeInfoDb::genome(neighborSubset))
        ## Scaffold Name - using the provided name. Should always exist.
        fauxIMG$Scaffold.Name <- as.character(GenomeInfoDb::seqnames(neighborSubset))
        ## Scaffold Length (bp) - should always be present.
        fauxIMG$Scaffold.Length..bp. <- as.numeric(GenomeInfoDb::seqlengths(neighborSubset))
        ## Scaffold GC % - unlikely to be present
        fauxIMG$Scaffold.GC.. <- ""
        ## Scaffold Read Depth - unlikely to be present
        fauxIMG$Scaffold.Read.Depth <- ""
        ## COG - could possibly add at least to some EMBL files, on the same terms as InterPRo or TIGRfam?
        fauxIMG$COG <- ""
        ## Pfam. N.b. requiring digits to avoid "prediction via PFAM" or whatever hits, and making the name format more IMG-like
        fauxIMG$Pfam <- ""
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="inference"))!=0) {
            if (any(grepl("PF[[:digit:]]", neighborSubset$inference)) || any(grepl("PFAM[[:digit:]]", neighborSubset$inference)) || any(grepl("pfam[[:digit:]]", neighborSubset$inference)))  {
                for (j in length(fauxIMG$Locus.Tag)) {
                    if (any(grepl("PF[[:digit:]]", neighborSubset$inference[j])) || any(grepl("PFAM[[:digit:]]", neighborSubset$inference[j])) || any(grepl("pfam[[:digit:]]", neighborSubset$inference[j])))  {
                        fauxIMG$Pfam[j] <- neighborSubset$inference[j]
                        fauxIMG$Pfam[j] <- gsub("PF"," pfam", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("PFAM"," pfam", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]][[:digit:]]", " ", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]]", " ", fauxIMG$Pfam[j])
                    } else  {
                         fauxIMG$Pfam[j] <- ""
                    }
                }
            }            
        } else if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="note"))!=0) {
            if (any(grepl("PF[[:digit:]]", neighborSubset$note)) || any(grepl("PFAM[[:digit:]]", neighborSubset$note)) || any(grepl("pfam[[:digit:]]", neighborSubset$note)))  {
                for (j in length(fauxIMG$Locus.Tag)) {
                    if (any(grepl("PF[[:digit:]]", neighborSubset$note[j])) || any(grepl("PFAM[[:digit:]]", neighborSubset$note[j])) || any(grepl("pfam[[:digit:]]", neighborSubset$note[j])))  {
                        fauxIMG$Pfam[j] <- neighborSubset$note[j]
                        fauxIMG$Pfam[j] <- gsub("PF"," pfam", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("PFAM"," pfam", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]][[:digit:]]", " ", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]]", " ", fauxIMG$Pfam[j])
                    } else  {
                         fauxIMG$Pfam[j] <- ""
                    }
                }
            }    
        } else if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="db_xref"))!=0) {
            if (any(grepl("PF[[:digit:]]", neighborSubset$db_xref)) || any(grepl("PFAM[[:digit:]]", neighborSubset$db_xref)) || any(grepl("pfam[[:digit:]]", neighborSubset$db_xref)))  {
                for (j in length(fauxIMG$Locus.Tag)) {
                    if (any(grepl("PF[[:digit:]]", neighborSubset$db_xref[j])) || any(grepl("PFAM[[:digit:]]", neighborSubset$db_xref[j])) || any(grepl("pfam[[:digit:]]", neighborSubset$db_xref[j])))  {
                        fauxIMG$Pfam[j] <- neighborSubset$db_xref[j]
                        fauxIMG$Pfam[j] <- gsub("PF"," pfam", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("PFAM"," pfam", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]][[:digit:]]", " ", fauxIMG$Pfam[j])
                        fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]]", " ", fauxIMG$Pfam[j])
                    } else  {
                         fauxIMG$Pfam[j] <- ""
                    }
                }
            }  
        } else {
            fauxIMG$Pfam <- ""
        }
        ## Tigrfam - these generally have a standard ID format
        fauxIMG$TIGRfam <- ""
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="inference"))!=0 && any(grepl("TIGR[[:digit:]]", neighborSubset$inference)))  {
            for (j in length(fauxIMG$Locus.Tag)) {
                if (any(grepl("TIGR[[:digit:]]", neighborSubset$inference[j]))) {
                    fauxIMG$TIGRfam[i] <- neighborSubset$inference[j]
                }
            }       
        } else if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="note"))!=0 && any(grepl("TIGR[[:digit:]]", neighborSubset$note))) {
            for (j in length(fauxIMG$Locus.Tag)) {
                if (any(grepl("TIGR[[:digit:]]", neighborSubset$note[j]))) {
                    fauxIMG$TIGRfam[i] <- neighborSubset$note[j]
                }
            }  
        } else if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="db_xref"))!=0 && any(grepl("TIGR[[:digit:]]", neighborSubset$db_xref))) {
            for (j in length(fauxIMG$Locus.Tag)) {
                if (any(grepl("TIGR[[:digit:]]", neighborSubset$db_xref[j]))) {
                    fauxIMG$TIGRfam[j] <- neighborSubset$db_xref[j]
                }
            }   
        } else {
            fauxIMG$TIGRfam <- ""
        }
        ## SMART.ID, SUPERFam.ID and CATH.FunFam.ID just started showing up in IMG data
        ## haven't seen them much in GenBank files so sticking with blank for now as a placeholder
        fauxIMG$SMART.ID <- ""
        fauxIMG$SUPERFam.ID <- ""
        fauxIMG$CATH.FunFam.ID <- ""
        ## Enzyme - aka EC number - shows up very unpredictably and with variable format so I'm inclined to stick with a blank space for now
        ## When it does show it is sometimes separately annotated as EC_number
        fauxIMG$Enzyme <- ""
        ## KO, IMG Term - the latter is unlikely to be present, and I haven't seen the former frequently, so sticking with a blank for now
        fauxIMG$KO <- ""
        fauxIMG$IMG.Term <- ""
        ## InterPro - these almost always have the right IPR-based code but just in case
        fauxIMG$InterPro <- ""
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="inference"))!=0) {
            if (any(grepl("InterPro", neighborSubset$inference)) || any(grepl("IPR[[:digit:]]", neighborSubset$inference)))  {
                for (j in length(fauxIMG$Locus.Tag)) {
                    if (any(grepl("InterPro", neighborSubset$inference[j])) || any(grepl("IPR[[:digit:]]", neighborSubset$inference[j])))  {
                        fauxIMG$InterPro[j] <- neighborSubset$inference[j]
                    }       
                }
            }  
        } else  if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="note"))!=0) {
            if (any(grepl("InterPro", neighborSubset$note)) || any(grepl("IPR[[:digit:]]", neighborSubset$note)))  {
                for (j in length(fauxIMG$Locus.Tag)) {
                    if (any(grepl("InterPro", neighborSubset$note[j])) || any(grepl("IPR[[:digit:]]", neighborSubset$note[j])))  {
                        fauxIMG$InterPro[j] <- neighborSubset$note[j]
                    }       
                }
            }  
        } else   if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="db_xref"))!=0) {
            if (any(grepl("InterPro", neighborSubset$db_xref)) || any(grepl("IPR[[:digit:]]", neighborSubset$db_xref)))  {
                for (j in length(fauxIMG$Locus.Tag)) {
                    if (any(grepl("InterPro", neighborSubset$db_xref[j])) || any(grepl("IPR[[:digit:]]", neighborSubset$db_xref[j])))  {
                        fauxIMG$InterPro[j] <- neighborSubset$db_xref[j]
                    }       
                }
            }  
        } else {
            fauxIMG$InterPro <- ""
        }
        ## data assembly for larger metadata file
        if (exists(x="fauxNeighborData")==FALSE) {
            fauxNeighborData <- fauxIMG
        } else {
            fauxNeighborData <- rbind(fauxNeighborData, fauxIMG)
        }
        ## data assembly for the context file
        if (exists(x="fullContext")==FALSE)  {
            fullContext <- contextTable
         } else {
            fullContext <- rbind(fullContext, contextTable)
        }          
        ## putting together the protein sequences, if applicable
        if (length(which(colnames(GenomicRanges::mcols(neighborSubset))=="translation"))!=0)  {
            fauxSeqIMG <-   neighborSubset$translation
            fauxSeqIMG <- as.data.frame(fauxSeqIMG)
            colnames(fauxSeqIMG) <- c("sequence")
            fauxSeqIMG$gene_oid <- as.character(fauxIMG$gene_oid)
            fauxSeqIMG <- fauxSeqIMG %>% dplyr::select(gene_oid, everything())
            if (exists(x="fauxNeighborSeqs")==FALSE) {
                fauxNeighborSeqs <- fauxSeqIMG
            } else {
                fauxNeighborSeqs <- rbind(fauxNeighborSeqs, fauxSeqIMG)
            }
        }
        rm(neighborSubset)
    }
    ## outputting the faux-IMG-formatted metadata file for the neighborhoods
    if (removeDupes == TRUE) {
        fauxNeighborData <- fauxNeighborData[!duplicated(fauxNeighborData$Locus.Tag), ]
    }
    fauxNeighborData$Pfam <- as.character(fauxNeighborData$Pfam)
    fauxNeighborData$TIGRfam <- as.character(fauxNeighborData$TIGRfam)
    fauxNeighborData$InterPro <- as.character(fauxNeighborData$InterPro)
    
    write.table(fauxNeighborData, file=fauxNeighborDataFile, row.names=FALSE, col.names = TRUE, sep="\t")
    ## writing the context file tooo
    write.table(fullContext, file=fauxContextFile, row.names=FALSE, col.names = TRUE, sep="\t")
    ## and the .fa file for the neighbor protein sequences
    fauxNeighborSeqList <- list()
    for (j in 1:length(fauxNeighborSeqs[,2])) {
        fauxNeighborSeqList[[fauxNeighborSeqs[j,1]]] <- fauxNeighborSeqs[j,2]
    }
    seqinr::write.fasta(fauxNeighborSeqList, names=names(fauxNeighborSeqList),file.out=fauxNeighborSeqsFile)
    ## and outputting the genes of interest separately
    fauxGeneData <- fauxNeighborData %>% dplyr::filter(Locus.Tag %in% goiList$locus_tag)
    if (removeDupes == TRUE) {
        fauxGeneData <- fauxGeneData[!duplicated(fauxGeneData$Locus.Tag), ]
    }
    write.table(fauxGeneData, file=fauxGeneDataFile, row.names=F, col.names = T, sep="\t")
    ## and their sequences too
    fauxGeneList <- fauxGeneData$gene_oid
    fauxGeneSeqs <- fauxNeighborSeqs %>% dplyr::filter(gene_oid %in% fauxGeneList)
    fauxGeneSeqList <- list()
    for (j in 1:length(fauxGeneSeqs[,2])) {
        fauxGeneSeqList[[fauxGeneSeqs[j,1]]] <- fauxGeneSeqs[j,2]
    }
    seqinr::write.fasta(fauxGeneSeqList, names=names(fauxGeneSeqList),file.out=fauxGeneSeqsFile)
    print("GenBank-formatted files in folder processed.")
    return(list(geneData=fauxGeneData, neighborData=fauxNeighborData))
} 

