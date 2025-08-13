#' Making a fake IMG metadata file from a GenBank file
#'
#' This function takes a folder of GenBank files as input and outputs fake IMG-style metadata files for genes of interest and for their neighbors, along with sequences etc.  Replaces generateNeighbors for GenBank input.
#' @param dataFolder Location of folder containing GenBank files to be processed. Character string, required.
#' @param neighborNum Number of neighbors surrounding the gene of interest. Integer, defaults to 10.
#' @param goiListInput Filename of file containing list of locus tags for gene of interest.  Character string, required.
#' @param geneName Name of family of genes of interest (for filenaming purposes).  Character string, required.
#' @param removeDupes Remove duplicate gene entries if they show up.  Boolean, defaults to TRUE.
#' @param scaffoldGenBase Starting faux scaffold ID number. Integer, defaults to 3000000000.
#' @param genomeGenBase Starting faux genome ID number. Integer, defaults to 4000000000.
#' @param includeIPR Dictates whether or not InterPro annotations should be included. Boolean, defaults to FALSE.
#' @param seqExtract Dictates whether a new .faa file for protein sequences will be prepared from an import sharing the genome name. Boolean, defaults to TRUE.
#' @return List containing IMG-style metadata for genes of interest, their neighbors, and a generateNeighbors-like neighborsContext file
#' @export
#' @importFrom magrittr %>% 
#' @importFrom utils read.csv write.csv write.table read.table
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' gbToIMGOutput <- gbToIMG(dataFolder="/data/here", 
#'                          neighborNum=5, 
#'                          goiListInput = "goiList.txt", 
#'                          geneName="genE")
#' }
gbToIMG <- function(dataFolder=dataFolder, 
                    neighborNum = 10, 
                    goiListInput = goiListInput, 
                    geneName=geneName, 
                    removeDupes=TRUE, 
                    scaffoldGenBase =30000000000, 
                    genomeGenBase=40000000000, 
                    includeIPR = FALSE,
                    seqExtract = TRUE)  {
  ## getting the list of .gb files in the folder
  ## note: will accept .gbk, .gbf as file suffixes
  fileList <- list.files(path=dataFolder, pattern="*\\.gb", full.names=TRUE, recursive=FALSE)
  scaffNum <- length(fileList)
  ## just say no to scientific notation
  options(scipen = 999)
  goiList <- read.csv(goiListInput, header=TRUE, sep = "\t",stringsAsFactors=FALSE)
  colnames(goiList) <- c("locus_tag")
  fileDate <- format(Sys.Date(),format="%Y%m%d")
  fileName <- paste(fileDate,"_gb2img_",geneName,sep="")
  fauxNeighborDataFile <- paste(fileName, "_neighborData.txt",sep="")
  fauxGeneDataFile <- paste(fileName, "_geneData.txt",sep="")
  fauxNeighborSeqsFile <- paste(fileName, "_neighborSeqs.fa",sep="")
  fauxGeneSeqsFile <- paste(fileName, "_geneSeqs.fa",sep="")
  fauxContextFile <- paste(fileName, "_neighborContext.txt",sep="")
  noAnnotList <- list()
  noGoiList <- list()
  asFmtList <- list()
  processedList <- list()
  noAnnotFile <- paste(fileName, "_noAnnotList.txt",sep="")
  noGoiFile <- paste(fileName, "_noGoiList.txt",sep="")
  asFmtFile <- paste(fileName, "_asFmtList.txt",sep="")
  processedFile <- paste(fileName, "_processedList.txt",sep="")
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
    antismash <- FALSE
    ## file input
##    gbTemp <- genbankr::readGenBank(fileList[i], ret.seq = FALSE)
    
    gbTemp <- gggenomes::read_gbk(fileList[i])
    ## skipping files that have no gene annotations (using genes, exons, and cds as a proxy for "gene annotations").
    ## will probably fail at someone's stupid non-standard file where all the info is in mRNA type entries or something
    if (length(which(gbTemp$type == "CDS")) == 0 && length(which(gbTemp$type == "gene"))==0 && length(which(gbTemp$type == "exon"))==0) {
      noAnnotList <- append(noAnnotList, fileList[i])
      tempWarning <- paste("There appear to be no annotations for ", fileList[i],". Moving on.", sep="")
      print(tempWarning)
      next
    }
    ## getting temporary gene info, privileging CDS (likely to have translations) over exons and genes
    ## but going with whatever has the most columns of metadata
    if (length(which(gbTemp$type == "CDS")) > length(which(gbTemp$type == "gene")) && any(grepl("a_s_domain",colnames(gbTemp))==TRUE)) {
      ## antiSmash GenBank files are Special
        asFmtList <- append(asFmtList, fileList[i])
        tempWarning <- paste("This is probably an antiSmash-annotated GenBank file - names and annotations may differ from the original version: ", fileList[i],".", sep="")
        print(tempWarning)
        antismash <-TRUE
        genesTemp <- gbTemp[which(gbTemp$type == "CDS"),]
    } else if (length(which(gbTemp$type == "CDS")) > 0 && length(which(gbTemp$type == "CDS")) > length(which(gbTemp$type == "exon")) ) {
        genesTemp <- gbTemp[which(gbTemp$type == "CDS"),]
      } else if (length(which(gbTemp$type == "exon")) > length(which(gbTemp$type == "gene"))) {
      genesTemp <- gbTemp[which(gbTemp$type == "exon"),]
    } else {
      genesTemp <- gbTemp[which(gbTemp$type == "gene"),]
    }
    ## locate GOI and subset the genomic neighborhood around it
    ## note that this assumes we're working with a locus tag to identify the genes
    ## but will back up to name (antismash), gene, protein_id, and transcript_id just in case
    if (length(which(colnames(genesTemp)=="locus_tag"))!=0 && length(which(genesTemp$locus_tag != "NA")) != 0) {
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
            tempSpan <- endIdx - startIdx
          } else {
            neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
            tempSpan <- append(tempSpan, endIdx - startIdx)
          }
        }
        noHits <- FALSE
        fauxIMG <- neighborSubset$locus_tag
      } else {
        ## we didn't see a goi hit despite having locus_tags
        noHits <- TRUE
      }
    } else { 
      ## there wasn't a locus_tag column
      noHits <- TRUE
    }
    ## if it is an antismash file which may use "name" not locus_tag
    if (noHits == TRUE && antismash == TRUE && length(which(colnames(genesTemp)=="name"))!=0)  {
      goiHits <- which(genesTemp$name %in% goiList$locus_tag)
      if (length(goiHits) >= 1) {
        for (j in 1:length(goiHits))  {
          if (goiHits[j] - neighborNum < 1) {
            startIdx <- 1
          } else {
            startIdx <- goiHits[j] - neighborNum
          }
          if (goiHits[j] + neighborNum > length(genesTemp$name)) {
            endIdx <- length(genesTemp$name)
          } else {
            endIdx <- goiHits[j] + neighborNum
          }
          if (j == 1) {
            neighborSubset <- genesTemp[startIdx:endIdx,]
            tempSpan <- endIdx - startIdx
          } else {
            neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
            tempSpan <- append(tempSpan, endIdx - startIdx)
          }
        }
        noHits <- FALSE
        fauxIMG <- neighborSubset$name
      } else if (exists(x="neighborSubset")==FALSE) {
        noHits <- TRUE
      }
      ## implicitly if there are no antismash "name" entries and yet we got here, noHits is still TRUE
    }
    if (noHits == TRUE && length(which(colnames(genesTemp)=="gene_id"))!=0)  {
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
            tempSpan <- endIdx - startIdx
          } else {
            neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
            tempSpan <- append(tempSpan, endIdx - startIdx)
          }
        }
        noHits <- FALSE
        fauxIMG <- neighborSubset$gene_id
      } else if (exists(x="neighborSubset")==FALSE) {
        noHits <- TRUE
      }
      ## implicitly if there are no gene_ids but we got here, noHits is still TRUE
    }
    if (noHits == TRUE && length(which(colnames(genesTemp)=="transcript_id"))!=0)  {
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
            tempSpan <- endIdx - startIdx
          } else {
            neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
            tempSpan <- append(tempSpan, endIdx - startIdx)
          }
        }
        noHits <- FALSE
        fauxIMG <- neighborSubset$transcript_id
      } else if (exists(x="neighborSubset")==FALSE) {
        noHits <- TRUE
      }
    }
    if (noHits == TRUE && length(which(colnames(genesTemp)=="protein_id"))!=0)  {
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
            tempSpan <- endIdx - startIdx
          } else {
            neighborSubset <- c(neighborSubset, genesTemp[startIdx:endIdx,])
            tempSpan <- append(tempSpan, endIdx - startIdx)
          }
        }
        noHits <- FALSE
        fauxIMG <- neighborSubset$protein_id
      } else if (exists(x="neighborSubset")==FALSE) {
        noHits <- TRUE
      }
    }
    if (exists(x="neighborSubset")==FALSE && noHits == TRUE) {
      noGoiList <- append(noGoiList, fileList[i])
      tempWarning <- paste("No gene of interest found in locus_tag, name, gene_id, protein_id, or transcript_id for: ", fileList[i],". Moving on.", sep="")
      print(tempWarning)
      next
    }
    ## assuming we got something from the above, we'll use that fauxIMG to seed the Locus.Tag field
    fauxIMG <- as.data.frame(fauxIMG)
    colnames(fauxIMG) <- "Locus.Tag"

    ## making the fake gene_oid
    ## numbering is based on the neighborhood subset so it's not going to be consistent if you redo this with different size neighborhoods
    fauxIMG$gene_oid <-  as.numeric(rownames(fauxIMG)) + scaffoldGenID
    ## while we're here, let's quickly make the neighborsContext file
    contextTable <-  as.character(fauxIMG$gene_oid)
    contextTable <- as.data.frame(contextTable)
    colnames(contextTable) <- c("gene_oid")
    contextTable$source_scaffold_id <- as.character(scaffoldGenID)
    locGOI <- which(fauxIMG$Locus.Tag %in% goiList$locus_tag)
    contextTable$source_gene_oid <- ""
    if (length(locGOI) == 1)  {
      contextTable$source_gene_oid <- fauxIMG$gene_oid[locGOI]
    } else {
      ## in case we have multiple genes of interest, le sigh
      for (j in 1:length(locGOI)) {
        minTemp <- fauxIMG$gene_oid[locGOI[j]] - neighborNum
        maxTemp <- fauxIMG$gene_oid[locGOI[j]] + neighborNum
        contextTable <- contextTable  %>% dplyr::mutate(source_gene_oid = ifelse(.data$gene_oid >= minTemp & .data$gene_oid <= maxTemp,
                                                                                 fauxIMG$gene_oid[locGOI[j]],
                                                                                 .data$source_gene_oid))
      }
    }
    ## Gene Product
    if (length(which(colnames(neighborSubset)=="product"))!=0 && antismash == FALSE) {
      fauxIMG$Gene.Product.Name <- as.character(neighborSubset$product)
    } else if (length(which(colnames(neighborSubset)=="label"))!=0 && antismash == TRUE) {
      fauxIMG$Gene.Product.Name <- as.character(neighborSubset$label)
    } else {
      fauxIMG$Gene.Product.Name <- ""
    }
    fauxIMG$Gene.Product.Name[which(is.na(fauxIMG$Gene.Product.Name))] <- ""    
    ## Genome ID
    fauxIMG$Genome.ID <- as.character(genomeGenID)
    ## Genome Name - absent an organism this will default to the filename/accession
    if (gbTemp$type[1] == "region" && gbTemp$organism[1] != "NA") {
      fauxIMG$Genome.Name <- as.character(gbTemp$organism[1])
    } else {
     fauxIMG$Genome.Name <- as.character(gbTemp$seq_id[1]) 
    }

##    fauxIMG$Genome.Name <-  as.character(genbankr::sources(gbTemp)$organism)
    ## Gene Symbol - going through less-useless failure modes before returning empty values.  Starts with name, then gene, then tries protein_id and gene_id.
    if (length(which(colnames(neighborSubset)=="name"))!=0 && length(which(neighborSubset$name != "NA")) != 0)  {
      fauxIMG$Gene.Symbol <- neighborSubset$name
    } else if  (length(which(colnames(neighborSubset)=="gene"))!=0 && length(which(neighborSubset$gene != "NA")) != 0)  {
      fauxIMG$Gene.Symbol <- neighborSubset$gene
    } else if  (length(which(colnames(neighborSubset)=="protein_id"))!=0 && length(which(neighborSubset$protein_id != "NA")) != 0)  {
      fauxIMG$Gene.Symbol <- neighborSubset$protein_id
    } else if  (length(which(colnames(neighborSubset)=="gene_id"))!=0 && length(which(neighborSubset$gene_id != "NA")) != 0)  {
      fauxIMG$Gene.Symbol <- neighborSubset$gene_id
          } else {
      fauxIMG$Gene.Symbol <- ""
          }
    fauxIMG$Gene.Symbol[which(is.na(fauxIMG$Gene.Symbol))] <- ""
    ## GenBank Accession - using protein IDs here because they fit more poorly elsewhere
    ## (Gene.Symbol is more usefully populated by the gene value, if there is one.)
    if (length(which(colnames(neighborSubset)=="protein_id"))!=0 && length(which(neighborSubset$protein_id != "NA")) != 0)  {
      fauxIMG$GenBank.Accession <- neighborSubset$protein_id
    } else {
      fauxIMG$GenBank.Accession <- ""
    }
    ## Chromosome - not likely to be useful and generally empty in IMG anyway
    fauxIMG$Chromosome <- ""
    ## Start Coord - should exist as long as genes do
    fauxIMG$Start.Coord <- neighborSubset$start
    ## End Coord - should exist as long as genes do
    fauxIMG$End.Coord <- neighborSubset$end
    ## Strand - should exist as long as genes do
    fauxIMG$Strand <- neighborSubset$strand
    ## DNA Sequence Length (bp) - should exist as long as genes do
    fauxIMG$DNA.Sequence.Length..bp. <- neighborSubset$end - neighborSubset$start + 1
    ## Amino Acid Sequence Length (aa) - should exist as long as genes do
    fauxIMG$Amino.Acid.Sequence.Length..aa. <- (fauxIMG$DNA.Sequence.Length..bp. - 3)/3
    ## Locus Type
    if (length(which(colnames(neighborSubset)=="type"))!=0 && length(which(neighborSubset$type != "NA")) != 0)  {
      fauxIMG$Locus.Type <- neighborSubset$type
    } else {
      fauxIMG$Locus.Type <- ""
    }
    ## Is Pseudogene
    if (length(which(colnames(neighborSubset)=="pseudo"))!=0 && length(which(neighborSubset$pseudo != "NA")) != 0)  {
      fauxIMG$Is.Pseudogene <- neighborSubset$pseudo
    } else {
      fauxIMG$Is.Pseudogene <- ""
    }
    ## Is Obsolete, Is Partial Gene, Add Date, Is Public - unlikely to auto-gen in any meaningful way
    fauxIMG$Is.Obsolete <- "No"
    fauxIMG$Is.Partial.Gene <- ""
    fauxIMG$Add.Date <- ""
    fauxIMG$Is.Public <- ""
    fauxIMG$Sequencing.Status <- ""
    ## Scaffold ID
    fauxIMG$Scaffold.ID <- as.character(scaffoldGenID)
    ## Scaffold External Accession - using the GenBank accession. Should always exist.
    fauxIMG$Scaffold.External.Accession <- as.character(neighborSubset$seq_id)
    ## Scaffold Name - using the provided name. Should always exist, but may be something dumb (e.g. ATCC number or strain number).
#    fauxIMG$Scaffold.Name <- as.character(genbankr::sources(gbTemp)$organism)
    if (gbTemp$type[1] == "region" && gbTemp$organism[1] != "NA") {
      fauxIMG$Scaffold.Name <- as.character(gbTemp$organism[1])
    } else {
      fauxIMG$Scaffold.Name <- as.character(gbTemp$seq_id[1]) 
    }
    
        ## Scaffold Length (bp) - should always be present.
##    fauxIMG$Scaffold.Length..bp. <- as.numeric(GenomeInfoDb::seqlengths(neighborSubset))
    
    if (gbTemp$type[1] == "region" && gbTemp$end[1] != "NA") {
      fauxIMG$Scaffold.Length..bp. <- as.numeric(gbTemp$end[1])
    } else {
      lastnum <- length(gbTemp$type)
      fauxIMG$Scaffold.Length..bp. <- as.character(gbTemp$end[lastnum]) 
    }
    ## Scaffold GC % - unlikely to be present, not a priority to calculate anyway
    fauxIMG$Scaffold.GC.. <- ""
    ## Scaffold Read Depth - unlikely to be present
    fauxIMG$Scaffold.Read.Depth <- ""
    ## Most annotation is honestly better added via InterProScan later
    ## COG - could possibly add at least to some EMBL files, on the same terms as InterPRo or TIGRfam?
    ## For now, defaulting to blank.
    fauxIMG$COG <- ""
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
    ##
    ## Searching the inference, note, and db_xref fields (if they exist) for more annotations
    ##
    ## (one of the more important but less homogenous sorts of data unfortunately)
    numGenes <- length(fauxIMG$Locus.Tag)
    ## using a for loop here is gonna increase processing time
    ## BUT it also prevents certain kinds of failure modes when these complex entries are weird list and vector formats
    ## so it is probably the lesser evil
    fauxIMG$Transmembrane.Helices <- ""
    fauxIMG$Signal.Peptides <- ""
    fauxIMG$Pfam <- ""
    fauxIMG$Tigrfam <- ""
    fauxIMG$InterPro <- ""
    
    ## antismash note: because of the complex nested features, protein families are generally associated with things other than the CDS, and
    ## deconvolution is actually a pain.

    for (j in 1:numGenes)  {
      tmAnno <- FALSE
      spAnno <- FALSE
      pfAnno <- FALSE
      tfAnno <- FALSE
      ipAnno <- FALSE
      if (antismash == TRUE) {
        ## is there a pfam_domain under dbxref
        pfdLoc <- which(gbTemp$type == "PFAM_domain" & gbTemp$start >= neighborSubset[j,]$start & gbTemp$end <= neighborSubset[j,]$end)
      if (length(pfdLoc) > 0 && any(grepl("pfa?m?[[:digit:]]",gbTemp$dbxref[[pfdLoc]], ignore.case=TRUE))) {
#          pfAsIdx <- grep("pfa?m?[[:digit:]]",neighborSubset$dbxref[[pfdLoc]], ignore.case=TRUE)
#         pfIdx <- grep("pfa?m?[[:digit:]]",neighborSubset$note[[j]], ignore.case=TRUE)
        ## this nasty fellow looks for pf (capitalization-indifferent), up to 2 additional characters, and 5 numbers
        ## it then collapses multiple hits (if any) into a space-separated string of accessions
          tmpPfam <- unlist(regmatches(gbTemp$dbxref[[pfdLoc]], gregexpr("(?i)P(?i)F\\p{L}*\\p{L}*[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", gbTemp$dbxref[[pfdLoc]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$Pfam[j] <- paste(unlist(tmpPfam), collapse=" ")
 #          fauxIMG$Pfam[j] <- paste(as.character(neighborSubset$inference[[j]])[pfIdx], collapse=" ")
          fauxIMG$Pfam[j] <- gsub("PF"," pfam", fauxIMG$Pfam[j])
          fauxIMG$Pfam[j] <- gsub("PFAM"," pfam", fauxIMG$Pfam[j])
 #          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]][[:digit:]]", " ", fauxIMG$Pfam[j])
#          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]]", " ", fauxIMG$Pfam[j])
#          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]]", " ", fauxIMG$Pfam[j])
          pfAnno <- TRUE
        }
        ## is there a tigrfam under a_s_domain
        tigrLoc <- which(gbTemp$type == "aSDomain" & gbTemp$start >= neighborSubset[j,]$start & gbTemp$end <= neighborSubset[j,]$end)
        if (length(tigrLoc) > 0 && any(grepl("TIGR[[:digit:]]", gbTemp$identifier[[tigrLoc]], ignore.case=TRUE))) {
 #         tfAsIdx <- grep("TIGR[[:digit:]]",neighborSubset$a_s_domain[[tigrLoc]], ignore.case=TRUE)
          #          tfIdx <- grep("TIGR[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          tmpTigr <- unlist(regmatches(gbTemp$identifier[[tigrLoc]], gregexpr("(?i)t(?i)i(?i)g(?i)r[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", gbTemp$identifier[[tigrLoc]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$Tigrfam[j]  <- paste(unlist(tmpTigr), collapse=" ")
          #          fauxIMG$Tigrfam[j] <- paste(as.character(neighborSubset$inference[[j]])[tfIdx], collapse=" ")
          #          fauxIMG$Tigrfam[j] <- gsub("TIGRfam:","",fauxIMG$Tigrfam[j], ignore.case=TRUE)               
          tfAnno <- TRUE
        } 
          ## currently no IPR but I'll add it under aSDomain in case it shows up
          if (length(tigrLoc) > 0 && any(grepl("IPR[[:digit:]]", gbTemp$identifier[[tigrLoc]], ignore.case=TRUE))) {          
          #          ipIdx <- grep("IPR[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          tmpIpr <- unlist(regmatches(gbTemp$identifier[[tigrLoc]], gregexpr("(?i)i(?i)p(?i)r[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", gbTemp$identifier[[tigrLoc]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$InterPro[j] <- paste(unlist(tmpIpr), collapse=" ")
          #          fauxIMG$InterPro[j] <- paste(as.character(neighborSubset$inference[[j]])[ipIdx], collapse=" ")
          #          fauxIMG$InterPro[j] <- gsub("InterPro:","",fauxIMG$InterPro[j], ignore.case=TRUE)
          ipAnno <- TRUE
        } 
      } else {
      if (length(which(colnames(neighborSubset)=="inference"))!=0) {
        ## check for a note field in inference
        ## tm helices
        ## probably via TMHMM if they exist, so...
        if (any(grepl("TMHMM", neighborSubset$inference[[j]], ignore.case=TRUE))) {
          ## who knows what horrible way these might manifest, better just to say something was detected
          fauxIMG$Transmembrane.Helices[j] <- "yes"
          tmAnno <- TRUE
        }
        ## signal peptides
        ## probably via SignalP if they are here
        if (any(grepl("SignalP", neighborSubset$inference[[j]], ignore.case=TRUE))) {
          fauxIMG$Signal.Peptides[j] <- "yes"
          spAnno <- TRUE
        }              
        ## pfam
        ## n.b. requiring digits to avoid "prediction via PFAM" or whatever hits
        ## also keeps the name format more IMG-like
        if (any(grepl("pfa?m?[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE))) {
#          pfIdx <- grep("pfa?m?[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          ## hideous regex matches case-insensitve pf, 0-2 additional characters, and 5 numbers
          ## the results'll get concatenated into a space-separated string
          tmpPfam <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)P(?i)F\\p{L}*\\p{L}*[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$Pfam[j] <- paste(unlist(tmpPfam), collapse=" ")
#          fauxIMG$Pfam[j] <- paste(as.character(neighborSubset$inference[[j]])[pfIdx], collapse=" ")
          fauxIMG$Pfam[j] <- gsub("PF"," pfam", fauxIMG$Pfam[j])
          fauxIMG$Pfam[j] <- gsub("PFAM","pfam", fauxIMG$Pfam[j])
#          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]][[:digit:]]", " ", fauxIMG$Pfam[j])
#          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]]", " ", fauxIMG$Pfam[j])
          pfAnno <- TRUE
                  }                    
        ## tigrfam
        ## ID format is pretty standard for these, but still requiring digits for the same reason
        if (any(grepl("TIGR[[:digit:]]", neighborSubset$inference[[j]], ignore.case=TRUE))) {
#          tfIdx <- grep("TIGR[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          tmpTigr <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)t(?i)i(?i)g(?i)r[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$Tigrfam[j]  <- paste(unlist(tmpTigr), collapse=" ")
#          fauxIMG$Tigrfam[j] <- paste(as.character(neighborSubset$inference[[j]])[tfIdx], collapse=" ")
#          fauxIMG$Tigrfam[j] <- gsub("TIGRfam:","",fauxIMG$Tigrfam[j], ignore.case=TRUE)               
          tfAnno <- TRUE
        } 
        ## interpro
        ## as with pfam and tigrfam
        if (any(grepl("IPR[[:digit:]]", neighborSubset$inference[[j]], ignore.case=TRUE)))  {
#          ipIdx <- grep("IPR[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          tmpIpr <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)i(?i)p(?i)r[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$InterPro[j] <- paste(unlist(tmpIpr), collapse=" ")
#          fauxIMG$InterPro[j] <- paste(as.character(neighborSubset$inference[[j]])[ipIdx], collapse=" ")
#          fauxIMG$InterPro[j] <- gsub("InterPro:","",fauxIMG$InterPro[j], ignore.case=TRUE)
          ipAnno <- TRUE
        }
      }
      if (length(which(colnames(neighborSubset)=="note"))!=0) {
        ## check for a note field under note
        ## tm helices
        ## probably via TMHMM if they exist, so...
        if (any(grepl("TMHMM", neighborSubset$note[[j]], ignore.case=TRUE)) && tmAnno == FALSE) {
          ## who knows what horrible way these might manifest, better just to say something was detected
          fauxIMG$Transmembrane.Helices[j] <- "yes"
          tmAnno <- TRUE
        }
        ## signal peptides
        ## probably via SignalP if they are here
        if (any(grepl("SignalP", neighborSubset$note[[j]], ignore.case=TRUE)) && spAnno == FALSE) {
          fauxIMG$Signal.Peptides[j] <- "yes"
          spAnno <- TRUE
        }              
        ## pfam
        ## n.b. requiring digits to avoid "prediction via PFAM" or whatever hits
        ## also keeps the name format more IMG-like
        if (any(grepl("pfa?m?[[:digit:]]",neighborSubset$note[[j]], ignore.case=TRUE)) && pfAnno == FALSE) {
 #         pfIdx <- grep("pfa?m?[[:digit:]]",neighborSubset$note[[j]], ignore.case=TRUE)
          tmpPfam <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)P(?i)F\\p{L}*\\p{L}*[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$Pfam[j] <- paste(unlist(tmpPfam), collapse=" ")
          #          fauxIMG$Pfam[j] <- paste(as.character(neighborSubset$inference[[j]])[pfIdx], collapse=" ")
          fauxIMG$Pfam[j] <- gsub("PF"," pfam", fauxIMG$Pfam[j])
          fauxIMG$Pfam[j] <- gsub("PFAM","pfam", fauxIMG$Pfam[j])
          #          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]][[:digit:]]", " ", fauxIMG$Pfam[j])
          #          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]]", " ", fauxIMG$Pfam[j])
          pfAnno <- TRUE
        }                    
        ## tigrfam
        ## ID format is pretty standard for these, but still requiring digits for the same reason
        if (any(grepl("TIGR[[:digit:]]", neighborSubset$note[[j]], ignore.case=TRUE)) && tfAnno == FALSE) {
          #          tfIdx <- grep("TIGR[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          tmpTigr <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)t(?i)i(?i)g(?i)r[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$Tigrfam[j]  <- paste(unlist(tmpPfam), collapse=" ")
          #          fauxIMG$Tigrfam[j] <- paste(as.character(neighborSubset$inference[[j]])[tfIdx], collapse=" ")
          #          fauxIMG$Tigrfam[j] <- gsub("TIGRfam:","",fauxIMG$Tigrfam[j], ignore.case=TRUE)               
          tfAnno <- TRUE
        } 
        ## interpro
        ## as with pfam and tigrfam
        if (any(grepl("IPR[[:digit:]]", neighborSubset$note[[j]], ignore.case=TRUE)) && ipAnno == FALSE)  {
          #          ipIdx <- grep("IPR[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          tmpIpr <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)i(?i)p(?i)r[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$InterPro[j] <- paste(unlist(tmpIpr), collapse=" ")
          #          fauxIMG$InterPro[j] <- paste(as.character(neighborSubset$inference[[j]])[ipIdx], collapse=" ")
          #          fauxIMG$InterPro[j] <- gsub("InterPro:","",fauxIMG$InterPro[j], ignore.case=TRUE)
          ipAnno <- TRUE
        }
      }
      if (length(which(colnames(neighborSubset)=="db_xref"))!=0) {
        ## check for a note field
        ## tm helices
        ## probably via TMHMM if they exist, so...
        if (any(grepl("TMHMM", neighborSubset$db_xref[[j]], ignore.case=TRUE)) && tmAnno == FALSE) {
          ## who knows what horrible way these might manifest, better just to say something was detected
          fauxIMG$Transmembrane.Helices[j] <- "yes"
          tmAnno <- TRUE
        }
        ## signal peptides
        ## probably via SignalP if they are here
        if (any(grepl("SignalP", neighborSubset$db_xref[[j]], ignore.case=TRUE)) && spAnno == FALSE) {
          fauxIMG$Signal.Peptides[j] <- "yes"
          spAnno <- TRUE
        }              
        ## pfam
        ## n.b. requiring digits to avoid "prediction via PFAM" or whatever hits
        ## also keeps the name format more IMG-like
        if (any(grepl("pfa?m?[[:digit:]]",neighborSubset$db_xref[[j]], ignore.case=TRUE)) && pfAnno == FALSE) {
          #         pfIdx <- grep("pfa?m?[[:digit:]]",neighborSubset$note[[j]], ignore.case=TRUE)
          tmpPfam <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)P(?i)F\\p{L}*\\p{L}*[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$Pfam[j] <- paste(unlist(tmpPfam), collapse=" ")
          #          fauxIMG$Pfam[j] <- paste(as.character(neighborSubset$inference[[j]])[pfIdx], collapse=" ")
          fauxIMG$Pfam[j] <- gsub("PF"," pfam", fauxIMG$Pfam[j])
          fauxIMG$Pfam[j] <- gsub("PFAM","pfam", fauxIMG$Pfam[j])
          #          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]][[:digit:]]", " ", fauxIMG$Pfam[j])
          #          fauxIMG$Pfam[j] <- gsub("\\.[[:digit:]]", " ", fauxIMG$Pfam[j])
          pfAnno <- TRUE
        }                    
        ## tigrfam
        ## ID format is pretty standard for these, but still requiring digits for the same reason
        if (any(grepl("TIGR[[:digit:]]", neighborSubset$db_xref[[j]]), ignore.case=TRUE) && tfAnno == FALSE) {
          #          tfIdx <- grep("TIGR[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          tmpTigr <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)t(?i)i(?i)g(?i)r[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$Tigrfam[j]  <- paste(unlist(tmpTigr), collapse=" ")
          #          fauxIMG$Tigrfam[j] <- paste(as.character(neighborSubset$inference[[j]])[tfIdx], collapse=" ")
          #          fauxIMG$Tigrfam[j] <- gsub("TIGRfam:","",fauxIMG$Tigrfam[j], ignore.case=TRUE)               
          tfAnno <- TRUE
        } 
        ## interpro
        ## as with pfam and tigrfam
        if (any(grepl("IPR[[:digit:]]", neighborSubset$db_xref[[j]], ignore.case=TRUE)) && ipAnno == FALSE)  {
          #          ipIdx <- grep("IPR[[:digit:]]",neighborSubset$inference[[j]], ignore.case=TRUE)
          tmpIpr <- unlist(regmatches(neighborSubset$inference[[j]], gregexpr("(?i)i(?i)p(?i)r[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]", neighborSubset$inference[[j]], perl=TRUE, ignore.case=TRUE)))
          fauxIMG$InterPro[j] <- paste(unlist(tmpIpr), collapse=" ")
          #          fauxIMG$InterPro[j] <- paste(as.character(neighborSubset$inference[[j]])[ipIdx], collapse=" ")
          #          fauxIMG$InterPro[j] <- gsub("InterPro:","",fauxIMG$InterPro[j], ignore.case=TRUE)
          ipAnno <- TRUE
        }
      }
      }
      
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
    ## probably defunct with a gggenomes switch
    if (seqExtract == TRUE && length(which(colnames(neighborSubset)=="translation"))!=0)  {
      fauxSeqIMG <-   neighborSubset$translation
      fauxSeqIMG <- as.data.frame(fauxSeqIMG)
      colnames(fauxSeqIMG) <- c("sequence")
      fauxSeqIMG$gene_oid <- as.character(fauxIMG$gene_oid)
      fauxSeqIMG <- fauxSeqIMG %>% dplyr::select("gene_oid", dplyr::everything())
      if (exists(x="fauxNeighborSeqs")==FALSE) {
        fauxNeighborSeqs <- fauxSeqIMG
      } else {
        fauxNeighborSeqs <- rbind(fauxNeighborSeqs, fauxSeqIMG)
      }
    } else if (seqExtract == TRUE) {
      ## gggenomes does not import translations currently
      ## so we need to grab them via seqinr.
      ## this means we have to have generated the .faa previously
      seqFile <- basename(fileList[i])
      seqFile <- gsub("\\.gbk", ".faa", seqFile)
      seqFile <- gsub("\\.gb", ".faa", seqFile)
      ## this'll look for whatever went into Locus.Tag (originally locus_tag, protein_id, or whatever) in fasta seq IDs
      ## and add the seq to the genes/neighbors fasta under the new gene_oid
      extSeqs <- seqinr::read.fasta(file=seqFile, seqtyp="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
      fauxSeqIMG <- fauxIMG$gene_oid
      fauxSeqIMG <- as.data.frame(fauxSeqIMG)
      colnames(fauxSeqIMG) <- c("gene_oid")
      fauxSeqIMG$sequence <- ""
                  for (j in 1:length(fauxIMG$gene_oid)) {
              
              seqIdx <-  grep(fauxIMG$Locus.Tag[j],names(extSeqs))
              fauxSeqIMG$sequence[j] <- extSeqs[[seqIdx]]
      }
      ## same initialization check
      if (exists(x="fauxNeighborSeqs")==FALSE) {
        fauxNeighborSeqs <- fauxSeqIMG
      } else {
        fauxNeighborSeqs <- rbind(fauxNeighborSeqs, fauxSeqIMG)
      }
    }
    processedList <- append(processedList, fileList[i])
    tempWarning <- paste("Processed: ", fileList[i],". Moving on.", sep="")
    print(tempWarning)
    rm(neighborSubset)
  }
  ## outputting the faux-IMG-formatted metadata file for the neighborhoods
  if (removeDupes == TRUE) {
    fauxNeighborData <- fauxNeighborData[!duplicated(fauxNeighborData$Locus.Tag), ]
  }
  fauxNeighborData$Pfam <- as.character(fauxNeighborData$Pfam)
  fauxNeighborData$Tigrfam <- as.character(fauxNeighborData$Tigrfam)
  fauxNeighborData$InterPro <- as.character(fauxNeighborData$InterPro)
  fauxNeighborData$gene_oid <- as.character(fauxNeighborData$gene_oid)
  fauxNeighborData$Scaffold.ID <- as.character(fauxNeighborData$Scaffold.ID)
  fauxNeighborData$Genome.ID <- as.character(fauxNeighborData$Genome.ID)
  ## let's make sure the column order of the metadata matches the IMG order
  if (includeIPR == FALSE) {
    colOrder = c("gene_oid",
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
                 "Locus.Type",
                 "Is.Pseudogene",
                 "Is.Partial.Gene",
                 "Add.Date",
                 "Is.Public",
                 "Sequencing.Status",
                 "Transmembrane.Helices",
                 "Signal.Peptides",
                 "Scaffold.ID",
                 "Scaffold.External.Accession",
                 "Scaffold.Name",
                 "Scaffold.Length..bp.",
                 "Scaffold.GC..",
                 "Scaffold.Read.Depth",
                 "COG",
                 "Pfam",
                 "Tigrfam",
                 "SMART.ID",
                 "SUPERFam.ID",
                 "CATH.FunFam.ID",
                 "Enzyme",
                 "KO",
                 "IMG.Term")
    fauxNeighborData <- fauxNeighborData %>% dplyr::select(dplyr::all_of(colOrder))
  } else {   
    colOrder = c("gene_oid",
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
                 "Locus.Type",
                 "Is.Pseudogene",
                 "Is.Partial.Gene",
                 "Add.Date",
                 "Is.Public",
                 "Sequencing.Status",
                 "Transmembrane.Helices",
                 "Signal.Peptides",
                 "Scaffold.ID",
                 "Scaffold.External.Accession",
                 "Scaffold.Name",
                 "Scaffold.Length..bp.",
                 "Scaffold.GC..",
                 "Scaffold.Read.Depth",
                 "COG",
                 "Pfam",
                 "Tigrfam",
                 "SMART.ID",
                 "SUPERFam.ID",
                 "CATH.FunFam.ID",
                 "Enzyme",
                 "KO",
                 "IMG.Term",
                 "InterPro")
    fauxNeighborData <- fauxNeighborData %>% dplyr::select(dplyr::all_of(colOrder))
  }
  write.table(fauxNeighborData, file=fauxNeighborDataFile, row.names=FALSE, col.names = TRUE, sep="\t")
  ## writing the context file too, with the correct column order
  ctxtOrd <- c("gene_oid", "source_gene_oid", "source_scaffold_id")
  fullContext <- fullContext %>% dplyr::select(dplyr::all_of(ctxtOrd))
  write.table(fullContext, file=fauxContextFile, row.names=FALSE, col.names = TRUE, sep="\t")
  ## and the .fa file for the neighbor protein sequences
if (seqExtract == TRUE) {
    fauxNeighborSeqList <- list()
  emptyNeighbors <- which(fauxNeighborSeqs[[2]]=="")
  fauxNeighborSeqs[[1]][emptyNeighbors] <- "No sequence found"
  for (k in 1:length(fauxNeighborSeqs[,2])) {
    fauxNeighborSeqList[[fauxNeighborSeqs[k,1]]] <- fauxNeighborSeqs[k,2]
  }
  seqinr::write.fasta(fauxNeighborSeqList, names=names(fauxNeighborSeqList),file.out=fauxNeighborSeqsFile)
  }
  ## and outputting the genes of interest separately
  fauxGeneData <- fauxNeighborData %>% dplyr::filter(.data$Locus.Tag %in% goiList$locus_tag)
  if (removeDupes == TRUE) {
    fauxGeneData <- fauxGeneData[!duplicated(fauxGeneData$Locus.Tag), ]
  }
  write.table(fauxGeneData, file=fauxGeneDataFile, row.names=F, col.names = T, sep="\t")
  ## and their sequences too
  fauxGeneList <- fauxGeneData$gene_oid
  if (seqExtract == TRUE) {
  fauxGeneSeqs <- fauxNeighborSeqs %>% dplyr::filter(.data$gene_oid %in% fauxGeneList)
  fauxGeneSeqList <- list()
  for (m in 1:length(fauxGeneSeqs[,2])) {
    fauxGeneSeqList[[fauxGeneSeqs[m,1]]] <- fauxGeneSeqs[m,2]
  }
  seqinr::write.fasta(fauxGeneSeqList, names=names(fauxGeneSeqList),file.out=fauxGeneSeqsFile)
  }
  ## make files for the lists of what did and didn't process correctly (what lacked annotations, a GoI, was antiSMASH-formatted, etc.)
  if (length(noAnnotList) != 0) {
    noAnnotList <- as.data.frame(unlist(noAnnotList))
    colnames(noAnnotList) <- c("unannotated")
    write.table(noAnnotList, file=noAnnotFile, row.names=F, col.names = T, sep="\t")
  }
  if (length(noGoiList) != 0) {
    noGoiList <- as.data.frame(unlist(noGoiList))
    colnames(noGoiList) <- c("no-GoI")
    write.table(noGoiList, file=noGoiFile, row.names=F, col.names = T, sep="\t")
  }
  if (length(asFmtList) != 0) {
    asFmtList <- as.data.frame(unlist(asFmtList))
    colnames(asFmtList) <- c("antismash")
    write.table(asFmtList, file=asFmtFile, row.names=F, col.names = T, sep="\t")
  }
  if (length(processedList) != 0) {
    processedList <- as.data.frame(unlist(processedList))
    colnames(processedList) <- c("processed")
    write.table(processedList, file=processedFile, row.names=F, col.names = T, sep="\t")
  }
  print("GenBank-formatted files in folder processed.")
  return(list(geneData=fauxGeneData, neighborData=fauxNeighborData, neighborsContext=fullContext))
} 