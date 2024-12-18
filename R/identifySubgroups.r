#' Standalone function to identify subgroups of proteins in larger IMG datasets, based on sequence
#'
#' Given a desired genome neighborhood and lists of genes of interest and their neighbors, along with their protein sequences, this program uses all-by-all-blast & a network-based clustering method to identify subsets of proteins that may represent meaningful subgroups.
#' @param geneList File location for list of subset of genes to analyze (single column with gene_oid header). String.
#' @param imgNeighbors File location for file with metadata for neighbors of your genes of interest.  String.
#' @param imgNeighborSeqs File location for fasta-formatted file with sequence information for neighbors of your gene of interest.  String.
#' @param geneName Character string with your gene of interest's name - single word preferred.  Character.
#' @param subgroupDesc Character string with the description of gene subset - single word (incl. hyphens) preferred.  Character.
#' @param cutoffType The cutoff type being used - "identity" or "evalue".  Character.
#' @param cutoffValue The percent ID or evalue cutoff after which edges are not formed when using tidygraph.  Number.
#' @param sysTerm The type of terminal (wsl vs linux/unix/macos) from which blastp and other commands will be run.  Character.
#' @param numThreads The number of processor threads to be devoted to certain steps. Number, defaults to 1.
#' @param defFamNum Starting number for families. Useful when running several consecutive protein types. Number, defaults to 0.
#' @param lightExport Indicates whether a simplified export format - gene name and family only - should be used. Boolean, defaults to FALSE.
#' @param screenPep Should peptide-friendly settings be used? T/F, defaults to FALSE.
#' @param alnClust Should we make MAFFT alignments for clusters?  T/F, defaults to FALSE.
#' @param hmmClust Should we make HMM models for clusters? T/F, defaults to FALSE.
#' @return Updated metadata for neighboring genes (additional files generated en route)
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @importFrom utils capture.output combn read.csv write.csv write.table
#' @importFrom stats dist end hclust na.omit start
#' @importFrom grDevices cairo_pdf cairo_ps col2rgb colorRampPalette dev.off hsv pdf rgb2hsv
#' @examples 
#' \dontrun{
#' identifySubgroupsOut <- identifySubgroups(geneList = "geneList.txt", 
#'                                           imgNeighbors = "imgNeighbors.txt", 
#'                                           imgNeighborSeqs = "imgneighborSeqs.fa", 
#'                                           geneName = "genE", 
#'                                           subgroupDesc="enzymes", 
#'                                           cutoffType = "identity", 
#'                                           cutoffValue = 40, 
#'                                           sysTerm = "nix", 
#'                                           numThreads = 8)
#' }
#' 
identifySubgroups <- function(geneList = geneList,
                              imgNeighbors = imgNeighbors,
                              imgNeighborSeqs = imgNeighborSeqs,
                              geneName = geneName,
                              subgroupDesc = subgroupDesc,
                              cutoffType = cutoffType,
                              cutoffValue = cutoffValue,
                              sysTerm = sysTerm,
                              numThreads = 1,
                              alnClust = FALSE,
                              hmmClust = FALSE,
                              defFamNum = 0,
                              lightExport = FALSE,
                              screenPep = FALSE) { 
    ## if key things are missing
    if(exists(x="imgNeighbors") == FALSE | exists(x="imgNeighborSeqs") == FALSE | exists(x="geneName") == FALSE | exists(x="geneList") == FALSE | exists(x="subgroupDesc") == FALSE | exists(x="cutoffType") == FALSE | exists(x="cutoffValue") == FALSE | exists(x="sysTerm") == FALSE) {
        print("Missing a required term")
        return(0)
    }
    ## stage-setting
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    fileName <- paste(fileDate,"_identifySubgroups_",geneName,"_",subgroupDesc,sep="")
    subSeqsFile <- paste(fileName, "_subgroupSeqs.fa",sep="")
    blastFile <- paste(fileName, "_Blast.txt",sep="")
    blastError <- paste(fileName, "_BlastError.txt",sep="")
    clustListFile <- paste(fileName, "_clusterList.txt",sep="")
    networkFileName <- paste(fileName, "_networkFull_cutoff_",cutoffValue,".pdf",sep="")
    gmlFileName <- paste(fileName, "_networkFull_cutoff_",cutoffValue,".gml",sep="")
    outFileName <- paste(fileName,"_metadata_cutoff_",cutoffValue,".txt",sep="")
    subIndex <- list()
    counter <- 1
    sSeqTable <- data.frame()
    subSeqs <- data.frame()
    ssReList <- list()
    tidyMonoBlast <- data.frame()
    clustListings <- data.table::data.table(stringsAsFactors = FALSE)
    ## tool check: do we have blast and mafft?
    ## note: assuming use of linux, macos, or wsl on windows
    ## wsl is easy to setup and native windows installs are harder to deal with, so
    ## gonna be lazy
    ## ok, checking on blast
    if(exists(x="blastp") == FALSE || exists(x="makeblastdb") == FALSE) {
        if (sysTerm == "wsl") {
            blastp <- system2(command = "wsl",
                              args = c("which",
                                       "blastp"),
                              stdout=TRUE)
            makeblastdb <-  system2(command = "wsl",
                                    args = c("which",
                                             "makeblastdb"),
                                    stdout=TRUE)
        } else if (sysTerm == "nix") {
            blastp <- system2(command = "which",
                              args = c("blastp"),
                              stdout = TRUE)
            makeblastdb <-  system2(command = "which",
                                    args = c("makeblastdb"),
                                    stdout = TRUE)
        } else {
            print("blast path absent or unrecognized")
            return(0)
        }
    }
    ## checking on mafft IF we've decided to align
    if (alnClust == TRUE) {
        if(exists(x="mafft") == FALSE) {
            if (sysTerm == "wsl") {
                mafft <- system2(command = "wsl",
                                 args = c("which",
                                          "mafft"),
                                 stdout = TRUE)
            } else if (sysTerm == "nix") {
                mafft <- system2(command = "which",
                                 args = c("mafft"),
                                 stdout=TRUE)
            } else {
                print("mafft path absent or unrecognized")
                return(0)
            }
        }
    }
    ## checking on hmmer tools IF we've decided to make and use HMMs
    if (hmmClust == TRUE) {
        if(exists(x="hmmalign") == FALSE) {
            if (sysTerm == "wsl") {
                mafft <- system2(command = "wsl",
                                 args = c("which",
                                          "hmmalign"),
                                 stdout = TRUE)
            } else if (sysTerm == "nix") {
                mafft <- system2(command = "which",
                                 args = c("hmmalign"),
                                 stdout=TRUE)
            } else {
                print("hmmalign path absent or unrecognized")
                return(0)
            }
        }
        if(exists(x="hmmbuild") == FALSE) {
            if (sysTerm == "wsl") {
                mafft <- system2(command = "wsl",
                                 args = c("which",
                                          "hmmbuild"),
                                 stdout = TRUE)
            } else if (sysTerm == "nix") {
                mafft <- system2(command = "which",
                                 args = c("hmmbuild"),
                                 stdout=TRUE)
            } else {
                print("hmmbuild path absent or unrecognized")
                return(0)
            }
        }
    }
    ## no to scientific notation
    options(scipen = 999)
    ## file input, starting with the expected IMG columns
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
    ## if interpro and hypofam exist, add them
    ## might want to add a "select families" version of this later
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
    ## also, just in case this is post-analyzeNeighbors and/or repnodeTrim
    ## these aren't needed for anything, necessarily, but it might be handy to have 'em in the output
    if (any(grepl("clustNum", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$clustNum <- imgNeighborsTemp$clustNum
    } 
    if (any(grepl("clustOrd", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$clustOrd <- imgNeighborsTemp$clustOrd
    }
    if (any(grepl("source_gene_oid", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$source_gene_oid <- imgNeighborsTemp$source_gene_oid
    }
    if (any(grepl("source_scaffold_id", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$source_scaffold_id <- imgNeighborsTemp$source_scaffold_id
    }
    if (any(grepl("efi_oid", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$efi_oid <- imgNeighborsTemp$efi_oid
    }
    if (any(grepl("isRepnode", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$isRepnode <- imgNeighborsTemp$isRepnode
    }
    if (any(grepl("repnodeIs", colnames(imgNeighborsTemp)))) {
        imgNeighborsTrimmed$repnodeIs <- imgNeighborsTemp$repnodeIs
    }
    ## trimming to the subgroup of genes of interest
    subGeneList <- read.csv(file=geneList, header=TRUE, sep=",", stringsAsFactors=FALSE)
    subNeighborsData <- imgNeighborsTrimmed[imgNeighborsTrimmed$gene_oid %in% subGeneList$gene_oid, ]
    ## taking care of the sequences themselves
    neighborSeqs <- seqinr::read.fasta(file=imgNeighborSeqs, seqtyp="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
    sSeqNames <- names(neighborSeqs)
    for (i in 1:length(sSeqNames)) {
        sSeqTable[i,1] <- sSeqNames[i]
        sSeqTable[i,2] <- neighborSeqs[[i]]
    }
    colnames(sSeqTable) <- c("gene_oid", "sequence")
    ## sometimes neighbors are rRNA or tRNA or whatever and don't have protein sequences
    ## and if this didn't go through analyzeNeighbors they might still be around
    badSeq <- "No sequence found"
    goodSeqs <- sSeqTable %>% dplyr::filter(.data$sequence != badSeq)
    subSeqs <- goodSeqs %>% dplyr::filter(.data$gene_oid %in% subGeneList$gene_oid)
    ssReList <- list()
    ## export the subgroup sequences separately after re-fasta-ifying
    for (i in 1:length(subSeqs[,1])) {
        ssReList[[subSeqs[i,1]]] <- subSeqs[i,2] 
    }
    subSeqs <- ssReList
    seqinr::write.fasta(ssReList,names=names(ssReList), file.out=subSeqsFile)
    ## preparing for all-by-all blast to ID similar (sub)sets of proteins
    if(dir.exists("subtemporary") == FALSE) {
        dir.create("subtemporary")
    }
    if(dir.exists("subtemporary/dbblast") == FALSE) {
        dir.create("subtemporary/dbblast")
    }
    dbtype <- "prot"
    blast_db <- "subSeqsDb"
    input <- subSeqsFile
    ## easier to cut down later than rerun - note altered settings for peptides
    if (screenPep == TRUE) {evalue <- "10"} else {evalue <- "1"}
    format <- "6"
    output <- "subSeqsDb"
    query <- subSeqsFile
    blastColNames <- c("qseqid", "sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
    print("Beginning all-by-all blast of proteins in your chosen subgroup.")
    ## when running under a wsl or unix/linux environment
    if (sysTerm == "wsl") {
        system2(command = "wsl",
                args = c(makeblastdb,
                         "-dbtype", "prot",
                         "-in", input,
                         "-out",  "subtemporary/dbblast/subSeqsDb"),
                stdout=FALSE)
        ## times out while waiting for stdout sometimes on very large datasets - let's see if stderr works better?
        blast_err <- system2(command = "wsl", 
                             args = c(blastp, 
                                      "-db", "subtemporary/dbblast/subSeqsDb",
                                      "-query", query,
                                      "-outfmt", "6",
                                      "-evalue", evalue,
                                      "-num_threads", numThreads,
                                      "-out", blastFile),
                             wait = TRUE,
                             stdout = FALSE,
                             stderr = TRUE)
        capture.output(blast_err, blastError)
        blast_out <- as.data.frame(read.csv(blastFile, header=FALSE, sep="\t" , stringsAsFactors=FALSE))
        colnames(blast_out) <- blastColNames
        tidyBlast <- tibble::as_tibble(blast_out)
    } else if ( sysTerm == "nix") {
        system2(command = makeblastdb,
                args = c("-dbtype", "prot",
                         "-in", input,
                         "-out", "subtemporary/dbblast/subSeqsDb"),
                stdout=FALSE)
        ## possibly implement stderr variant here too?  need to repeat more and see whether this is a more general issue.
        ## or see if the stdout file switch works better on wsl too...?
        system2(command = blastp, 
                args = c("-db", "subtemporary/dbblast/subSeqsDb",
                         "-query", query,
                         "-outfmt", "6",
                         "-evalue", evalue,
                         "-num_threads", numThreads,
                         "-out", blastFile),
                stdout = blastFile)
        blast_out <- as.data.frame(read.csv(blastFile, header=FALSE, sep="\t" , stringsAsFactors=FALSE))
        colnames(blast_out) <- blastColNames
        tidyBlast <- tibble::as_tibble(blast_out)
    }
    ## moving the pairwise results into a matrix form
    ## and removing things other than evalue and gene_oid interactions
    print("All-by-all blast of proteins in your subgroup complete.")
    rm(blast_out)
    tidyBlast <- as.data.frame(tidyBlast)
    ## fixing occasional weird issues with the blast info types
    tidyBlast$qseqid <- as.character(tidyBlast$qseqid)
    tidyBlast$sseqid <- as.character(tidyBlast$sseqid)
    tidyBlast$length <- as.numeric(tidyBlast$length)
    tidyBlast$pident <- as.numeric(tidyBlast$pident)
    tidyBlast$bitscore <- as.numeric(tidyBlast$bitscore)
    tidyBlast$evalue <- as.numeric(tidyBlast$evalue)
    ##  dealing with NA issues that're sometimes gonna happen - population should generally be small enough to be deletable
    tidyBlast <- na.omit(tidyBlast)
    ## making the dataset skinnier - only IDs, percentID, length, eval, and bitscore are maintained later
    tidyBlast <- tidyBlast %>% dplyr::select("qseqid", "sseqid", "length", "pident", "bitscore", "evalue")
    ## remove full-on dupes - whyyyyy blast why
    tidyBlast <- tidyBlast %>% dplyr::distinct()
    ## even in blastfmt 6, you sometimes get stupid-ass duplicates (multiple hits in one gene and so on)
    ## here, we only pass on the unique qseqid-sseqid pairs with the highest bitscore (as a proxy for %ID/length combo)
    tidyBlast <- tidyBlast %>% dplyr::group_by(.data$qseqid, .data$sseqid) %>% dplyr::arrange(dplyr::desc(.data$bitscore), .by_group=TRUE)
    yetTidierBlast <- dplyr::distinct_at(tidyBlast, dplyr::vars(.data$qseqid, .data$sseqid), .keep_all=TRUE)
    yetTidierBlast <- yetTidierBlast %>% dplyr::group_by(.data$qseqid, .data$sseqid) %>% dplyr::arrange(dplyr::desc(.data$bitscore), .by_group=TRUE)
    ## a somewhat stringent length requirement - on the one hand, this may miss a few fusion proteins.  
    ## On the other, this avoids, like, FeS or heme motifs binding a tiny region really well.  
    ## Bitscore is OK for the relative calls, but not great for single matches with crappy length.
    ## anyway arbitrarily using 2/3s the seq length as a rule-out point
    ## and pegging it to the larger protein to prevent 30 aa peptides from nevertheless matching 3000 aa NRPS monsters
    yetTidierBlast$qseqLength <- as.numeric(nchar(subSeqs[as.character(yetTidierBlast$qseqid)]))
    yetTidierBlast$sseqLength <- as.numeric(nchar(subSeqs[as.character(yetTidierBlast$sseqid)]))
    ## need to implement this better - is it fast enough with mutate?
    yetTidierBlast <- yetTidierBlast %>% dplyr::rowwise() %>% dplyr::mutate(maxSeqLength = max(.data$qseqLength, .data$sseqLength))
    ## imposing a length cutoff (match is % of whichever is larger) for matches
    ## to avoid high-identity but short-length matches
    ## note modded values for peptides
    if (screenPep == TRUE) {
        tidyMonoBlast <- yetTidierBlast %>% dplyr::filter(.data$length>=.45*.data$maxSeqLength)   
    } else {
        tidyMonoBlast <- yetTidierBlast %>% dplyr::filter(.data$length>=.65*.data$maxSeqLength)   
    }
    ## keeping only names and percent id or evalue now that we are done with length and bitscore, then reshaping it
    rm(yetTidierBlast)
    if (cutoffType == "identity")  {
        tidyMinBlast <- tidyMonoBlast %>% dplyr::select("qseqid", "sseqid", "pident")
        tidyMinBlast <- as.data.frame(tidyMinBlast)
        tidyMinBlastTrimmed <- tidyMinBlast %>% dplyr::filter(.data$pident>=cutoffValue)
    } else if (cutoffType == "evalue") {
        tidyMinBlast <- tidyMonoBlast %>% dplyr::select("qseqid", "sseqid", "evalue")
        tidyMinBlastTrimmed <- tidyMinBlast %>% dplyr::filter(.data$evalue<=cutoffValue)
        tidyMinBlast <- as.data.frame(tidyMinBlast)
    } 
    rm(tidyMonoBlast)
    ## now onto tidygraph for the actual clustering
    ## make the object using the cutoff, and do so as an undirected network (deleting self-linking edges)
    tidyBlastNetwork <- tidygraph::as_tbl_graph(tidyMinBlastTrimmed, directed=FALSE)
    tidyBlastNetworkTrimmed <- tidyBlastNetwork %>% tidygraph::activate(edges) %>% dplyr::filter(.data$from != .data$to)  
    ## calculating some stuff - dunno if i'll keep all of these in the end
    tidyBlastNetworkTrimmed <- tidyBlastNetworkTrimmed %>% tidygraph::activate(nodes) %>% dplyr::mutate(group = tidygraph::group_infomap())
    tidyBlastNetworkTrimmed <- tidyBlastNetworkTrimmed %>% tidygraph::activate(nodes) %>% dplyr::mutate(cluster = tidygraph::group_components())
    tidyBlastNetworkTrimmed <- tidyBlastNetworkTrimmed %>% tidygraph::activate(edges) %>% dplyr::mutate(betweenness = tidygraph::centrality_edge_betweenness())
    # tidyBlastNetworkTrimmed <- tidyBlastNetworkTrimmed %>% tidygraph::activate(nodes) %>% dplyr::mutate(center_dist = tidygraph::node_distance_to(tidygraph::node_is_center()))
    nodeNum <- length(unique(tidyMinBlast$qseqid))
    kNum <- ceiling(nodeNum/20)
    tidyBlastNetworkTrimmed <- tidyBlastNetworkTrimmed %>% tidygraph::activate(nodes) %>% dplyr::mutate(keyplayer = tidygraph::node_is_keyplayer(k = kNum))
    ## removing singletons
    keepCluster <- which(table(igraph::V(tidyBlastNetworkTrimmed)$cluster) > 1)
    tidyBlastClusters <- tidyBlastNetworkTrimmed %>% tidygraph::activate(nodes) %>% dplyr::filter(.data$cluster %in% keepCluster)
    print("Singletons not shown.")
    ## a few things need to be done differently for e-value vs. %id data presentation
    if (screenPep == TRUE) {matchLength <- "45"} else {matchLength <- "65"}
    if (cutoffType == "identity")  {
        titleText <- paste("Cluster-based identification of protein subgroups, (",
                           subgroupDesc,
                           ") among the neighbors of ",
                           geneName,
                           ".",
                           sep="")
        subtitleText <- paste("Edge values represent %ID, with a cutoff of ",
                              cutoffValue,
                              " and a ",
                              matchLength,
                              "% match length used for edges. Analyzed on ",
                              fileDate,
                              ".",
                              sep="")
        tidyBlastNetworkPic <- ggraph::ggraph(tidyBlastClusters, layout="stress") + 
            ggraph::geom_edge_link0(ggplot2::aes(edge_alpha=.data$pident/100),
                                    edge_colour="black",
                                    edge_width=.5,
                                    show.legend=FALSE) +    
            ggraph::geom_node_point(ggplot2::aes(color=factor(.data$cluster)),
                                    size=7)
        ##note: could add size of node as a centrality or keyplayer measure
    } else if (cutoffType == "evalue") {
        titleText <- paste("Cluster-based identification of protein subgroups, (",
                           subgroupDesc,
                           ") among the neighbors of ",
                           geneName,
                           ".",
                           sep="")
        subtitleText <- paste("Edge values represent expectation value, with a cutoff of ",
                              cutoffValue,
                              " and a ",
                              matchLength,
                              "% match length used for edges. Analyzed on ",
                              fileDate,
                              ".",
                              sep="")
        minEval <- min(tidyMinBlastTrimmed$evalue)
        tMinEval <- log10(1/minEval)
        tidyBlastNetworkPic <- ggraph::ggraph(tidyBlastClusters, layout="stress") + 
            ggraph::geom_edge_link0(ggplot2::aes(edge_alpha=log10(1/.data$evalue)/tMinEval),
                                    edge_colour="black",
                                    edge_width=.5,
                                    show.legend=FALSE) +    
            ggraph::geom_node_point(ggplot2::aes(color=factor(.data$cluster)),
                                    size=7)
        ##note: could add size of node as a centrality or keyplayer measure
    }
    ##
    tidyBlastNetworkPic <- tidyBlastNetworkPic + 
        ggplot2::labs(title=titleText,
                      subtitle=subtitleText) +
        ggplot2::theme(plot.title=ggplot2::element_text(color="black",
                                                        size=12,
                                                        margin=ggplot2::margin(10,0,10,0)),
                       plot.subtitle=ggplot2::element_text(color="grey66",
                                                           size=10,                  
                                                           margin=ggplot2::margin(10,0,10,0))) +
        ggraph::theme_graph(base_family="sans") +
        ggplot2::theme(legend.position="bottom", plot.margin=ggplot2::unit(c(.2,.2,.2,.2), "cm")) +
        ggplot2::scale_color_viridis_d()
    ## pic of the full-sized cluster
    ggplot2::ggsave(filename=networkFileName, tidyBlastNetworkPic, height=10, width=20, dpi=75, units="in", device="pdf")
    ##let's export that as a Cytoscape-friendly file too
    ## igraph - could try saveXGMML in GeneNetworkBuilder?
    igraph::write_graph(tidyBlastNetworkTrimmed, file = gmlFileName, format="gml")
    ## and now to get that data exported
    clustSeqs <- list()
    clusterOutput <- data.frame(igraph::V(tidyBlastNetworkTrimmed)$name, igraph::V(tidyBlastNetworkTrimmed)$cluster)
    colnames(clusterOutput) <- list("gene_oid", "cluster")
    clusterNum <- sort(unique(clusterOutput$cluster))
    ## in case the defined family starting number pushes us into a new order of magnitude above the cluster numbers
    totClust <- length(clusterNum) + defFamNum
    padNum <- nchar(totClust)
    for (i in 1:length(clusterNum)) {
        k <- i + defFamNum
        clustSubset <- clusterOutput[which(clusterOutput$cluster == clusterNum[i]),]
        inClustLength <- length(clustSubset$gene_oid)
        ## note: this should fix hypofam numbering (hypofam_01 not hypofam_1)
        inClustNum <- paste("subfam_",stringr::str_pad(k, padNum, "left", "0"),sep="")
        ## for each protein in the cluster, getting the sequences
        for (j in 1:inClustLength) {
            inClustGene <- clustSubset$gene_oid[j]
            clustSeqs[[as.character(inClustGene)]] <- subSeqs[[inClustGene]]
            ##                   clustSeqs <- append(clustSeqs,hypoSeqs[[inClustGene]])
            if (i == 1 && j == 1) {
                clustListings <- data.table::data.table(stringsAsFactors = FALSE)
                clustListings$gene_oid[1] <- inClustGene
                clustListings$Hypofam[1] <- inClustNum
            } else {
                clustListings <- rbind(clustListings, data.frame(inClustGene, inClustNum), use.names=FALSE)
            }
        }
        dirname <- paste(fileName,"_subClusters_",geneName,sep="")
        if(dir.exists(dirname) == FALSE) {
            dir.create(dirname)
        }
        colnames(clustListings) <- c("gene_oid","subClust")
                                        # write plain .fa files and mafft-aligned files for cluster members
        clusterFile <- paste(dirname,"/",fileName,"_tgCluster_",stringr::str_pad(k, padNum, "left", "0"),".fa",sep="")
        seqinr::write.fasta(clustSeqs, names=names(clustSeqs),file.out=clusterFile)
        ## let's align members of a given cluster wheeee
        if (alnClust == TRUE) {
            mafftInput <- clusterFile
            clusterAlnFile <- paste(dirname,"/",fileName,"_tgCluster_",stringr::str_pad(k, padNum, "left", "0"),"_mafft.fa",sep="")
            mafftOutput <- clusterAlnFile
            ## see previous note re: running with --quiet
            if (sysTerm == "wsl") {
                clustSeqsMafftOut <- system2(command = "wsl", 
                                             args = c("mafft", 
                                                      "--quiet",
                                                      "--auto",
                                                      mafftInput,
                                                      ">",
                                                      mafftOutput),
                                             wait = TRUE,
                                             stdout = TRUE)
            } else { 
                clustSeqsMafftOut <- system2(command = "mafft", 
                                             args = c("--auto",
                                                      "--quiet",
                                                      mafftInput,
                                                      ">",
                                                      mafftOutput),
                                             wait = TRUE,
                                             stdout = TRUE)
                ## stdout change here too?
                                        #        system2(command = "mafft", 
                                        #          args = c("--auto",
                                        #            "--quiet",
                                        #            mafftInput,
                                        #            ">",
                                        #            mafftOutput),
                                        #          stdout = mafftOutput)
            }
            ## NOTE:  add HMM making here, if desired???
            ## or maybe within the sysTerm if/else bit.
            ## if hmmClust == TRUE
            ## get 5% of the seqs or 5 seqs, whichever is more
            ## hmmbuild a model
            ## hmmalign all seqs against it
        }
        inClustLength <- 0
        clustSeqs <- list()
    }
                                        # write out the cluster info separately, just in case
    write.table(clustListings, file=clustListFile, row.names=FALSE, sep="\t")
    ## add the subgroup cluster IDs to the neighbor metadata table, under the column name Hypofam.
    ##
    if (lightExport == TRUE) {
        exportData <- imgNeighborsTrimmed %>% dplyr::select("gene_oid", "Hypofam")
    } else {     
        exportData <- imgNeighborsTrimmed
    }
    for (i in 1:length(exportData$gene_oid)) {
        findSub <- grepl(exportData$gene_oid[i], clustListings$gene_oid)
        if(any(TRUE %in% findSub) == TRUE) {
            clustIdx <- which(findSub == TRUE)
            if (exportData$Hypofam[i] == "") {
                exportData$Hypofam[i] <- as.character(clustListings$subClust[clustIdx])
            } else {
                exportData$Hypofam[i] <- paste(exportData$Hypofam[i]," ",as.character(clustListings$subClust[clustIdx]),sep="")
            }
        }
    }
    ## export the data sorted for easier integration?
    write.table(exportData[order(exportData$gene_oid),], file=outFileName, row.names=FALSE, col.names = TRUE, sep="\t", quote=FALSE)
    print("Subgroups of proteins have been clustered, and members of clusters have been aligned and annotated.")
    return(imgNeighborsTrimmed)
}
