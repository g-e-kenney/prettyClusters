#' Subfunction to identify subgroups of hypothetical proteins in larger IMG datasets
#'
#' Given a desired genome neighborhood and lists of genes of interest and their neighbors, along with their protein sequences, this program uses all-by-all-blast & a clustering method of choice to identify subsets of hypothetical proteins that may represent meaningfully associated proteins.
#' @param imgGenesData Data frame from neighborPrep with metadata for your genes of interest.
#' @param imgNeighborsData Data frame from neighborPrep with metadata for neighbors of your genes of interest.
#' @param imgNeighborSeqs Data frame from neighborPrep with sequence information for neighbors of your gene of interest.
#' @param geneName Character string with your gene of interest's name
#' @param clustMethod Character string with the name of the desired clustering tool
#' @param pidCutoff The percent ID cutoff below which edges are not formed when using tidygraph
#' @param alphaVal The alpha value cutoff used for pvclust
#' @param bootStrap The bootstrap number used for pvclust
#' @param sysTerm The type of terminal (wsl vs linux/unix/macos) from which blastp and other commands will be run
#' @param numThreads The number of processor threads to be devoted to certain steps
#' @param screenPep Should we use peptide-friendly defaults in this run? 
#' @param alnClust Should we make MAFFT alignments of all clusters?
#' @param hmmClust Should we make HMM models of all clusters?
#' @return Updated metadata for neighboring genes (additional files generated en route)
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @importFrom utils capture.output combn read.csv write.csv write.table
#' @importFrom stats dist end hclust na.omit start
#' @importFrom grDevices cairo_pdf cairo_ps col2rgb colorRampPalette dev.off hsv pdf rgb2hsv
#' @examples 
#' \dontrun{
#' imgNeighborsData <- neighborHypothetical(imgGenesData = imgGenesData, 
#'                                          imgNeighborsData = imgNeighborsData, 
#'                                          imgNeighborSeqs=imgNeighborSeqs, 
#'                                          geneName = "genE", 
#'                                          clustMethod = "tidygraph", 
#'                                          pidCutoff = 35, 
#'                                          sysTerm = "nix", 
#'                                          numThreads = 5) 
#' }
#' 
neighborHypothetical <- function(imgGenesData = imgGenesData,
                                 imgNeighborsData = imgNeighborsData,
                                 imgNeighborSeqs=imgNeighborSeqs,
                                 geneName = geneName,
                                 clustMethod = clustMethod,
                                 pidCutoff = pidCutoff,
                                 alphaVal = alphaVal,
                                 bootStrap = bootStrap,
                                 sysTerm = sysTerm,
                                 numThreads = numThreads,
                                 screenPep = screenPep,
                                 alnClust = alnClust,
                                 hmmClust = hmmClust) { 
                                        # first step: flag hypothetical proteins  in the neighbordata file
    ## have some variables and stuff
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    if (screenPep == TRUE) {
        fileName <- paste(fileDate,"_neighborHypothetical_",geneName,sep="")
        hypoSeqsFile <- paste(fileName, "_pepSeqs.fa",sep="")
        blastFile <- paste(fileName, "_pepBlast.txt",sep="")
        blastError <- paste(fileName, "_pepBlastError.txt",sep="")
        pvClustFile <- paste(fileName,"_pepPvClust.txt",sep="")
        pvClustTree <- paste(fileName,"_pepPvClustTree.pdf",sep="")
        clustListFile <- paste(fileName, "_pepClusterList.txt",sep="")
        networkFileName <- paste(fileName, "_pepNetworkFull_pid_",pidCutoff,".pdf",sep="")
        clustNetworkFileName <- paste(fileName, "_pepNetworkClusters_pid_",pidCutoff,".pdf",sep="")
    } else {
        fileName <- paste(fileDate,"_neighborHypothetical_",geneName,sep="")
        hypoSeqsFile <- paste(fileName, "_hypoSeqs.fa",sep="")
        blastFile <- paste(fileName, "_hypoBlast.txt",sep="")
        blastError <- paste(fileName, "_hypoBlastError.txt",sep="")
        pvClustFile <- paste(fileName,"_hypoPvClust.txt",sep="")
        pvClustTree <- paste(fileName,"_hypoPvClustTree.pdf",sep="")
        clustListFile <- paste(fileName, "_hypoClusterList.txt",sep="")
        networkFileName <- paste(fileName, "_hypoNetworkFull_pid_",pidCutoff,".pdf",sep="")
        clustNetworkFileName <- paste(fileName, "_hypoNetworkClusters_pid_",pidCutoff,".pdf",sep="")
    }
    nLength <- length(imgNeighborsData$gene_oid)
    hypoIndex <- list()
    counter <- 1
    hSeqTable <- data.frame()
    hypoSeqs <- data.frame()
    hsReList <- list()
    tidyMonoBlast <- data.frame()
    clustListings <- data.table::data.table(stringsAsFactors = FALSE)
                                        # making sure we have the external tools for this
    ## note:  on windows, the easiest way to Deal is to install blast and mafft via Windows Subsystem for Linux
    ## FIIK how to deal with them as full windows installs.
    ## let's assume users went with the wsl-based installation options as suggested
    print("Analyzing hypothetical proteins.")
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
                                        # looking for things with no pfam or tigrfam IDs
    ## this is flagging probable hypothetical proteins.  
    ## COG, EC, etc. are all too vague/sketchy to count as helpful annotation for these purposes
    ## gene products and symbols are worse
    ## and interpro is sadly not yet a default in IMG data exports
    imgNeighborsData$Pfam <- as.character(imgNeighborsData$Pfam)
    imgNeighborsData$Tigrfam <- as.character(imgNeighborsData$Tigrfam)
    ## intentionally EXCLUDING InterPro listings here, since those include some things that are vague pretty vague (superfamilies)
    ## and for those we may still want to run hypothetical protein analysis.
    ## but there's no good way to sort vague from specific interpro families by number (?)
    ## IF we are doing peptide screening, we are assuming that annotations might be extra-crappy
    if (screenPep == TRUE) {
        hypoTemp <- imgNeighborsData %>% dplyr::filter(.data$Amino.Acid.Sequence.Length..aa. <= 150)
        hypoIndex <- hypoTemp$gene_oid
        rm(hypoTemp)
    } else if (screenPep == FALSE) {     
        ## we only want metadata for hypothetical proteins only
        for (i in 1:nLength) {
            if(imgNeighborsData$Pfam[i] == "" && imgNeighborsData$Tigrfam[i] == "") {
                hypoIndex[[counter]] <-  imgNeighborsData$gene_oid[i]    
                counter <- counter + 1    
            } else {
                hypoIndex <- hypoIndex
            }
        }
    }
    hypoNeighbors <- imgNeighborsData[imgNeighborsData$gene_oid %in% hypoIndex, ]
    ## make a data frame with amino acid sequences for all the hypothetical proteins
    hSeqNames <- names(imgNeighborSeqs)
    for (i in 1:length(hSeqNames)) {
        hSeqTable[i,1] <- hSeqNames[i]
        hSeqTable[i,2] <- imgNeighborSeqs[[i]]
    }
    colnames(hSeqTable) <- c("gene_oid", "sequence")
    ## sometimes neighbors are rRNA or tRNA or whatever
    ## and this is what IMG will have as a "sequence" for them
    ## NOTE: gb2img & incorpIprScan should fake this adequately, but might need a quick visual check...
    badSeq <- "No sequence found"
    ## NOTE: as of may 2022 IMG seems to have been tweaking what happens with RNA; this is adjusted to suit
    badSeq2 <- "not found"
    badSeq3 <- ""
    goodSeqs <- hSeqTable %>% dplyr::filter(.data$sequence != badSeq)
    goodSeqs <- goodSeqs %>% dplyr::filter(.data$sequence != badSeq2)
    goodSeqs <- goodSeqs %>% dplyr::filter(.data$sequence != badSeq3)
    hypoSeqs <- goodSeqs %>% dplyr::filter(.data$gene_oid %in% hypoIndex)
    hsReList <- list()
    ## export the hypothetical sequences separately after re-fasta-ifying
    for (i in 1:length(hypoSeqs[,1])) {
        hsReList[[hypoSeqs[i,1]]] <- hypoSeqs[i,2] 
    }
    hypoSeqs <- hsReList
    seqinr::write.fasta(hsReList,names=names(hsReList), file.out=hypoSeqsFile)
                                        # preparing for all-by-all blast to ID similar sets of proteins
    if(dir.exists("temporary") == FALSE) {
        dir.create("temporary")
    }
    if(dir.exists("temporary/dbblast") == FALSE) {
        dir.create("temporary/dbblast")
    }
    dbtype <- "prot"
    blast_db <- "hypoSeqsDb"
    input <- hypoSeqsFile
    ## easier to cut latter than rerun
    ## note that for peptide screening we want to be relatively lax (particularly for RiPPs) - %ID will be more useful than e-value
    if (screenPep == TRUE) {
        evalue <- "10"
    } else {
        evalue <- "1"
    }
    format <- "6"
    output <- "hypoSeqsDb"
    query <- hypoSeqsFile
    blastColNames <- c("qseqid", "sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
    if (screenPep == TRUE) {
        print("Beginning all-by-all blast of peptides under 150 aa.")
    } else {
        print("Beginning all-by-all blast of hypothetical proteins.")
    }
    ## when running under a wsl or unix/linux environment
    if (sysTerm == "wsl") {
        system2(command = "wsl",
                args = c(makeblastdb,
                         "-dbtype", "prot",
                         "-in", input,
                         "-out",  "temporary/dbblast/hypoSeqsDb"),
                stdout=FALSE)
        ## times out while waiting for stdout sometimes on very large datasets - let's see if stderr works better?
        ## ...it seems to.
        blast_err <- system2(command = "wsl", 
                             args = c(blastp, 
                                      "-db", "temporary/dbblast/hypoSeqsDb",
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
                         "-out", "temporary/dbblast/hypoSeqsDb"),
                stdout=FALSE)
        ## possibly implement stderr variant here too?  need to repeat more and see whether this is a more general issue.
        ## or see if the stdout file switch works better on wsl too...?
        system2(command = blastp, 
                args = c("-db", "temporary/dbblast/hypoSeqsDb",
                         "-query", query,
                         "-outfmt", "6",
                         "-evalue", evalue,
                         "-num_threads", numThreads,
                         "-out", blastFile),
                stdout = blastFile)
                                        #   capture.output(blast_err, blastError)
        blast_out <- as.data.frame(read.csv(blastFile, header=FALSE, sep="\t" , stringsAsFactors=FALSE))
        colnames(blast_out) <- blastColNames
        tidyBlast <- tibble::as_tibble(blast_out)
        ##   write.table(blast_out, blastFile, sep="\t", quote=FALSE)
        ##   tidyBlast <- blast_out %>% tibble::as_tibble() %>% tidyr::separate(col = value, into = blastColNames, sep = "\t", convert = TRUE)
    }
    ## moving the pairwise results into a matrix form
    ## and removing things other than evalue and gene_oid interactions
    print("All-by-all blast of proteins complete.")
    rm(blast_out)
    tidyBlast <- as.data.frame(tidyBlast)
    ## for huge datasets, gotta remove the ginormous objects
    tidyBlast$qseqid <- as.character(tidyBlast$qseqid)
    tidyBlast$sseqid <- as.character(tidyBlast$sseqid)
    tidyBlast$length <- as.numeric(tidyBlast$length)
    tidyBlast$pident <- as.numeric(tidyBlast$pident)
    tidyBlast$bitscore <- as.numeric(tidyBlast$bitscore)
    ##  dealing with NA issues that're sometimes gonna happen - population should generally be small enoguh to be deletable
    tidyBlast <- na.omit(tidyBlast)
    ## making the dataset skinnier - only IDs, percentID, length, and bitscore are maintained later
    tidyBlast <- tidyBlast[,-c(5:11)]
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
    ## anyway arbitrarily using 2/3s the seq length as a rule-out point for proteins
    ## and pegging it to the larger protein to prevent 30 aa peptides from nevertheless matching 3000 aa NRPS monsters
    ## caveat: for peptides we'll be less stringent, since alternate start sites or divergent termini hurt unduly
    yetTidierBlast$qseqLength <- as.numeric(nchar(hypoSeqs[as.character(yetTidierBlast$qseqid)]))
    yetTidierBlast$sseqLength <- as.numeric(nchar(hypoSeqs[as.character(yetTidierBlast$sseqid)]))
    ## need to implement this better - is it fast enough with mutate?
    yetTidierBlast <- yetTidierBlast %>% dplyr::rowwise() %>% dplyr::mutate(maxSeqLength = max(.data$qseqLength, .data$sseqLength))
                                        # for (i in 1:length(yetTidierBlast$qseqLength)) {
                                        #   yetTidierBlast$maxSeqLength <- max(c(yetTidierBlast$qseqLength[i], yetTidierBlast$sseqLength[i]))
                                        # }
    if (screenPep == TRUE) {
        tidyMonoBlast <- yetTidierBlast %>% dplyr::filter(.data$length>=.45*.data$maxSeqLength)   
    } else {
        tidyMonoBlast <- yetTidierBlast %>% dplyr::filter(.data$length>=.65*.data$maxSeqLength)   
    }
    ## keeping only names and percent id now that we are done with length and bitscore, then reshaping it
    rm(yetTidierBlast)
    tidyMinBlast <- as.data.frame(tidyMonoBlast[,c(1:3)])
    rm(tidyMonoBlast)
    ## note that if using eval I'd do a log transform here to avoid crazy results
    ## not necessary for %ID though
    if (clustMethod == "pvClust") {
        ## to consider - make the cluster and distance methods easy for the user to alter?
        ## bootstrap > 10 would be better but it does not scale linearly with sequence number
        tidyBlastMatrix <- tidyMinBlast %>% tidyr::spread(.data$sseqid, .data$pident)
                                        # alternate phrasing for some errors i have occasionally hit?  keeping these here in case they are widespread
        ## tidyBlastMatrix <- tidyMinBlast %>% group_by_at(vars(-pident)) %>% mutate(row_id=1:n()) %>% ungroup() %>% spread(key=sseqid, value=pident) %>% select(-row_id)
        ##  tidyBlastMatrix <- tidyMinBlast %>% tidyr::pivot_wider(names_from = sseqid, values_from = pident)
        tbRnames <- tidyBlastMatrix[,1]
        matLength <- length(tidyBlastMatrix[1,])
                                        # using the max eval (or 0% ID) to fill in NAs (things that didn't even have matches at that value)
        tidyBlastMatrix[is.na(tidyBlastMatrix)] <- 0
        tidyBlastMatrixL2 <- as.matrix(tidyBlastMatrix[,2:matLength])
        rownames(tidyBlastMatrixL2) <- tbRnames
        pvBlast <- pvclust::pvclust(tidyBlastMatrixL2, method.dist="euclidean", method.hclust="ward.D2",nboot=bootStrap, parallel=TRUE)
                                        # output info as a .txt file and output the dendrogram with the appropriate cutoff
        ## given the slow speed of the initial pvclust, consider making the alphaVal finalization interactive
        ##if(interactive()) {#
        ##  plot(pvBlast)
        ##  pvclust::pvrect(pvBlast, alphaVal, pv="au",type="geq")
        ##  aValOK <- readline(prompt="Accept alpha value cutoff for cluster (y/n): ")
        ##  while (aValOK !="y") {
        ##      alphaVal <- readline(prompt="New alpha value cutoff (.00-1.00):  ")
        ##    plot(pvBlast)
        ##      pvclust::pvrect(pvBlast, alphaVal, pv="au",type="geq")
        ##      aValOK <- readline(prompt="Accept alpha value cutoff for cluster (y/n): ")
        ##  }  
        ##}
        pdf(file=pvClustTree)
        plot(pvBlast)
        pvclust::pvrect(pvBlast, alphaVal, pv="au",type="geq")
        dev.off()
        print("hypothetical proteins clustered via pvclust")
        pvBlastPick <- pvclust::pvpick(pvBlast, alphaVal)
        clusterNum <- length(pvBlastPick$clusters)
        ## post-clustering, generates aligned and unaligned .fa files and a gene_oid/cluster concordance
        ## adding a floor so that if a gene is present in <1% of neighborhoods we don't report the cluster
        ## (2-gene clusters and so on are... very rarely gonna be useful)
        ## once clusters have been picked, sets up .fa files for the sequences, as-is and aligned
        ## it also creates a table matching gene_oids and clusters for input to the later metadata file
        atLeastGenes <- floor(.005*length(tbRnames))
        clustSeqs <- list()
        for (i in 1:clusterNum) {
            inClustLength <- length(pvBlastPick$clusters[[i]])
            if (inClustLength >= atLeastGenes) {

                inClustNum <- paste("hypofam_",i,sep="")
                                        # for each gene in the cluster, getting the sequences
                for (j in 1:inClustLength) {
                    inClustGene <- pvBlastPick$clusters[[i]][j]
                    clustSeqs[[as.character(inClustGene)]] <- hypoSeqs[[inClustGene]]
                    ##      clustSeqs <- append(clustSeqs,hypoSeqs[[inClustGene]])
                    if (i == 1 && j == 1) {
                        clustListings <- data.table::data.table(stringsAsFactors = FALSE)
                        clustListings$gene_oid[1] <- inClustGene
                        clustListings$Hypofam[1] <- inClustNum
                    } else {
                        clustListings <- rbind(clustListings, data.frame(inClustGene, inClustNum), use.names=FALSE)
                    }
                }
                dirname <- paste(fileName,"_hypoClusters_",geneName,sep="")
                if(dir.exists(dirname) == FALSE) {
                    dir.create(dirname)
                }
                colnames(clustListings) <- c("gene_oid","hypoClust")
                                        # write plain .fa files and mafft-aligned files for cluster members
                clusterFile <- paste(dirname,"/",fileName,"_pvCluster_",i,".fa",sep="")
                seqinr::write.fasta(clustSeqs, names=names(clustSeqs),file.out=clusterFile)
                                        # let's align members of a given cluster wheeee
                ## this is valuable since it lets you evaluate pretty quickly whether or not the groups are sane
                ## even looking at the filesize - did a large increase in gaps 
                mafftInput <- clusterFile
                clusterAlnFile <- paste(dirname,"/",fileName,"_pvCluster_",i,"_mafft.fa",sep="")
                mafftOutput <- clusterAlnFile
                ## running in --quiet because mafft has a bunch of non-deadly errors (e.g. one shows up if mafft run was run after su)
                if (sysTerm == "wsl") {
                    clustSeqsMafftOut <- system2(command = "wsl", 
                                                 args = c(mafft, 
                                                          "--quiet",
                                                          "--auto",
                                                          mafftInput,
                                                          ">",
                                                          mafftOutput),
                                                 wait = TRUE,
                                                 stdout = TRUE)
                } else { 
                    clustSeqsMafftOut <- system2(command = mafft, 
                                                 args = c("--auto",
                                                          "--quiet",       
                                                          mafftInput,
                                                          ">",
                                                          mafftOutput),
                                                 wait = TRUE,
                                                 stdout = TRUE)
                }
            } else {
                if (i == 1) {
                    clustListings <- data.table::data.table(stringsAsFactors = FALSE)
                } else {
                    clustListings <- clustListings
                }
            }
            inClustLength <- 0
            clustSeqs <- list()
        }
    } else if (clustMethod == "tidygraph") {
        ## the tidygraph option
        ## make the object using the percent ID cutoff, and do so as an undirected network (deleting self-linking edges)
        tidyMinBlastTrimmed <- tidyMinBlast %>% dplyr::filter(.data$pident>=pidCutoff)
        tidyBlastNetwork <- tidygraph::as_tbl_graph(tidyMinBlastTrimmed, directed=FALSE)
        tidyBlastNetworkTrimmed <- tidyBlastNetwork %>% tidygraph::activate(edges) %>% dplyr::filter(.data$from != .data$to)
        ## highlighting the gene of interest, in case it is a hypothetical protein
        ## not useful if it isn't, but if it is, it provides an easy way to track whether cluster-making is crap
        tidyBlastNetworkTrimmed <- tidyBlastNetworkTrimmed %>% tidygraph::activate(nodes) %>% dplyr::mutate(isGene = ifelse(.data$name %in% imgGenesData$gene_oid, "yes", "no"))
        if (screenPep == TRUE) {
            titleText <- paste("Cluster-based identification of peptide subgroups in neighbors of ",geneName, ".",sep="")
        } else {
            titleText <- paste("Cluster-based identification of hypothetical protein subgroups in neighbors of ",geneName, ".",sep="")
        }
        subtitleText <- paste(" %ID cutoff of ", pidCutoff, " and a 65% match length used for edges. Analyzed on ", fileDate,".",sep="")
        ## pic of the full-sized cluster
        tidyBlastNetworkPic <- ggraph::ggraph(tidyBlastNetworkTrimmed, layout="stress") + 
            ggraph::geom_edge_link0(ggplot2::aes(edge_alpha=.data$pident/100), edge_colour="black", edge_width=.5) + 
            ggraph::geom_node_point(ggplot2::aes(color=.data$isGene), size=3) + 
            ggplot2::labs(title=titleText, subtitle=subtitleText) + 
            ggplot2::theme(plot.title=ggplot2::element_text(color="black", size=18, margin=ggplot2::margin(10,0,10,0)),
                           plot.subtitle=ggplot2::element_text(color="grey", size=12, margin=ggplot2::margin(10,0,10,0))) +
            ggraph::theme_graph(base_family="sans") + 
            ggplot2::scale_color_manual(values=c("yes"="firebrick","no"="steelblue")) + 
            ggplot2::theme(legend.position="none", plot.margin=ggplot2::unit(c(.2,.2,.2,.2), "cm"))  
        ggplot2::ggsave(filename=networkFileName, tidyBlastNetworkPic, height=10, width=20, dpi=75, units="in", device="pdf")
        ## now getting the individual clusters("components") and keeping the ones above our atLastGenes cutoff
        tidyBlastNetworkTrimmed <- tidyBlastNetworkTrimmed %>% tidygraph::activate(nodes) %>% dplyr::mutate(cluster=tidygraph::group_components())
        atLeastGenes <- floor(.005*length(unique(tidyMinBlast$qseqid)))
        keepCluster <- which(table(igraph::V(tidyBlastNetworkTrimmed)$cluster) >= atLeastGenes)
        tidyBlastClusters <- tidyBlastNetworkTrimmed %>% tidygraph::activate(nodes) %>% dplyr::filter(.data$cluster %in% keepCluster)
        ## making a graph for just the larger clusters
        tidyBlastClusterPic <- ggraph::ggraph(tidyBlastClusters, layout="stress") + 
            ggraph::geom_edge_link0(ggplot2::aes(edge_alpha=.data$pident/100), edge_colour="black", edge_width=.5) + 
            ggraph::geom_node_point(ggplot2::aes(color=.data$isGene), size=3) + 
            ggplot2::labs(title=titleText, subtitle=subtitleText) + 
            ggplot2::theme(plot.title=ggplot2::element_text(color="black", size=18), plot.subtitle=ggplot2::element_text(color="grey", size=12)) +
            ggraph::theme_graph(base_family="sans") + 
            ggplot2::scale_color_manual(values=c("yes"="firebrick","no"="steelblue")) + 
            ggplot2::theme(legend.position="bottom", plot.margin=ggplot2::unit(c(.2,.2,.2,.2), "cm")) 
        ggplot2::ggsave(filename=clustNetworkFileName, tidyBlastClusterPic, height=10, width=20, dpi=75, units="in", device="pdf")
        ## and now to get that data exported
        ## once we get it out of some awkward list formats
        clustSeqs <- list()
        clusterOutput <- data.frame(igraph::V(tidyBlastClusters)$name, igraph::V(tidyBlastClusters)$cluster)
        colnames(clusterOutput) <- list("gene_oid", "cluster")
        clusterNum <- sort(unique(clusterOutput$cluster))
        totClust <- length(clusterNum)
        padNum <- nchar(totClust)
        for (i in 1:length(clusterNum)) {
            clustSubset <- clusterOutput[which(clusterOutput$cluster == clusterNum[i]),]
            inClustLength <- length(clustSubset$gene_oid)
            if (screenPep == TRUE) {
                inClustNum <- paste("pepfam_",stringr::str_pad(i, padNum, "left", "0"),sep="")
            } else {
                inClustNum <- paste("hypofam_",stringr::str_pad(i, padNum, "left", "0"),sep="")
            }
                                        # for each gene in the cluster, getting the sequences
            for (j in 1:inClustLength) {
                inClustGene <- clustSubset$gene_oid[j]
                clustSeqs[[as.character(inClustGene)]] <- hypoSeqs[[inClustGene]]
                ##                   clustSeqs <- append(clustSeqs,hypoSeqs[[inClustGene]])
                if (i == 1 && j == 1) {
                    clustListings <- data.table::data.table(stringsAsFactors = FALSE)
                    clustListings$gene_oid[1] <- inClustGene
                    clustListings$Hypofam[1] <- inClustNum
                } else {
                    clustListings <- rbind(clustListings, data.frame(inClustGene, inClustNum), use.names=FALSE)
                }
            }
            if (screenPep == TRUE) {
                dirname <- paste(fileName,"_pepClusters_",geneName,sep="")
                if(dir.exists(dirname) == FALSE) {
                    dir.create(dirname)
                }
                colnames(clustListings) <- c("gene_oid","pepClust")
            } else {
                dirname <- paste(fileName,"_hypoClusters_",geneName,sep="")
                if(dir.exists(dirname) == FALSE) {
                    dir.create(dirname)
                }
                colnames(clustListings) <- c("gene_oid","hypoClust")
            }
            ## writing plain .fa files and (in a sec) mafft-aligned files for cluster members
            clusterFile <- paste(dirname,"/",fileName,"_tgCluster_",stringr::str_pad(i, padNum, "left", "0"),".fa",sep="")
            seqinr::write.fasta(clustSeqs, names=names(clustSeqs),file.out=clusterFile)
            if (alnClust == TRUE) {
                ## let's align members of a given cluster wheeee
                mafftInput <- clusterFile
                clusterAlnFile <- paste(dirname,"/",fileName,"_tgCluster_",stringr::str_pad(i, padNum, "left", "0"),"_mafft.fa",sep="")
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
                    ## stdout change here too?  possibly not necessary.
                                        #        system2(command = "mafft", 
                                        #          args = c("--auto",
                                        #            "--quiet",
                                        #            mafftInput,
                                        #            ">",
                                        #            mafftOutput),
                                        #          stdout = mafftOutput)
                }
                ## NOTE:  add HMM making here, if desired???
                ## if hmmClust == TRUE
                ## get 5% of the seqs or 5 seqs, whichever is more
                ## hmmbuild a model
                ## hmmalign all seqs against it
            }
            inClustLength <- 0
            clustSeqs <- list()
        }
    }
    if (alnClust == TRUE) {
        print("Proteins have been clustered and aligned.")    
    } else {
        print("Proteins have been clustered.")    
    }  
    ## write out the cluster info separately, just in case
    write.table(clustListings, file=clustListFile, row.names=FALSE, sep="\t", quote=FALSE)
    ## add the hypo cluster IDs to the neighbor metadata table, under the column name Hypofam.
    colnames(clustListings) <- c("gene_oid", "subClust")
    for (i in 1:length(imgNeighborsData$gene_oid)) {
        findHypo <- which(clustListings$gene_oid == imgNeighborsData$gene_oid[i])
        if(any(findHypo) == TRUE) {
#        findHypo <- grepl(imgNeighborsData$gene_oid[i], clustListings$gene_oid)
#        if(any(TRUE %in% findHypo) == TRUE) {
#            clustIdx <- which(findHypo == TRUE)
            if (imgNeighborsData$Hypofam[i] == "") {
                ## if we had no annotations, we can just write this one in
                imgNeighborsData$Hypofam[i] <- as.character(clustListings$subClust[findHypo])
            } else {
                ## since this could be a peptide run OR postdate some other annotation
                ## let's make sure that we don't overwrite any existing annotations
                imgNeighborsData$Hypofam[i] <- paste(imgNeighborsData$Hypofam[i],as.character(clustListings$subClust[findHypo]),sep=" ")
            }
        }
    }
    if (screenPep == TRUE) {
        print("Members of peptide subgroups have been annotated.")    
    } else {
        print("Hypothetical protein subgroups have been annotated.")    
    }  
    return(imgNeighborsData)
    ## everything else is output as .txt files i guess
}
