#' Subfunction to cluster genomic neighborhoods based on the similarity of their protein family sets
#'
#' This function cluster-based methods to identify subgroups of gene clusters containing genes of interest that have similar genomic neighborhoods.
#' @param imgGenesTrimmed Data frame from analyzeNeighbors containing metadata for genes of interest. Character string, required.
#' @param imgNeighborsTrimmed Data frame from analyzeNeighbors containing metadata for neighbors of genes of interest. Character string, required.
#' @param geneName Name of gene of interest as string.  Required.
#' @param neighborMatrixData Matrix from analyzeNeighbors/neighborMatrix containing binary data re: protein family presence in gene clusters. Required.
#' @param autoClust Should clusters be automatically assigned? Boolean, defaults to TRUE in parent function.
#' @param clustMethod String specifying what tool (tidygraph, pvclust) will be used to identify clusters. Defaults to "tidygraph" in parent function.
#' @param bootStrap Bootstrap rounds for pvclust. Integer, defined in parent function.
#' @param alphaVal Alpha value cutoff for pvclust. Number, defined in parent function.
#' @param coreGeneName Name of gene of interest as string, without any suffix. String, defined in parent function.
#' @param tgCutoff Number (0-1) representing edge similarity to keep for tidygraph. Default 0.6.
#' @return List with matrix object and updated metadata, with misc. files output along the way
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @importFrom utils capture.output combn read.csv write.csv write.table
#' @importFrom stats dist end hclust na.omit start
#' @importFrom grDevices cairo_pdf cairo_ps col2rgb colorRampPalette dev.off hsv pdf rgb2hsv
#' @examples
#' \dontrun{
#' neighborClustersOut <- neighborClusters(imgGenesTrimmed = imgGenesTrimmed, 
#'                                         imgNeighborsTrimmed = imgNeighborsTrimmed, 
#'                                         geneName = "genE", 
#'                                         neighborMatrixData = neighborMatrixData, 
#'                                         autoClust = TRUE, 
#'                                         clustMethod = "tidygraph", 
#'                                         coreGeneName = "genE", 
#'                                         tgCutoff=.6)
#' }
#'
neighborClusters <- function(imgGenesTrimmed = imgGenesTrimmed,
                             imgNeighborsTrimmed = imgNeighborsTrimmed,
                             geneName = geneName,
                             neighborMatrixData = neighborMatrixData,
                             autoClust = autoClust,
                             clustMethod = clustMethod,
                             alphaVal = alphaVal,
                             bootStrap = bootStrap,
                             coreGeneName = coreGeneName,
                             tgCutoff=tgCutoff) { 
    ## make alignments for all family members?
    ## calulate abundance of all families
    ## (input .fasta list for all, extract and mafft-align) 
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    fileName <- paste(fileDate,"_neighborClusters_",geneName,sep="")
    finaltxtname <- paste(fileName,"_clusterList.txt",sep="")
    finalMetaGenes <- paste(fileName,"_geneMetadata.txt",sep="")
    finalMetaNeighbors <- paste(fileName,"_neighborMetadata.txt",sep="")
    finalpdfname <- paste(fileName,"_heatmap.pdf",sep="")
    finalpsname <- paste(fileName,"_heatmap.ps",sep="")
    finalpngname <- paste(fileName,"_heatmap.png",sep="")
    pvMatrixFile <- paste(fileName,"_pvClust.txt",sep="")
    pvMatrixTree <- paste(fileName,"_pvClustTree.pdf",sep="")
    networkFile <- paste(fileName,"_networkFile.txt",sep="")
    networkFigureFile <- paste(fileName,"_networkFigure.pdf",sep="")
    clustColorFile <- paste(fileName,"_clustColors.txt",sep="")
                                        # data import
    geneNames <- neighborMatrixData[,1]
    matrixData <- data.matrix(neighborMatrixData[,2:ncol(neighborMatrixData)])
    rownames(matrixData) <- geneNames
    tMatrixData <- t(matrixData)
    ## clustering based on trait-set similarity and gene-set similarity
    ## if automatic cluster determination is off, everything'll output in order but clusters will not be automatically assigned
    ## this method uses pvclust to auto-identify the groups of bgcs
    ## note that it is not doing anything to the genes or changing the final display or anything
    if (autoClust == FALSE)  {
        rowDistance  <- dist(matrixData, method = "euclidean")
        rowCluster  <- hclust(rowDistance, method = "ward.D2")
        colDistance  <- dist(tMatrixData, method = "euclidean")
        colCluster  <- hclust(colDistance, method = "ward.D2")
                                        # output in cluster order
        initialOrder <- matrixData[rev(rowCluster$labels[rowCluster$order]), colCluster$labels[colCluster$order]]
                                        # export of clustered table and heatmap.
        orderedGenes <- rownames(initialOrder)
        orderedMatrix <- data.table::as.data.table(initialOrder)
        orderedMatrix$clustNum <- ""
        orderedMatrix$gene_oid <- orderedGenes
        ## no clusters defined - the user can do so manually later if they so desire
        imgGenesTrimmed$clustNum <- ""
        imgNeighborsTrimmed$clustNum <- ""
    } else {
        if (clustMethod == "pvClust") {
            pvMatrix<- pvclust::pvclust(tMatrixData, method.dist="manhattan", method.hclust="ward.D2",nboot=bootStrap, parallel=TRUE)
            ## might need to experiment with alpha here to get a good cluster starting point
            pvMatrixPick <- pvclust::pvpick(pvMatrix, alpha=alphaVal)
            utils::capture.output(pvMatrixPick, file=pvMatrixFile)
            cairo_pdf(filename=pvMatrixTree)
            plot(pvMatrix)
            pvclust::pvrect(pvMatrix, alphaVal, pv="au",type="geq")
            dev.off()
            clusterNum <- length(pvMatrixPick$clusters)
            clusterGroups <- data.frame(matrix("", nrow=1,ncol=2),stringsAsFactors=FALSE)
            ## tying gene and cluster together non-annoyingly    
            for (i in 1:clusterNum) {
                inClustNum <- length(pvMatrixPick$clusters[[i]])
                for (j in 1:inClustNum) {
                    inClustGene <- pvMatrixPick$clusters[[i]][j]
                    if (i == 1 && j == 1) {
                        clusterGroups[1,1] <- inClustGene
                        clusterGroups[1,2] <- i
                    } else {
                        clusterGroups <- rbind(clusterGroups, c(inClustGene, i))
                    }
                    next
                }
            }
            colnames(clusterGroups) <- c("clustGene","clustNum")
            ## getting the family clustering done for the figure 
            colDistance  <- dist(tMatrixData, method = "euclidean")
            colCluster  <- hclust(colDistance, method = "ward.D2")
            ## note that here we're sticking with the pvClust order
            rowCluster <- pvMatrix$hclust
            initialOrder <- matrixData[rowCluster$labels[rowCluster$order], rev(colCluster$labels[colCluster$order])]
            orderedGenes <- rownames(initialOrder)
            orderedMatrix <- data.table::as.data.table(initialOrder)
            orderedMatrix$clustNum <- ""
            orderedMatrix$gene_oid <- orderedGenes
            imgGenesTrimmed$clustNum <- ""
            imgNeighborsTrimmed$clustNum <- ""
            ## getting the cluster number tied to the gene in the general metadata doc
            for (k in 1:length(orderedMatrix$gene_oid)) {
                ## not all genes are in pvclust-IDed clusters)
                if (any(grepl(orderedMatrix$gene_oid[k], clusterGroups$clustGene)) == FALSE) {
                    next
                } else {
                    clustGeneLoc <- grep(orderedMatrix$gene_oid[k], clusterGroups$clustGene)
                    orderedMatrix$clustNum[k] <- clusterGroups$clustNum[clustGeneLoc]
                    originLoc <- grep(orderedMatrix$gene_oid[k], imgGenesTrimmed$gene_oid)
                    imgGenesTrimmed$clustNum[originLoc] <- clusterGroups$clustNum[clustGeneLoc]
                    ## this should add the info to the neighbors file for pretty clustering
                    neighborLoc <- grep(orderedMatrix$gene_oid[k], imgNeighborsTrimmed$source_gene_oid)
                    for (j in 1:length(neighborLoc)) {
                        imgNeighborsTrimmed$clustNum[neighborLoc[j]] <- clusterGroups$clustNum[clustGeneLoc]
                    }
                }
            }
        } else if (clustMethod == "tidygraph") {
            ## the tidygraph option
            ## calculating euclidean distance between gene cluster varieties 
            eucDistBGC <- dist(matrixData, method="euclidean", diag=FALSE, upper=TRUE, p=2)
            eucDistPairs <- data.frame(t(combn(rownames(matrixData),2)), as.numeric(eucDistBGC))
            names(eucDistPairs) <- c("node1", "node2", "eucDist")
            ## dealing with the identical-neighborhood issue, do not want 0-weighted edges
            eucDistPairs$eucDist <- eucDistPairs$eucDist + 1
            ## default using something handwavily similar to the gene stuff - let's toss the worst 40% of edges
            eucCutoff <- max(eucDistPairs$eucDist)*tgCutoff           
            eucPairsTrimmed <- eucDistPairs %>% dplyr::filter(.data$eucDist <= eucCutoff)
            if (length(eucPairsTrimmed$node1) < .1*length(eucDistPairs$node1)) {
                networkError <- paste("At your euclidian distance cutoff (tgCutoff), ",tgCutoff,", more than 90% of network edges are removed - try increasing. Proceeding without clustering.",sep="")
                print(networkError)
                rowDistance  <- dist(matrixData, method = "euclidean")
                rowCluster  <- hclust(rowDistance, method = "ward.D2")
                colDistance  <- dist(tMatrixData, method = "euclidean")
                colCluster  <- hclust(colDistance, method = "ward.D2")
                                        # output in cluster order
                initialOrder <- matrixData[rev(rowCluster$labels[rowCluster$order]), colCluster$labels[colCluster$order]]
                                        # export of clustered table and heatmap.
                orderedGenes <- rownames(initialOrder)
                orderedMatrix <- data.table::as.data.table(initialOrder)
                orderedMatrix$clustNum <- ""
                orderedMatrix$gene_oid <- orderedGenes
                ## no clusters defined - the user can do so manually later if they so desire
                imgGenesTrimmed$clustNum <- ""
                imgNeighborsTrimmed$clustNum <- ""
                autoClust <- FALSE
            } else  {
                write.table(eucPairsTrimmed, networkFile, row.names=FALSE,sep="\t", quote=FALSE) 
                ## zomg network
                eucNetwork <- tidygraph::as_tbl_graph(eucPairsTrimmed, directed=FALSE)    
                ## calculating some stuff - dunno if i'll keep all of these in the end
                eucNetwork <- eucNetwork %>% tidygraph::activate(nodes) %>% dplyr::mutate(group = tidygraph::group_infomap())
                eucNetwork <- eucNetwork %>% tidygraph::activate(edges) %>% dplyr::mutate(betweenness = tidygraph::centrality_edge_betweenness())
                eucNetwork <- eucNetwork %>% tidygraph::activate(nodes) %>% dplyr::mutate(center_dist = tidygraph::node_distance_to(tidygraph::node_is_center()))
                nodeNum <- length(unique(eucDistPairs$node1))
                kNum <- ceiling(nodeNum/20)
                eucNetwork <- eucNetwork %>% tidygraph::activate(nodes) %>% dplyr::mutate(center = tidygraph::node_is_center(), keyplayer = tidygraph::node_is_keyplayer(k = kNum, maxsec = 1200))
                titleText <- paste("Identifying genome neighborhood clusters for ", coreGeneName," via tidygraph",sep="")
                subtitleText <- paste("Edge values represent euclidean distance between neighborhoods, with a cutoff of ",
                                      tgCutoff,
                                      "% of the maximum distance. Unconnected nodes (neighborhoods) not shown. Analyzed ",
                                      fileDate,
                                      ".",
                                      sep="")
                eucNetworkPic <- ggraph::ggraph(eucNetwork , layout="stress") + 
                    ggraph::geom_edge_link0(edge_colour="grey66", edge_width=.25, show.legend=FALSE) + 
                    ggraph::geom_node_point(ggplot2::aes(color=factor(.data$group)), size=7) + 
                    ggplot2::labs(title=titleText, subtitle=subtitleText) + 
                    ggplot2::theme(plot.title=ggplot2::element_text(color="black",
                                                                    size=12,
                                                                    margin=ggplot2::margin(10,0,10,0)),
                                   plot.subtitle=ggplot2::element_text(color="grey66",
                                                                       size=10,
                                                                       margin=ggplot2::margin(10,0,10,0))) +
                    ggraph::theme_graph(base_family="sans") + 
                    ggplot2::theme(legend.position="bottom", plot.margin=ggplot2::unit(c(.2,.2,.2,.2), "cm"))   + 
                    ggplot2::scale_color_viridis_d() 
                ggplot2::ggsave(filename=networkFigureFile, eucNetworkPic, height=10, width=20, dpi=75, units="in", device="pdf")
                clusterNum <- length(unique(igraph::V(eucNetwork)$group))
                clusterGroups <- data.frame(list("clustGene"=igraph::V(eucNetwork)$name))
                clusterGroups$clustNum <- igraph::V(eucNetwork)$group
                ## we still want a hierarchical clustering object though
                ## this uses the same defaults as we use in the non-pvclust version
                rowDistance  <- dist(matrixData, method = "euclidean")
                rowCluster  <- hclust(rowDistance, method = "ward.D2")
                colDistance  <- dist(tMatrixData, method = "euclidean")
                colCluster  <- hclust(colDistance, method = "ward.D2")
                                        # output in cluster order
                initialOrder <- matrixData[rev(rowCluster$labels[rowCluster$order]), colCluster$labels[colCluster$order]]
                                        # export of clustered table and heatmap.
                orderedGenes <- rownames(initialOrder)
                orderedMatrix <- data.table::as.data.table(initialOrder)
                orderedMatrix$gene_oid <- orderedGenes
                orderedMatrix$clustNum <- ""
                imgGenesTrimmed$clustNum <- ""
                imgNeighborsTrimmed$clustNum <- ""
                ## getting the cluster number tied to the gene in the general metadata doc
                for (k in 1:length(orderedMatrix$gene_oid)) {
                    ## in case the cluster ID method leaves some neighborhoods out in the proverbial cold
                    if (any(grepl(orderedMatrix$gene_oid[k], clusterGroups$clustGene)) == FALSE) {
                        next
                    } else {
                        clustGeneLoc <- grep(orderedMatrix$gene_oid[k], clusterGroups$clustGene)
                        orderedMatrix$clustNum[k] <- clusterGroups$clustNum[clustGeneLoc]
                        originLoc <- grep(orderedMatrix$gene_oid[k], imgGenesTrimmed$gene_oid)
                        imgGenesTrimmed$clustNum[originLoc] <- clusterGroups$clustNum[clustGeneLoc]
                        ## this should add the info to the neighbors file for pretty clustering
                        neighborLoc <- grep(orderedMatrix$gene_oid[k], imgNeighborsTrimmed$source_gene_oid)
                        for (j in 1:length(neighborLoc)) {
                            imgNeighborsTrimmed$clustNum[neighborLoc[j]] <- clusterGroups$clustNum[clustGeneLoc]
                        }
                    }
                }
            }
        }
    }        
    ## we also want to add the order of the cluster output just in case
    rowidx <- order(orderedMatrix$gene_oid)
    imgGenesTrimmed$clustOrd <- ""
    imgNeighborsTrimmed$clustOrd <- ""
    for (i in 1:length(rowidx)) {
        geneIdx <- grep(orderedMatrix$gene_oid[i], imgGenesTrimmed$gene_oid)
        imgGenesTrimmed$clustOrd[geneIdx] <- i
        ## this should add the order info to the neighbors file for pretty clustering
        neighborIdx  <- grep(orderedMatrix$gene_oid[i], imgNeighborsTrimmed$source_gene_oid)
        for (j in 1:length(neighborIdx)) {
            imgNeighborsTrimmed$clustOrd[neighborIdx[j]] <- i
        }
    }
    ## adding this  facilitates merging back into the Cytoscape network
    imgGenesTrimmed$name <- imgGenesTrimmed$efi_oid
    write.table(imgGenesTrimmed, finalMetaGenes, row.names=FALSE, sep="\t", quote=FALSE) 
    write.table(imgNeighborsTrimmed, finalMetaNeighbors, row.names=FALSE, sep="\t", quote=FALSE) 
    orderedMatrix$clustOrd <- rownames(orderedMatrix)
    write.table(orderedMatrix, finaltxtname, row.names=FALSE, sep="\t", quote=FALSE)
                                        # drawing the heatmap
    heatmapName <- paste("genome neighborhood clustering for ",coreGeneName,sep="")
    ## depending on whether or not pvClust has been used
    if (autoClust == TRUE) {
        heatAnnoRow <- as.data.frame(imgGenesTrimmed$clustNum, stringsAsFactors=FALSE)
        rownames(heatAnnoRow) <- imgGenesTrimmed$gene_oid
        colnames(heatAnnoRow) <- "clustNum"
        clustNumTemp <- unique(heatAnnoRow$clustNum)
        clustNumTemp[clustNumTemp==""]<-"0"
        clustNumTemp[clustNumTemp=="none"]<-"0"
        clustNumTemp <- sort(as.numeric(clustNumTemp))
        ## to deal with genes that did not get assigned a cluster
        if ("0" %in% as.character(clustNumTemp)) {
            clustColorsTemp <- viridis::viridis(length(unique(imgGenesTrimmed$clustNum)) -1)
            clustColorsTemp <- append("#FFFFFFFF", clustColorsTemp)
            clustNumWords <- numbers2words(clustNumTemp[-1])
            clustNumWords <- append("none",clustNumWords)
                                        #      clustConcord <- data.frame(list(clustNum=clustN))
        } else {
            clustColorsTemp <- viridis::viridis(length(unique(imgGenesTrimmed$clustNum)))
            clustNumWords <- numbers2words(sort(clustNumTemp))
        }
        clustColors <- data.frame(list(clustNum=clustNumTemp, clustColor=clustColorsTemp, clustName=clustNumWords))
        write.table(clustColors, clustColorFile, col.names=TRUE, sep="\t", quote=FALSE)    
        ## this all deals with genes that pvClust did  not feel confident in assigning a cluster
        ## those show up as 'none' and have a white cluster color
        ## i.e. it looks like they have none
        heatAnnoRow[,1][heatAnnoRow[,1]==""]<-"0"
        heatRows <- rownames(heatAnnoRow)
        clustColors$clustNum <- as.character(clustColors$clustNum)
        heatAnnoRow <- heatAnnoRow %>% dplyr::left_join(clustColors, by="clustNum")
        heatAnnoRow <- as.data.frame(heatAnnoRow[,3])
        rownames(heatAnnoRow) <- heatRows
        colnames(heatAnnoRow) <- "clustNum"
        clustNumListTemp <- clustColorsTemp
        names(clustNumListTemp) <- clustNumWords
        clustColorList <- list(clustNum=clustNumListTemp)
        ## this adds colors at the left (gene) dendrogram corresponding to the cluster number
        prettyHeatmap <- pheatmap::pheatmap(matrixData,
                                            annotation_row=heatAnnoRow,
                                            cluster_cols = colCluster,
                                            cluster_rows = rowCluster,
                                            color=c("#FFFFFF", "#000000"),
                                            fontsize=8,
                                            fontsize_row=5,
                                            fontsize_col=5,
                                            cellwidth=5,
                                            cellheight=5,
                                            border=FALSE,
                                            legend=FALSE,
                                            annotation_names_row=FALSE,
                                            annotation_colors=clustColorList,
                                            main=heatmapName,
                                            angle_col=45,
                                            filename=finalpdfname)
    } else if (autoClust == FALSE) {
        ## this just has the raw heatmap without cluster info
        prettyHeatmap <- pheatmap::pheatmap(matrixData,
                                            cluster_cols = colCluster,
                                            cluster_rows = rowCluster,
                                            color=c("#FFFFFF", "#000000"),
                                            fontsize=8,
                                            fontsize_row=5,
                                            fontsize_col=5,
                                            cellwidth=5,
                                            cellheight=5,
                                            border=FALSE,
                                            legend=FALSE,
                                            main=heatmapName,
                                            angle_col=45,
                                            filename=finalpdfname)
    }    
    picHeight <- length(matrixData[,1])/10
    picWidth <-  length(matrixData[1,])/10 
                                        #  cairo_pdf(filename=finalpdfname, width=picWidth, height=picHeight)
                                        #    prettyHeatmap
                                        #  dev.off()
                                        #  cairo_ps(filename=finalpsname,  width=picWidth, height=picHeight)
                                        #    prettyHeatmap
                                        #  dev.off()
                                        #  png(file=finalpngname, width=picWidth, height=picHeight, units="in")
                                        #    prettyHeatmap
                                        #  dev.off()
    print("Genes of interest clustered by genome neighborhood similarity.")
    neighborClustersOut <- list(imgGenesTrimmed = imgGenesTrimmed, imgNeighborsTrimmed = imgNeighborsTrimmed, orderedMatrix = orderedMatrix)
    return(neighborClustersOut)
}
