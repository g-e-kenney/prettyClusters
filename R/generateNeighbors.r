#' Standalone Function to Generate Lists of Neighboring Genes from IMG gene_oids
#'
#' This function uses IMG metadatasets to generate relatively attractive gene cluster/genomic neighborhood diagrams that are scaled and that generate vector graphics
#' @param imgGenes What is the file with the metadata for your gene family of interest? Filename as string ("filename.txt")
#' @param imgGeneSeqs What is the file with sequences of your gene family of interest? Filename as string ("filename.fa")
#' @param neighborNumber How many neighbors do you want to look at on each side of the gene? Integer.
#' @param includeGene Do you want your genes of interest in the diagrams? T/F value, defaults to TRUE.
#' @param geneName What is the name of your gene? Gene name as string ("genE")
#' @return Lists of neighboring genes to upload to your IMG Gene Cart
#' @export
#' @examples
#' generateNeighborsOutput <- generateNeighbors(imgGenes="geneFile.txt", imgGeneSeqs="geneSeqs.fa", neighborNumber=10, geneName="genE")
#'
generateNeighbors <- function(imgGenes = imgGenes, imgGeneSeqs = imgGeneSeqs, neighborNumber = neighborNumber, includeGene = TRUE, geneName = geneName) {
                                        # boring setup
    fileDate <- format(Sys.Date(),format="%Y%m%d")
    fileName <- paste(fileDate,"_generateNeighbors_",geneName,sep="")
    fileNameContext <- paste(fileName, "_context.txt", sep="")
    fileNameSeqs <- paste(fileName,".fa",sep="")
    neighbors <- list()
    gene_oid <- list()
    source_gene_oid <- list()
    source_scaffold_id <- list()
    ## input is a standard IMG file for your genes of interest
    ## can be generated via blast or families or whatever
    inputGenes <- read.csv(imgGenes, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
    imgGeneSeqs <- seqinr::read.fasta(file=imgGeneSeqs, seqtype="AA", whole.header=FALSE, as.string=TRUE, set.attributes=FALSE)
    ## let's quickly just save that fasta file with simplified headers for EFI use
    write.fasta(imgGeneSeqs, names=names(imgGeneSeqs),file.out=fileNameSeqs)
    ## and onwards
    geneList <- inputGenes$gene_oid
    scaffList <- inputGenes$Scaffold.ID
    geneNum <- length(geneList)
    nnum <- neighborNumber*2+1
    nstart <- neighborNumber*-1
    ## this assumes your gene is at the middle with the same range on both sides; you could change it manually
    ## a symmetric setup is proably safer though since the gene could be on either strand...
    yourgene <- (nnum/2) + .5
                                        # iterating through the gene list and making neighbors    
    for (i in 1:geneNum)  {
        for (j in 1:nnum) {
            if(j != yourgene) {
                neighbors <- rbind(neighbors, geneList[i] + nstart + j - 1)
				gene_oid <- rbind(gene_oid, geneList[i] + nstart + j - 1)
				source_gene_oid <- rbind(source_gene_oid, geneList[i])
				source_scaffold_id <- rbind(source_scaffold_id,scaffList[i])
            } else {
                if(includeGene == TRUE) {
                    neighbors <- rbind(neighbors, geneList[i] + nstart + j - 1)
                    gene_oid <- rbind(gene_oid, geneList[i] + nstart + j - 1)
                    source_gene_oid <- rbind(source_gene_oid, geneList[i])
                    source_scaffold_id <- rbind(source_scaffold_id,scaffList[i])
                } else {
                    print("Your gene of interest is not included in the neighbors metadata.")
                }
            }
        }
    }
    # exporting things
    ## this neighbor list contains the gene_oid for each neighbor, the source_gene_oid (i.e. the gene of of interest it neighbors), and the source scaffold
    ## this can be useful later for tracking stuff like "neighbors" that are not on the same scaffold
    neighborsContext <- data.frame(gene_oid, source_gene_oid, source_scaffold_id)
    colnames(neighborsContext) <- c("gene_oid", "source_gene_oid", "source_scaffold_id")
    neighborsContext <- as.matrix(neighborsContext)
    write.table(neighborsContext, file=fileNameContext, row.names=FALSE, col.names = TRUE, quote=FALSE, sep="\t")
    neighborLength <- length(neighbors)
    neighbors <- as.matrix(neighbors)
    ## now making the versions that'll be useful for uploading to IMG
    ## particularly given the 20k gene cart limit
    ## (splits into a bunch of files that can be uploaded separately)
    if(neighborLength < 20000)  {
        fileNameOut <- paste(fileName, "_neighbors.txt", sep="")
        colnames(neighbors) <- "gene_oid"
        write.table(neighbors, file=fileNameOut, row.names=FALSE, col.names = TRUE, quote=FALSE, sep="\t")
    } else {
        num20k <- ceiling(neighborLength/20000)
        for (i in 1:num20k) {
            neighborsStart <- 1 + (i-1)*20000
            neighborsEnd <- 20000*i
            neighborsPart <- neighbors[neighborsStart:neighborsEnd]
            colnames(neighborsPart) <- "gene_oid"
            fileNameOut <- paste(fileName, "_neighbors_",i,".txt", sep="")
            write.table(neighbors, file=fileNameOut, row.names=FALSE, col.names = TRUE, quote=FALSE, sep="\t")
        }
    }
    print("Neighbor list generated.")  
	return(list(neighbors = neighbors, neighborsContext = neighborsContext))
}
