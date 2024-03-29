% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepNeighbors.r
\name{prepNeighbors}
\alias{prepNeighbors}
\title{Function for preliminary analysis of neighboring genes from IMG}
\usage{
prepNeighbors(
  imgGenes = imgGenes,
  imgNeighbors = imgNeighbors,
  geneSeqs = geneSeqs,
  neighborSeqs = neighborSeqs,
  neighborsContext = neighborsContext,
  geneName = geneName,
  neighborNumber = neighborNumber,
  efiRepnodes = FALSE,
  trimShortClusters = FALSE,
  hypoAnalysis = FALSE,
  sysTerm = "nix",
  clustMethod = "tidygraph",
  numThreads = 1,
  neighborThreshold = 0.025,
  alphaVal = 0.95,
  bootStrap = 10,
  pidCutoff = 35,
  pepScreen = FALSE,
  pepMax = 150,
  alnClust = FALSE,
  hmmClust = FALSE
)
}
\arguments{
\item{imgGenes}{What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")}

\item{imgNeighbors}{What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")}

\item{geneSeqs}{What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")}

\item{neighborSeqs}{What is the file with the metadata for your gene of interest? Filename as string ("filename.txt")}

\item{neighborsContext}{File from neighborGenerate tying genes of interest to their neighbors. Filename as string ("filename.txt")}

\item{geneName}{What is the name of your gene? String ("genE")}

\item{neighborNumber}{How many neighbors do you want to look at on each side of the gene? Integer.}

\item{efiRepnodes}{Does the input dataset consist of EFI repnodes?  T/F, defaults to FALSE.}

\item{trimShortClusters}{Should gene clusters with fewer than the minimum neighbor number be removed? T/F value, defaults to FALSE.}

\item{hypoAnalysis}{Should hypothetical proteins be clustered and analyzeD? T/F value, defaults to FALSE.}

\item{sysTerm}{If running hypoAnalysis, what terminal are you using? String value ("wsl", "nix"), defaults to "nix".}

\item{clustMethod}{What clustering method should be used? String value ("tidygraph", "pvclust"), defaults to "tidygraph".}

\item{numThreads}{How many threads should processes use?  Number depends on your processor, defaults to 1 to be safe.}

\item{neighborThreshold}{What percent of gene clusters should a protein family occur in to be of interest? Number, defaults to 0.05}

\item{alphaVal}{What alpha value cutoff should be used for pvclust? Number from 0-1, defaults to 0.95}

\item{bootStrap}{How many bootstrap rounds does pvclust get? Integer, defaults to 10.}

\item{pidCutoff}{Below what percent ID should edges be deleted? Number from 1-100, defaults to 35.}

\item{pepScreen}{Should subgroups of peptides be identified (annotated or not)?  T/F, defaults to FALSE.}

\item{pepMax}{Maximum size (in aa) for peptides in pepScreen. Number, defaults to 150.}

\item{alnClust}{Should MAFFT alignments be made of members of a hypothetical protein cluster? T/F, defaults to FALSE.}

\item{hmmClust}{Should HMM models be made for a given hypothetical protein cluster? T/F, defaults to FALSE.}
}
\value{
List with trimmed metadata sets for both genes of interest and neighboring genes (additional files generated en route)
}
\description{
This function uses IMG metadatasets for genes of interest and their genomic neighbors to do some QC and identify subfamilies of hypothetical proteins.
}
\examples{
\dontrun{
prepNeighborsOutput <- prepNeighbors(imgGenes="geneFile.txt", 
                                     imgNeighbors="neighborsFile.txt", 
                                     geneSeqs="geneSeqs.fa", 
                                     neighborSeqs="neighborSeqs.fa", 
                                     neighborsContext = "context.txt", 
                                     geneName="genE",  
                                     neighborNumber=10, 
                                     sysTerm="wsl")
}

}
