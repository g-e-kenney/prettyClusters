% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neighborHypothetical.r
\name{neighborHypothetical}
\alias{neighborHypothetical}
\title{Subfunction to identify subgroups of hypothetical proteins in larger IMG datasets}
\usage{
neighborHypothetical(
  imgGenesData = imgGenesData,
  imgNeighborsData = imgNeighborsData,
  imgNeighborSeqs = imgNeighborSeqs,
  geneName = geneName,
  clustMethod = clustMethod,
  pidCutoff = pidCutoff,
  alphaVal = alphaVal,
  bootStrap = bootStrap,
  sysTerm = sysTerm,
  numThreads = numThreads,
  screenPep = screenPep,
  pepMax = pepMax,
  alnClust = alnClust,
  hmmClust = hmmClust
)
}
\arguments{
\item{imgGenesData}{Data frame from neighborPrep with metadata for your genes of interest.}

\item{imgNeighborsData}{Data frame from neighborPrep with metadata for neighbors of your genes of interest.}

\item{imgNeighborSeqs}{Data frame from neighborPrep with sequence information for neighbors of your gene of interest.}

\item{geneName}{Character string with your gene of interest's name}

\item{clustMethod}{Character string with the name of the desired clustering tool}

\item{pidCutoff}{The percent ID cutoff below which edges are not formed when using tidygraph}

\item{alphaVal}{The alpha value cutoff used for pvclust}

\item{bootStrap}{The bootstrap number used for pvclust}

\item{sysTerm}{The type of terminal (wsl vs linux/unix/macos) from which blastp and other commands will be run}

\item{numThreads}{The number of processor threads to be devoted to certain steps}

\item{screenPep}{Should we use peptide-friendly defaults in this run?}

\item{pepMax}{What's the maximum length to use for the "peptide" defaults?}

\item{alnClust}{Should we make MAFFT alignments of all clusters?}

\item{hmmClust}{Should we make HMM models of all clusters?}
}
\value{
Updated metadata for neighboring genes (additional files generated en route)
}
\description{
Given a desired genome neighborhood and lists of genes of interest and their neighbors, along with their protein sequences, this program uses all-by-all-blast & a clustering method of choice to identify subsets of hypothetical proteins that may represent meaningfully associated proteins.
}
\examples{
\dontrun{
imgNeighborsData <- neighborHypothetical(imgGenesData = imgGenesData, 
                                         imgNeighborsData = imgNeighborsData, 
                                         imgNeighborSeqs=imgNeighborSeqs, 
                                         geneName = "genE", 
                                         clustMethod = "tidygraph", 
                                         pidCutoff = 35, 
                                         sysTerm = "nix", 
                                         numThreads = 5) 
}

}
