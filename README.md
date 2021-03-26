# prettyClusters
A set of tools analyze and make non-hideous publication-friendly diagrams of genomic neighborhoods in microbial genomes.

## Important note
This is very much a work in progress, and I'm a biochemist doing terrible things to code. There will be bugs. I'll do what I can to address them; if you've come up with a fix, I'm happy to try to incorporate it!

## Why?
I spend a lot of time working with bacterial gene clusters.  I wanted some sort of tool that could:
- Export diagrams of gene clusters as a vector (not bitmap!) file suitable for figure layout.
- Export diagrams with multiple gene clusters in the same relative scale.
- Import gene metadata from the JGI's [IMG database](https://img.jgi.doe.gov/index.html), which has many genomes, metagenomes, and so on that are absent from NCBI's NR database (and UniProt), or that are present but poorly annotated in databases like NCBI/WGS.  
- Handle hypothetical and predicted proteins helpfully (i.e. by identifying new groups of hypothetical proteins that frequently appear in the genomic neighborhoods of interest).
- Integrate into workflows using sequence similarity networks generated via the [EFI-EST toolset](https://efi.igb.illinois.edu/efi-est/).
- Interrogate similarity of genomic neighborhoods without relying on the sequence similarity of genes of interest - I did not want to have to make the assumption that sequence similarity and gene cluster similarity necessarily track, since that's not always a sound assumption.
I haven't encountered anything that quite handles all of those things, so...

## The `prettyClusters` toolset
The Wiki entries contain a more detailed description of the use of specific functions.
### The core toolset
- [`generateNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/generateNeighbors)
- [`prepNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/prepNeighbors)
- [`analyzeNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/analyzeNeighbors)
- [`prettyClusterDiagrams`](https://github.com/g-e-kenney/prettyClusters/wiki/prettyClusterDiagrams)

### Accessory components
- [`gbToIMG`](https://github.com/g-e-kenney/prettyClusters/wiki/gbToIMG)
- [`incorpIprScan`](https://github.com/g-e-kenney/prettyClusters/wiki/incorpIprScan)
- [`repnodeTrim`](https://github.com/g-e-kenney/prettyClusters/wiki/repnodeTrim)

### External requirements
Beyond the (many) R packages used, `analyzeNeighbors` uses local installs of [mafft](https://cytoscape.org/) and [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) in its analysis of hypothetical proteins. It's highly likely I'll be adding steps that use [hmmer](http://hmmer.org) as well. Installation instructions are going to be system-specific, but if using Windows, these tools are currently only set up to deal with installations handled via the Windows Subsystem for Linux.  It is anticipated that these tools will often be paired with use of the [EFI-EST toolset](https://efi.igb.illinois.edu/efi-est/).  Those are online, but [Cytoscape](https://cytoscape.org/) is used for visualizing the sequence-similarity networks generated there.  Direct integration with [InterProScan](https://github.com/ebi-pf-team/interproscan) is not planned, since it is Linux-only, but access is required to produce the input for `incorpIprScan`.  

### Getting started
This package is nowhere near ready for legit repositories, so you are stuck with the development version.  Probably easiest to install via:
```
devtools::install_github("g-e-kenney/prettyClusters")
```
A sample session:
```
# load the package (a number of others should autoload)
# core packages include data.table, dplyr, gggenes, ggplot2, ggraph, pheatmap, pvclust, scales, seqinr, stringr, tibble, tidygraph, tidyr, and utils
# palette-focused packages include fishualize, ghibli, lisa, nord, rtist, scico, viridis, and wesanderson
# working with genbank files will add in genbankr and GenomicRanges
library(prettyClusters)

# run the initial neighbor generation 
generateNeighborsOut <- generateNeighbors(imgGenes = "imgGeneMetadata.txt", imgGeneSeqs = "imgGeneSeqs.fa", neighborNumber = 10, includeGene = TRUE, geneName = "genE") 

# after the neighbor list has been uploaded to the IMG database and the neighbor metadata files have been downloaded
prepNeighborsOut <- prepNeighbors(imgGenes = generateNeighborsOut$gene_oid, imgNeighbors = "imgNeighborMetadata.txt", geneSeqs = "imgGeneSeqs.fa", neighborSeqs = "imgNeighborSeqs.fa", neighborsContext = generateNeighborsOut$neighborsContext, geneName = "genE", neighborNumber = 10, sysTerm = "wsl", efiRepnodes = FALSE, neighborThreshold = 0.025, hypoAnalysis = TRUE, clustMethod = "tidygraph", numThreads = 7, alphaVal = 0.95, bootStrap = 10, pidCutoff = 35, trimShortClusters = TRUE)

# after generating an EFI-EST SSN and looking at the repnode options
repnodeTrimOut <- repnodeTrim(imgGenes = "imgGenesTrimmed.txt", imgNeighbors = "imgNeighborsTrimmed.txt", imgGeneSeqs = "imgGeneSeqs.fa", imgNeighborSeqs = "imgNeighborSeqs.fa", geneName = "genE", efiFullMetadata = "efiMetadataFull.csv", efiFinalMetadata = "efiMetadataRepnodes95.csv")

# passing on the trimmed data for neighborhood-based cluster analysis
analyzeNeighborsOut <- analyzeNeighbors(imgGenes =  repnodeTrimOut$repGenesTrimmed, imgNeighbors = repnodeTrimOut$repNeighborsTrimmed, efiRepnodes = TRUE, neighborThreshold = 0.025, geneName = "genE", autoClust = TRUE, clustMethod = "tidygraph", alphaVal = 0.95, bootStrap = 10, tgCutoff = 0.65)	

# generating the actual diagrams
prettyClusterDiagrams(imgGenes = neighborClustersOut$imgGenesTrimmed, imgNeighbors = neighborClustersOut$imgNeighborsTrimmed, geneFormat = "geneFormat.txt", geneName = "genE", efiRepnodes = TRUE, neighborNumber = 10, annotateGenes = TRUE, standAlone = FALSE, markClusters = TRUE, autoColor = TRUE, colorType = "fishualize", paletteInput = "Scarus_hoefleri", showScaffold = FALSE, alignToCore=TRUE, labelGenes = FALSE)
```
(Again, see the Wiki entries for per-function details, and soon for a vignette).  Illustrating the output of some of the components (or, in the case of the cluster diagrams themselves, just under 10% of the output):
<img src="https://github.com/g-e-kenney/prettyClusters/raw/master/20210106_pretty-cluster-general-01.png " width="100%" alt="genome neighborhood diagram output">
Notably, sequence similarity and genome neighborhood similarity are not always tightly coupled.  The analyses in `prettyClusters` make it possible to investigate a protein family along both axes.

## Development
### Planned additions
- A vignette!  
- Import from UniProt and GFF/GFF-3 formatted files.  As with `gbToIMG` this will likely be a separate function that can replace `generateNeighbors`, and it will likely suffer from the same data heterogeneity (and lack of gene family annotations) that that tool does.  Supplementation with `incorpIprScan` is likely still going to be advisable.
- Automatic generation of genome neighborhood diagrams for specific clusters (or for random representatives of a cluster) in `prettyClusterDiagrams`, since the diagrams get unwieldy with very large datasets...
- Single scale bar in `prettyClusterDiagrams`.  Probably will try generating a final fake "gene cluster" with single-nt "genes" every kb or something?
- Generation of HMMs for hypothetical protein families identified in `analyzeNeighbors`
- User-supplied HMMs for annotation of predefined custom protein (sub)families in `analyzeNeighbors`
- Options to let the user specify distance and clustering methods in `prepNeighbors` and `analyzeNeighbors`
### Known issues
- Auto-annotation in `prettyClusterDiagrams` may miss genes if their initial family categorization (or ORF-finding) was poor!  This is doubly the case for GenBank-derived files (use of `incorpIprScan` can improve annotation, but not ORF-finding).
- Generation of hypothetical protein and genome neighborhood clusters is approximate and sensitive to user-supplied cutoffs and distance/clustering methods. No way around this.
- Forward- and reverse-facing genes are on the same vertical level in `prettyClusterDiagrams`. I personally find it visually clearer to have forward genes above the line and reverse genes below, but will need to probably do a bunch more digging into [gggenes](https://github.com/wilkox/gggenes) and [ggplot2](https://github.com/tidyverse/ggplot2) to figure out if/how I can make it happen.
- Distance between genes of interest and their neighbors is not taken into account.  I have done a few initial analyses using weighted distances rather than binary present/absent values; it is not clear they've added much more info than the binary value-based analyses, and they're more complicated to run.  Could revisit?
- Use of non-WSL `mafft` and `blast` installs on Windows isn't built-in at the moment; it's a low priority, given how easy WSL is to get up and running.
