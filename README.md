# prettyClusters
A set of tools analyze and make non-hideous publication-friendly diagrams of bacterial gene clusters.

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

## The core `prettyClusters` tool set
### Components of the prettyClusters toolset
- `generateNeighbors`
- `prepNeighbors`
- `analyzeNeighbors`
- `prettyClusterDiagrams`

### Accessory components
- `gbToIMG`
- `incorpIprScan`
- `repnodeTrim`

### External requirements
Beyond the (many) R packages used, `analyzeNeighbors` uses local installs of [mafft](https://cytoscape.org/) and [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) in its analysis of hypothetical proteins. It's highly likely I'll be adding steps that use [hmmer](http://hmmer.org) as well. Installation instructions are going to be system-specific, but if using Windows, these tools are currently only set up to deal with installations handled via the Windows Subsystem for Linux.  It is anticipated that these tools will often be paired with use of the [EFI-EST toolset](https://efi.igb.illinois.edu/efi-est/).  Those are online, but [Cytoscape](https://cytoscape.org/) is used for visualizing the sequence-similarity networks generated there.  When 

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
Illustrating the output of some of the components (or, in the case of the cluster diagrams themselves, just under 10% of the output):
<img src="https://github.com/g-e-kenney/prettyClusters/raw/master/20210106_pretty-cluster-general-01.png " width="100%" alt="genome neighborhood diagram output">
Notably, sequence similarity and genome neighborhood similarity are not always tightly coupled.  The analyses in `prettyClusters` make it possible to investigate a protein family along both axes.

### `generateNeighbors`
This tool takes advantage of the fully numeric and contiguous nature of gene_oids in the IMG database.  Given an IMG metadata table for a set of genes of interest (identified via BLAST, protein family-based filtering, or other methods), this tool generates IDs for genes that ought to be in the same neighborhood. These gene lists can be used to download the data-rich IMG metadata files for all all of those genes from the IMG database.  Note that if you're starting from GenBank files, you can skip this and use the `gbToIMG` accessory tool instead (possibly accompanied by `incorpIprScan` to improve annotations), although like all tools that involve parsing GenBank files, it's a bit finicky.
#### Use of `generateNeighbors`
```
generateNeighborsOut <- generateNeighbors(imgGenes = "imgGeneMetadata.txt", imgGeneSeqs = "imgGeneSeqs.fa", neighborNumber = 10, includeGene = TRUE, geneName = "genE") 
```
#### Options
- `imgGenes` Filename. IMG-formatted metadata file for genes of interest. Required.
- `imgGeneSeqs` Filename. Fasta-formatted sequence file for genes of interest. Required.
- `neighborNumber` Integer. The number of neighbors upstream and downstream of the gene of interest to be analyzed. Required.
- `includeGene` T/F. Signals whether to include the genes of interest in the lists of "neighbors."  Currently this is suggested, since it makes diagram generation later easier. Defaults to TRUE.
- `geneName` Text. This is to the gene name in autogenerated filenames.  Required.
#### Output
- `20210101_generateNeighbors_genE.fa` File. Fasta-formatted file for the input sequences with simplified headers containing only the gene_oid.
- `20210101_generateNeighbors_genE_context.txt` File. Tab-delimited table with three columns: gene_oid (the neighbor gene_oid), source_gene_oid (the gene_oid for which the neighbor was generated) and scaffold_id (the scaffold on which the original gene_oid was found.)
- `20210101_generateNeighbors_genE_neighbors.txt` File. Tab-delimited table with a single column (gene_oid).  There may be multiple numbered derivatives ('20210101_generateNeighbors_genE_neighbors_1.txt') if there are 20k+ neighbors.
- `generateNeighborsOut` List. Contains `generateNeighborsOut$gene_oid` (a data frame with a single column listing the neighboring genes) and `generateNeighborsOut$neighborsContext` (a data frame with three columns:  gene_oid (the neighbor gene_oid), source_gene_oid (the gene_oid for which the neighbor was generated) and scaffold_id (the scaffold on which the original gene_oid was found.))
      
### `prepNeighbors`
This is a wrapper for two subfunctions:
- `neighborHypothetical`
This tool identifies sets of hypothetical proteins that are overrepresented in the genomic neighborhood of the genes of interest.  Beyond R packages, this tool uses the NCBI blast+ package and MAFFT, which need to be pre-installed.  Note that for identifying subgroups of hypothetical proteins, an all-by-all blast step is used, which can lead to long computation time on large datasets!  [Tidygraph](https://github.com/thomasp85/tidygraph) (with [ggraph](https://github.com/thomasp85/ggraph)) and/or [pvclust](https://github.com/shimo-lab/pvclust) are used identify and visualize possible cutoff points for defining those subgroups of hypothetical proteins. For computational reasons, tidygraph is likely to be quicker for large datasets.
- `neighborTrim`
This tool identifies genomic neighborhoods where contig ends result in small neighborhoods.  Since truncated neighborhoods can mis-weight genome neighborhood clustering, they can be removed from further analyses.
#### Use of `prepNeighbors`
```
prepNeighborsOut <- prepNeighbors(imgGenes = "imgGeneMetadata.txt", imgNeighbors = "imgNeighborMetadata.txt", geneSeqs = "imgGeneSeqs.fa", neighborSeqs = "imgNeighborSeqs.fa", neighborsContext = "neighborsContext.txt", geneName = "genE", neighborNumber = 10, sysTerm = "wsl", efiRepnodes = FALSE, neighborThreshold = 0.025, hypoAnalysis = TRUE, clustMethod = "tidygraph", numThreads = 7, alphaVal = 0.95, bootStrap = 10, pidCutoff = 35, trimShortClusters = TRUE) 
```
#### Options
- `imgGenes` Filename. IMG-formated metadata file for genes of interest. Required.
- `imgNeighbors` Filename. IMG-formated metadata file for neighbors of genes of interest. Required.
- `geneSeqs` Filename. FASTA-formatted sequence file for genes of interest. Required.
- `neighborSeqs` Filename. FASTA-formatted sequence file for neighbors of genes of interest. Required.
- `neighborsContext` Filename. File generated in `generateNeighbors` that ties the metadata of the neighbors to the gene_oids and scaffold IDs of the original genes of interest.  Required.
- `geneName` Text. the name of your gene family of interest - used for filenames and figure labels. Required.
- `neighborNum` Integer. The number of neighboring genes on either side of the gene of interest that will be investigated. Required.
- `sysTerm` Text. Used to identify the sort of terminal (current options: "wsl" for the Linux subsystem on Windows and "nix" for generic Linux, MacOS, etc.)  Used to run blast and mafft commands. Required.
- `efiRepnodes` T/F. Just used for filename generation. Defaults to FALSE.
- `neighborThreshold` Number. (0-1.00). Fraction of gene clusters in which a protein family must occur to be included in further analyses.  Broadly, used to determine how widely the net for subsets of hypothetical proteins is cast. Defaults to 0.025 (2.5%).
- `hypoAnalysis` T/F. Signals whether or not subsets of similar hypothetical proteins will be identified. Defaults to TRUE.
- `clustMethod` Text. The name of the package to be used to identify clusters of similar hypothetical proteins. [Tidygraph](https://github.com/thomasp85/tidygraph)  and pvclust are the current options, with the former suggested. Defaults to "tidygraph".
- `numThreads` Integer. Number of threads to use during analysis. Maximum value depends on your computer. Defaults to 1.
- `alphaVal` Number. (0.0-1.0). Alpha value cutoff for analyzing pvClust data.  Defaults to 0.95.
- `bootStrap` Number. Bootstrap number for pvClust.  Can be fairly computationally intensive: 100 is a good starting point for small datasets (e.g. 100 genes of interest, 5 genes on either size), and 10 for larger (500 genes of interest or more, 10+ genes in either direction.)  Defaults to 10.
- `pidCutoff` Number. (0-100). %ID cutoff to be used for edge generation when analyzing hypothetical proteins. Defaults to 35 (35% ID).
- `trimShortClusters` T/F. Signals whether or not to remove gene clusters where a truncated contig interrupts the cluter and a sub-sized cluster remains. (Recommended when running further analyses and not just generating cluster diagrams.)  Defaults to TRUE.
#### Output
- `20210101_neighborHypothetical_genE_hypoSeqs.fa` File. FASTA file containing amino acid sequences of all putative hypothetical proteins
- `20210101_neighborHypothetical_genE_blastFile.txt` File. All-by-all BLAST results for hypothetical sequences in a tab-delimited format
- `20210101_neighborHypothetical_genE_pvClustFile.txt` File. Present only if pvclust is used. Cluster information in text file format.
- `20210101_neighborHypothetical_genE_pvClustTree.pdf` File. Present only if pvclust is used. Image of pvclust tree at final alpha value cutoff, as a .pdf file.
- `20210101_neighborHypothetical_genE_pvCluster_x.fa` File. Present only if pvclust is used. In subfolder with other equivalent files. FASTA file containing amino acid sequences of members of a specific hypothetical protein cluster. "x" is the cluster number. 
- `20210101_neighborHypothetical_genE_pvCluster_x_mafft.fa` File. Present only if pvclust is used. In subfolder with other equivalent files. FASTA file containing amino acid sequences of members of a specific hypothetical protein cluster, MAFFT-aligned. "x" is the cluster number. 
- `20210101_neighborHypothetical_genE_networkFull_xx.pdf` File. Present only if tidygraph is used for clustering.  "xx" is the percent ID cutoff. PDF of full network, with the gene of interest highlighted (if applicable)
- `20210101_neighborHypothetical_genE_networkClusters_xx.pdf` File. Present only if tidygraph is used. Image of tidygraph network at the percent ID cutoff of "xx" with small clusters removed.
- `20210101_neighborHypothetical_genE_tgCluster_x.fa` File. Present only if tidygraph is used. In subfolder with other equivalent files. FASTA file containing amino acid sequences of members of a specific hypothetical protein cluster.  "x" is the cluster number. 
- `20210101_neighborHypothetical_genE_tgCluster_x_mafft.fa` File. Present only if pvclust is used. In subfolder with other equivalent files.  FASTA file containing amino acid sequences of members of a specific hypothetical protein cluster, MAFFT-aligned.  "x" is the cluster number. 
- `20210101_neighborHypothetical_genE_clusterList.txt` File. Tab-delimited text file with two columns listing the gene_oid of every hypothetical protein and its cluster group (entered in the annotation as a "Hypofam" family.)
- `20210101_neighborTrim_genE_imgGenesTrimmed.txt` File. IMG-formatted tab-delimited metadata text file of genes of interest with genes from short contigs removed.
- `20210101_neighborTrim_genE_imgNeighborsTrimmed.txt` File. IMG-formatted tab-delimited metadata text file of neighbors of genes of interest with genes from short contigs or from incorrect scaffolds. Contains Hypofam annotations if neighborHypothetical was run.
- `20210101_neighborTrim_genE_imgNeighborsContextTrimmed.txt` File. Tab-delimited list of IMG gene_oids for neighbors, their source genes, and the scaffolds of the source genes, with neighbors from short contigs or on incorrect scaffolds removed.
- `20210101_neighborTrim_genE_imgGeneSequencesTrimmed.fa` File. FASTA-formatted set of amino acid sequences, trimmed to match the updated gene list.
- `20210101_neighborTrim_genE_imgNeighborsSequencesTrimmed.fa` File. FASTA-formatted set of amino acid sequences, trimmed to match the updated neighbor list.  

### `analyzeNeighbors`
This is separate from `prepNeighbors` because most of the analyses make sense to do _after_ submitting data to EFI for SSN generation, so that they can be run on representative nodes and their neighbors only, or full datasets.  However, EFI processing is not required to proceed with data analysis.  This is a wrapper for four subfunctions:
- `neighborCatalog` Identifies what protein families are present in all genomic neighborhoods, and which are above the neighborThreshold cutoff and will be used in further analyses.
- `neighborHere` Identifies which of those families are present in the genomic neighborhoods of genes of interest. Currently a binary present/absent value; it's not clear that weighting 
- `neighborMatrix` Generates a horrible matrix based on that binary data.
- `neighborCluster` Clusters genomic neighborhoods on the basis of their similarity, using that data.  Automatic assignment to clusters is an option, using tidygraph (currently recommended) or pvClust.  The exported metadata can - among other things - be re-imported to Cytoscape while viewing the EFI-derived SSN.
#### Use of `analyzeNeighbors`
```
analyzeNeighborsOut <- analyzeNeighbors(imgGenes =  "repnodeGeneMetadata.txt", imgNeighbors = "repnodeNeighborMetadata.txt", efiRepnodes = TRUE, neighborThreshold = 0.025, geneName = geneName, autoClust = TRUE, clustMethod = "tidygraph", alphaVal = 0.95, bootStrap= 10)	
```
#### Options
- `imgGenes` Filename. IMG-formatted metadata file for genes of interest. Can also be the name of the equivalent object from the previous suite member (i.e. neighborPrep, if no EFI repnodes are used). Required.
- `imgNeighbors` Filename. IMG-formatted metadata file for neighbors of genes of interest. Can also be the name of the equivalent object from the previous suite member (i.e. neighborPrep, if no EFI repnodes are used). Required.
- `efiRepnodes` T/F. Just for labeling filenames (will add "repnodes" to distinguish files from full-network analyses). Defaults to FALSE.
- `neighborThreshold` Number. (0.0-1.0) The fraction of genomic neighborhoods that a protein family needs to show up in to be highlighted in further analyses. Required.
- `geneName` Text. The name of your gene family of interest - used for filenames and figure labels. Required.
- `autoClust` T/F. Signals whether or not to automatically assign cluster numbers for subsets of gene clusters, based on their similarity.  If used, requires a choice between tidygraph (currently recommended) and pvclust.  If not chosen, sets of genome neighborhoods will still be clustered according to their similarity via hierarchical clustering, as will sets of protein families that co-occur in similar sets of genomic neighborhoods. However, cutoffs for cluster identity will need to be drawn by the user, on their own, and added manually to the metadata file for use with `prettyClusterDiagrams`. Defaults to TRUE.
- `clustMethod` Text. Required if autoClust is TRUE.  "pvclust" or "tidygraph" are the two options I currently have working.  I'd suggest [tidygraph](https://github.com/thomasp85/tidygraph) (paired with [ggraph](https://github.com/thomasp85/ggraph)) based on the datasets I've looked at so far, but [pvclust](https://github.com/shimo-lab/pvclust)  is still an option. Defaults to "tidygraph".
- `alphaVal` Number. (0.0-1.0). Alpha value cutoff for use with [pvclust](https://github.com/shimo-lab/pvclust).  Defaults to 0.95.
- `bootStrap` Number. Number of bootstrap rounds for use with [pvclust](https://github.com/shimo-lab/pvclust). Defaults to 10.
#### Output
- `20210101_neighborCatalog_genE_family-abundance.txt` File. Tab-delimited table with protein families over the 2.5 % abundance cutoff, and their actutal abundances.
- `20210101_neighborHere_genE_neighborBinary.csv` File. CSV table with binary values representing assignment of protein families of interest to neighboring genes.
- `20210101_neighborMatrix_genE_neighborMatrix.csv` File. CSV table containing a matrix of binary values representing the presence or absence of protein families in the neighborhoods of specific genes of interest.
- `20210101_neighborClusters_genE_clusterList` File. Text file listing clusters.
- `20210101_neighborClusters_genE_clusterGeneMetadata.txt` File. Tab-delimited table of IMG metadata for genes of interest, with genome neighborhood cluster info added.
- `20210101_neighborClusters_genE_clusterNeighborMetadata` File. Tab-delimited table of IMG metadata for neighbors of genes of interest, with genome neighborhood cluster info added.
- `20210101_neighborClusters_genE_heatmap.pdf` File. PDF of heatmap for hierarchically clustered genome neighborhoods and protein families, with cluster info.
- `20210101_neighborClusters_genE_heatmap.png` File. PNG of heatmap for hierarchically clustered genome neighborhoods and protein families, with cluster info.
- `20210101_neighborClusters_genE_heatmap.ps` File. PS of heatmap for hierarchically clustered genome neighborhoods and protein families, with cluster info.  (I've occasionally had some PDF export issues, so this is just a backup?  I dunno.)
- `20210101_neighborClusters_genE_pvClust.txt` File. Text file with pvclust assignments, if applicable.
- `20210101_neighborClusters_genE_pvClustTree.pdf` File. PDF with pvClust tree and identified clusters, if applicable.
- `20210101_neighborClusters_genE_networkFile.txt` File. PDF with node/edge information for tidygraph network (can be imported into Cytoscape), if applicable.
- `20210101_neighborClusters_genE_networkFigureFile.pdf` File. PDF with tidygraph-generated network representing genome neighborhood similarity, if applicable.
- `20210101_neighborClusters_genE_clustColorFile.txt` File. Key to coloring of clusters, just in case.
- `neighborClustersOut` List.  Contains `neighborClustersOut$imgGenesTrimmed` (data frame with IMG metadata for genes of interest, with genome neighborhood cluster info added), `neighborClustersOut$imgNeighborsTrimmed` (data frame with IMG metadata for neighbors of genes of interest, with genome neighborhood cluster info added), and `neighborClustersOut$orderedMatrix` (matrix of binary values representing the presence or absence of protein families in the neighborhoods of specific genes of interest)

### `prettyClusterDiagrams`
This will make what are hopefully visually passable gene cluster diagrams, exported in vector and bitmap formats. Color schemes are  generated from chosen palettes.  Gene types to be highlighted are specified by the user, either manually or via auto-assignment based on membership in specific protein families (Pfam, TIGRfam, and IMG Term, along with hypothetical protein sets identified in `prepNeighbors` - InterPro to be added as soon as IMG adds it to their metadata export).  This is designed to be run with other components of this suite, but if not run in tandem with `prepNeighbors` and `analyzeNeighbors`, additional QC components will confirm that genes are co-localized to the same scaffold and so on.  Similarly, manual assignment of cluster numbers for visualization is also possible. Note that visualization depends heavily on [gggenes](https://github.com/wilkox/gggenes), with default settings biased towards use cases that involve visualization of many gene clusters. 
#### Use of `prettyClusterDiagrams`
```
prettyClusterDiagrams(imgGenes = "repnodeGeneMetadata.txt", imgNeighbors = "repnodeNeighborMetadata.txt", geneFormat = geneFormat, geneName = "genE", efiRepnodes = TRUE, neighborNumber = 10, annotateGenes = TRUE, standAlone = FALSE, markClusters = TRUE, autoColor = TRUE, colorType = "nord", paletteInput = "aurora", showScaffold = FALSE, alignToCore=TRUE, labelGenes = FALSE)
```
#### Options
- `imgGenes` Filename. IMG-formatted metadata file for genes of interest. Can also be the name of the equivalent object or object subset from the previous suite member. Required.
- `imgNeighbors` Filename. IMG-formatted metadata file for neighboring genes. Can also be the name of the equivalent object or object subset from the previous suite member. Required.
- `geneFormat` Filename. Comma-delimited text file that provides key for annotation. Columns: geneSymbol (genE, annotation gene name to be assigned), fusion (yes/no, is this a fused version of other genes in the list), Requirement (all/any, which families must match to assign annotation - all/any), Pfam (pfam00001, with additional families separated by spaces), Tigrfam (TIGR00001, additional families separated by spaces), Hypofam (hypofam1, ditto), IMG.Term (0001, ditto), InterPro (IPR000001, ditto), Color (#000000, if manual color assignment is required, otherwise leave empty).  Required.
- `geneName` Text. The name of your gene family of interest - used for filenames and figure labels. Required.
- `efiRepnodes`	T/F. Signals when you are working with a repnode subset - just for labeling filenames. Defaults to FALSE.
- `neighborNumber` Integer. The number of neighboring genes that will be shown on either side. Required.
- `annotateGenes` T/F. Uses existing annotation and user-provided family-based rules assign color to genes in diagrams. If true, will require the `geneFormat` file.  Recommended so that you don't end up with a million different colors for spuriously annotated genes, unless you already have manually curated annotation.  (Currently somewhat slow, though.)  Defaults to TRUE.
- `standAlone` T/F. If true, the assumption is that the data file did not go through some or all of the rest of the pipeline, and that QC steps from `analyzeNeighbors` should be performed. Defaults to FALSE.
- `markClusters` T/F. Signals whether or not cluster labels should be included in the diagram. Note that this assumes you've run `analyzeNeighbors` and have a clustNum and a clustOrd column in the neighbor metadata file.  You can alternately manually add clustNum values to the input neighbor metadata file. Defaults to FALSE.
- `autoColor` T/F. Signals whether or not to use a user-provided palette of colors for genes (in `geneFormat`) or whether to auto-generate one. Defaults to TRUE.
- `colorType` Text. Identifies what color library (nord, etc.) to use.  Options include [fishualize](https://github.com/nschiett/fishualize), [ghibli](https://github.com/ewenme/ghibli), [lisa](https://github.com/tyluRp/lisa), [nord](https://github.com/cran/nord), [rtist](https://github.com/tomasokal/rtist), [scico](https://github.com/thomasp85/scico), [viridis](https://github.com/sjmgarnier/viridis), and [wesanderson](https://github.com/karthik/wesanderson) because finding color sets that work well for a given gene number and layout can be hard, and I'd rather like this to be able to live up to its name. Defaults to "nord".
- `paletteInput` Text. Identify what palette within that color library to choose (e.g. aurora within nord). Defaults to "aurora".
- `showScaffold` T/F. Signals whether or not to include the scaffold ID in labels.  Not recommended as a default since it makes the labels really long. Defaults to FALSE.
- `alignToCore`	T/F. Signals whether or not to center the gene family of interest in the figure. (Generally recommended). Defaults to TRUE.
- `labelGenes` T/F. Signals whether or not individual genes should be labeled in diagrams (not recommended for larger datasets!) Defaults to FALSE.
#### Output 
- `20210101_prettyClusters_with-axes.png` File. PNG file for genome neighborhood diagrams, containing coordinates for each cluster.
- `20210101_prettyClusters_with-axes.pdf` File. PDF file for genome neighborhood diagrams, containing coordinates for each cluster.
- `20210101_prettyClusters_no-axes.png` File. PNG file for genome neighborhood diagrams, without coordinates for each cluster.
- `20210101_prettyClusters_no-axes.pdf` File. PDF file for genome neighborhood diagrams, without coordinates for each cluster.
- `20210101_prettyClusters_annotation.txt` File. Tab-delimited table of minimal IMG-derived metadata with added metadata from this function (BGC IDs, gene assignments, etc.) added.

## Accessory functions
### `gbToIMG`
This takes a directory of GenBank files and makes fake IMG-formatted metadata files for them, with the goal of making non-IMG files more compatible with the rest of this toolset.  This includes faux IMG-style gene_oids.  GenBank files from the ENA or the NCBI databases often contain less information on the annotation sources used to assign genes, so expectations should be tempered accordingly (antiSmash-annotated .gb files contain a great deal but differ too much from standard .gb files to be parsed by this tool.)  Relies on the [genbankr](https://bioconductor.org/packages/release/bioc/html/genbankr.html) and [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) packages heavily.  
#### Use of `gbToIMG`
```
gbToIMGOutput <- gbToIMG(dataFolder="/user/data/gbfiles", goiListInput = "goiList.txt", neighborNum=5, geneName="genE", removeDupes=TRUE)
```
#### Options
- `dataFolder` Folder path. For folder containing all .gb/.gbk files to be analyzed. Required.
- `goiListInput` Filename. For a text file containing a list of genes of interest by the names they are likely to be identified by in their GenBank files (probably locus_tag). Names are provided on single lines, with a header.  Required. 
- `neighborNum` Integer. Number of neighbors to be provided for each gene of interest. Required. 
- `geneName` Character string. Name of gene family of interest (purely for file naming). Required. 
- `removeDupes` Boolean.  Removes duplicated entries (probably an OK default, unless you have copies of genomes with different annotations that you want evaluated independently.  Required.
- `scaffoldGenBase` Integer. Initial value for generating IMG-style faux scaffold IDs (new IDs are generated for each run currently). Defaults to 30000000000.   - `genomeGenBase` Integer. Initial value for generating IMG-style faux genome IDs (new IDs are generated for each run currently). Defaults to 40000000000.    
#### Output
- `20210101_gb2img_genE_geneSeqs.fa` File. Fasta-formatted file for the protein sequences of genes of interest with simplified headers containing only the IMG-style gene_oids.
- `20210101_gb2img_genE_neighborSeqs.fa` File. Fasta-formatted file for the protein sequences of neighbors of genes of interest with simplified headers containing only the IMG-style gene_oids.
- `20210101_gb2img_genE_neighborcontext.txt` File. Tab-delimited table with three columns: gene_oid (the neighbor gene_oid), source_gene_oid (the gene_oid for which the neighbor was generated) and scaffold_id (the scaffold on which the original gene_oid was found.)
- `20210101_gb2img_genE_geneData.txt` File. Tab-delimited metadata table, IMG-styled, for genes of interest.
- `20210101_gb2img_genE_neighborData.txt` File. Tab-delimited metadata table, IMG-styled, for neighbors of genes of interest.
- `gbToIMGOutput` List. Contains `gbToIMGOutput$geneData` (a data frame of IMG-styled metadata for genes of interest), `gbToIMGOutput$neighborData` (a data frame of IMG-styled metadata for neighbors of genes of interest), and `gbToIMGOutput$neighborsContext` (a data frame connecting neighbors to the genes of interest they are associated with, along with their scaffolds).

### `incorpIprScan`
NCBI GenBank-derived metadata often does not list the families from which gene annotations are derived, which is a problem for this workflow.  This function takes a tab-delimited InterProScan output file and an IMG-formatted metadata file for the same genes and adds Pfam, TIGRfam, and/or InterPro families to the metadata file.  Note that you may need to run [InterProScan](https://github.com/ebi-pf-team/interproscan) on your local cluster or install it on a linux machine or VM; there is not a Windows or MacOS-compatible version, and actually running it is thus not included in this workflow.  That said, a basic command outputting the required tab-delimited output file looks like this:  
```
./interproscan.sh -i neighborSeqs.fa -t p -o neighborSeqs_iprScan.txt -f TSV -iprlookup
```
#### Use of `incorpIprScan`
```
incorpIprScanOutput <- incorpIprScan(iprScanSource = "iprScanfile.txt", imgNeighborsSource = "imgNeighborsData.txt", geneName = "genE", addPfam = TRUE, addTigrfam = TRUE, addIPRfam = TRUE)
```
#### Options
- `iprScanSource` Folder path. For folder containing all .gb/.gbk files to be analyzed. Required.
- `imgNeighborsSource` Filename. For a text file IMG-formatted metadata for neighbors of genes of interest.  Required. 
- `geneName` Character string. Name of gene family of interest (purely for file naming). Required. 
- `addPfam` Boolean.  Specifies whether or not Pfam annotations should be added to the metadata sheet. 
- `addTigrfam` Boolean.  Specifies whether or not Pfam annotations should be added to the metadata sheet. 
- `addIprfam` Boolean.  Specifies whether or not Pfam annotations should be added to the metadata sheet. 
#### Output
- `20210101_incorpIprScan_genE_neighborData.txt` File. Fasta-formatted file for the protein sequences of genes of interest with simplified headers containing only the IMG-style gene_oids.
- `incorpIprScanOutput` List. Contains `imgNeighborsData` (a data frame of IMG-styled metadata for neighbors of genes of interest, now with more annotation info.)

### `repnodeTrim`
Generally, after running `prepNeighbors`, I submit the trimmed protein sequences for my genes of interest to the [EFI-EST server](https://efi.igb.illinois.edu/efi-est/) via option C (a user-uploaded fasta file).  The EFI-EST toolsets are designed to work with UniProt data, and so genes are assigned faux-UniProt IDs, which can complicate tying the output back to analyses that use the IMG gene_oid values.  This tool connects the two IDs.  Additionally, during the EFI-EST SSN setup, sequences below/above length cutoffs are often trimmed, and a %ID cutoff is established, above which highly similar sequences are grouped and a representative node is chosen for SSN visualization. These "repnodes" are also useful for avoiding over-weighting of a network towards data from very highly sequenced species and genera (e.g. pathogens) and for working with very large and computationally intensive networks.  Thus, it's often helpful to go ahead using only the representative nodes.  This function provides trimmed versions of gene and neighbor metadata and sequence for this purpose (but use of this function is not required for the . 
#### Use of `repnodeTrim`
```
repnodeTrimOut <- repnodeTrim(imgGenes = "imgGeneMetadata.txt", imgNeighbors = "imgNeighborMetadata.txt", imgGeneSeqs = "imgGeneSeqs.fa", imgNeighborSeqs = "imgNeighborSeqs.fa", geneName = "genE", efiFullMetadata = "efiMetadataFull.csv", efiFinalMetadata = "efiMetadataRepnodes95.csv")
```
#### Options
- `imgGenes` Filename. IMG-formatted metadata file for genes of interest. Required.
- `imgNeighbors` Filename. IMG-formatted metadata file for neighbors of genes of interest. Required.
- `geneSeqs` Filename. FASTA-formatted sequence file for genes of interest. Required.
- `neighborSeqs` Filename. FASTA-formatted sequence file for neighbors of genes of interest. Required.
- `geneName` Text. the name of your gene family of interest - used for filenames and figure labels. Required.
- `efiFullMetadata` Filename. The node metadata file, exported for the full network from Cytoscape. Note for large networks: Layout and visualization of the network is not necessary to export this file!  Required.
- `efiFinalMetadata` Filename.  A second node metadata file, exported from the network at the appropriate repnode cutoff. Required.
#### Output
- `20210101_repnodeTrim_genE_imgGenes.txt` File. Tab-delimited table, contains metadata for all input genes, but with EFI and repnode IDs added.
- `20210101_repnodeTrim_genE_imgNeighbors.txt` File. Tab-delimited table, contains metadata for all input neighbors, but with EFI and repnode ID for the associated genes of interest added.
- `20210101_repnodeTrim_genE_repnodeGenes.txt` File. Tab-delimited table, contains metadata (including EFI IDs) for repnodes only.
- `20210101_repnodeTrim_genE_repnodeNeighbors.txt` File. Tab-delimited table, contains metadata for neighbors of repnodes only.
- `20210101_repnodeTrim_genE_repnodeGeneSeqs.fa` File. FASTA-formatted file with amino acid sequences for repnodes only.
- `20210101_repnodeTrim_genE_repnodeNeighborSeqs.fa` File. FASTA-formatted file with amino acid sequences for neighbors of repnodes only.
- `repnodeTrimOut` List. Contains `repnodeTrimOut$repGenesTrimmed` (data frame containing the metadata for the repnodes) and `repnodeTrimOut$repNeighborsTrimmed` (data frame containing the metadata for the neighbors of repnodes)

## Development
### Planned additions
- A vignette!  
- Import from UniProt and GFF/GFF-3 formatted files.  As with `gbToIMG` this will likely be a separate function that can replace `generateNeighbors`, and it will likely suffer from the same data heterogeneity that that tool does.  
- Automatic generation of genome neighborhood diagrams for specific clusters (or for random representatives of a cluster) in `prettyClusterDiagrams`.
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
