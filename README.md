# prettyClusters
A set of tools analyze and make non-hideous figure-friendly diagrams of bacterial gene clusters.

## Important note
This is very much a work in progress, and I'm a biochemist doing terrible things to code.  There will be bugs.  I'll do what I can to address them; if you've come up with a fix, I'm happy to try to incorporate it!

## Why?
I spend a lot of time looking at or making bacterial gene clusters.  I wanted some sort of tool that could:
- Export diagrams of gene clusters as a vector (not bitmap!) file suitable for figure layout.
- Export diagrams with multiple gene clusters in the same relative scale.
- Accept data from the JGI's [IMG database](https://img.jgi.doe.gov/index.html), which has many genomes, metagenomes, and so on that are absent from NCBI's NR database (and UniProt), or that are present but poorly annotated in databases like NCBI/WGS.  
- Handle hypothetical and predicted proteins helpfully (i.e. by identifying new groups of hypothetical proteins that frequently appear in the genomic neighborhoods of interset).
- Integrate into workflows using sequence similarity networks generated via the [EFI-EST toolset](https://efi.igb.illinois.edu/efi-est/).
- Interrogate similarity of genomic neighborhoods without relying on the sequence similarity of genes of interest - I did not want to have to make the assumption that sequence similarity and gene cluster similarity necessarily track, since that's not always a sound assumption.
I haven't encountered anything that quite handles all of those things, so...

## Components of the `prettyClusters` tool set
- `neighborGen`
- `neighborPrep`
- `repnodeTrim`
- `neighboranalysis`
- `prettyClusters`
- `additional`

### `neighborGen`
This tool takes advantage of the fully numeric and contiguous nature of gene_oids in the IMG database.  Given an IMG metadata table for a set of genes of interest (identified via BLAST, protein family-based filtering, or other methods), this tool generates IDs for genes that ought to be in the same neighborhood. These gene lists can be used to download the data-rich IMG metadata files for all all of those genes from the IMG database.
#### Use of `neighborGen`
```
neighborGen(imgGenes = "geneFile.txt", neighborNumber = 10, includeGene = TRUE, geneName = "genE") 
```
Required:
- `imgGenes`
			Filename:  IMG-formated metadata file for genes of interest. 
- `neighborNumber`
			Number:  The number of neighbors upstream and downstream of the gene of interest to be analyzed
- `includeGene`
      T/F: Include the genes of interest in the lists of "neighbors."  Currently this is suggested, since it makes diagram generation later easier.
- `geneName`      
      Text: Include the gene name in autogenerated filenames.
      
### `neighborPrep`
- `hypoProtein`
This tool identifies sets of hypothetical proteins that are overrepresented in the genomic neighborhood of the genes of interest.  Beyond R packages, this tool uses the NCBI blast+ package and MAFFT, which need to be pre-installed.  [Tidygraph](https://github.com/thomasp85/tidygraph) (with [ggraph](https://github.com/thomasp85/ggraph)) and/or [pvclust](https://github.com/shimo-lab/pvclust) are used identify and visualize possible cutoff points.
- `neighborTrim`
This tool identifies genomic neighborhoods where contig ends result in small neighborhoods.  These can - if desired - be ignored in further analyses.
#### Use of `neighborPrep`
```
neighborPrep(imgGenesTemp = "geneFile.txt", imgNeighborsTemp = "neighborFile.txt", geneSeqs = "genes.fa", neighborSeqs = "neighbors.fa", geneName = "genE", efiRepnodes = FALSE, hypoAnalysis = TRUE, neighborNum = 10, neighborThreshold = 0.05, trimShortClusters = TRUE,  alphaVal = .99, imgNeighborsContext = imgNeighborsContext, bootStrap = bootStrap, sysTerm = sysTerm, numThreads = numThreads, clustMethod = clustMethod, pidCutoff = pidCutoff, ...) 
```
Required:
- `imgGenesTemp`
			Filename:  IMG-formated metadata file for genes of interest. 
- `imgNeighborsTemp`
      Filename:  IMG-formated metadata file for neighbors of genes of interest. 
- `geneSeqs`
      Filename:  FASTA-formatted sequence file for genes of interest. 
- `neighborSeqsSeqs`
      Filename:  FASTA-formatted sequence file for neighbors of genes of interest.     
- `geneName`
			Text: the name of your gene family of interest - used for filenames and figure labels.
- `efiRepnodes`
			T/F: just for labeling filenames.  Generally apt to be false at this stage.
- `neighborNum`
			Num: the number of neighboring genes on either side of the gene of interest that will be investigated.
- `efiRepnodes`
			T/F: just for labeling filenames.  Generally apt to be false at this stage.
- `hypoAnalysis`
			T/F: signals whether or not subsets of similar hypothetical proteins will be identified.
- `neighborThreshold`
			Number: 0-1.00.  Used to determine how widely the net for subsets of hypothetical proteins is cast.
- `sysTerm`
      Text: system terminal types (current options: "wsl" for the Linux subsystem on Windows and "nix" for generic Linux, MacOS, etc.)  Used to run blast and mafft commands.
- `numThreads`
			Num: Number of threads to use during analysis.
- `clustMethod`
			Text: the name of the package to be used to identify clusters of similar hypothetical proteins. [Tidygraph](https://github.com/thomasp85/tidygraph)  and pvclust are the current options, with the former suggested.
- `pidCutoff`
			Num: 0-100. %ID cutoff to be used for edge generation when analyzing hypothetical proteins.
- `alphaVal`
			Num: 0-1.  Alpha value cutoff for analyzing pvClust data.
- `bootStrap`
			Num: 0+  Bootstrap number for pvClust.  Can be computationally intensive, 100 is a good starting point for small datasets (e.g. 100 genes of interest, 5 genes on either size), and 10 for larger (500 genes of interest or more, 10+ genes in either direction.)
- `trimShortClusters`
			T/F: Toggles whether or not to remove gene clusters where a truncated contig interrupts the cluter and a sub-sized cluster remains.

### `repnodeTrim`
Generally, at this point I submit the trimmed protein sequences for my genes of interest to the EFI-EST server via option C (a user-uploaded fasta file).  The EFI-EST toolsets are designed to work with UniProt data, and so genes are assigned faux-UniProt IDs, which can complicate tying the output back to analyses that use the IMG gene_oid values.  This tool connects the two IDs.  Additionally, during the EFI-EST SSN setup, sequences below/above length cutoffs are often trimmed, and a %ID cutoff is established, above which a representative node is chosen for SSN visualization. These "repnodes" are also useful for avoiding over-weighting of a network towards data from very highly sequenced species and genera (e.g. pathogens).  Thus, it's often helpful to go ahead using only the representative nodes.  This function provides trimmed versionf of gene and neighbor metadata and sequence for this purpose. 
#### Use of `repnodeTrim`
```
repnodeTrim(imgGenesTrimmed = "genesTrimmed.txt", imgNeighborDataTrimmed="neighborsTrimmed.txt", imgGeneSeqsTrimmed = "geneSeqsTrimmed.fa", imgNeighborSeqsTrimmed = "neighborSeqsTrimmed.fa", geneName = "genE", efiFullMetadata = "efiMetadataFull.txt", efiFinalMetadata = "efiMetadataR95.txt") 
```
### `neighborAnalysis`
Separate from imgNeighborPrep because most of the analyses make sense to do _after_ submitting data to EFI for SSN generation.
- `neighborCatalog`
  This identifies what protein families are present in all genomic neighborhoods, and which are above the neighborThreshold cutoff and will be used in further analyses.
- `neighborHere`
  This identifies which of those families are present in the genomic neighborhoods of genes of interest. Currently a binary present/absent value; it's not clear that weighting 
- `neighborMatrix`
  This generates a horrible matrix based on that binary data.
- `neighborCluster`
  This clusters genomic neighborhoods on the basis of their similarity, using that data.  Automatic assignment to clusters is an option, using tidygraph (currently recommended) or pvClust.  The exported metadata can - among other things - be re-imported to Cytoscape while viewing the EFI-derived SSN.
#### Use of `neighborAnalysis`
```
neighborAnalysis(imgGenesTrimmed = "repnodeGenes.txt", imgNeighborsTrimmed = "repnodeNeighbors.txt", efiRepnodes = TRUE, neighborThreshold = 0.05, geneName = "genE", clusterOrder = TRUE, alphaVal = alphaVal, autoClust = TRUE, clustMethod = "tidygraph",...)	
```
Required:
- `imgGenesTrimmed`
			Filename:  IMG-formated metadata file for genes of interest. Can also be the name of the equivalent object from the previous suite member (i.e. neighborPrep, if no EFI repnodes are used).
- `imgNeighborsTrimmed`
      Filename:  IMG-formated metadata file for neighbors of genes of interest. Can also be the name of the equivalent object from the previous suite member (i.e. neighborPrep, if no EFI repnodes are used).
- `geneName`
			Text: the name of your gene family of interest - used for filenames and figure labels.
- `efiRepnodes`
			T/F: just for labeling filenames.
- `neighborThreshold`
      Num: 0-100.  families of neighbors need to show up in at least this percent of genome neighborhoods to be highlighted in further analyses
- `autoClust`
      T/F: automatically assign cluster numbers for subsets of gene clusters, based on their similarity.  If used, requires a choice between tidygraph (currently recommended) and pvclust.  If not chosen, sets of genome neighborhoods will still be clustered according to their similarity via hierarchical clustering, as will sets of protein families that co-occur in similar sets of genomic neighborhoods. However, cutoffs for cluster identity will need to be drawn by the user, on their own.
Optional:
- `clustMethod`
			Text: required if autoClust is TRUE.  "pvclust" or "tidygraph" are the two options I currently have working.  I'd suggest [tidygraph](https://github.com/thomasp85/tidygraph) (paired with [ggraph](https://github.com/thomasp85/ggraph)) based on the datasets I've looked at so far, but [pvclust](https://github.com/shimo-lab/pvclust)  is still an option.
- `alphaVal`
			Num: 0.00-1.00 alpha value cutoff for use with [pvclust](https://github.com/shimo-lab/pvclust) .
- `pidCutoff`
			Num: 0-100.  Percent ID cutoff for use in tidygraph cluster analysis.  (Edges representing under this percent identity will be deleted.)  Somewhere in the neighborhood of 35% ID is a good starting point.

### prettyClusters
This will make what are hopefully visually passable gene cluster diagrams, exported in vector and bitmap formats. Color schemes are  generated from chosen palettes.  Gene types to be highlighted are specified by the user, either manually or via auto-assignment based on membership in specific protein families (Pfam, TIGRfam, and IMG Term - InterPro to be added as soon as IMG adds it to their metadata export).  This is designed to be run with other components of this suite, but if not run in tandem with neighborPrep and neighborAnalysis, additional QC components will confirm that genes are co-localized to the same scaffold and so on.  Similarly, manual assignment of cluster numbers for visualization is also possible. Note that visualization depends heavily on [gggenes](https://github.com/wilkox/gggenes), with default settings biased towards use cases that involve visualization of many gene clusters. 
#### Use of `prettyClusters`
```
prettyClusters(imgGenesTrimmed = "geneFile.txt", imgNeighborsTrimmed = "neighborFile.txt", geneName = "genE", neighborNumber = 10, checkScaffold = FALSE, efiRepnodes = TRUE, markClusters = TRUE, alignToCore=TRUE, colorType = "nord", annotateGenes = annotateGenes, labelGenes = labelGenes, paletteInput = paletteInput, showScaffold = FALSE, checkScaffold = FALSE, efiRepnodes=TRUE, neighborNumber = 10)
```
- `imgGenesTrimmed`
			Filename:  IMG-formated metadata file for gene. Can also be the name of the equivalent object from the previous suite member.
- `imgNeighborsTrimmed`
			Filename:  IMG-formated metadata file for neighboring genes.  Can also be the name of the equivalent object from the previous suite member.
- `geneName`
			Text: the name of your gene family of interest - used for filenames and figure labels.
- `neighborNumber`
			Num: the number of neighboring genes that will be shown on either side.
- `checkScaffold`
			T/F: if T, assumes that neighboring genes may not have been checked for presence on the same scaffold (recommended if you've skipped imgNeighborsPrep)
- `efiRepnodes`
			T/F: just for labeling filenames.
- `markClusters`
			T/F: note that this assumes you've run imgNeighborhoodAnalysis and have a clustNum column in the neighbor metadata file.  you can alternately manually add clustNum values to imgNeighborsTrimmed.
- `alignToCore`
			T/F: center the gene family of interest in the figure (generally recommended)
- `colorType`
			Text:  what color library (nord, etc.) to use.  Options include [fishualize](https://github.com/nschiett/fishualize), [ghibli](https://github.com/ewenme/ghibli), [lisa](https://github.com/tyluRp/lisa), [nord](https://github.com/cran/nord), [rtist](https://github.com/tomasokal/rtist), [scico](https://github.com/thomasp85/scico), [viridis](https://github.com/sjmgarnier/viridis), and [wesanderson](https://github.com/karthik/wesanderson) because finding color sets that work well for a given gene number and layout can be hard, and I'd rather like this to be able to live up to its name.
- `paletteInput`
			Text: what palette within that color library to choose (e.g. aurora within nord).
- `annotateGenes`
			T/F: uses family-based rules to choose what genes to color. if true, will require the geneFormat file.  recommended so that you don't end up with a million different colors for spuriously annotated genes.  (Currently somewhat slow, though.)
- `geneFormat`
			Filename: text file that provides key for annotation.
- `labelGenes`
			T/F: label genes in figure.  Not recommended for larger figures.
- `showScaffold`
			T/F: include the scaffold ID in labels.  Not recommended as a default since it makes the labels really long.

### Accessory functions
#### `numbers2words`
Derived from the function found [here](http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html).
#### `gffToIMG`
A function to make faux-IMG-formatted metadata tables for genes of interest and their neighbors from a set of input .gff files and a list of genes of interest.  Not yet finished.
#### `uniProtToIMG`
A function to make faux-IMG-formatted metadata tables for genes of interest and their neighbors from a set of input UniProt metadata tables.  Not yet finished.

## Known issues
Partly just to keep track of these for myself!
### Input from UniProt or GFF/GFF3-formatted files does not work
Still working on the best way to get data from these formats.  GFF3 and GenBank files are much more variable in content, complicating some of the analyses in this toolset annoyingly.  UniProt is more consistent but is also protein- rather than nucleotide-based, and some of the analyses on this toolset are nucleotide focused, so I need to figure out how to handle the getting-neighbors bits.
### No option for single scale bar in prettyClusters
Working on it. Probably via generating a final fake "gene cluster" with single-NT "genes" every kb or something?
### No option to indicate the end of contigs in prettyClusters
Considering doing something similar for this, but it's a lower priority.
### Forward and reverse-facing genes are on the same line.
I personally find this aesthetically clearer and want to implement it, but will need to probably do a bunch more digging into  [gggenes](https://github.com/wilkox/gggenes) and [ggplot2](https://github.com/tidyverse/ggplot2) to figure out if/how I can make it happen.
### Clusters generated in hypoProteins or neighborClusters are not totally accurate
As with all similar analyses: this is easily affected by missing data in input files, and there's often no right answer anyway.  But 
### Auto-annotation in prettyClusters misses genes
This currently relies on existing annotation regarding which families a given protein has been predicted to belong to.  It currently fails at things that match only broader superfamilies, or that lack membership to characterized protein families.  (And, for that matter, older or poorly annotated genes.)  Manual annotation may be necessary for this!  I may eventually add a bit where one uses a custom HMM (or set of HMMs) to add better annotation for these.
### Distance between genes of interest and their neighbors is not taken into account
I have done a few initial analyses using weighted distances rather than binary present/absent values; it is not clear they've added much more info than the binary value-based analyses, and they're more complicated to run.  Could revisit?
