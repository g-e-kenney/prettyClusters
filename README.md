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
- [`generateNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-generateNeighbors)
- [`prepNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-prepNeighbors)
- [`analyzeNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-analyzeNeighbors)
- [`prettyClusterDiagrams`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-prettyClusterDiagrams)

### Accessory components
- [`gbToIMG`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-gbToIMG)
- [`incorpIprScan`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-incorpIprScan)
- [`repnodeTrim`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-repnodeTrim)

### Using `prettyClusters`
- A [basic installation guide](https://github.com/g-e-kenney/prettyClusters/wiki/Installation-guide) is the best starting point.
- I've got a simple walkthrough for a [standard run](https://github.com/g-e-kenney/prettyClusters/wiki/a-standard-prettyClusters-run) using the `prettyClusters` toolset.
- This is very much a work in progress, and doing this is using `prettyClusters` in difficult mode, but I've got a [rough workflow](https://github.com/g-e-kenney/prettyClusters/wiki/Running-prettyClusters-with-GenBank-files) for this as well.

### `prettyClusters` output
Illustrating the output of some of the components (or, in the case of the cluster diagrams themselves, just under 10% of the output):
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
