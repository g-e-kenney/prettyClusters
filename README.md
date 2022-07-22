# prettyClusters
A set of tools analyze and make non-hideous publication-friendly diagrams of genomic neighborhoods in microbial genomes.

## Table of Contents
- Important [introductory note](https://github.com/g-e-kenney/prettyClusters/#important-note)
- [Why](https://github.com/g-e-kenney/prettyClusters/#why) make `prettyClusters`?
- [Components](https://github.com/g-e-kenney/prettyClusters/#the-prettyclusters-toolset) of `prettyClusters`, and how you use them
- [Development](https://github.com/g-e-kenney/prettyClusters/#development) information, including updates (most recent: 20220719) and a to-do list

## Important note
This is very much a work in progress, and I'm a biochemist doing terrible things to code: there will be bugs. I'll do what I can to address them, and if you've come up with a fix, I'm happy to try to incorporate it!  Also, I want to emphasize that this package is utterly reliant on some excellent R packages (particularly [gggenes](https://github.com/wilkox/gggenes) and [tidygraph](https://github.com/thomasp85/tidygraph)).

## Why?
I spend a lot of time working with bacterial gene clusters.  I wanted some sort of tool that could:
- Export diagrams of gene clusters as vector files for figures.
- Export diagrams with multiple gene clusters in the same relative scale.
- Import gene metadata from the JGI's [IMG database](https://img.jgi.doe.gov/index.html), which has many genomes, metagenomes, and so on that are absent from NCBI's NR database (and UniProt), or that are present but poorly annotated in databases like NCBI/WGS.  Also, the metadata is more detailed and less obnoxious than GenBank files.
- Handle hypothetical and predicted proteins helpfully (i.e. by identifying new groups of hypothetical proteins that frequently appear in the genomic neighborhoods of interest).
- Integrate into workflows using sequence similarity networks generated via the [EFI-EST toolset](https://efi.igb.illinois.edu/efi-est/).
- Interrogate similarity of genomic neighborhoods without relying on the sequence similarity of genes of interest - I did not want to have to make the assumption that sequence similarity and gene cluster similarity necessarily track, since that's not always a sound assumption.

I haven't encountered anything that quite handles all of those things in one go, so...

## The `prettyClusters` toolset
The Wiki entries contain a more detailed description of the use of specific functions.
### The core toolset
- [`generateNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-generateNeighbors)
- [`prepNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-prepNeighbors)
- [`analyzeNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-analyzeNeighbors)
- [`prettyClusterDiagrams`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-prettyClusterDiagrams)

### Accessory components
- [`gbToIMG`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-gbToIMG)
- [`identifySubgroups`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-identifySubgroups)
- [`incorpIprScan`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-incorpIprScan)
- [`repnodeTrim`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-repnodeTrim)
- [`trimFasta`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-trimFasta)

### Using `prettyClusters`
- A [basic installation guide](https://github.com/g-e-kenney/prettyClusters/wiki/Installation-guide) is the best starting point.
- I've got a simple walkthrough for a [standard run](https://github.com/g-e-kenney/prettyClusters/wiki/a-standard-prettyClusters-run) using the `prettyClusters` toolset.
- This is very much a work in progress, and doing this is using `prettyClusters` in difficult mode, but I've got a [rough workflow](https://github.com/g-e-kenney/prettyClusters/wiki/Preparing-data-from-non-IMG-sources-for-prettyClusters) for getting data from non-IMG sources, and [another workflow](https://github.com/g-e-kenney/prettyClusters/wiki/Running-prettyClusters-with-GenBank-files) for standardizing it for use in `prettyClusters`.

### `prettyClusters` output
Illustrating the output of some of the components (or, in the case of the cluster diagrams themselves, just under 10% of the output):
<img src="https://github.com/g-e-kenney/prettyClusters/raw/master/data/20210106_pretty-cluster-general-01.png " width="100%" alt="genome neighborhood diagram output">
Notably, sequence similarity and genome neighborhood similarity are not always tightly coupled.  The analyses in `prettyClusters` make it possible to investigate a protein family along both axes.

## Development
### Recent updates
- `prettyClusterDiagrams` got a bunch of small bugfixes (20220719) that had to do with manually specifying colors or with people visualizing a very small number of gene clusters.
- `prepNeighbors` and its subfunction `neighborHypothetical` as well as `identifySubgroups` got some updates (20220510) that make identification of peptides (<150 aa) work better; in `prepNeighbors` this means that peptide families can be flagged independently of annotation (since they fail to get annotated at relatively high rates, particularly for things like RiPP precursor peptides).
- `prettyClusterDiagrams` got an update (20220510) that finally adds in a single scale bar for the entire diagram.
- Recent changes in IMG export formats are now (20220510) taken into account (a few modules got tweaks for this).
- Similarly affecting several modules, manual addition of genes is now (20220510) supported (i.e. if a RiPP precursor peptide between genes 2000000 and 2000001 is unannotated, you can manually add it as 2000000.5 in the neighbor metadata, neighbor fasta sequence list, and gene/neighbor context file, and everything downstream'll work.)  This obviously scales poorly in large datasets with systemic annotation issues, but it can be helpful if you've got just a handful of missing sequences. 
- `prettyClusterDiagrams` got several updates (20220108) that improve handling for less common datasets - very small datasets, datasets where there are no hypothetical or unannotated proteins, and (if you are making diagrams for subgroups of clusters) figure generation for very small groups.
- `prepNeighbors` and its subcomponents (along with `generateNeighbors`) got some additional updates (20220106) that correct data import errors when certain characters are present in gene metadata.
- `prepNeighbors` got some updates (20211214) that correct handling of smaller gene neighborhoods.
- `gbToIMG` got some big fixes (20211129) that improve stability when it encounters problems (a GenBank file with no gene of interest, an AntiSmash-formatted GenBank file, a GenBank file with no annotations, etc.) and that improve output annotations, including both the metadata format and the content (particularly organism and scaffold info.)
### Up next
- Some updates to the suggested workflow when starting with user-annotated genomes (improved scripts and annotation recommendations.)
- Also in `prepNeighbors`, generation of HMMs for hypothetical protein families identified in `prepNeighbors` and `identifySubgroups` (and with it the ability to turn on and off MSA and HMM generation in both tools.)
- Probably `averageCluster` and `compareCluster` as described below.
### Planned additions
#### Writeup as a paper?
- I may try to write this up as a paper (despite the shame of publicizing my hideous code).  This toolset has been useful in enough projects that it might be nice to have something legit to cite?
#### New modules
- Maybe another sort of per-cluster diagram - one with a "typical" genome neighborhood, illustrating order and abundance of neighbors?  Still working on how to automate this well.  Will probably be `averageCluster`.
- Similarly, possibly a third sort of per-cluster diagram with highlighting of %ID between homologs - more along the lines of [clinker](https://github.com/gamcil/clinker), but without the GenBank input.  Or like [gggenomes](https://github.com/thackl/gggenomes), but protein-only.  (Possibly employing one of those tools if I can figure out a way to easily do so!)  Likely to be `compareCluster`.
- For use when working with really large families, `repnodePreTrim`- using the [EFI-EST toolset](https://efi.igb.illinois.edu/efi-est/) before even generating neighborhoods, as a way of getting a more manageable dataset.  
- A vignette?  Might be redundant given [the wiki](https://github.com/g-e-kenney/prettyClusters/wiki).  
#### Smaller tweaks
- User-supplied HMMs for annotation of predefined custom protein (sub)families as a standalone subfunction.
- Options to let the user specify distance and clustering methods in `prepNeighbors` and `analyzeNeighbors`.
- Additional tools or instructions for import from UniProt metadata and from GFF/GFF-3 formatted files - haven't decided whether to keep guiding people towards GenBank as an input format (making `gbToIMG` the entry for non-IMG data), or whether to develop dedicated tools for one or two other common formats.  If so, this will likely be a separate function that can also replace `generateNeighbors`, and it will likely suffer from the same data heterogeneity (and lack of gene family annotations) that that `gbToIMG` does.  Supplementation with `incorpIprScan` is likely still going to be advisable.
### Known issues
- There may be complications for people on newer M1 Macs.  I've highlighted a few in the [installation guide](https://github.com/g-e-kenney/prettyClusters/wiki/Installation-guide) as people have reported them, but I don't have a Mac new enough to test this myself or to validate workarounds.  I don't anticipate quite the same level of problems for Windows 10 vs. Windows 11, but I also don't have a PC new enough to test any Windows 11-specific issues.
- Auto-annotation in `prettyClusterDiagrams` may miss genes if their initial family categorization (or ORF-finding, particularly for small genes like RiPP precursors) was poor!  This is doubly the case for GenBank-derived files (use of `incorpIprScan` can improve annotation, but not ORF-finding).  If this is a big problem for your genome neighborhoods, you can consider manually updating problem genes (adding metadata lines and amino acid sequences, which I know is terrible) or re-annotating your sequences ([reannotated data](https://github.com/g-e-kenney/prettyClusters/wiki/Preparing-data-from-non-IMG-sources-for-prettyClusters#re-annotating-nucleic-acid-sequences) can be reintegrated into a `prettyClusters` workflow.)
- Generation of hypothetical protein and genome neighborhood clusters is approximate and sensitive to user-supplied cutoffs, to distance/clustering methods, to overrepresentation of closely related gene clusters, and to the conservation of genes outside of the gene cluster limits. There are limited ways around these problems, and they come with their own compromises.  (Overrepresentation at least can be dealt with using representative sequences chosen via [EFI-EST](https://efi.igb.illinois.edu/efi-est/) (repnodes) or [CD-HIT](https://github.com/weizhongli/cdhit).)
- Forward- and reverse-facing genes are on the same vertical level in `prettyClusterDiagrams`. I personally find it visually clearer to have forward genes above the line and reverse genes below, but will need to probably do a bunch more digging into [gggenes](https://github.com/wilkox/gggenes) and [ggplot2](https://github.com/tidyverse/ggplot2) to figure out if/how I can make it happen.
- Distance between genes of interest and their neighbors is not taken into account.  I have done a few initial analyses using weighted distances rather than binary present/absent values; it is not clear they've added much more info than the binary value-based analyses, and they're more complicated to run.  Could revisit?
- Use of non-WSL `mafft` and `blast` installs on Windows isn't built-in at the moment; it's a low priority, given how easy WSL is to get up and running.
