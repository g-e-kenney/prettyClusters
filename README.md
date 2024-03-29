# prettyClusters
A set of tools analyze and make non-hideous publication-friendly diagrams of genomic neighborhoods in microbial genomes.

## Table of Contents
- [Introductory note](https://github.com/g-e-kenney/prettyClusters/#important-note)
- [Why](https://github.com/g-e-kenney/prettyClusters/#why) make or use `prettyClusters`?
- [Components](https://github.com/g-e-kenney/prettyClusters/#the-prettyclusters-toolset) of `prettyClusters`, and how you use them
- [Wiki](https://github.com/g-e-kenney/prettyClusters/wiki) for `prettyClusters`, with a lot of detail
- [Development](https://github.com/g-e-kenney/prettyClusters/#development) information, including updates (most recent: 20230805) and a to-do list

## Important note
This is very much a work in progress, and I'm a chemical biologist who has wandered over to the dark side: there will be bugs. I'll do what I can to address them, and if you've come up with a fix, I'm happy to try to incorporate it!  Also, I want to emphasize that this package is utterly reliant on some excellent R packages (particularly [gggenes](https://github.com/wilkox/gggenes) and [tidygraph](https://github.com/thomasp85/tidygraph)).

## Why?
I spend a lot of time staring at microbial gene clusters, trying to tie new or understudied enzyme (super)families to pathways.  On a high level, I wanted to be able to start with a gene (super)family of interest and:

- Examine similarity of genomic neighborhoods without relying on the sequence similarity of proteins (a point of differentiation between this tool and [EFI-GNT](https://efi.igb.illinois.edu/efi-gnt/) - I did not want to have to make the assumption that sequence similarity and gene cluster similarity necessarily track, since that's not a sound assumption for all protein families involved in secondary metabolite production.)
- Interrogate similarity of genomic neighborhoods without relying on similarity to known gene clusters or gene cluster components (a point of differentiation between this tool and [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE), since gene clusters for novel types of natural products (and for pathways that aren't involved in secondary metabolite production) are missed by pipelines that rely on [antiSMASH](https://antismash.secondarymetabolites.org)-like approaches.)
- Interrogate similarity of genomic neighborhoods without relying solely on Pfam to define gene cluster content (and without stringent limits on data visualization) - [IMG-ABC](https://img.jgi.doe.gov/cgi-bin/abc/main.cgi) has recently started offering prettyClusters-like heatmaps and neighborhood networks, but visualization is limited to a fairly small number of BGCs and you can only search for clusters using Pfams as hooks (which, as a connoisseur of hypothetical proteins and poorly annotated peptides, I find limiting). [CAGECAT](https://cagecat.bioinformatics.nl/) - which also popped up after I started hacking at this - is also nifty and uses broadly similar approaches, but it has not-dissimilar limitations.
- Explore the diversity _within_ a class of gene clusters, anchored around a gene type of interest.  Tools like [deepBGC](https://github.com/Merck/deepbgc), [SanntiS](https://github.com/Finn-Lab/SanntiS), and [GECCO](https://github.com/zellerlab/GECCO) are intriguing ways to try to identify classes of gene clusters that are missed by traditional tools, but I wanted something better suited to examining patterns within a (new) family. 

On a more granular level, I wanted to be able to do a few additional things:
- Handle hypothetical and predicted proteins helpfully (i.e. by identifying new groups of hypothetical or sporadically annotated proteins that frequently appear in the genomic neighborhoods of interest).  Hypothetical proteins remain a major stumbling point for all annotation-dependent approaches to gene cluster exploration, and I didn't want to lose useful information because of it.
- Import gene metadata from the JGI's [IMG database](https://img.jgi.doe.gov/index.html), which has many genomes, metagenomes, and so on that are absent from NCBI's NR database (and UniProt), or that are present but poorly annotated in databases like NCBI/WGS. (I also wanted to be able to handle InterProScan-supplemented GenBank files to address that poor annotation.) Unfortunately few tools play well with IMG data by default - unfortunate, since the annotation metadata is of relatively reliable quality and is available in formats that aren't GFF3 or GenBank.
- Handle subfamilies helpfully (i.e. annotations for subtypes of huge families can be added - generally automatically - and taken into account when visualizing datasets)
- Integrate into workflows using sequence similarity networks generated via the [EFI-EST toolset](https://efi.igb.illinois.edu/efi-est/).

Put it all together, and it's still a niche that doesn't seem to be quite filled by other prettier or better-written tools... so here we are.

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
- Also check out the [troubleshooting](https://github.com/g-e-kenney/prettyClusters/wiki/Troubleshooting-common-issues) list for cryptic-sounding errors that I've seen pop up enough times to be worth noting.

### `prettyClusters` output
Illustrating the output of some of the components (or, in the case of the cluster diagrams themselves, just under 10% of the output):
<img src="https://github.com/g-e-kenney/prettyClusters/raw/master/inst/20210106_pretty-cluster-general-01.png " width="100%" alt="genome neighborhood diagram output">
Notably, sequence similarity and genome neighborhood similarity are not always tightly coupled.  The analyses in `prettyClusters` make it possible to investigate a protein family along both axes.

### Citing `prettyClusters`
Until I get a proper paper out, you can use `citation()` to generate a basic citation for the version of `prettyclusters` that you are currently using.
```R
> citation("prettyClusters")

In order to cite the package ‘prettyClusters’ in publications, please use:

  Kenney GE (2023). _prettyClusters: Exploring and Classifying Genomic
  Neighborhoods Using IMG-Like Data_. R package version 0.2.3.

A BibTeX-formatted citation for LaTeX users is:
  @Manual{,
    title = {prettyClusters: Exploring and Classifying Genomic Neighborhoods Using IMG-Like Data},
    author = {G. E. Kenney},
    year = {2023},
    note = {R package version 0.2.3},
  }
```

## Development
### Current Version
- Version 0.2.3.0  See the [release notes](https://github.com/g-e-kenney/prettyClusters/blob/master/NEWS.md).
### Recent updates
- Several additional patches (20230804), most relevantly an update to Pfam and InterPro family listings, plus some tweaks that'll help with handling messy and possibly duplicated genomes from mixed sources in `prettyClusters`.
- Quick patch (20230802): `gbToIMG` fixes to handle certain kinds of iffy .gbk files.
- Quick patch to 0.2.1 (20230714), because I got around to updating packages and it turns out the spring CRAN release for [gggenes](https://github.com/wilkox/gggenes) is pickier about designating gene directions. 
- Version increment to 0.2.0.0 (20230712), reflecting accumulated small changes and hotfixes.
- Small tweak (20230711) to `prepNeighbors` and its subfunction `neighborHypothetical`- you can now assign all proteins below a given size a hypothetical family regardless of annotation. Increase the aa cutoff with caution on large datasets, since this uses an all-by-all blast.  Also added a whole bunch of wiki updates for clarity (including updated default runs & clarification of default vs. required input).
- Minor bugfixes (20230208) to `prepNeighbors` and `incorpIprScan`
- Some additions to the [troubleshooting](https://github.com/g-e-kenney/prettyClusters/wiki/Troubleshooting-common-issues) list (20230109).
- Small `prepNeighbors` bugfix (20221004) for people who are running without trimming truncated gene clusters!
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
- Some quick fixes to handle some JGI gene naming updates
- User-supplied HMMs for annotation of predefined custom protein (sub)families as a standalone subfunction.
- Looking into setting up a GUI with an eye towards online deployment.  The activation energy for getting a new user up and running on R is unfortunately real!
- Some updates to the suggested workflow when starting with user-annotated genomes (improved scripts and annotation recommendations.)
- Also in `prepNeighbors`, generation of HMMs for hypothetical protein families identified in `prepNeighbors` and `identifySubgroups` (and with it the ability to turn on and off MSA and HMM generation in both tools.)
### Longer term plans
#### Write up as a paper? 
- I may try to write this up as a paper (despite the shame of publicizing my hideous code).  This toolset has been useful in enough projects that it might be nice to have something legit to cite?
#### New modules & major tweaks
- A way to handle non-gene things as neighborhood "anchors" - regulator binding sites, riboswitches, etc.?  Annotated elements like tRNAs are more straightforward; user-supplied coordinates and pseudo-genes (analogous to approaches used for unannotated peptides) might be the way to go.
- A way to handle datasets where there is no one gene of interest to "anchor" the cluster (i.e. a cluster has 2+ of a larger set of gene families within a constrained genomic region, or within the genome.)  Basically loosening some of the `analyzeNeighbors` requirements.
- May try to auto-generate a "typical" genome neighborhood, illustrating order and abundance of neighbors in a given gene cluster family?  Still working on how to automate this well.  Will probably be `averageCluster`.
- Similarly, possibly a third sort of per-cluster diagram with highlighting of %ID between homologs - more along the lines of [clinker](https://github.com/gamcil/clinker), but without the GenBank input.  Or like [gggenomes](https://github.com/thackl/gggenomes), but protein-only.  (Possibly employing one of those tools if I can figure out a way to easily do so!)  Likely to be `compareCluster`.
- For use when working with really large families, `repnodePreTrim`- using the [EFI-EST toolset](https://efi.igb.illinois.edu/efi-est/) before even generating neighborhoods, as a way of getting a more manageable dataset.  "If you need more than 128 GB of RAM to open the SSN, you may find this helpful..."
#### Smaller tweaks

- Options to let the user specify distance and clustering methods in `prepNeighbors` and `analyzeNeighbors`. (Perhaps you prefer Jaccard distance, for example?)
- Additional tools or instructions for import from UniProt metadata and from GFF/GFF-3 formatted files - haven't decided whether to keep guiding people towards GenBank as an input format (making `gbToIMG` the entry for non-IMG data), or whether to develop dedicated tools for one or two other common formats that suck less.  If so, this will likely be a separate function that can also replace `generateNeighbors` - though any GFF/GFF3-based option will likely suffer from the same data heterogeneity (and lack of gene family annotations) that that `gbToIMG` does.  Supplementation with `incorpIprScan` is likely still going to be advisable in those situations.
### Known issues
- I am keeping track of [workarounds](https://github.com/g-e-kenney/prettyClusters/wiki/Troubleshooting-common-issues) for some common issues that are (for various reasons) hard to fix.
- There may be installation complications for people on newer Apple Silicon Macs.  I've highlighted a few in the [installation guide](https://github.com/g-e-kenney/prettyClusters/wiki/Installation-guide) as people have reported them, and as I've gotten a chance to try things out myself.  I don't anticipate quite the same level of problems for Windows 10 vs. Windows 11, but I also don't have a PC new enough to test any Windows 11-specific issues.
- Auto-annotation in `prettyClusterDiagrams` may miss genes if their initial family categorization (or ORF-finding, particularly for small genes like RiPP precursors) was poor!  This is doubly the case for GenBank-derived files (use of `incorpIprScan` can improve annotation, but not ORF-finding).  If this is a big problem for your genome neighborhoods, you can consider manually updating problem genes (adding metadata lines and amino acid sequences, which I know is terrible) or re-annotating your sequences ([reannotated data](https://github.com/g-e-kenney/prettyClusters/wiki/Preparing-data-from-non-IMG-sources-for-prettyClusters#re-annotating-nucleic-acid-sequences) can be reintegrated into a `prettyClusters` workflow.)
- Generation of hypothetical protein and genome neighborhood clusters is approximate and sensitive to user-supplied cutoffs, to distance/clustering methods, to overrepresentation of closely related gene clusters, and to the conservation of genes outside of the gene cluster limits. There are limited ways around these problems, and they come with their own compromises.  (Overrepresentation at least can be dealt with using representative sequences chosen via [EFI-EST](https://efi.igb.illinois.edu/efi-est/) (repnodes) or [CD-HIT](https://github.com/weizhongli/cdhit).)
- Forward- and reverse-facing genes are on the same vertical level in `prettyClusterDiagrams`. I personally find it visually clearer to have forward genes above the line and reverse genes below, but will need to probably do a bunch more digging into [gggenes](https://github.com/wilkox/gggenes) and [ggplot2](https://github.com/tidyverse/ggplot2) to figure out if/how I can make it happen.
- Distance between genes of interest and their neighbors is not taken into account.  I have done some analyses using weighted distances rather than binary present/absent values; it is not clear they've added much more info than the binary value-based analyses, and they're more complicated to run.  Could revisit?
- Similarly, number of genes from a given family is not taken into account. Some families are much more specific while others are more like superfamilies, and weighting abundance of a given family by copy number alone can thus be misleading (I've seen this in some BGC families.) But it's something I should probably revisit and benchmark with some bigger datasets.
- Use of non-WSL `mafft` and `blast` installs on Windows isn't built-in at the moment; it's a low priority, given how easy WSL is to get up and running.
