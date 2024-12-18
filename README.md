
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![DOI](https://zenodo.org/badge/326540314.svg)](https://zenodo.org/badge/latestdoi/326540314)
[![](https://img.shields.io/badge/devel%20version-0.2.5-blue.svg)](https://github.com/g-e-kenney/prettyClusters)
[![](https://img.shields.io/github/last-commit/g-e-kenney/prettyClusters.svg)](https://github.com/g-e-kenney/prettyClusters/commits/master)
<!-- badges: end -->

# prettyClusters

A set of tools analyze and make publication-friendly diagrams of genome
neighborhoods in microbial genomes.

Particular focuses: investigation and classification of genome
neighborhoods that do not encode well-studied types of natural product
pathways, integration of hypothetical proteins into analyses, and
compatibility with IMG datasets

## Table of Contents

- [Introductory
  note](https://github.com/g-e-kenney/prettyClusters/#important-notes)
- [Components](https://github.com/g-e-kenney/prettyClusters/#the-prettyclusters-toolset)
  of `prettyClusters`, and how you use them
- [Development](https://github.com/g-e-kenney/prettyClusters/#development)
  information, including updates (most recent: 20241216) and a to-do
  list

## Important notes

There are a number of other excellent tools for looking at genome
neighborhoods in various ways, including
[antiSMASH](https://antismash.secondarymetabolites.org),
[EFI-GNT](https://efi.igb.illinois.edu/efi-gnt/),
[BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE) and
[BiG-SCAPE-CORASON](https://bigscape-corason.secondarymetabolites.org/),
[CAGECAT](https://cagecat.bioinformatics.nl/) (a server wrapper for
[cblaster](https://github.com/gamcil/cblaster) and
[clinker](https://github.com/gamcil/clinker)),
[zol/fai](https://github.com/Kalan-Lab/zol), and
[socialGene](https://socialgene.github.io/). However, you may want to
use `prettyClusters` if:

- you’re starting out with a gene/protein of interest and don’t know
  much about the pathways it’s part of yet.
- you are working on pathways that don’t look like the type of
  well-studied pathways found in
  [MiBIG](https://mibig.secondarymetabolites.org/) and detected by
  [antiSMASH](https://antismash.secondarymetabolites.org).
- you are working on pathways that have a lot of hypothetical proteins,
  or that have a lot of proteins with vague/blanket annotations.
- you want to visualize both genome neighborhood content and
  relationships, and potentially to correlate it to sequence similarity
  for your protein of interest.
- you use a lot of datasets from the JGI’s [IMG
  database](https://img.jgi.doe.gov/index.html).

Additionally: this is a work in progress, and I’m a chemical biologist:
there will be bugs and inefficient code, but I do my best to fix them as
they’re identified.

## The `prettyClusters` toolset

The Wiki entries for each function contain a detailed description of
theuse of specific functions.

### The core toolset

- [`generateNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-generateNeighbors).
  Sets up the list of neighboring genes for your gene of interest.
- [`prepNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-prepNeighbors).
  QC, and provisional assignment of hypothetical protein families.
- [`analyzeNeighbors`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-analyzeNeighbors)
  Quantification of common types of neighboring genes, genome
  neighborhood classification.
- [`prettyClusterDiagrams`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-prettyClusterDiagrams)
  Visualization of genome neighborhoods.

### Accessory components

- [`gbToIMG`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-gbToIMG).
  Processing GenBank files for a `prettyClusters` workflow.
- [`incorpIprScan`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-incorpIprScan).
  Incorporate InterProScan annotations (primarily needed when using
  GenBank input.)
- [`repnodeTrim`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-repnodeTrim).
  When integrating EFI-EST SSNs, trim the active dataset to
  representative nodes.
- [`identifySubgroups`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-identifySubgroups).
  Assign provisional subgroup annotations within a large group of
  related proteins.
- [`trimFasta`](https://github.com/g-e-kenney/prettyClusters/wiki/Function:-trimFasta).
  Provide a trimmed multiFASTA file for a subset of proteins.

### Installation

- A [basic installation
  guide](https://github.com/g-e-kenney/prettyClusters/wiki/Installation-guide)
  is the best starting point. Includes cross-platform installation
  instructions for everything needed to make this package work.

### Use

- For users working with IMG data, I’ve got a walkthrough for a
  [standard
  run](https://github.com/g-e-kenney/prettyClusters/wiki/a-standard-prettyClusters-run).
- For users working with data from other sources, I’ve got a [rough
  workflow](https://github.com/g-e-kenney/prettyClusters/wiki/Preparing-data-from-non-IMG-sources-for-prettyClusters)
  for getting data from non-IMG sources, and [another
  workflow](https://github.com/g-e-kenney/prettyClusters/wiki/Running-prettyClusters-with-GenBank-files)
  for standardizing it for use in `prettyClusters`.
- Also check out the
  [troubleshooting](https://github.com/g-e-kenney/prettyClusters/wiki/Troubleshooting-common-issues)
  list for cryptic-sounding errors that I’ve seen pop up enough times to
  be worth noting, as well as some additional limitations.

### Output

Illustrating the output of some of the components (or, in the case of
the cluster diagrams themselves, just under 10% of the output for this
example):
<img src="https://github.com/g-e-kenney/prettyClusters/raw/master/inst/20210106_pretty-cluster-general-01.png " alt="genome neighborhood diagram output" width="100%"/>
Notably, sequence similarity and genome neighborhood similarity are not
always tightly coupled. The analyses in `prettyClusters` make it
possible to investigate a protein family along both axes.

### Citing

Until I get a proper paper out, you can use `citation()` to generate a
basic citation for the version of `prettyClusters` that you are
currently using.

``` r
citation("prettyClusters")
#> To cite package 'prettyClusters' in publications use:
#> 
#>   Kenney G (2024). _prettyClusters: Exploring and Classifying Genomic
#>   Neighborhoods Using IMG-Like Data_. R package version 0.2.5.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {prettyClusters: Exploring and Classifying Genomic Neighborhoods Using IMG-Like Data},
#>     author = {G. E. Kenney},
#>     year = {2024},
#>     note = {R package version 0.2.5},
#>   }
```

Additionally, this package makes use of other tools, including:

[BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/references.html#references)

    Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., Madden T.L.  BMC Bioinformatics (2008) 10:421 

[MAFFT](https://mafft.cbrc.jp/alignment/software/)

    Nakamura T., Yamada K.D., Tomii K., Katoh K.  Bioinformatics (2018) 34:2490–2492 

[HMMER](https://www.hmmer.org) - note that citing the website is
preferred.

    HMMER 3.4 (Aug 2023); http://hmmer.org/
    Eddy S. R. PLOS Comp. Biol. (2011) 7:e1002195

[gggenes](https://github.com/wilkox/gggenes)

``` r
citation("gggenes")
#> To cite package 'gggenes' in publications use:
#> 
#>   Wilkins D (2023). _gggenes: Draw Gene Arrow Maps in 'ggplot2'_. R
#>   package version 0.5.1, <https://CRAN.R-project.org/package=gggenes>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {gggenes: Draw Gene Arrow Maps in 'ggplot2'},
#>     author = {David Wilkins},
#>     year = {2023},
#>     note = {R package version 0.5.1},
#>     url = {https://CRAN.R-project.org/package=gggenes},
#>   }
```

[tidygraph](https://tidygraph.data-imaginist.com/)

``` r
citation("tidygraph")
#> To cite package 'tidygraph' in publications use:
#> 
#>   Pedersen T (2024). _tidygraph: A Tidy API for Graph Manipulation_. R
#>   package version 1.3.1,
#>   <https://CRAN.R-project.org/package=tidygraph>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {tidygraph: A Tidy API for Graph Manipulation},
#>     author = {Thomas Lin Pedersen},
#>     year = {2024},
#>     note = {R package version 1.3.1},
#>     url = {https://CRAN.R-project.org/package=tidygraph},
#>   }
```

## Development

### Current Version

- Version 0.2.5 See the [release
  notes](https://github.com/g-e-kenney/prettyClusters/blob/master/NEWS.md).

### Recent updates

- Improved handling of coloring edge cases in `prettyClusterDiagrams`
- Minor tweaks in peptide handling and NA handling in `prepNeighbors`
  and `identifySubgroups` (and in metadata/graph export for the latter)
- `incorpIprScan` deals gracefully with some InterProScan updates
- An unhelpful and non-fatal blast error is provisionally silenced in
  `prepNeighbors` and `identifySubgroups`

### Longer term plans

#### New modules & major tweaks

- Looking into setting up a shiny and potentially interactive GUI, if
  not full online deployment. The activation energy for getting a new
  user up and running on R is unfortunately real!
- A way to handle non-gene things as neighborhood “anchors” - regulator
  binding sites, riboswitches, etc.? Annotated elements like tRNAs are
  more straightforward; user-supplied coordinates and pseudo-genes
  (analogous to approaches used for unannotated peptides) might be the
  way to go.
- A gene co-occurrence tool to handle datasets where there is no one
  gene of interest to “anchor” the cluster (i.e. a cluster has 2+ of a
  larger set of gene families within a constrained genomic region, or
  within the genome.) Basically loosening some of the `analyzeNeighbors`
  requirements. Will need a modified approach for diagram generation
  without a core GoI. `coGenes` or something similar.
- May try to auto-generate a “typical” genome neighborhood, illustrating
  order and abundance of neighbors in a given gene cluster family? Still
  working on how to automate this well. Would be `averageCluster`.
- Similarly, possibly a third sort of per-cluster diagram with
  highlighting of %ID between homologs - more along the lines of
  [clinker](https://github.com/gamcil/clinker), but without the GenBank
  input. Or like [gggenomes](https://github.com/thackl/gggenomes), but
  protein-only. (Possibly employing one of those tools if I can figure
  out a way to easily do so!) Likely to be `compareCluster`.

#### Smaller tweaks

- For use when working with really large families, tweak `repnodeTrim`
  and add an option for using the [EFI-EST
  toolset](https://efi.igb.illinois.edu/efi-est/) before even generating
  neighborhoods, as a way of getting a more manageable dataset. “If you
  need more than 128 GB of RAM to open the SSN, you may find this
  helpful…”
- Generation of HMMs for hypothetical protein families identified in
  `prepNeighbors` and `identifySubgroups` (and with it the ability to
  turn on and off MSA and HMM generation in both tools.)
- User-supplied HMMs for annotation of predefined custom protein
  (sub)families as a standalone subfunction.
- Options to let the user specify distance and clustering methods in
  `prepNeighbors` and `analyzeNeighbors`. (Perhaps you prefer Jaccard to
  Euclidean distance, for example?)
