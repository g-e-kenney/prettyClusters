# prettyClusters 0.3.0
Boring but helpful dependency updates:
* GenBank import (gbToIMG) now uses gggenomes, since genbankr is no longer in recent versions of Bioconductor.  This means some changes to the workflow (most importantly, the requirement that you extract amino acid sequences from the GenBank file in advance, though I'll look into less-annoying ways to do this) - check the wiki!
* This actually also means that antiSMASH-annotated GenBank files are now handled acceptably.


# prettyClusters 0.2.5
* Small issue with NA handling in prepNeighbors and identifySubgroups fixed
* Issue with identifySubgroups metadata output fixed
* Gene coloring in prettyClusterDiagrams better handles annotating a very small number of genes and some other edge cases
* Now analzyeNeighbors defaults to auto-categorizing genome neighborhood types
* Assigning a trial DOI for this one

# prettyClusters 0.2.4
Minipatches:
* incorpIprScan now handles some updated InterProScan formatting (hello NCBIfam)
* Unhelpful blastdb warning that doesn't affect function silenced in neighborHypothetical and identifySubgroups

# prettyClusters 0.2.3
A quartet of small patches:
* Dealing with the deprecation of a bit of tidyselect syntax
* Handling an annoying scientific notation issue for gene/scaffold/genome IDs beyond a certain length (currently only an issue when importing files via gbToIMG)
* Cleaning up handling of complex situations where multiple genomes with multiple gene clusters share the same genome name
* Updating Pfam & InterPro family info.

# prettyClusters 0.2.2
* gbToIMG patch: a few small errors fixed that'll improve handling of moderately mangled Genbank files, keep IDs from clashing with IMG IDs, and otherwise improve integration of Genbank data.

# prettyClusters 0.2.1
* Patch to address more stringent handling of gene direction in the recent CRAN release of gggenes.  Got behind on my package updates...

# prettyClusters 0.2.0
* Stable enough that I feel comfortable telling other people to try it!  (But still in development, obvs.)

## Major changes
* Version number increment reflecting the accumulation of small patches and updated documentation (upping my package game)
* Expansion of the custom protein family identification tools in prepNeighbors

## Minor improvements
* Added a `NEWS.md` file to track changes to the package.
* Updated documentation in the main functions and in the wiki
