# prettyClusters 0.2.4

# # prettyClusters 0.2.4;
Minipatches:
* incorpIprScan now handles some updated InterProScan formatting (hello NCBIfam)
* Unhelpful blastdb warning silenced in neighborHypothetical

# prettyClusters 0.2.3
A quartet of small patches:
* Dealing with the deprecation of a bit of tidyselect syntax
* Handling an annoying scientific notation issue for gene/scaffold/genome IDs beyond a certain length (currently only an issue when importing files via gbToIMG)
* Cleaning up handling of complex situations where multiple genomes with multiple gene clusters share the same genome name
* Updating Pfam & InterPro family info.

# prettyClusters 0.2.2
* gbToIMG patch: a few small errors fixed that'll improve handling of moderately mangled Genbank files, keep IDs from clashing with IMG IDs, and otherwise improve integration of Genbank data.

# prettyClusters 0.2.1.0
* Patch to address more stringent handling of gene direction in the recent CRAN release of gggenes.  Got behind on my package updates...

# prettyClusters 0.2.0.0
* Stable enough that I feel comfortable telling other people to try it!  (But still in development.)

## Major changes
* Version number increment reflecting the accumulation of small patches and updated documentation (upping my package game)
* Expansion of the custom protein family identification tools in prepNeighbors

## Minor improvements
* Added a `NEWS.md` file to track changes to the package.
* Updated documentation in the main functions and in the wiki
