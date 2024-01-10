
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Familias 2.6.1

Familias has existed since 2000 as a Windows program and can be
downloaded from <https://familias.no>. The program calculates
likelihoods and probabilities used to infer family relationships based
on autosomal marker data. The core code of the Windows program
programmed in C++ has been available since 2012 in R. The mentioned
website also contains references describing the implementation and
validation.

No further development of this package is planned. For forensic pedigree
analyses and visualisations we rather recommend the `pedsuite` packages,
see <https://magnusdv.github.io/pedsuite/>. In particular, the package
`pedFamilias` facilitates conversion of .fam files into the `pedsuite`
format.

## Installation

Install from GitHub by as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install Familias from GitHub
devtools::install_github("thoree/Familias")
```
