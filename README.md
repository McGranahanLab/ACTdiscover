# ACT-Discover
Aneuploidy in Circulating Tumour DNA (ACT-Discover) is an R package to detect somatic copy number aberrations (SCNAs) in circulating tumour DNA (ctDNA) by leveraging haplotype phasing of paired tumour biopsies or patient-derived xenograft (PDX) samples. ACT-Discover can be used to infer karyotype heterogeneity during tumour evolution. For more details on the method and use of ACT-Discover please read our publication [*ACT-Discover*: identifying karyotype heterogeneity in pancreatic cancer evolution using ctDNA. Huebner, Black et al. *Genome Med* (2023)](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01171-w)

## Installation guide

```
# install.packages("devtools")
library(devtools)
devtools::install_github("McGranahanLab/ACTdiscover")
``` 

## Troubleshooting

This package depends on additional packages being installed. If the installation of *ACT-Discover* fails, try installing the following packages:

```
# install.packages("devtools")
# install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
install.packages("tidyverse")
install.packages("cowplot")
```

## Reference
Huebner, A., Black, J.R.M., Sarno, F. et al. ACT-Discover: identifying karyotype heterogeneity in pancreatic cancer evolution using ctDNA. Genome Med 15, 27 (2023). https://doi.org/10.1186/s13073-023-01171-w

