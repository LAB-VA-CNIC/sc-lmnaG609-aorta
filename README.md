# Single Cell Analysis of LmnaG609 and WT mice aortas

## Mice:

Mouse studies were conducted with male and female 14-week-old LmnaG609G/G609G mice with ubiquitous progerin expression9 (provided by C. López-Otín, Universidad de Oviedo, Spain) and wild-type Lmna+/+ littermates as controls. Atherosclerosis studies were conducted with atheroprone progeroid Apoe–/– LmnaG609G/G609G mice[[1]](#1).

## Single-cell RNA-seq:

Aortas (including aortic arch and thoracic aorta) from Lmna+/+ and LmnaG609G/G609G mice were dissected, cleaned of perivascular fat, and opened longitudinally. Viable single cell suspensions were obtained using a previously described protocol with minor modifications[[2]](#2). We analyzed two replicate samples per genotype, and each sample contained pooled cells from the aortas of two male and two female animals to avoid possible sex-related bias. Cells were 126 loaded onto a Chromium Next GEM Chip G (10x Genomics), and libraries were created using the Next GEM Single cell 3’Library preparation kit v3.1 (10x Genomics) and indexed using the Chromium i7 Multiplex kit (10x Genomics). Libraries were paired-end sequenced in a HiSeq 4000 system (Illumina), and single cell transcriptomes were obtained using the 10x Genomics Cell Ranger 3.1.0 pipeline and analyzed with the Scater[[3]](#3) and Seurat[[4]](#4) R packages. We analyzed a total of 34,152 cells after removing predicted doublets and low-quality cells.



McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi:10.1093/bioinformatics/btw777 (https://doi.org/10.1093/bioinformatics/btw777).

Yuhan Hao and Stephanie Hao and Erica Andersen-Nissen and William M. Mauck III and Shiwei Zheng and Andrew Butler and Maddie J. Lee and Aaron J. Wilk and Charlotte Darby and Michael Zagar and Paul Hoffman and Marlon Stoeckius and Efthymia Papalexi and Eleni P. Mimitou and Jaison Jain and Avi Srivastava and Tim Stuart and Lamar B. Fleming and Bertrand Yeung and Angela J. Rogers and Juliana M. McElrath and Catherine A. Blish and Raphael Gottardo and Peter Smibert and Rahul Satija (2021). "Integrated analysis of multimodal single-cell data." Cell, Vo.l. 184, Issue 13,,
  title = {Integrated analysis of multimodal single-cell data},
  journal = {Cell},
  year = {2021},
  doi = {10.1016/j.cell.2021.04.048},
  url = {https://doi.org/10.1016/j.cell.2021.04.048},
}
