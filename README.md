# NHD_snRNA

This repository contains the code used for snRNA-seq analysis by Zhou et al. (2022) of Nasu-Hakola Disease (NHD) patient brains. 

## Code
Included are the codes necessary to replicate the analyses.

 - `cell type clustering`: clustering analyses of all cells and subclustering of each cell type, using `Seurat` and `Harmony`.
 - `nichenet`: ligand-receptor interaction analyses between other brain cells and microglia, using `nichenet`.
 - `UCell_plot`: geneset score calculation and plotting, using `UCell`.
 
## Data
Raw and processed data are available at the Gene Expression Omnibus (GEO) database under accession number GSE190013. 

Processed data, UMAP coordinates and annotations are freely available to download and visualize at UCSC Cell Browser (https://nhd-brain.cells.ucsc.edu).


## Acknowledgements
This work was supoorted by the NIH (RF1 AG051485, R21 AG059176, and RF1 AG059082) (M. Colonna), Centene (M. Colonna), the Strategic Research Program for Brain Sciences from Japan Agency for Medical Research and Development (AMED, JP21wm0425019) (A. Kakita), the Collaborative Research Project of the Brain Research Institute, Niigata University, Japan (A. Kakita),  AMED (JP21wm0425019) (M. Takao), and Grants-in-Aids from the Ministry of Education, Culture, Sports, Science and Technology (MEXT) for Scientific Research on Innovative Areas (KAKENHI: 22H04798) (M. Takao).

We would like to thank Mari Tada (Niigata University, Japan), Akiyoshi Kakita (Niigata University, Japan), and Masaki Takao (National Center of Neurology and Psychiatry, Japan) for their contributions of postmortem tissue. In addition, we would like to thank Vincent Peng and Marina Terekhova for advice on bioinformatic analyses.
