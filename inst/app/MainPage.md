
# MoSpaTA (Version 0.3A)

Welcome to the Mouse Spatial Omic Testis Atlas (MoSpaTA), an interactive tool for exploring spatial and single-cell transcriptomic data of mouse testis.

This early alpha release of MoSpaTA aims to foster discussion and interest in developing a unified spatial omic atlas of the testis.

Conrad Lab 2024

Eisa Mahyari (@eisamahyari)


See: github.com/eisascience/MoSpaTA


## Underlying Data

MoSpaTA integrates the following key resources:

- **STseq MT1**: This dataset represents a spatial transcriptomic atlas of a healthy adult wild-type (WT) mouse. The cell-level data are derived from algorithms that define cell boundaries using a DAPI image, covering a total of approximately 35,000 cells to optimize boundary selection and enhance the performance of this interactive atlas.

- **MoDSTA (Mouse Developmental Single-cell Testis Atlas)**: Comprising 177,891 testis cells from 20 publicly available 10X Genomics samples, this dataset includes cells from 2 embryonic, 4 postnatal, and 14 WT adult samples. For performance optimization, MoDSTA has been randomly downsampled to 30,000 cells. Future versions will include integrated analyses combining STseq MT1 and MoDSTA datasets.

- **Machine Learning on MoDSTA**: Using SDA (soft clustering), our machine-learning models learn transcriptional patterns within high-dimensional transcriptomic spaces. Several SDA models have been trained on subsets of MoDSTA to enable specialized analysis.

## User Manual 

Navigate MoSpaTA using the left-side navigation bar to access various features.

### Load Data Tab

Utilize the five buttons available for data management, with the "Load All" button recommended for general use. If data loading issues arise, try loading each dataset individually for troubleshooting.

### Gene Expression Tab

Explore gene expression by typing the names of single or multiple genes. View gene expression maps derived from UMAP analyses of MoDSTA and STseq MT1 datasets, along with their spatial representations in the testis.

### MoDSTA MetaData Tab

Access visual representations and metadata related to cell type labels from the MoDSTA dataset.

### STseq MT1 MetaData Tab

Explore metadata and related visual content for the STseq MT1 dataset.

### SDA Component Browser

After loading the data (see Load Data tab), you can explore various SDA runs. 

There are several SDA runs that were performed on MoDSTA; some are runs on subsets of MoDSTA, others are simply different seeds of the combo objects. The aim is to have machine-learning train on every part of this data to find signatures of potential interest. These signatures (also called components) are effectively multi-geneic weighted signatures that fine-tune to a signature; providing surgical dissection of the signature. 

First, an SDA run needs to be selected, e.g., sda_FullCombo. 

Second, click on the two Apply Buttons, to project the SDA run in to each dataset. This produces a score, for each cell of each data for each component in each selected run. A run may have 20-100 components, which were selected a priori heuristically. 

Third, navigate the components either using the drop down or the next/prev buttons. 


### SDA Runs Table:

| ID | Description | Name |
|----|--------|-------------|
| 569507 | Other Minor Clusters of the Testis | Minor Clusters |
| 569509 | Leydig Cell Clusters | Leydig Clusters |
| 569516 | Germ Cell Clusters | GC Clusters |
| 569520 | Sertoli Cell Clusters | SC Clusters |
| 572405 | WT Adults only Combo | WT8Week Full Combo |
| 573391 | Juviniles and Prenatal only Combo | JuvPreNatal Combo |
| 574419 | Full Combo | Full Combo DS-seed1 |
| 574426 | Full Combo | Full Combo DS-seed2 |
| 574430 | Full Combo | Full Combo DS-seed3 |
| WTCombo | Adults Only | Adults Combo DS-seed1 |
| WTCombo2 | Adults Only | Adults Combo DS-seed2 |
| 568747 | Full Combo | Full Combo DS-seed4 |

* DS = Down sample



The genes and their weights are shown in the first figure. The top weighted genes are shown. The weights have direction, so the positive genes score towards the positive direction and the negative loaded genes, score towards the negative direction. The key point is the grouping of genes in each direction, as a 'soft' cluster of transcriptional correlation between those identified genes. In other words, the component has found a transcriptional signature between them. 

Next, How this scores cells in MoDSTA and STseqMT1 is visualized.



