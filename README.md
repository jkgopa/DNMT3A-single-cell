# DNMT3A-single-cell

scRNAseq and single cell MIBI analysis in R for Rauch & Gopakumar paper

**Computational Analysis of scRNA-seq sequencing data**<br>
The BCL files were demultiplexed using 8 base-pair 10X sample indexes and “cellranger mkfastq” to generate paired-end FASTQ. We ran “cellranger count” to align the reads to the mouse UCSC mm10 reference genome using STAR aligner as well as perform filtering, barcode counting, and UMI counting. The alignment results were used to quantify the expression level of mouse genes and generation of gene-barcode matrix. Each sample's cellranger matrix was then loaded in a SeuratObject_4.1.0 using Seurat_4.1.1 (https://github.com/satijalab/seurat). The R code used for further analysis is available on GitHub, and summarized below.
QC and filtering parameters. Low quality cells, doublets, and potential dead cells were removed according to the percentage of mitochondrial genes and number of genes expressed in each cell (nFeature_RNA > 200 & nFeature_RNA < 4400 & percent.mt < 15). Filtering parameters were based on visualized QC metrics using the percent.mt/nCount_RNA and nFeature_RNA/nCount_RNA plots.
Normalization, variable feature selection, scaling and linear dimensionality reduction. Samples were then normalized using a global-scaling normalization method called Log Normalize. Using the FindVariableFeatures function, the top 2000 most highly variable genes were identified and used for downstream analysis to highlight biological signal. Cell cycle scores were assigned using Seurat's CellCycleScoring function. Linear transformation scaling was performed to prepare the data for dimensional reduction techniques. Cell cycle scores were regressed out during data scaling, followed by principal component analysis. 
Integration and clustering. The different objects were then integrated using Harmony (https://github.com/immunogenomics/harmony). Based on the ElbowPlot, we selected 8 significant dimensions to cluster the cells. Dimensionality reduction via UMAP embedding was performed on the integrated dataset, followed by the FindNeighbors and FindClusters functions. 
Doublet detection. Doublet detection and removal was performed on the 30-week dataset using doubletFinder_v3 (https://github.com/chris-mcginnis-ucsf/DoubletFinder).
Downstream analysis and data visualization. Identities of the cell clusters were determined using canonical cell type markers. For identification of cluster signature genes, we used the FindAllMarkers function in the Seurat package using the Wilcoxon Rank Sum test, testing genes that were detected in at least 25% of cells in either of the two groups with at least a 0.25 difference (log-scale) between the two groups, and at least 3 cells expressing the feature in at least one of the two groups. Differential gene expression within the Resident-like macrophage cluster was calculated using the FindMarkers function, using a negative binomial generalized linear model. We tested genes that were detected in at least 10% of cells in WT or mutant populations, with at least a 0.25 difference (log-scale) between the two groups of cells, and at least 3 cells expressing the feature in at least one of the two groups. The difference between the proportion of celltypes across different samples was calculated by a proportion test using the Single Cell Proportion Test R package (https://github.com/rpolicastro/scProportionTest).
Analysis of the replicate 24-week dataset. For the replicate (24-week) dataset, an analogous analysis was performed with three modifications: (1) filtering parameters were nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 15, again based on visualizing QC metrics; (2) we used sctransform-based normalization; (3) we selected 11 dimensions to cluster the cells, again based on the ElbowPlot.

**Multiplexed ion beam imaging time of flight (MIBI-TOF)**<br>
MIBI-TOF Imaging:<br>
Imaging was performed using a MIBI-TOF instrument with a Hyperion ion source. Xe+primary ions were used to sequentially sputter pixels for a given FOV. The following imaging parameters were used:<br>
• Pulse setting: Medium<br>
• Field size: 400 μm2 at 1024 x 1024 pixels<br>
• Dwell time: 0.5 ms<br>
• Median gun current on tissue: 4.11 nA Xe+<br>

Low-level Image Processing:<br>
Multiplexed image sets were extracted, slide background-subtracted, denoised, and aggregate filtered as previously described19. All parameters for these steps can be found in Extended Data Table 3.

Single Cell Segmentation:<br>
Nuclear segmentation and localization was performed using DeepCell on aggregate filtered images (https://deepcell.org).

Comments describing steps of analysis and generated figures are included in 3 .R files.
