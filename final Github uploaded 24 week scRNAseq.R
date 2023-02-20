#install.packages("cellpypes")

#set working directory to current location of script
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

library(Seurat)
library(data.table)
library(tidyverse)
library(cellpypes)
library(harmony)
library(dplyr)

# Basic function to convert human to mouse gene names
convert_human_to_mouse <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(mouse_gene in mouse_genes){
        output = append(output,mouse_gene)
      }
    }
  }
  
  return (output)
}

#read file that describes samples
samplekey <- fread("SampleKey.csv")
list_of_file_names <- samplekey$Sample
list_of_sample_names <- samplekey$SampleShort

#read in each seurat object
for (i in 1:length(list_of_file_names)) {
  file_name <- list_of_file_names[[i]]
  data_path <- paste0("path to file")
  print(data_path)
  data <- Read10X(data.dir = data_path)
  assign(paste0("object_", i), CreateSeuratObject(counts = data, project = list_of_sample_names[[i]]))
}

#merge objects and label objects
merged_object <- merge(object_1, y = c(object_2,object_3,object_4,object_5,object_6), 
                       project = "merged")
merged_object@meta.data

current.sample.ids <- unique(merged_object@meta.data$orig.ident)
new.sample.ids <- c("wt", "d3a", "tet2",
                     "wt", "wt", "wt")
merged_object@meta.data$genotype <- plyr::mapvalues(x = merged_object@meta.data$orig.ident, 
                                                    from = current.sample.ids, to = new.sample.ids)

#subset object
DefaultAssay(merged_object) <- "RNA"
merged_object <- PercentageFeatureSet(merged_object, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(merged_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
merged_object_subset <- subset(merged_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)

merged_object_subset <- readRDS("merged_object_subset.rds")
merged_object_subset@meta.data

#normalize data
DefaultAssay(merged_object_subset) <- "RNA"
merged_object_subset <- SCTransform(merged_object_subset, assay = 'RNA',
                                    new.assay.name = 'SCT',
                                    vars.to.regress = c('percent.mt', 'orig.ident'))

#regress out cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
m.s.genes <- convert_human_to_mouse(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convert_human_to_mouse(cc.genes.updated.2019$g2m.genes)

merged_object_subset <- CellCycleScoring(
  merged_object_subset,
  s.features = m.s.genes,
  g2m.features = m.g2m.genes,
  assay = 'SCT',
  set.ident = TRUE
)
merged_object_subset$CC.Difference <- merged_object_subset$S.Score - merged_object_subset$G2M.Score

merged_object_subset <- SCTransform(
  merged_object_subset,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'orig.ident', 'CC.Difference'))

#dimensionality reduction
merged_object_subset <- RunPCA(merged_object_subset)
merged_object_subset <- RunUMAP(merged_object_subset, reduction='pca', dims = 1:30)
merged_object_subset <- FindNeighbors(merged_object_subset, reduction='pca', dims = 1:30)
merged_object_subset <- FindClusters(merged_object_subset)

#save merged object
saveRDS(merged_object_subset, file = "jk_2022_merged_object_subset_cellcycle_regressed.rds")
###################################################################################################################################
merged_object_subset <- readRDS("jk_2022_merged_object_subset_cellcycle_regressed.rds")

#use harmony to integrate datasets (set seed = 2)
set.seed(2)
seurat_obj <- merged_object_subset %>%
  RunHarmony("orig.ident", assay.use="SCT")

#dimensionality reduction of integrated datasets
DefaultAssay(seurat_obj) <- "SCT"
#Extended Figure 5i
ElbowPlot(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, reduction='harmony', dims = 1:11)
seurat_obj <- FindNeighbors(seurat_obj, reduction='harmony', dims = 1:11)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.8)

GeneralPlot <- DimPlot(seurat_obj, label = TRUE)

new_ids_24wk <- RenameIdents(object = seurat_obj, 
                             `0` = "B-cells", `1` = "B-cells", 
                             `2` = "T-cells", `3` = "T-cells",
                             `4` = "Trem2hi Macrophages", `5` = "Tissue-resident Macrophages",
                             `6` = "B-cells", `7` = "T-cells",
                             `8` = "T-cells", `9` = "B-cells",
                             `10` = "Inflammatory Macrophages", `11` = "Inflammatory Macrophages",
                             `12` = "T-cells", `13` = "T-cells",
                             `14` = "Dendritic Cells", `15` = "Cxcr6 T-cells",
                             `16` = "Cxcr6 T-cells", `17` = "Cycling T-cells",
                             `18` = "Cycling T-cells", `19` = "T-cells",
                             `20` = "Monocytes", `21` = "Other",
                             `22` = "T-cells", `23` = "T-cells",
                             `24` = "Neutrophils", `25` = "Dendritic Cells",
                             `26` = "Neutrophils", `27` = "Cycling T-cells",
                             `28` = "Other", `29` = "B-cells",
                             `30` = "T-cells", `31` = "Other",
                             `32` = "T-cells", `33` = "Dendritic Cells",
                             `34` = "B-cells", `35` = "Cxcr6 T-cells",
                             `36` = "T-cells")

paired = c("B-cells"="#A6CDE2","Cxcr6 T-cells"="#1E78B4","Inflammatory Macrophages"="#74C476",
           "Monocytes"="#B15928","T-cells"="#F59899","Tissue-resident Macrophages"="#E11E26",
           "Trem2hi Macrophages"="#FCBF6E","Dendritic Cells"="#F47E1F",
           "Cycling T-cells"="#6A3E98","Neutrophils" = "#FAF39B","Other" = "#808080")
my_cols2 <- paired[order(as.integer(names(paired)))]
scales::show_col(my_cols2)
DimPlot(new_ids_24wk, cols = my_cols2)

new_ids_24wk@meta.data$genotype <- factor(new_ids_24wk@meta.data$genotype, levels=c('wt','d3a','tet2'))
table(new_ids_24wk@meta.data$orig.ident)
new_ids_24wk@meta.data$orig.ident <- factor(new_ids_24wk@meta.data$orig.ident, levels=c('wt_1','wt_2','wt','tet2','wt_d3a','d3a_d3a'))
DimPlot(new_ids_24wk, cols = my_cols2)

#Extended Data 7a
DimPlot(new_ids_24wk, cols = my_cols2, split.by = 'orig.ident', ncol = 2)

#QC
new_ids_24wk@meta.data$GenesPerUMI <- new_ids_24wk@meta.data$nFeature_RNA/new_ids_24wk@meta.data$nCount_RNA

#Extended Table 5 and Extended Figure 5h
QC_input_2022_sample <- as_tibble(subset(new_ids_24wk@meta.data, orig.ident == "d3a_d3a")) %>% 
  dplyr::select(percent.mt, nCount_RNA, nFeature_RNA, GenesPerUMI)
cluster_QC_2022_sample <- as.data.frame(do.call(cbind, lapply(QC_input_2022_sample, summary)))
cluster_QC_2022_sample$sample <- "d3a_d3a"
remaining_sample_type_list <- list("wt_d3a","tet2","wt","wt_1","wt_2")
for (x in remaining_sample_type_list) {
  QC_input_2022_sample <- subset(new_ids_24wk@meta.data, orig.ident == x) %>% dplyr::select(percent.mt, nCount_RNA, nFeature_RNA, GenesPerUMI)
  cluster_QC_2022_next_sample <- as.data.frame(do.call(cbind, lapply(QC_input_2022_sample, summary)))
  cluster_QC_2022_next_sample$sample <- x
  cluster_QC_2022_sample <- cluster_QC_2022_sample %>% 
    bind_rows(cluster_QC_2022_next_sample)
}


cell_counts <- table(Idents(new_ids_24wk), new_ids_24wk@meta.data$genotype)

cell_counts_orig_ident <- table(Idents(new_ids_24wk), new_ids_24wk@meta.data$orig.ident)

CHIP_cluster_type <- new_ids_24wk
FMC8C1_cells <- WhichCells(CHIP_cluster_type, expression = Folr2 > 0 & Mrc1 > 0 & Ccl8 > 0 & Cxcl1 > 0)
CHIP_cluster_type@meta.data$genotype <- factor(CHIP_cluster_type@meta.data$genotype, levels=c('wt','d3a','tet2'))

CHIP_cluster_type@meta.data$orig.ident <- factor(CHIP_cluster_type@meta.data$orig.ident, levels=c('wt_1','wt_2','wt','tet2','wt_d3a','d3a_d3a'))

#Extended Figure 7b
DimPlot(CHIP_cluster_type, cells.highlight = FMC8C1_cells, sizes.highlight = 2, split.by = 'orig.ident', ncol = 2)

#Figure 4d
Idents(object = CHIP_cluster_type, cells = FMC8C1_cells) <- 'Folr2+ Mrc1+ Ccl8+ Cxcl1+'
chip_cluster_counts <- data.frame(unclass(table(CHIP_cluster_type@meta.data$genotype, Idents(CHIP_cluster_type))))
chip_cluster_counts_percentage <- data.frame(apply(chip_cluster_counts, 1, function(x){x*100/sum(x,na.rm=T)}))
subset <- chip_cluster_counts_percentage[1,] %>% as_tibble()
subset

#Figure 4e
myeloid_clusters <- subset(x = seurat_obj, 
                           idents = c(1,6,9))

myeloid_count_table <- data.frame(unclass(table(Idents(myeloid_clusters), myeloid_clusters@meta.data$genotype)))
myeloid_count_table_percentage <- apply(myeloid_count_table, 1, function(x){x*100/sum(x,na.rm=T)})
myeloid_count_tibble <- myeloid_count_table_percentage %>% as_tibble()
myeloid_count_tibble$genotype <- c("d3a", "tet2", "wt")


