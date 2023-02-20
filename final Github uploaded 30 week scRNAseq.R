# Load libraries
library(Seurat)
library(data.table)
library(tidyverse)
library(cellpypes)
library(harmony)
library(dplyr)
library(gdata)
library(scProportionTest)
library(DoubletFinder)
library(parallel)
library(ggtree)
library(EnhancedVolcano)

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

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
m.s.genes <- convert_human_to_mouse(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convert_human_to_mouse(cc.genes.updated.2019$g2m.genes)

# Read in matrices
wt_451.data <- Read10X(data.dir = "wt_451")
dnmt3a_451.data <- Read10X(data.dir = "dnmt3a_451")
tet2_451.data <- Read10X(data.dir = "tet2_451")
dnmt3a_452.data <- Read10X(data.dir = "dnmt3a_452")
tet2_452.data <- Read10X(data.dir = "tet2_452")

# Create Seurat object
plaque <- CreateSeuratObject(counts = cbind(wt_451.data, dnmt3a_451.data, tet2_451.data, dnmt3a_452.data, tet2_452.data), project = "PLAQUE", min.cells = 5)
plaque@meta.data$genotype <- c(rep("wt_451", ncol(wt_451.data)), rep("dnmt3a_451", ncol(dnmt3a_451.data)), rep("tet2_451", ncol(tet2_451.data)), rep("dnmt3a_452", ncol(dnmt3a_452.data)), rep("tet2_452", ncol(tet2_452.data)))

# Filter, normalize, scale, linear dimensionality reduction
plaque[["percent.mt"]] <- PercentageFeatureSet(plaque, pattern = "^mt-")
plot1 <- FeatureScatter(plaque, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(plaque, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plaque <- subset(plaque, subset = nFeature_RNA > 200 & nFeature_RNA < 4400 & percent.mt < 15)
plaque <- Seurat::NormalizeData(plaque, verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%   CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) %>% ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA(pc.genes = plaque@var.genes, npcs = 20, verbose = FALSE)
ElbowPlot(plaque)

# Run harmony integration
plaque_harmony <- plaque %>%
  RunHarmony("genotype")
ElbowPlot(plaque_harmony)
plaque_harmony <- RunUMAP(plaque_harmony, reduction='harmony', dims = 1:8)
plaque_harmony <- FindNeighbors(plaque_harmony, reduction='harmony', dims = 1:8)
plaque_harmony <- FindClusters(plaque_harmony, resolution = 0.7)

# Doublet detection and removal
sweep.list <- paramSweep_v3(plaque_harmony, PCs = 1:min.pc, num.cores = detectCores() - 1)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)
nExp <- round(0.03*nrow(plaque_harmony@meta.data)) 
plaque_doublets <- doubletFinder_v3(plaque_harmony, PCs = 1:8, pN = 0.25, pK = 0.09, nExp = nExp, reuse.pANN = FALSE, sct = FALSE)
plaque_harmony = plaque_doublets[, plaque_doublets@meta.data[, DF.name] == "Singlet"]

# Find markers to enable cluster annotation (supplementary file)
markers <- FindAllMarkers(plaque_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
write.csv(top10, file = "top10.csv")

# Name and color clusters
new_ids_30wk <- RenameIdents(object = plaque_harmony, 
                             `0` = "B cells", `1` = "Inflammatory macrophages", 
                             `2` = "T cells", `4` = "Trem2hi macrophages",
                             `3` = "Resident-like macrophages", `5` = "T cells",
                             `6` = "Neutrophils", `9` = "Dendritic cells",
                             `7` = "B cells", `10` = "Cxcr6+ T cells",
                             `8` = "Monocytes", `11` = "Dendritic cells",
                             `12` = "B cells", `16` = "Mixed lymphocytes",
                             `14` = "Proliferating T cells", `13` = "Proliferating T cells", '17' = "Vascular fibroblasts", '15' = "Mixed lymphocytes")
new_ids_30wk@meta.data$orig.ident <- new_ids_30wk@meta.data$genotype
current.sample.ids <- unique(new_ids_30wk@meta.data$orig.ident)
new.sample.ids <- c("wt", "wt", "wt", "dnmt3a", "tet2")
new_ids_30wk@meta.data$genotype <- plyr::mapvalues(x = new_ids_30wk@meta.data$orig.ident, from = current.sample.ids, to = new.sample.ids)
paired = c("B cells"="#A6CDE2","Cxcr6+ T cells"="#1E78B4","Inflammatory macrophages"="#74C476",
           "Monocytes"="#B15928","T cells"="#F59899","Resident-like macrophages"="#E11E26",
           "Trem2hi macrophages"="#FCBF6E","Dendritic cells"="#F47E1F",
           "Proliferating T cells"="#6A3E98","Neutrophils" = "#FAF39B", "Vascular fibroblasts" = "#808080", "Mixed lymphocytes" = "#006600")
my_cols2 <- paired[order(as.integer(names(paired)))]
scales::show_col(my_cols2)
new_ids_30wk@meta.data$genotype <- factor(new_ids_30wk@meta.data$genotype, levels=c('wt','dnmt3a','tet2'))
new_ids_30wk@meta.data$orig.ident <- factor(new_ids_30wk@meta.data$orig.ident, levels=c('wt_451','dnmt3a_451', 'dnmt3a_452', 'tet2_451', 'tet2_452'))

# UMAP Fig. 3b
pdf(file = "30 wk UMAP split by sample -dbl recolor.pdf",
    width = 11, 
    height = 3)
DimPlot(new_ids_30wk, cols = my_cols2, split.by = "orig.ident")
dev.off()

# Feature plot Ext Data Fig. 6a 
CHIP_features_feature_plot <- c("Folr2", "F13a1", "Lyve1", "Cxcl1", "Pf4", "Ccl8", "Mrc1", "Ccl2", "Ccl7", "Jun", "Egr1", "Fos")
pdf(file = "ext6a 30 wk UMAP CHIP Features.pdf",
    width = 10, 
    height = 30)
FeaturePlot(new_ids_30wk, features = CHIP_features_feature_plot, split.by = 'genotype', min.cutoff = 'q25', label.size = 2) & 
  theme(axis.title.x = element_blank()) & theme(axis.title.y = element_blank()) & theme(legend.position = "right")
dev.off()

# Mononuclear phagocyte heatmap Fig. 3c
celltype_up_markers <- FindAllMarkers(new_ids_30wk, min.pct = 0.25, only.pos = TRUE)
celltype_up_markers_paired_down <- celltype_up_markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
activefilter = c('Monocytes', 'Dendritic cells', 'Inflammatory macrophages',
                 'Trem2hi macrophages','Resident-like macrophages')
myeloid_FindALLMarkers <- filter(celltype_up_markers_paired_down,cluster %in% activefilter)%>%
  mutate(cluster =  factor(cluster, levels = activefilter)) %>%
  arrange(cluster)
myeloid_FindALLMarkers_heatmap_genes <- myeloid_FindALLMarkers$gene
all_myeloid <- subset(x = new_ids_30wk, idents = c("Resident-like macrophages", "Inflammatory macrophages", "Trem2hi macrophages", "Monocytes", "Dendritic cells"))
Idents(all_myeloid) <- factor(Idents(all_myeloid), levels=c('Monocytes', 'Dendritic cells', 'Inflammatory macrophages', 'Trem2hi macrophages','Resident-like macrophages'))
pdf(file = "3c_30wk_allmyeloid_FindALLMarkers_heatmap.pdf",
    width = 12, 
    height = 13)
DoHeatmap(all_myeloid, features = myeloid_FindALLMarkers_heatmap_genes, raster = TRUE)
dev.off()

# Mononuclear phagocyte proportions-of-clusters analysis (Fig. 3d)
Idents(all_myeloid) <- all_myeloid@meta.data$genotype
all_myeloid_prop_test <- sc_utils(all_myeloid)
all_myeloid_D3A_prop_test <- permutation_test(
  all_myeloid_prop_test, cluster_identity = "celltype",
  sample_1 = "wt", sample_2 = "dnmt3a",
  sample_identity = "genotype"
)
permutation_plot(all_myeloid_D3A_prop_test)
all_myeloid_D3A_plot_data <- copy(all_myeloid_D3A_prop_test@results$permutation)
all_myeloid_D3A_plot_data[, clusters := fct_reorder(factor(clusters), obs_log2FD)]
all_myeloid_TET_prop_test <- permutation_test(
  all_myeloid_prop_test, cluster_identity = "celltype",
  sample_1 = "wt", sample_2 = "tet2",
  sample_identity = "genotype"
)
permutation_plot(all_myeloid_TET_prop_test)
all_myeloid_TET_plot_data <- copy(all_myeloid_TET_prop_test@results$permutation)
all_myeloid_TET_plot_data[, clusters := fct_reorder(factor(clusters), obs_log2FD)]
mac_mono_TET_plot_data$Genotype <- "Tet2 KO"
mac_mono_D3A_plot_data$Genotype <- "Dnmt3a KO"
all_myeloid_combo_prop_plot <- all_myeloid_TET_plot_data %>% 
  bind_rows(all_myeloid_D3A_plot_data)
write.csv(all_myeloid_combo_prop_plot, "all_myeloid_scprop_data.csv")

# Cell cycle analysis (Ext Data Fig. 5f)
macrophages <- subset(new_ids_30wk, idents = c("Resident-like macrophages", "Trem2hi macrophages", "Inflammatory macrophages"))
macrophages@meta.data$type_geno <- paste0(macrophages@meta.data$genotype," ",macrophages@meta.data$celltype)
phase_count_table2 <- data.frame(unclass(table(macrophages@meta.data$type_geno, macrophages@meta.data$Phase)))
phase_count_table_percentage2 <- apply(phase_count_table2, 1, function(x){x*100/sum(x,na.rm=T)})
phase_count_tibble2 <- phase_count_table_percentage2 %>% as_tibble()
phase_count_tibble2$Phase <- c("G1", "G2M", "S")
phase_data_long2 <- gather(phase_count_tibble2, CellType, Percent, "dnmt3a Inflammatory macrophages":"wt Trem2hi macrophages", factor_key=TRUE)
macrophage_geno_cell_cycle <- ggplot(phase_data_long2, aes(fill=Phase, y=Percent, x=CellType)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Cell Cluster")
write.csv(phase_count_tibble2, "ext5f_30wk_10X_mac_genotype_cellcycle_phases.csv")

# Volcano plots Fig. 4a
d3aTR.de.markers_all_1 <- FindMarkers(new_ids_30wk, ident.1 = "dnmt3a Resident-like macrophages", ident.2 = "wt Resident-like macrophages", logfc.threshold = 0.25, test.use = "negbinom",
                                      min.pct = 0.1, only.pos = FALSE) %>% rownames_to_column(var = 'gene') %>% arrange(desc(avg_log2FC))
t2TR.de.markers_all_1 <- FindMarkers(new_ids_30wk, ident.1 = "tet2 Resident-like macrophages", ident.2 = "wt Resident-like macrophages", logfc.threshold = 0.25, test.use = "negbinom",
                                     min.pct = 0.1, only.pos = FALSE) %>% rownames_to_column(var = 'gene') %>% arrange(desc(avg_log2FC))
pdf(file = "volcano_TR_d3a.pdf",
    width = 5, 
    height = 4.5)
ggplot(d3aTR.de.markers_all_1, aes(x = avg_log2FC, y = -log10(p_val_adj))) + geom_vline(xintercept=c(-0.5,0.5),linetype=3) + geom_hline(yintercept=12,linetype=3) + geom_point(alpha = 0.4, shape = 16, size = 0.6, col = "darkgray") + theme_classic(base_size = 14) + geom_point(data = filter(d3aTR.de.markers_all_1, -log10(p_val_adj)>12 & (avg_log2FC>0.5 | avg_log2FC<(-0.5))), alpha = 1, shape = 21, size = 1.5, fill = "#ba3434") + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) + geom_text_repel(
  data = filter(d3aTR.de.markers_all_1, -log10(p_val_adj)>12 & (avg_log2FC>0.5 | avg_log2FC<(-0.5))),
  aes(label = gene),
  size = 4,
  box.padding = unit(0.3, "lines"),
  point.padding = unit(0.3, "lines"),
  max.overlaps = 25) 
dev.off()
pdf(file = "volcano_TR_t2.pdf",
    width = 5, 
    height = 4.5)
ggplot(t2TR.de.markers_all_1, aes(x = avg_log2FC, y = -log10(p_val_adj))) + geom_vline(xintercept=c(-0.5,0.5),linetype=3) + geom_hline(yintercept=12,linetype=3) + geom_point(alpha = 0.4, shape = 16, size = 0.6, col = "darkgray") + theme_classic(base_size = 14) + geom_point(data = filter(t2TR.de.markers_all_1, -log10(p_val_adj)>12 & (avg_log2FC>0.5 | avg_log2FC<(-0.5))), alpha = 1, shape = 21, size = 1.5, fill = "#ba3434") + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) + geom_text_repel(
  data = filter(t2TR.de.markers_all_1, -log10(p_val_adj)>12 & (avg_log2FC>0.5 | avg_log2FC<(-0.5))),
  aes(label = gene),
  size = 4,
  box.padding = unit(0.3, "lines"),
  point.padding = unit(0.3, "lines"),
  max.overlaps = 25) 
dev.off()

# Violin plots Fig. 4b
TR <- subset(x = new_ids_30wk, idents = c("Resident-like macrophages"))
Idents(TR) <- TR@meta.data$genotype
CHIP_features_2 <- c("Folr2","Mrc1","Ccl8","Cxcl1","Jun","Egr1")
pdf(file = "4b 30 wk violin plots orig_ident 112722-only myeloid.pdf",
    width = 4, 
    height = 10)
VlnPlot(object = TR, features = CHIP_features_2, ncol = 2, cols = c("grey", "salmon", "steelblue1")) & theme(axis.title.x = element_blank()) & 
  theme(axis.title.y = element_blank()) & theme(axis.text.x = element_blank())
dev.off()

# CHIP cluster identification Fig. 4c
CHIP_cluster_type <- new_ids_30wk
FMC8C1_cells <- WhichCells(CHIP_cluster_type, expression = Folr2 > 0 & Mrc1 > 0 & Ccl8 > 0 & Cxcl1 > 0)
CHIP_cluster_type@meta.data$genotype <- factor(CHIP_cluster_type@meta.data$genotype, levels=c('wt','dnmt3a','tet2'))
CHIP_cluster_type@meta.data$orig.ident <- factor(CHIP_cluster_type@meta.data$orig.ident, levels=c('wt_451','dnmt3a_451','tet_451','dnmt3a_452','tet2_452'))
pdf(file = "4c 30 wk CHIP cluster orig_ident 112722.pdf",
    width = 12, 
    height = 5)
DimPlot(CHIP_cluster_type, cells.highlight = FMC8C1_cells, cols.highlight = 'cyan4', sizes.highlight = 2, split.by = 'genotype', ncol = 3)
dev.off()

# Enumerate CHIP cluster Fig. 4d
Idents(object = CHIP_cluster_type, cells = FMC8C1_cells) <- 'Folr2+ Mrc1+ Ccl8+ Cxcl1+'
chip_cluster_counts <- data.frame(unclass(table(CHIP_cluster_type@meta.data$genotype, Idents(CHIP_cluster_type))))
write.csv(chip_cluster_counts, "4d 30 wk chip cluster counts.csv")

# Monocyte subset analysis (Ext Data Fig. 5d)
new_ids_30wk_sub <- FindSubCluster(
  new_ids_30wk,
  cluster="Monocytes",
  "RNA_snn",
  subcluster.name = "subsets",
  resolution = 0.12,
  algorithm = 1
)
DimPlot(new_ids_30wk_sub, reduction = "umap", group.by = "subsets", label = TRUE, label.size = 3)
new_ids_30wk_sub <- SetIdent(new_ids_30wk_sub, value = new_ids_30wk_sub@meta.data$subsets)
mono <- subset(x = new_ids_30wk_sub, idents = c("Monocytes_0","Monocytes_1"))
sub_up_markers <- FindAllMarkers(mono, min.pct = 0.25, only.pos = TRUE)
sub_up_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> sub_top20
cluster_counts <- data.frame(unclass(table(mono@meta.data$genotype, mono@meta.data$subsets)))
write.csv(cluster_counts, "mono sub counts by genotype.csv")

#Lymphoid proportions-of-clusters analysis (Ext Data Fig. 5e)
lymphoid <- subset(x = new_ids_30wk, idents = c("Cxcr6+ T cells", "T cells", 
                                                "B cells", "Proliferating T cells"))
Idents(lymphoid) <- TR@meta.data$genotype
lymphoid_prop_test <- sc_utils(lymphoid)
lymphoid_D3A_prop_test <- permutation_test(
  lymphoid_prop_test, cluster_identity = "celltype",
  sample_1 = "wt", sample_2 = "dnmt3a",
  sample_identity = "genotype"
)
permutation_plot(lymphoid_D3A_prop_test)
lymphoid_D3A_plot_data <- copy(lymphoid_D3A_prop_test@results$permutation)
lymphoid_D3A_plot_data[, clusters := fct_reorder(factor(clusters), obs_log2FD)]
lymphoid_TET_prop_test <- permutation_test(
  lymphoid_prop_test, cluster_identity = "celltype",
  sample_1 = "wt", sample_2 = "tet2",
  sample_identity = "genotype"
)
permutation_plot(lymphoid_TET_prop_test)
lymphoid_TET_plot_data <- copy(lymphoid_TET_prop_test@results$permutation)
lymphoid_TET_plot_data[, clusters := fct_reorder(factor(clusters), obs_log2FD)]
lymphoid_TET_plot_data$Genotype <- "Tet2 KO"
lymphoid_D3A_plot_data$Genotype <- "Dnmt3a KO"
lymphoid_combo_prop_plot <- lymphoid_TET_plot_data %>% 
  bind_rows(lymphoid_D3A_plot_data)
write.csv(lymphoid_combo_prop_plot, "ext5e_30wk_lymphoid_scprop_data.csv")