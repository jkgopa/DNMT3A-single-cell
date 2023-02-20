library(tidyverse)
library(dplyr)
library(data.table)
library(ggpubr)

#read in single cell data
root1_4 <- read.csv("...root1_4/single_cell_output/cell_table_size_normalized.csv") %>% 
  as.tibble()
root5_8 <- read.csv(".../root5_8/single_cell_output/cell_table_size_normalized.csv") %>% 
  as.tibble()
root9 <- read.csv (".../root9/single_cell_output/cell_table_size_normalized.csv") %>% 
  as.tibble()

#merge data from all 9 roots and define cell types
eight_roots <- full_join(root1_4, root5_8)
all_roots <- full_join(eight_roots, root9)
all_roots$cell_type <- ifelse(all_roots$SMA>0,"Media",
                              ifelse(all_roots$CD206>0 & (all_roots$CD68>0 | all_roots$F4.80>0) & all_roots$CD452>0, "CD206_Macs",
                                     ifelse(all_roots$CD11c>0 & all_roots$CD68>0 & all_roots$CD452>0 & all_roots$MHCII>0, "Lesional_Inflamm_Macs",
                                            ifelse(all_roots$TREM2>0 & all_roots$CD68>0 & all_roots$CD452>0 & all_roots$MHCII<0.5, "Lesional_Lipid_Macs",
                                                   ifelse(all_roots$F4.80>0 & all_roots$CD452>0 & all_roots$MHCII>0, "Other_Macrophages",
                                                          ifelse(all_roots$F4.80==0 & all_roots$CD452>0 & all_roots$MHCII>0, "Myeloid_Cells",
                                                                 ifelse(all_roots$CD45.1>0, "Host_Hematopoietic_Cells",
                                                                        ifelse(all_roots$CD3>0 & all_roots$CD452>0, "T_cells",
                                                                               ifelse(all_roots$B220>0 & all_roots$CD452>0, "B_cells",
                                                                                      ifelse((all_roots$CD31>0 | all_roots$VWF>0)& all_roots$CD105>0, "Activated_Endothelium",
                                                                                             ifelse((all_roots$CD31>0 | all_roots$VWF>0)& all_roots$CD105==0, "Non_Activated_Endothelium","Other")))))))))))
all_roots$ID <- paste(all_roots$label, all_roots$cell_type, sep="-")

#read in cell id's for each cell found in adventitia of 9 roots
i = 1
adv_filedirectory = ".../CellsByRegion/Adventitia/Point"
adv_pointname = paste0(adv_filedirectory,i,".csv")
point1_adv_csv <- read.csv(adv_pointname, header = FALSE)
point1_all <- filter(all_roots, fov == "Point1")

adventitia_cells <- point1_all[point1_all$label %in% point1_adv_csv[,1],]
j <- 2

while (j < 10) {
  adv_pointname = paste0(adv_filedirectory,j,".csv")
  point_adv_csv <- read.csv(adv_pointname, header = FALSE)
  point_all_cells <- filter(all_roots, fov == paste0("Point",j))
  point_adv_cells <- point_all_cells[point_all_cells$label %in% point_adv_csv[,1],]
  adventitia_cells <- full_join(adventitia_cells,point_adv_cells)
  j = j+1
}

#calculate density of cell types in each root
adventitia_celltypes_per_point <- as.data.frame.matrix(table(adventitia_cells$fov, adventitia_cells$cell_type))
point_density <- c(318645.964, 394205.74, 322254.173, 351044.195, 360103.462, 341590.796, 310608.637, 442344.588, 274843.496)
density_celltypes_per_point <- adventitia_celltypes_per_point/point_density
density_celltypes_per_point_tibble <- rownames_to_column(density_celltypes_per_point, var = "Point") %>% as_tibble()

#label genotypes of each root
density_celltypes_per_point_tibble$genotype <- ifelse(density_celltypes_per_point_tibble$Point == "Point1", "Tet2",
                                                      ifelse(density_celltypes_per_point_tibble$Point == "Point2", "D3a",
                                                             ifelse(density_celltypes_per_point_tibble$Point == "Point3", "D3a",
                                                                    ifelse(density_celltypes_per_point_tibble$Point == "Point4", "WT",
                                                                           ifelse(density_celltypes_per_point_tibble$Point == "Point5", "Tet2",
                                                                                  ifelse(density_celltypes_per_point_tibble$Point == "Point6", "D3a",
                                                                                         ifelse(density_celltypes_per_point_tibble$Point == "Point7", "WT",
                                                                                                ifelse(density_celltypes_per_point_tibble$Point == "Point8", "Tet2",
                                                                                                       ifelse(density_celltypes_per_point_tibble$Point == "Point9", "WT",
                                                                                                              NA)))))))))

density_celltypes_per_point_long <- gather(density_celltypes_per_point_tibble, cell_type, density, Activated_Endothelium:T_cells, factor_key=TRUE)
density_celltypes_per_point_long$genotype <- factor(density_celltypes_per_point_long$genotype, levels = c("Tet2", "D3a", "WT"))

#Figure 5b
density_CD206_per_point <- density_celltypes_per_point_long %>% filter(cell_type == "CD206_Macs")
write.csv(density_CD206_per_point, "CD206DensityPerPoint.csv")

ggbarplot(density_CD206_per_point, x = "cell_type", y = "density",
          add = c("mean_se", "point"),
          color = "genotype", fill = "genotype", alpha = 0.5,
          position = position_dodge(0.9))


#####################################################################################################################################################################################
#Distance Matrix Analysis
#####################################################################################################################################################################################
point_dist_home <- "/Users/Jk/Desktop/JG103_roots.nosync/CellsByRegion/Full_Root/Dist_matrixes/Point"

point = 1
combined = NULL

#creating matrix of distances of cells from CD206 cells
while (point < 10) {
  point_dist_filepath <- paste0(point_dist_home,point,".csv")
  point_dist_matrix <- read.csv(point_dist_filepath, header = FALSE) %>% 
    as.tibble()
  
  fov_num <- paste0("Point",point)
  
  root_cells <- dplyr::filter(all_roots, fov == fov_num)
  root_cells$cell_type <- ifelse(root_cells$SMA>0,"Media",
                                 ifelse(root_cells$CD206>0 & (root_cells$CD68>0 | root_cells$F4.80>0) & root_cells$CD452>0, "CD206_Macs",
                                        ifelse(root_cells$CD11c>0 & root_cells$CD68>0 & root_cells$CD452>0 & root_cells$MHCII>0, "Lesional_Inflamm_Macs",
                                               ifelse(root_cells$TREM2>0 & root_cells$CD68>0 & root_cells$CD452>0 & root_cells$MHCII<0.5, "Lesional_Lipid_Macs",
                                                      ifelse(root_cells$F4.80>0 & root_cells$CD452>0 & root_cells$MHCII>0, "Other_Myeloid_Cells",
                                                             ifelse(root_cells$F4.80==0 & root_cells$CD452>0 & root_cells$MHCII>0, "Other_Myeloid_Cells",
                                                                    ifelse(root_cells$CD45.1>0, "Host_Hematopoietic_Cells",
                                                                           ifelse(root_cells$CD3>0 & root_cells$CD452>0, "T_cells",
                                                                                  ifelse(root_cells$B220>0 & root_cells$CD452>0, "B_cells",
                                                                                         ifelse((root_cells$CD31>0 | root_cells$VWF>0)& root_cells$CD105>0, "Activated_Endothelium",
                                                                                                ifelse((root_cells$CD31>0 | root_cells$VWF>0)& root_cells$CD105==0, "Non_Activated_Endothelium","Other")))))))))))
  root_cells$ID <- paste(root_cells$label, root_cells$cell_type, sep="-")
  
  names(point_dist_matrix) <- root_cells$ID
  point_dist_matrix$cell_ID <- root_cells$ID
  head(point_dist_matrix)
  
  x <- nrow(point_dist_matrix)
  df.long <- pivot_longer(point_dist_matrix, cols=1:x, names_to = "cell_type", values_to = "distance")
  cells_minus_self <- subset(df.long, df.long[ , 3] != 0)
  cells_minus_self$point_num <- point
  cells_minus_self$parent <- sub(".*-", "", cells_minus_self$cell_ID)
  cells_minus_self$child <- sub(".*-", "", cells_minus_self$cell_type)
  cells_CD206_parent <- dplyr::filter(cells_minus_self, parent == "CD206_Macs")
  combined <- bind_rows(combined, cells_CD206_parent)
  point = point+1
}

table(combined$point_num)  
write.csv(combined, "...CD206_Cells_FullRoots_All_Points_newPopula.csv")
combined <- read_csv(".../CD206_Cells_FullRoots_All_Points_newPopula.csv")

#label genotypes
combined$genotype <- ifelse(combined$point_num == "1", "Tet2",
                            ifelse(combined$point_num == "2", "D3a",
                                   ifelse(combined$point_num == "3", "D3a",
                                          ifelse(combined$point_num == "4", "WT",
                                                 ifelse(combined$point_num == "5", "Tet2",
                                                        ifelse(combined$point_num == "6", "D3a",
                                                               ifelse(combined$point_num == "7", "WT",
                                                                      ifelse(combined$point_num == "8", "Tet2",
                                                                             ifelse(combined$point_num == "9", "WT",
                                                                                    NA)))))))))

table(combined$genotype, combined$point_num)

#creating list of adventitial cells
advent_cell_list_home <- "/Users/Jk/Desktop/JG103_roots.nosync/CellsByRegion/Adventitia/Point"

adv_point = 1
combined_adv_list = NULL

while (adv_point < 10) {
  adv_point_filepath <- paste0(advent_cell_list_home,adv_point,".csv")
  adv_point_list <- read.csv(adv_point_filepath) %>% 
    as.tibble() %>% 
    rename(label = V1, region = V2)
  adv_point_list$point_num <- adv_point
  combined_adv_list <- bind_rows(combined_adv_list, adv_point_list)
  adv_point = adv_point+1
  print(adv_point)
}

table(combined_adv_list$point_num)

#subset combined distance cell matrix to just adventitia; first select for parents in adventitia, then children in adventitia
#cd206 cells in adventitia
cd206_adventitia_dist <- combined %>% separate(cell_ID, c("label", "parent_type"), sep = "-") %>%
  mutate(label = as.integer(label))
adventitia_CD206_cells <- cd206_adventitia_dist %>%
  semi_join(combined_adv_list, by = c("label", "point_num"))

table(adventitia_CD206_cells$point_num)
table(adventitia_CD206_cells$genotype) 

#children in adventitia
child_combined_adv_list <- combined_adv_list %>% rename(child_label = label)
child_cd206_adventitia_dist <- adventitia_CD206_cells %>% separate(cell_type, c("child_label", "child_type"), sep = "-") %>%
  mutate(child_label = as.integer(child_label))
adventitia_cells_only <- child_cd206_adventitia_dist %>%
  semi_join(child_combined_adv_list, by = c("child_label", "point_num")) %>% select (-c(X1))

table(adventitia_cells_only$genotype, adventitia_cells_only$child_type)


#number of adventitial CD206 macs with at least 3 cells of interest within an expanding area: Figure 5c,d
j <- NULL
cell_type_list_3 <- c("Activated_Endothelium","Other_Myeloid_Cells")
cluster_size <- 3

for (j in cell_type_list_3) {
  dist3 <- 25
  combined_num_of_clustered_cells_by_genotype <- NULL
  
  while (dist3 < 1000) {
    cells_within_dist3 <- subset(adventitia_cells_only, adventitia_cells_only[ , 5] < dist3)
    cells_types_within_dist3 <- filter(cells_within_dist3, cells_within_dist3$child == j)
    
    number_of_cells_within_dist3 <- cells_types_within_dist3 %>% count(label, point_num)
    
    clusters <- number_of_cells_within_dist3 %>% 
      filter(n > cluster_size)
    
    unique_cells_counted3 <- inner_join(clusters, cells_types_within_dist3, by = c("label","point_num"))
    
    remove_repeats3 <- unique_cells_counted3 %>%
      distinct(label, point_num, .keep_all = TRUE)
    
    num_of_clusters3 <- remove_repeats3 %>% count(genotype)
    num_of_clusters3$distance <- dist3
    
    combined_num_of_clustered_cells_by_genotype <- bind_rows(combined_num_of_clustered_cells_by_genotype, num_of_clusters3)
    dist3 <- dist3 + 25
  }
  
  print(j)
  
  Endo_clusters <- ggplot(combined_num_of_clustered_cells_by_genotype, aes(x = distance, y = n, color = genotype)) +
    geom_point() +
    geom_line() +
    ggtitle(paste(j, "- Number of Cells Per Cluster: ", cluster_size))
  
  Endo_clusters
  
  ggsave(Endo_clusters, file=paste0(".../ClustersWithinRadius/", j,"_clustersize_new_",cluster_size,".png"), width = 14, height = 10, units = "cm")  
}

write.csv(combined_num_of_clustered_cells_by_genotype, ".../Other_Myeloid_Table.csv")


