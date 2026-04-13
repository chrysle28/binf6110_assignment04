### ===============================
## BINF110 Assignment 4
## Analysis of mice scRNA sequencing dataset

## Kenneth Gamueda
## 2026-04-13

### === PACKAGES USED ===========
library(tidyverse)
library(dplyr)
library(ggplot2)
library(Seurat)
library(glmGamPoi)
library(SingleR)
library(celldex)
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)

### === 1 | QUALITY CONTROL =======
seurat <- LoadSeuratRds("seurat_ass4.rds")
colnames(seurat@meta.data)
seurat # 156572 samples

# Add mitochondrial percent and visualize quality metrics
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
plot_QC_pre <- VlnPlot(seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
) &
  xlab(NULL) &
  plot_annotation(
    title = "Quality Metrics Pre-Filtering",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )

plot_QC_pre[[1]] <- plot_QC_pre[[1]] + labs(
  title = "Unique Genes",
  y = "Number of Unique Genes"
)
plot_QC_pre[[2]] <- plot_QC_pre[[2]] + labs(
  title = "Total Molecules",
  y = "Number of Molecules"
)
plot_QC_pre[[3]] <- plot_QC_pre[[3]] + labs(
  title = "Mitochondrial Percent",
  y = "% mitochondria"
)
plot_QC_pre

# Visualize correlation of features
FeatureScatter(seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  raster = F
) +
  ggtitle("nCount_RNA vs Mitochondrial %") # uncorrelated

FeatureScatter(seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

FeatureScatter(seurat,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt",
  raster = F
) +
  ggtitle("nFeature_RNA vs Mitochondrial %") # uncorrelated

# Subset based on visualizations
seurat <- subset(seurat, subset = nFeature_RNA > 750 & percent.mt < 10)
plot_QC_post <- VlnPlot(seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
) &
  xlab(NULL) &
  plot_annotation(
    title = "Quality Metrics Post-Filtering",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )

plot_QC_post[[1]] <- plot_QC_post[[1]] + labs(
  title = "Unique Genes",
  y = "Number of Unique Genes"
)
plot_QC_post[[2]] <- plot_QC_post[[2]] + labs(
  title = "Total Molecules",
  y = "Number of Molecules"
)
plot_QC_post[[3]] <- plot_QC_post[[3]] + labs(
  title = "Mitochondrial Percent",
  y = "% mitochondria"
)
plot_QC_post
seurat # 156572 samples --> 148381 samples

# Combined QC plots
wrap_elements(plot_QC_pre) / wrap_elements(plot_QC_post)

# Run SCTransform for normalization and scaling
seurat <- SCTransform(seurat, ncells = 2000, conserve.memory = T) # default of 3000 variable features
seurat

# Run PCA on transformed Seurat object and determine optimal k
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
ElbowPlot(seurat, ndims = 50) +
  plot_annotation(
    title = "Elbow Plot of Principal Components",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )
# k = 33

seurat <- FindNeighbors(seurat, dims = 1:33)
seurat <- FindClusters(seurat, resolution = 0.6)

# Create UMAP cluster plot
seurat <- RunUMAP(seurat, dims = 1:33)
plot_umap <- DimPlot(seurat,
  reduction = "umap",
  label = T,
  label.size = 4.5,
  label.box = T,
  repel = T,
  raster = F
) +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  plot_annotation(
    title = "Unannotated UMAP Clustering (Resolution = 0.6)",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )
plot_umap


# Save Seurat object with SCT assay and UMAP reduction
saveRDS(seurat, file = "seurat_UMAP.rds")

# Dimensional plots to check if any variable needs to be regressed out
DimPlot(seurat, group.by = "orig.ident")
DimPlot(seurat, group.by = "biosample_id")

### === 2.1 | AUTOMATIC ANNOTATION =======
seurat <- readRDS("seurat_UMAP.rds")

# Prepare celldex reference datasets
ref_mrs <- celldex::MouseRNAseqData()
ref_imm <- celldex::ImmGenData()

# Automatic annotation with SingleR and celldex references
ann_mrs <- SingleR(
  test = as.SingleCellExperiment(seurat),
  ref = ref_mrs,
  labels = ref_mrs$label.fine,
  de.method = "wilcox"
)

ann_imm <- SingleR(
  test = as.SingleCellExperiment(seurat),
  ref = ref_imm,
  labels = ref_imm$label.main,
  de.method = "wilcox"
)

table(ann_mrs$labels)
table(ann_imm$labels)

# Diagnostic plots
plotScoreHeatmap(ann_mrs,
  main = "Diagnostic Heatmap for Annotation with MouseRNAseq Data"
)

plotScoreHeatmap(ann_imm,
  main = "Diagnostic Heatmap for Annotation with ImmGen Data"
)


### === 2.2 | MANUAL ANNOTATION =======
markers_all <- FindAllMarkers(seurat,
  only.pos = T,
  logfc.threshold = 0.25,
  min.pct = 0.25,
  test.use = "wilcox"
)

markers_top <- markers_all %>%
  group_by(cluster) %>%
  filter(!grepl("^Gm", gene)) %>% # uncharacterized gene models
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 5)
markers_top

write.csv(markers_all, "markers_all.csv", row.names = F)
write.csv(markers_top, "markers_top.csv", row.names = F)

# Manual annotation of markers using supplementary material from Kazer et al. (2025)
markers_ppr <- c(
  "0"  = "Calb2+ Olfactory Sensory Neurons", # B830017H08Rik, Calb2, Kirrel3, Nqo1, Rims3
  "1"  = "Olfactory Sensory Neurons", # Gramd1c, Scgn
  "2"  = "Resting Basal Cells", # Acaa1b, Serpinb5, Krt15, Anxa8, Defb1
  "3"  = "Fcrls+Trem2+ Macrophages", # C1qa, C1qb, Ms4a7, C1qc, Slamf9
  "4"  = "Dlg2+ Olfactory Sensory Neurons", # 5730403I07Rik, Kcnmb3, Lrrc3b, Kirrel2, Chgb
  "5"  = "Dlg2+ Olfactory Sensory Neurons", # Hectd2os, S100a5, Kcnmb3, Lrrc3b
  "6"  = "B Cells", # Ighd
  "7"  = "Gzma++Cma1+ Natural Killer Cells", # Ncr1, Prf1, Klrb1c, Xcl1, Gzma
  "8"  = "Fabp4+ Capillary Cells", # Gpihbp1, Rbp7, Emcn, Btnl9
  "9"  = "IFN-Stim Monocytes", # Ifi205, Cd209a, Phf11a, Cxcl9, Ms4a4c
  "10" = "Differentiating Olfactory Sensory Neurons", # Gap43, Stk32a, Ly6h, Prph, Insm1
  "11" = "Non-Classical Monocytes", # Adgre4, Treml4, Ear2
  "12" = "Fibroblasts", # Dpep1, Apod, Ecrg4, Dpt, Hhip
  "13" = "Dlg2+ Olfactory Sensory Neurons", # Dlg2, S100a5, Pcp4l1, Nphs1
  "14" = "Igfbp2+Nrcam+ Basal Cells", # Car10, Abca4, Nrcam, Vit, Slc47a1
  "15" = "Mature Neutrophils", # Siglecf, Il1r2, Cxcr2, Csf3r, Acod1
  "16" = "Fcrls+Trem2+ Macrophages", # Fcrls, P2ry12, Pf4, Gpr34, Trem2
  "17" = "Dlg2+ Olfactory Sensory Neurons", # Dlg2, Hectd2os, S100a5, Pcp4l1, Nphs1
  "18" = "Gp2+Lyz2+ Goblet/Secretory Cells", # Gp2, Reg3g, Serpinb11, Gabrp, Kcnj16
  "19" = "Cxcl17+Ccl9+ Serous Cells", # Tac1, Wfdc18, Barx2, Cited4, Lcn11
  "20" = "Lyz2+Cpxm2+ Ionocytes", # Kl, Smbd1, Clcnka, Ascl3, Galnt13
  "21" = "B Cell Progenitors/Immature B Cells", # Rag1, Vpreb3, Iglc1, Fam129c, Spib
  "22" = "Neural Progenitors", # Neurod1, Neurog1, Grp, St18, Scg2
  "23" = "Pparg+Muc5b+ Goblet/Secretory Cells", # Chil6, Slc22a20, Vmo1, Ggt1, C2cd4b
  "24" = "Mural Cells (Smooth Muscle/Pericyte)", # Map3k7cl, Higd1b, Pln, Myh11, Des
  "25" = "Camp+Il1b+ Immature Neutrophils", # Retnlg, Mmp8, Ceacam10, Fpr1
  "26" = "Sustentacular Cells", # Muc2, Gldn, Sec14l3, Nipal1
  "27" = "Ciliated Cells", # Bpifa1, Tmem212, Dthd1, Sntn, Fam216b
  "28" = "Glandular Cells", # Bpifb9a, Bpifb9b, Bpifb5, 2310003L06Rik, Odam
  "29" = "Tuft Cells", # Mpz, Gpr37l1, Hmx3, Fabp7, Il25
  "30" = "Camp++Lta4h+ Immature Neutrophils", # Fcnb, Orm1, Camp, Ngp
  "31" = "Scgb-b27+Cck+ Nasal Epithelial Cells", # Car6, Scgb2b27, Scgb1b27, Csn3, Fbp1
  "32" = "Cycling Cells", # Hist1h1a, Hist1h1b, Hist1h2ab, Hist1h3c, Hist1h2bm
  "33" = "Bglap+ Osteoblasts", # Bglap, Bglap2, Ibsp, Col22a1, Car3
  "34" = "OBP-Related Cells", # 5430402E10Rik, Obp1b, Mup4, Obp1a, Lcn11
  "35" = "Plasmacytoid Dendritic Cells", # Ccr9, Cd300c, Cox6a2, Siglech, Cd209d
  "36" = "Dapl1+Pglyrp1+ Endothelial Cells", # Dnase2b, Paqr5, Mal, Mfsd2a, Abcb1a
  "37" = "Vomeronasal Neurons", # Pvalb, Ancv1r, Ankrd63, Cacng2, Calr4
  "38" = "Gpx6+Ces1a+ Goblet/Secretory Cells", # Bpifb6, Ces1a, Bpifb4, Ptgds, 5430419D17Rik
  "39" = "Krt13+Il1a+ Epithelial Cells", # Lce3a, Crct1, Clca4a, Csta1, Cysrt1
  "40" = "Meg3+MHC-II+ Epithelial Cells" # Fezf2, Ctf2, Otop1, Gulo
)

# Manual annotation of markers using CellMarker 2.0 and PanglaoDB
markers_web <- c(
  "0"  = "Spiral ganglion neurons", # Calb2, Kirrel3, Nqo1, Rims3
  "1"  = "Olfactory epithelial cells", # Gramd1c, Scgn
  "2"  = "Basal keratinocytes", # Acaa1b, Serpinb5, Krt15, Anxa8, Defb1
  "3"  = "Macrophages", # C1qa, C1qb, Ms4a7, C1qc, Slamf9
  "4"  = "Neurons", # 5730403I07Rik, Kcnmb3, Lrrc3b, Kirrel2, Chgb
  "5"  = "Neurons", # Hectd2os, S100a5, Kcnmb3, Lrrc3b
  "6"  = "B cells", # Ighd
  "7"  = "Natural killer cells", # Ncr1, Prf1, Klrb1c, Xcl1, Gzma
  "8"  = "Endothelial cells ", # Gpihbp1, Rbp7, Emcn, Btnl9
  "9"  = "Macrophages", # Ifi205, Cd209a, Phf11a, Cxcl9, Ms4a4c
  "10" = "Neurons", # Gap43, Stk32a, Ly6h, Prph, Insm1
  "11" = "Macrophages", # Adgre4, Treml4, Ear2
  "12" = "Fibroblasts", # Dpep1, Apod, Ecrg4, Dpt, Hhip
  "13" = "Epithelial cells", # Dlg2, S100a5, Pcp4l1, Nphs1
  "14" = "Neurons", # Car10, Abca4, Nrcam, Vit, Slc47a1
  "15" = "Neutrophils", # Siglecf, Il1r2, Cxcr2, Csf3r, Acod1
  "16" = "Microglial cells", # Fcrls, P2ry12, Pf4, Gpr34, Trem2
  "17" = "Epithelial cells", # Dlg2, Hectd2os, S100a5, Pcp4l1, Nphs1
  "18" = "Secretory epithelial cells", # Gp2, Reg3g, Serpinb11, Gabrp, Kcnj16
  "19" = "Epithelial cells", # Tac1, Wfdc18, Barx2, Cited4
  "20" = "Ionocytes", # Kl, Smbd1, Clcnka, Ascl3, Galnt13
  "21" = "B cells", # Rag1, Vpreb3, Iglc1, Fam129c, Spib
  "22" = "Neural precursor cells", # Neurod1, Neurog1, Grp, St18, Scg2
  "23" = "Epithelial cells", # Vmo1, Ggt1, C2cd4b
  "24" = "Mural cells", # Map3k7cl, Higd1b, Pln, Myh11, Des
  "25" = "Neutrophils", # Retnlg, Mmp8, Ceacam10, Fpr1
  "26" = "Goblet cells", # Muc2, Sec14l3
  "27" = "Secretory epithelial cells", # Bpifa1, Tmem212, Dthd1, Sntn, Fam216b
  "28" = "Glandular cells", # Pigr, Wfdc18, Qsox1, Aqp5
  "29" = "Glial cells", # Mpz, Gpr37l1, Hmx3, Fabp7,
  "30" = "Neutrophils", # Fcnb, Orm1, Camp, Ngp
  "31" = "Salivary mucous cells", # Car6, Scgb2b27, Scgb1b27
  "32" = "White blood cells", # Hist1h1a, Hist1h1b, Hist1h2ab, Hist1h3c, Hist1h2bm
  "33" = "Osteoblasts", # Bglap, Bglap2, Ibsp, Col22a1, Car3
  "34" = "OBP Cells",
  "35" = "Plasmacytoid dendritic cells", # Ccr9, Cd300c, Cox6a2, Siglech, Cd209d
  "36" = "Endothelial cells", # Paqr5, Mal, Mfsd2a, Abcb1a
  "37" = "Neurons", # Ankrd63, Cacng2, Calr4
  "38" = "Oligodendrocytes", # Ptgds
  "39" = "Keratinocytes", # Crct1, Csta1, Cysrt1
  "40" = "Neurons" # Fezf2
)

# Finding markers for unannotated clusters with web sources
Idents(seurat) <- "seurat_clusters"
markers_28 <- FindMarkers(seurat,
  ident.1 = 28,
  min.pct = 0.25
)
markers_28 %>%
  filter(!(row.names(.) %in% c("Bpifb9a", "Bpifb9b", "Bpifb5", "2310003L06Rik", "Odam"))) %>%
  head(n = 5)
# Cabs1: germ cell
# Pigr: glandular cells
# Wfdc18: glandular cells
# Qsox1: glandular cells
# Aqp5: glandular cells
# glandular cells = goblet, secretory, alveolar pneumocytes, Paneth

# Visualizing markers for cluster 34 (ambiguous cluster with markers for different cell types)
markers_34 <- FindMarkers(seurat,
  ident.1 = 34,
  min.pct = 0.25
)
markers_34 %>%
  arrange(desc(avg_log2FC)) %>%
  head(n = 5)

# Violin plot of top marker and feature plots for top 4 markers of cluster 34
VlnPlot(seurat, features = c("Gm14743")) +
  NoLegend()
FeaturePlot(seurat,
  features = c("Gm14743", "5430402E10Rik", "Obp1b", "Mup4"),
  raster = F
) +
  plot_annotation(
    title = "Feature Plots of Top Markers for Cluster 34",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )

# Finding markers for clusters annotated as the same cell type
seurat_OSN <- subset(seurat, idents = c("4", "5", "13", "17"))
markers_OSN <- FindAllMarkers(seurat_OSN,
  only.pos = T,
  logfc.threshold = 0.25,
  min.pct = 0.01,
  test.use = "wilcox"
)
head(markers_OSN)

markers_OSN_top <- markers_OSN %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  filter(!grepl("^Gm", gene)) %>% # uncharacterized gene models
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5)
markers_OSN_top
write.csv(markers_OSN_top, "markers_OSN_top.csv", row.names = F)

# Feature plots of known cluster markers
FeaturePlot(seurat,
  features = c(
    "Omp", # OSNs
    "Epcam", # epithelial
    "Ptprc", # immune
    "Flt1" # endothelial
  ),
  raster = F
) +
  plot_annotation(
    title = "Feature Plots of Lineage Markers",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )

# Feature plots of epithelial cell clusters
FeaturePlot(seurat,
  features = c(
    "Krt5", # basal
    "Foxj1", # ciliated
    "Ltf", # serous
    "Bpifb9b", # glandular
    "Reg3g", # goblet/secretory
    "Cftr", # ionocyte
    "Trpm5", # tuft
    "Sec14l3" # sustentacular
  ),
  raster = F,
  ncol = 4
) +
  plot_annotation(
    title = "Feature Plots of Markers of Epithelial Cell Subtypes",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )

# Adding all labels to Seurat object
seurat$mrs_labels <- ann_mrs$labels[match(rownames(seurat@meta.data), rownames(ann_mrs))]
seurat$imm_labels <- ann_imm$labels[match(rownames(seurat@meta.data), rownames(ann_imm))]
seurat$ppr_labels <- unname(markers_ppr[as.character(seurat$seurat_clusters)])
seurat$web_labels <- unname(markers_web[as.character(seurat$seurat_clusters)])

# UMAP plot (paper annotations)
plot_ppr <- DimPlot(seurat,
  reduction = "umap",
  group.by = "ppr_labels",
  label = T,
  label.size = 3.25,
  label.box = T,
  repel = T,
  raster = F
) +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  ggtitle("Annotated UMAP Clustering (Literature Source)") +
  theme(plot.title = element_text(
    size = 18,
    face = "bold",
    hjust = 0.5
  ))
plot_ppr

# UMAP plot (web annotations)
plot_web <- DimPlot(seurat,
  reduction = "umap",
  group.by = "web_labels",
  label = T,
  label.size = 3.25,
  label.box = T,
  repel = T,
  raster = F
) +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  ggtitle("Annotated UMAP Clustering (Web Source)") +
  theme(plot.title = element_text(
    size = 18,
    face = "bold",
    hjust = 0.5
  ))
plot_web

# Combined plots
plot_ppr + plot_web

# UMAP plot (automatic annotation with MouseRNASeq)
plot_mrs <- DimPlot(seurat,
  reduction = "umap",
  group.by = "mrs_labels",
  label = T,
  label.size = 3.25,
  label.box = T,
  repel = T,
  raster = F
) +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  ggtitle("Annotated UMAP Clustering (MouseRNAseqData)") +
  theme(plot.title = element_text(
    size = 18,
    face = "bold",
    hjust = 0.5
  ))
plot_mrs

# UMAP plot (automatic annotation with ImmGen)
plot_imm <- DimPlot(seurat,
  reduction = "umap",
  group.by = "imm_labels",
  label = T,
  label.size = 3.25,
  label.box = T,
  repel = T,
  raster = F
) +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  ggtitle("Annotated UMAP Clustering (ImmGenData)") +
  theme(plot.title = element_text(
    size = 18,
    face = "bold",
    hjust = 0.5
  ))
plot_imm

# Combined plots
plot_mrs + plot_imm

### === 3 | DE ANALYSIS =======
# Using paper annotations for cell identity
seurat <- RenameIdents(seurat, markers_ppr)
seurat$cell_type <- Idents(seurat)
table(seurat$cell_type)

# Pseudobulk
seurat_bulk <- AggregateExpression(seurat, assays = "RNA", return.seurat = T, group.by = c("mouse_id", "time", "cell_type"))
seurat_bulk$cell_time <- paste(seurat_bulk$cell_type, seurat_bulk$time, sep = "_")

Idents(seurat_bulk) <- "cell_time"
table(Idents(seurat_bulk))

seurat_DE <- FindMarkers(seurat_bulk,
  ident.1 = "IFN-Stim Monocytes_D05",
  ident.2 = "IFN-Stim Monocytes_Naive",
  test.use = "DESeq2"
)

seurat_DE$gene <- rownames(seurat_DE)
seurat_DE$expression <- case_when(
  seurat_DE$p_val_adj < 0.05 & seurat_DE$avg_log2FC > 1 ~ "Upregulated",
  seurat_DE$p_val_adj < 0.05 & seurat_DE$avg_log2FC < -1 ~ "Downregulated",
  TRUE ~ "NS"
)

# Key-value mapping for EnhancedVolcano
color_map <- c("Upregulated" = "darkred", "Downregulated" = "royalblue", "NS" = "gray")
keyvals <- color_map[seurat_DE$expression]
names(keyvals) <- names(color_map)[match(keyvals, color_map)]
keyvals[is.na(keyvals)] <- "gray"

plot_DE <- EnhancedVolcano(seurat_DE,
  lab = rownames(seurat_DE),
  title = "5 Days Post Infection vs Naive (IFN-Stim Monocytes)",
  x = "avg_log2FC",
  y = "p_val_adj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 4,
  cutoffLineType = "twodash",
  cutoffLineWidth = 0.8,
  colCustom = keyvals,
  legendPosition = "top",
  boxedLabels = T,
  drawConnectors = T,
  gridlines.major = F,
  gridlines.minor = T
)
plot_DE
# Up: Ly6c2, Gda, Slfn4
# Down: Rnd3, Klrd1, Cd209d

# Add cell_time column to original Seurat object
seurat$cell_time <- paste(seurat$cell_type, seurat$time, sep = "_")
Idents(seurat) <- "cell_time"
table(seurat$cell_time)

# Create violin plots of DE genes (pseudobulked)
plot_up_bulk <- VlnPlot(seurat_bulk, features <- c("Ly6c2", "Gda", "Slfn4"), idents = c("IFN-Stim Monocytes_D05", "IFN-Stim Monocytes_Naive"), group.by = "time") &
  xlab("Time") &
  plot_annotation(
    title = "Upregulated Genes 5 Days Post Infection (Pseudobulk)",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )
plot_up_bulk

plot_down_bulk <- VlnPlot(seurat_bulk, features <- c("Rnd3", "Klrd1", "Cd209d"), idents = c("IFN-Stim Monocytes_D05", "IFN-Stim Monocytes_Naive"), group.by = "time") &
  xlab("Time") &
  plot_annotation(
    title = "Downregulated Genes 5 Days Post Infection (Pseudobulk)",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )
plot_down_bulk

wrap_elements(plot_up_bulk) / wrap_elements(plot_down_bulk)

# Create violin plots of DE genes (non-pseudobulked)
plot_up <- VlnPlot(seurat, features <- c("Ly6c2", "Gda", "Slfn4"), idents = c("IFN-Stim Monocytes_D05", "IFN-Stim Monocytes_Naive"), group.by = "time") &
  xlab("Time") &
  plot_annotation(
    title = "Upregulated Genes 5 Days Post Infection",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )
plot_up

plot_down <- VlnPlot(seurat, features <- c("Rnd3", "Klrd1", "Cd209d"), idents = c("IFN-Stim Monocytes_D05", "IFN-Stim Monocytes_Naive"), group.by = "time") &
  xlab("Time") &
  plot_annotation(
    title = "Downregulated Genes 5 Days Post Infection",
    theme = theme(plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ))
  )
plot_down

wrap_elements(plot_up) / wrap_elements(plot_down)

### === 4.1 | GSEA =======
gene_conv <- bitr(seurat_DE$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)
seurat_DE <- seurat_DE %>%
  left_join(gene_conv, by = c("gene" = "SYMBOL"))
seurat_DE <- seurat_DE[!is.na(seurat_DE$ENTREZID), ]

# List of ranked genes (metric is log of p-value * log2FC)
seurat_DE$rank <- -log10(seurat_DE$p_val) * sign(seurat_DE$avg_log2FC)
ranked_genes <- seurat_DE$rank
names(ranked_genes) <- seurat_DE$ENTREZID
ranked_genes <- sort(ranked_genes, decreasing = T)

# Labels for plots
sign_labels <- c(
  "activated" = "Upregulated (Day 5)",
  "suppressed" = "Downregulated (Day 5)"
)

# GSEA with GO
gsea_GO <- gseGO(
  geneList = ranked_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  by = "fgsea"
)

plot_gsea_GO <- dotplot(gsea_GO,
  showCategory = 10,
  split = ".sign"
) +
  facet_grid(~.sign, labeller = as_labeller(sign_labels)) +
  ggtitle("GSEA of GO Biological Processes (Day 5 vs Naive)") +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    strip.text = element_text(size = 14)
  )
plot_gsea_GO

# GSEA with KEGG
gsea_KEGG <- gseKEGG(
  geneList = ranked_genes,
  organism = "mmu",
  keyType = "ncbi-geneid",
  by = "fgsea"
)

plot_gsea_KEGG <- dotplot(gsea_KEGG,
  showCategory = 10,
  split = ".sign"
) +
  facet_grid(~.sign, labeller = as_labeller(sign_labels)) +
  ggtitle("GSEA of KEGG Pathways (Day 5 vs Naive)") +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    strip.text = element_text(size = 14)
  )
plot_gsea_KEGG

wrap_elements(plot_gsea_GO) + wrap_elements(plot_gsea_KEGG)

### === 4.2 | ORA =======
sig_genes <- subset(seurat_DE, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
head(sig_genes)

# ORA with GO
ora_GO <- compareCluster(ENTREZID ~ expression,
  data = sig_genes,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  universe = seurat_DE,
  keyType = "ENTREZID",
  ont = "BP"
)

plot_ora_GO <- dotplot(ora_GO,
  showCategory = 10
) +
  ggtitle("ORA of GO Biological Processes (Day 5 vs Naive)") +
  theme(plot.title = element_text(
    size = 18,
    face = "bold",
    hjust = 0.5
  ))
plot_ora_GO

# ORA with KEGG
ora_KEGG <- compareCluster(ENTREZID ~ expression,
  data = sig_genes,
  fun = "enrichKEGG",
  organism = "mmu",
  keyType = "ncbi-geneid",
  universe = seurat_DE
)

plot_ora_KEGG <- dotplot(ora_KEGG,
  showCategory = 10
) +
  ggtitle("ORA of KEGG Pathways (Day 5 vs Naive)") +
  theme(plot.title = element_text(
    size = 18,
    face = "bold",
    hjust = 0.5
  ))
plot_ora_KEGG

wrap_elements(plot_ora_GO) + wrap_elements(plot_ora_KEGG)
