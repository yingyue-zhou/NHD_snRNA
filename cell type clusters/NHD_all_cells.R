### analysis of NHD, integration, all cells
### integration

rm(list = ls())

# load packages
library(tidyverse)
library(Seurat)
library(sctransform)
library(harmony)
library(gplots)
library(ggplot2)
library(dplyr)


load(file = "~/NHD_preprocess.Robj")

# SCTransform
merged_seurat <- SCTransform(merged_seurat, method = "glmGamPoi", vars.to.regress = c("mitoRatio","batch"), verbose = FALSE)

DefaultAssay(merged_seurat) <- "SCT"
# These are now standard steps in the Seurat workflow for visualization and clustering
merged_seurat <- RunPCA(merged_seurat, verbose = FALSE)
ElbowPlot(merged_seurat, ndims = 50)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30, verbose = FALSE)

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30, verbose = FALSE)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.8, verbose = FALSE)
DimPlot(merged_seurat, label = TRUE, split.by = "orig.ident")

# reorder
Idents(merged_seurat) <- "sample2"
levels(merged_seurat) <- c("Control1","Control2","Control3","Control4","Control5",
                           "Control6","Control7","Control8","Control9","Control10","Control11",
                           "NHD1","NHD2","NHD3")

# save UMAP before Harmony
library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(12,"Paired"))(14)
pdf("~/NHD__before_harmony_sample.pdf", width = 7, height = 5)
DimPlot(merged_seurat, cols = pal, group.by = "sample2")
dev.off()


# Harmony
merged_seurat <- merged_seurat %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, assay.use="SCT")

merged_seurat <- merged_seurat %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = c(0.8,1,2,3,4,5))

Idents(merged_seurat) <- "SCT_snn_res.5"
DimPlot(object = merged_seurat, reduction = "umap", label = T)

# clean data, remove doublets
# check cell type marker genes
FeaturePlot(merged_seurat, features = c("PVALB", "SST", "VIP", "RELN", "FGF13"))
FeaturePlot(merged_seurat, features = c("SYT1", "SLC17A7", "GAD1"))
FeaturePlot(merged_seurat, features = c("MOG", "CSF1R", "AQP4", "VCAN"))
FeaturePlot(merged_seurat, features = c("CEMIP", "CALD1", "MGP", "FLT1"))
FeaturePlot(merged_seurat, features = c("CD163", "F13A1", "NKG7", "FGF13"))
FeaturePlot(merged_seurat, features = c("TRPC4", "TAGLN", "NKG7", "F13A1"))
FeaturePlot(merged_seurat, features = c("FSTL5"), label = T)
FeaturePlot(merged_seurat, features = c("nFeature_RNA"))

# check cluster markers
sample.seurat <- merged_seurat[, sample(colnames(merged_seurat), size = 4000, replace=F)]
clu13 <- FindMarkers(sample.seurat, ident.1 = 13)

# remove doublets neurons，oligos (57,66,74,75,79,82,83,85,86)
Idents(merged_seurat) <- "SCT_snn_res.5"
merged_seurat <- subset(merged_seurat, idents = c(57,66,74,75,79,82,83,85,86), invert=T)

merged_seurat <- merged_seurat %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)

merged_seurat <- merged_seurat %>% FindClusters(resolution = c(1,1.2,1.6,2,2.4,3))
Idents(merged_seurat) <- "SCT_snn_res.3"
DimPlot(object = merged_seurat, reduction = "umap", label = T)

save(merged_seurat, file = "~/NHD_harmony.Robj")

# save umap after Harmony for QC
library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(12,"Paired"))(14)
pdf("~/umap_after_harmony_sample.pdf", width = 7, height = 5)
DimPlot(object = merged_seurat, reduction = "umap", cols = pal, group.by = "sample2", pt.size = 0.01, shuffle = T)
dev.off()

pal <- colorRampPalette(brewer.pal(8,"Set2"))(8)
pdf("~/umap_after_harmony_batch.pdf", width = 7, height = 5)
DimPlot(object = merged_seurat, reduction = "umap", cols = pal, group.by = "batch", pt.size = 0.01, shuffle = T)
dev.off()

pal <- colorRampPalette(brewer.pal(8,"Pastel1"))(8)
pdf("~/umap_after_harmony_sex.pdf", width = 7, height = 5)
DimPlot(object = merged_seurat, reduction = "umap", cols = pal, group.by = "sex", pt.size = 0.01, shuffle = T)
dev.off()


# cell type annotation
Idents(merged_seurat) <- "SCT_snn_res.3"
DimPlot(object = merged_seurat, reduction = "umap", label = T)
merged_seurat$celltype1 <- NA
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 %in% c(5,6,11,12,13,15,16,18,19,21,22,24,25,27,29,34,35,36,38,39,40,46,52,54))] <- "Ex_neuron"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 %in% c(2,10,20,26,32,42,45))] <- "In_neuron"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 %in% c(0,3,7,8,9,30,31,53))] <- "Oligo"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 %in% c(4,14,23,28,37,51))] <- "Astro"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 %in% c(1,44,48))] <- "Micro"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 == 33)] <- "Endo"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 %in% c(17,41))] <- "OPC"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 == 49)] <- "SMC"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 == 47)] <- "Peri"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 == 43)] <- "Fibro"
merged_seurat$celltype1[which(merged_seurat$SCT_snn_res.3 == 50)] <- "IMM"
DimPlot(object = merged_seurat, reduction = "umap", label = T, group.by = "celltype1")

# find cell type markers
DefaultAssay(merged_seurat) <- "RNA"
merged_seurat <- NormalizeData(merged_seurat)
Idents(merged_seurat) <- "celltype1"
markers <- FindAllMarkers(merged_seurat, test.use = "MAST", latent.vars = c("age","sex","PMI"))
write.table(markers,"~/markers_MAST.txt", sep= "\t")

# re-order cluster for umap plot
Idents(merged_seurat) <- "celltype1"
levels(merged_seurat) <- c("Ex_neuron","In_neuron","Oligo","OPC","Astro","Micro","Endo","SMC","Peri","Fibro","IMM")
pal <- c("#1f78b4","#a6cee3","#33a02b","#b2df8a","#e3211c","#fb9a99","#b15928","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#fddaec")

# save UMAP
pdf("~/umap_celltype1.pdf", width = 7, height = 5)
DimPlot(object = merged_seurat, reduction = "umap", label = T, cols = pal, pt.size = 0.01)
dev.off()

# cell type marker heatmap
Idents(merged_seurat) <- "celltype1"
levels(merged_seurat) <- c("Ex_neuron","In_neuron","Oligo","OPC","Astro","Micro","Endo","SMC","Peri","Fibro","IMM")

pdf("~/heatmap_marker.pdf", width = 7, height = 5)
DoHeatmap(merged_seurat, features = c("SYT1","RBFOX1","GAD1","SST","VIP","PVALB","LAMP5","MOBP","MBP","PLP1","PDGFRA","VCAN","SOX6",
                                      "AQP4","GJA1","GFAP","CSF1R","C3","FYB1","FLT1","CLDN5","TAGLN","ABCC9","CEMIP","NKG7","CD3E"),
          group.colors = c("#1f78b4","#a6cee3","#33a02b","#b2df8a","#e3211c","#fb9a99","#b15928","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#fddaec"))
dev.off()


# stats
nuclei.number <- table(merged_seurat$celltype1, merged_seurat$sample2)
write.table(nuclei.number,"~/nuclei_number_sample2.txt", sep= "\t")

# nuclei distribution
nuclei.number <- table(merged_seurat$sample2)
nuclei.number <- as.data.frame(nuclei.number)
names(nuclei.number) <- c("sample2","Freq")
nuclei.number$sample2 <- factor(nuclei.number$sample2, levels = c("Control1","Control2","Control3","Control4","Control5",
                                                                  "Control6","Control7","Control8","Control9","Control10","Control11",
                                                                  "NHD1","NHD2","NHD3") %>% rev()) 
ggplot(data=nuclei.number, aes(x=sample2, y=Freq)) +
  geom_bar(stat="identity", fill = "#8c6bb1") + 
  coord_flip() + 
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(filename = "~/nuclei_number.pdf",device = "pdf", height = 11, width = 11)

# median of gene and UMI counts
median <- merged_seurat@meta.data %>%
  dplyr::select(nCount_RNA,nFeature_RNA) %>%
  summarise(nUMI=median(nCount_RNA), nGene=median(nFeature_RNA))
### nUMI6943, nGene3150

# UMI distribution
UMI <- merged_seurat@meta.data %>%
  dplyr::select(sample2, nCount_RNA) %>%
  group_by(sample2) %>%
  summarise(nUMI=median(nCount_RNA))
UMI <- as.data.frame(UMI)
names(UMI) <- c("sample2","Freq")
UMI$sample2 <- factor(UMI$sample2, levels = c("Control1","Control2","Control3","Control4","Control5",
                                              "Control6","Control7","Control8","Control9","Control10","Control11",
                                              "NHD1","NHD2","NHD3") %>% rev()) 
ggplot(data=UMI, aes(x=sample2, y=Freq)) +
  geom_bar(stat="identity", fill = "#8c6bb1") + 
  coord_flip() + 
  ylab("UMI") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(filename = "~/UMI_number.pdf",device = "pdf", height = 11, width = 11)

# nGene distribution
nGene <- merged_seurat@meta.data %>%
  dplyr::select(sample2, nFeature_RNA) %>%
  group_by(sample2) %>%
  summarise(nUMI=median(nFeature_RNA))
nGene <- as.data.frame(nGene)
names(nGene) <- c("sample2","Freq")
nGene$sample2 <- factor(nGene$sample2, levels = c("Control1","Control2","Control3","Control4","Control5",
                                                  "Control6","Control7","Control8","Control9","Control10","Control11",
                                                  "NHD1","NHD2","NHD3") %>% rev()) 
ggplot(data=nGene, aes(x=sample2, y=Freq)) +
  geom_bar(stat="identity", fill = "#8c6bb1") + 
  coord_flip() + 
  ylab("nGene") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(filename = "~/Gene_number.pdf",device = "pdf", height = 11, width = 11)


save(merged_seurat, file = "~/NHD_harmony.Robj")





# violin plot for FGF2 expression
DefaultAssay(merged_seurat) <- "RNA"
merged_seurat <- NormalizeData(merged_seurat)
Idents(merged_seurat) <- "celltype1"
p <- VlnPlot(merged_seurat, features = "FGF2", cols = c("#a6cee3", "#fb9a99"), split.by = "disease", pt.size = 0) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", position = position_dodge(0.9), color = "#081d58") + 
  stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 3, method = "wilcox") + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
pdf("~/FGF2_RNA.pdf", width = 11, height = 5)
print(p)
dev.off()

# Feature plot for TYROBP expression
DefaultAssay(merged_seurat) <- "RNA"
pdf("~/TYROBP.pdf", width = 5, height = 5)
FeaturePlot(merged_seurat, features = "TYROBP", pt.size = 0.01)
dev.off()



