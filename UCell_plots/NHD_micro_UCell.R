### GSE199378

library(tidyverse)
library(Seurat)
library(sctransform)
library(harmony)
library(gplots)
library(ggplot2)
BiocManager::install("glmGamPoi")
library(glmGamPoi)

dir.create("~/GSE199378/analysis")
setwd("~/GSE199378/analysis")

### convert h5ad to seurat
GSM6038245 <- Read10X_h5("~/GSE199378_RAW/GSM6038245_794_002_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
S1 <- CreateSeuratObject(counts = GSM6038245$`Gene Expression`, project = "GSM6038245", min.cells = 3, min.features = 200)
S1 <- AddMetaData(S1, metadata = "LPS+IFNg", col.name = "treatment")

GSM6038247 <- Read10X_h5("~/GSE199378_RAW/GSM6038247_794_003_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
S2 <- CreateSeuratObject(counts = GSM6038247$`Gene Expression`, project = "GSM6038247", min.cells = 3, min.features = 200)
S2 <- AddMetaData(S2, metadata = "IL4+IL13", col.name = "treatment")

GSM6038249 <- Read10X_h5("~/GSE199378_RAW/GSM6038249_794_004_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
S3 <- CreateSeuratObject(counts = GSM6038249$`Gene Expression`, project = "GSM6038249", min.cells = 3, min.features = 200)
S3 <- AddMetaData(S3, metadata = "IL10", col.name = "treatment")

GSM6038251 <- Read10X_h5("~/GSE199378_RAW/GSM6038251_794_006_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
S4 <- CreateSeuratObject(counts = GSM6038251$`Gene Expression`, project = "GSM6038251", min.cells = 3, min.features = 200)
S4 <- AddMetaData(S4, metadata = "LPS+IFNg", col.name = "treatment")

GSM6038253 <- Read10X_h5("~/GSE199378_RAW/GSM6038253_794_007_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
S5 <- CreateSeuratObject(counts = GSM6038253$`Gene Expression`, project = "GSM6038253", min.cells = 3, min.features = 200)
S5 <- AddMetaData(S5, metadata = "IL4+IL13", col.name = "treatment")

GSM6038255 <- Read10X_h5("~/GSE199378_RAW/GSM6038255_794_008_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
S6 <- CreateSeuratObject(counts = GSM6038255$`Gene Expression`, project = "GSM6038255", min.cells = 3, min.features = 200)
S6 <- AddMetaData(S6, metadata = "IL10", col.name = "treatment")

merged_seurat <- merge(x = S1, y = c(S2,S3,S4,S5,S6), add.cell.ids = c("GSM6038245","GSM6038247","GSM6038249","GSM6038251","GSM6038253","GSM6038255"))

# calculate mitoRatio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# QC
merged_seurat@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

merged_seurat@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 1500)

merged_seurat@meta.data %>%   
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 1000)

merged_seurat@meta.data %>% 
  ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 0.2)

merged_seurat <- merged_seurat %>%
  subset(
    subset = (nCount_RNA > 1500 & nCount_RNA < 30000) &
      (nFeature_RNA > 1000) & 
      (mitoRatio < 0.2)
  )

save(merged_seurat, file = "~/GSE199378/analysis/GSE199378_preprocess.Robj")

# SCTransform
merged_seurat <- SCTransform(merged_seurat, method = "glmGamPoi", vars.to.regress = c("mitoRatio"), verbose = FALSE)

# clustering
merged_seurat <- RunPCA(merged_seurat, verbose = FALSE)
ElbowPlot(merged_seurat, ndims = 50)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20, verbose = FALSE)

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20, verbose = FALSE)
merged_seurat <- FindClusters(merged_seurat, resolution = c(0.6,0.8,1,2), verbose = FALSE)
Idents(merged_seurat) <- "SCT_snn_res.0.6"
DimPlot(merged_seurat, label = TRUE)

Idents(merged_seurat) <- "SCT_snn_res.0.6"
merged_seurat <- subset(merged_seurat, idents = c(11,12), invert = T)
merged_seurat <- RunPCA(merged_seurat, verbose = FALSE)
ElbowPlot(merged_seurat, ndims = 50)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20, verbose = FALSE)

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20, verbose = FALSE)
merged_seurat <- FindClusters(merged_seurat, resolution = c(0.6,0.8,1,2), verbose = FALSE)
Idents(merged_seurat) <- "SCT_snn_res.0.6"
DimPlot(merged_seurat, label = TRUE)
DimPlot(merged_seurat, label = TRUE, split.by = "treatment")

save(merged_seurat, file = "~/GSE199378/analysis/GSE199378_std.Robj")

DefaultAssay(merged_seurat) <-"RNA"
merged_seurat <- NormalizeData(merged_seurat)
Idents(merged_seurat) <- "treatment"
markers <- FindAllMarkers(merged_seurat)
write.table(markers,"~/GSE199378/analysis/markers_treatment_RNA.txt", sep= "\t")

markers <- read.table("~/GSE199378/analysis/markers_treatment_RNA.txt", sep= "\t")

# save UMAP
pdf("~/GSE199378/analysis/GSE199378_umap.pdf", width = 10, height = 5)
DimPlot(merged_seurat, label = TRUE, cols = c("#cab2d6", "#99c6fd", "#fdada2"), split.by = "treatment", pt.size = 0.01)
dev.off()

pdf("~/GSE199378/analysis/GSE199378_umap_treatment.pdf", width = 5, height = 5)
DimPlot(merged_seurat, label = TRUE, cols = c("#cab2d6", "#99c6fd", "#fdada2") %>% rev(), group.by = "treatment", pt.size = 0.01)
dev.off()

# dot plot
DefaultAssay(merged_seurat) <- "RNA"
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pdf("~/GSE199378/analysis/GSE199378_dotplot_top10.pdf", width = 15, height = 5)
DotPlot(merged_seurat, features = top10$gene) + RotatedAxis() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

for (i in unique(top10$cluster)) {
  genes <- top10 %>% filter(cluster == i)
  if (i=="LPS+IFNg") {
    pal <- "#9d66fa"
  } else if (i=="IL4+IL13") {
    pal <- "#6689fa"
  } else {
    pal <- "#fb9a99"
  }
  p <- DotPlot(merged_seurat, features = genes$gene) + RotatedAxis() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_color_gradient(low = "#d9d9d9", high = pal)
  pdf(paste0("~/GSE199378/analysis/GSE199378_dotplot_",i,".pdf"), width = 7, height = 5)
  print(p)
  dev.off()
}






### calculate geneset score using UCell
library(UCell)
### polarization signature on its own dataset
#read table of gene lists, fill empty cells with NA
M_sig <- read.table("~/GSE199378/analysis/markers_treatment_RNA.txt", sep= "\t")

LPS_IFNg <- M_sig %>% filter(cluster == "LPS+IFNg") %>% top_n(wt=avg_log2FC, n=100)
IL4_IL13 <- M_sig %>% filter(cluster == "IL4+IL13") %>% top_n(wt=avg_log2FC, n=100)
IL10 <- M_sig %>% filter(cluster == "IL10") %>% top_n(wt=avg_log2FC, n=100)

#take gene names from table, omit NA cells
markers <- list()
markers$LPS_IFNg <- na.omit(as.vector(LPS_IFNg$gene))
markers$IL4_IL13 <- na.omit(as.vector(IL4_IL13$gene))
markers$IL10 <- na.omit(as.vector(IL10$gene))

#add module score to micro subset
DefaultAssay(merged_seurat) <- "RNA"
merged_seurat <- AddModuleScore_UCell(obj = merged_seurat, features = markers)
signature.names <- paste0(names(markers), "_UCell")

VlnPlot(merged_seurat, features = signature.names, split.by = "treatment")

### adding p-value to violin plots
library(ggpubr)
# https://www.biostars.org/p/458261/
Idents(merged_seurat) <- "treatment"
levels(merged_seurat)
merged_seurat@meta.data$treatment <- rep(c("LPS+IFNg","IL4+IL13","IL10"), times = c(7000, 7000, 7863))

vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(merged_seurat, features = signature, 
            pt.size = 0, cols = c("#cab2d6", "#99c6fd", "#fdada2"),
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + geom_boxplot(width=0.3, fill="white") + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + stat_compare_means(comparisons = test_sign, label = "p.format")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 0.1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_p.pdf")
  ggsave(file_name, width = 7, height = 6)
}

gene_sig <- c("LPS_IFNg_UCell")
comparisons <- list(c("LPS+IFNg","IL4+IL13"),c("IL4+IL13","IL10"),c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/LPS_IFNg_UCell_RNA", test_sign = comparisons)

gene_sig <- c("LPS_IFNg_UCell")
comparisons <- list(c("IL4+IL13","IL10"),c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/LPS_IFNg_UCell1_RNA", test_sign = comparisons)

gene_sig <- c("LPS_IFNg_UCell")
comparisons <- list(c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/LPS_IFNg_UCell2_RNA", test_sign = comparisons)

gene_sig <- c("IL4_IL13_UCell")
comparisons <- list(c("LPS+IFNg","IL4+IL13"),c("IL4+IL13","IL10"),c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/IL4_IL13_UCell_RNA", test_sign = comparisons)

gene_sig <- c("IL4_IL13_UCell")
comparisons <- list(c("IL4+IL13","IL10"),c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/IL4_IL13_UCell1_RNA", test_sign = comparisons)

gene_sig <- c("IL4_IL13_UCell")
comparisons <- list(c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/IL4_IL13_UCell2_RNA", test_sign = comparisons)

gene_sig <- c("IL10_UCell")
comparisons <- list(c("LPS+IFNg","IL4+IL13"),c("IL4+IL13","IL10"),c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/IL10_UCell_RNA", test_sign = comparisons)

gene_sig <- c("IL10_UCell")
comparisons <- list(c("IL4+IL13","IL10"),c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/IL10_UCell1_RNA", test_sign = comparisons)

gene_sig <- c("IL10_UCell")
comparisons <- list(c("LPS+IFNg","IL10"))
vp_case1(gene_signature = gene_sig, file_name = "~/GSE199378/analysis/IL10_UCell2_RNA", test_sign = comparisons)





### VENN diagram of mac signature
M_sig <- read.table("~/GSE199378/analysis/markers_treatment.txt", sep= "\t")

LPS_IFNg <- M_sig %>% filter(cluster == "LPS+IFNg") %>% top_n(wt=avg_log2FC, n=100)
IL4_IL13 <- M_sig %>% filter(cluster == "IL4+IL13") %>% top_n(wt=avg_log2FC, n=100)
IL10 <- M_sig %>% filter(cluster == "IL10") %>% top_n(wt=avg_log2FC, n=100)

library(BioVenn)
setwd("~/GSE199378/analysis")
list_x <- na.omit(as.vector(LPS_IFNg$gene))
list_y <- na.omit(as.vector(IL4_IL13$gene))
list_z <- na.omit(as.vector(IL10$gene))
biovenn <- draw.venn(list_x, list_y, list_z, subtitle="Venn", xtitle = "LPS+IFNg", ytitle = "IL4+IL13", ztitle = "IL10",
                     x_c = "#cab2d6", y_c = "#99c6fd", z_c = "#fdada2",
                     nrtype="abs", output = "pdf", filename = "mac_sig_venn.pdf")
biovenn$xy_only
# character(0)
biovenn$yz_only
# [1] "F13A1"  "GPNMB"  "FUCA1"  "MS4A6A" "FTL"
