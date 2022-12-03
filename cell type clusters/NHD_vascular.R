### re-cluster of NHD vascular cells

rm(list = ls())

# load packages
library(tidyverse)
library(Seurat)
library(sctransform)
library(harmony)
library(gplots)
library(ggplot2)
library(dplyr)


# Re-cluster of vascular cells
load("~/NHD_harmony.Robj")
Idents(merged_seurat) <- "celltype1"
vascu <- subset(merged_seurat, idents = c("Endo","Fibro","SMC","Peri"))
vascu <- vascu %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, assay.use="SCT")
ElbowPlot(vascu, ndims = 50)
vascu <- vascu %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)
vascu <- vascu %>% FindClusters(resolution = c(0.4,0.6,1,2))
Idents(vascu) <- "SCT_snn_res.2"
DimPlot(object = vascu, reduction = "umap", label = T)

FeaturePlot(vascu, features = c("nFeature_RNA"))


# remove doublets by nFeature (8)
Idents(vascu) <- "SCT_snn_res.2"
vascu <- subset(vascu, idents = c(13), invert = T)
ElbowPlot(vascu, ndims = 50)
vascu <- vascu %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)
vascu <- vascu %>% FindClusters(resolution = c(0.4,0.6,1,1.2,1.6))
Idents(vascu) <- "SCT_snn_res.0.4"
DimPlot(object = vascu, reduction = "umap", label = T)
DimPlot(object = vascu, reduction = "umap", label = T, split.by = "orig.ident")

# find marker
DefaultAssay(vascu) <- "RNA"
vascu <- NormalizeData(vascu)
vascu.markers <- FindAllMarkers(vascu, test.use = "MAST", latent.vars = c("age","sex","PMI"))


Idents(vascu) <- "SCT_snn_res.0.4"
FeaturePlot(vascu, features = c("CEMIP","LAMA2","FLT1"))
FeaturePlot(vascu, features = c("GFAP"))
FeaturePlot(vascu, features = c("RGS5"), split.by = "orig.ident")
VlnPlot(vascu, features = c("COL4A1"), group.by = "orig.ident")
### SLC4A4 (M.FB,5),ABCA9 (P.FB,4),TAGLN (aSMC,3), TRPC4 (pericyte,2), TLL1 (ven, 7)
FeaturePlot(vascu, features = c("TRPC4","TAGLN","SLC4A4", "ABCA9", "TLL1"), label = T)
### ARL15 (ART), 
FeaturePlot(vascu, features = c("ABCB1","ANGPT2","APOD","ARL15"), label = T)

# cell type annotation
vascu$celltype1 <- NA
vascu$celltype1[which(vascu$SCT_snn_res.0.4 %in% c(0,3,6))] <- "Endo"
vascu$celltype1[which(vascu$SCT_snn_res.0.4 == 1)] <- "Peri"
vascu$celltype1[which(vascu$SCT_snn_res.0.4 == 5)] <- "aSMC"
vascu$celltype1[which(vascu$SCT_snn_res.0.4 == 4)] <- "P.FB"
vascu$celltype1[which(vascu$SCT_snn_res.0.4 == 2)] <- "M.FB"
DimPlot(object = vascu, reduction = "umap", label = T, group.by = "celltype1")

save(vascu, file = "~/vascular/NHD_vascu.Robj")


Idents(vascu) <- "celltype1"
vascu.markers <- FindAllMarkers(vascu, test.use = "MAST", latent.vars = c("age","sex","PMI"))
write.table(vascu.markers,"~/vascular/vascu_markers_celltype.txt", sep= "\t")


# save UMAP
library(RColorBrewer)
pal <- brewer.pal(11, "Paired")
pdf("~/vascular/UMAP_vascu.pdf", width = 5, height = 5)
DimPlot(vascu, cols = pal, group.by = "celltype1", pt.size = 0.1, label = T)
dev.off()


# DE analysis
Idents(vascu) <- "celltype1"
vascu$celltype.disease <- paste(Idents(vascu), vascu$disease, sep = "_")
Idents(vascu) <- "celltype.disease"
levels(vascu)
n = length(levels(vascu)) # define cluster number

DefaultAssay(vascu) <- "RNA"
vascu <- NormalizeData(vascu)
de <- data.frame()
for (i in unique(vascu@meta.data$celltype1)){
  nhd <- paste(i, "NHD", sep = "_")
  ctrl <- paste(i, "Control", sep = "_")
  de.genes <- FindMarkers(vascu, ident.1 = nhd, ident.2 = ctrl, test.use = "MAST", latent.vars = c("sex","age","PMI"))
  de.genes$cluster <- i
  de.genes$gene <- rownames(de.genes)
  de <- rbind(de, de.genes)
}
write.table(de,"~/vascular/nhd_de_vascu_MAST_RNA.txt", sep= "\t")


# marker dot plot
DefaultAssay(vascu) <- "RNA"
Idents(vascu) <- "celltype1"
levels(vascu) <- c("Endo","Peri","aSMC","M.FB","P.FB")
levels(vascu) <- c("Endo","Peri","aSMC","M.FB","P.FB") %>% rev()
pdf("~/vascular/dotplot_vascu_rev.pdf", width = 6, height = 5)
DotPlot(vascu, features = unique(c("FLT1", "ABCB1","DLC1","ABCC9",
                                   "TAGLN","ACTA2","SLC4A4","KCNMA1","ABCA9","LAMA2"), assay = "RNA") %>% rev()) + RotatedAxis()
dev.off()


# calculate number of DEGs per cluster
sig.number <- data.frame("up", "down", "cluster")
for (i in unique(vascu@meta.data$celltype1)){
  temp <- data.frame("up", "down", "cluster")
  sig.up <- de %>%
    filter(cluster == i) %>%
    filter(avg_log2FC>0.5 & p_val_adj<0.05)
  up.number <- nrow(sig.up)
  temp$X.up. <- up.number
  sig.down <- de %>%
    filter(cluster == i) %>%
    filter(avg_log2FC<0) %>%
    filter(abs(avg_log2FC)>0.5 & p_val_adj<0.05)
  down.number <- nrow(sig.down)
  temp$X.down. <- down.number
  temp$X.cluster. <- i
  sig.number <- rbind(sig.number, temp)
}
write.table(sig.number,"~vascular/siggene_number.txt", sep= "\t")



### volcano plot
de <- read.table("~/vascular/nhd_de_vascu_MAST_RNA.txt", sep= "\t")

### endo ###
endo.de <- de %>%
  filter(cluster == "Endo")
up <- endo.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("AP|/.",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("LINC",gene))
down <- endo.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("AP|/.",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("LINC",gene))
sig <- endo.de %>% 
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("AP|/.",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes)
anno_list <- sig$gene

endo.de$diffexpressed <- "NO"
endo.de$diffexpressed[endo.de$avg_log2FC > 0.5 & endo.de$p_val_adj < 0.05] <- "UP"
endo.de$diffexpressed[endo.de$avg_log2FC < -0.5 & endo.de$p_val_adj < 0.05] <- "DOWN"

library(ggrepel)
endo.de  %>% 
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj)), color = diffexpressed) +
  geom_point(size = 0.5) +
  geom_point(data = up, color = "#e42625") +
  geom_point(data = down, color = "#030a8f") +
  geom_text_repel(data = dplyr::filter(endo.de, gene %in% anno_list), size = 5,
                  aes(label=gene),box.padding = 0.4, color = "black", max.overlaps = 50) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("NHD vs. Control") +
  ylab("-log10(adj. p-value)") +
  xlab("log2(Fold Change)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "~/vascular/endo_de_RNA_PMI.pdf",device = "pdf", height = 9, width = 11)



### pericytes ###
de <- read.table("~/vascular/nhd_de_vascu_MAST_RNA_PMI.txt", sep= "\t")

cap <- de %>%
  filter(cluster == "Peri")
up <- cap %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene))
down <- cap %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene))
sig <- cap %>% 
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("AP|/.",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("LINC",gene)) 
anno_list <- sig$gene

cap$diffexpressed <- "NO"
cap$diffexpressed[cap$avg_log2FC > 0.5 & cap$p_val_adj < 0.05] <- "UP"
cap$diffexpressed[cap$avg_log2FC < -0.5 & cap$p_val_adj < 0.05] <- "DOWN"

library(ggrepel)
cap  %>% 
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj)), color = diffexpressed) +
  geom_point(size = 0.5) +
  geom_point(data = up, color = "#e42625") +
  geom_point(data = down, color = "#030a8f") +
  geom_text_repel(data = dplyr::filter(cap, gene %in% anno_list), size=5,
                  aes(label=gene),box.padding = 0.4, color = "black", max.overlaps = 50) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("NHD vs. Control") +
  ylab("-log10(adj. p-value)") +
  xlab("log2(Fold Change)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "~/vascular/peri_de_RNA_PMI.pdf",device = "pdf", height = 9, width = 9)



### violin plot
DefaultAssay(vascu) <- "RNA"
vascu <- NormalizeData(vascu)
Idents(vascu) <- "celltype1"
celltype.averages <- AverageExpression(vascu, add.ident = "sample2", assays = "RNA", slot = "data")
write.table(celltype.averages,"~/vascular/vascu_avg_sample_RNA.txt", sep= "\t")

avg <- read.table(file="~/vascular/vascu_avg_sample_RNA.txt", header = T)


### pericytes
avg <- avg %>%
  rownames_to_column("gene")
celltype = "RNA.Peri"
avg_peri <- avg[, c("gene",paste0(celltype,"_Control1"),paste0(celltype,"_Control2"),paste0(celltype,"_Control3"),
                    paste0(celltype,"_Control4"),paste0(celltype,"_Control5"),
                    paste0(celltype,"_Control6"),paste0(celltype,"_Control7"),paste0(celltype,"_Control8"),paste0(celltype,"_Control9"),
                    paste0(celltype,"_Control10"),paste0(celltype,"_Control11"),paste0(celltype,"_NHD1"),
                    paste0(celltype,"_NHD2"),paste0(celltype,"_NHD3"))]

de <- read.table("~/vascular/nhd_de_vascu_MAST.txt", sep= "\t")
celltype = "Peri"
peri <- de %>%
  filter(cluster == celltype)

merge <- merge(peri, avg_peri, by.x = "gene")

peri_vln <- merge %>% 
  pivot_longer(!c(colnames(.)[1:7]), names_to = "sample", values_to = "exprs") %>%
  mutate(disease = ifelse(
    grepl(pattern = "RNA.Peri_NHD1|RNA.Peri_NHD2|RNA.Peri_NHD3", sample),"NHD","Control"
  ))

peri_vln$disease <- as.factor(peri_vln$disease)


library(ggpubr)
genes <- c("PDE7B","COL4A2","SLC6A17","NAMPT","IL18R1","NEBL")
for (i in genes) {
  peri_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease)) + 
    scale_color_manual(values=c("#4168ff", "#fc7153")) +
    labs(title = i) +
    geom_boxplot(fill = "white") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + theme_classic() + theme(legend.position="none")
  ggsave(filename = paste0("~/vascular/peri_vlnplot/",i,".pdf"),device = "pdf", width = 3, height = 5)
}


peri.up <- peri %>%
  filter(avg_log2FC>0.5 & p_val_adj<0.05) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  top_n(wt = avg_log2FC, n = 10)
for (i in peri.up$gene) {
  peri_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease)) + 
    scale_color_manual(values=c("#4168ff", "#fc7153")) +
    labs(title = i) +
    geom_boxplot(fill = "white") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + theme_classic() + theme(legend.position="none")
  ggsave(filename = paste0("~/vascular/peri_vlnplot/",i,".pdf"),device = "pdf", width = 3, height = 5)
}

peri.down <- peri %>%
  filter(avg_log2FC<(-0.5) & p_val_adj<0.05) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  top_n(wt = avg_log2FC, n = -20)
for (i in peri.down$gene) {
  peri_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease)) + 
    scale_color_manual(values=c("#4168ff", "#fc7153")) +
    labs(title = i) +
    geom_boxplot(fill = "white") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + theme_classic() + theme(legend.position="none")
  ggsave(filename = paste0("~/vascular/peri_vlnplot/",i,".pdf"),device = "pdf", width = 3, height = 5)
}


### endo
avg <- read.table(file="~/vascular/vascu_avg_sample_RNA.txt", header = T)

avg <- avg %>%
  rownames_to_column("gene")
celltype = "RNA.Endo"
avg_endo <- avg[, c("gene",paste0(celltype,"_Control1"),paste0(celltype,"_Control2"),paste0(celltype,"_Control3"),
                    paste0(celltype,"_Control4"),paste0(celltype,"_Control5"),
                    paste0(celltype,"_Control6"),paste0(celltype,"_Control7"),paste0(celltype,"_Control8"),paste0(celltype,"_Control9"),
                    paste0(celltype,"_Control10"),paste0(celltype,"_Control11"),paste0(celltype,"_NHD1"),
                    paste0(celltype,"_NHD2"),paste0(celltype,"_NHD3"))]

de <- read.table("~/vascular/nhd_de_vascu_MAST_RNA.txt", sep= "\t")
celltype = "Endo"
endo <- de %>%
  filter(cluster == celltype)

merge <- merge(endo, avg_endo, by.x = "gene")

endo_vln <- merge %>% 
  pivot_longer(!c(colnames(.)[1:7]), names_to = "sample", values_to = "exprs") %>%
  mutate(disease = ifelse(
    grepl(pattern = "Endo_NHD1|Endo_NHD2|Endo_NHD3", sample),"NHD","Control"
  ))

endo_vln$disease <- as.factor(endo_vln$disease)



library(ggpubr)
genes <- c("PDE7B","COL4A2","SLC6A17","NAMPT","IL18R1","IL1RL1")
for (i in genes) {
  peri_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs), color = disease)) + 
    scale_color_manual(values=c("#4168ff", "#fc7153")) +
    labs(title = i) +
    geom_boxplot(fill = "white") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + theme_classic() + theme(legend.position="none")
  ggsave(filename = paste0("~/vascular/peri_vlnplot/",i,".pdf"),device = "pdf")
}

endo.up <- endo %>%
  filter(avg_log2FC>0.5 & p_val_adj<0.05) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  top_n(wt = avg_log2FC, n = 20)
for (i in endo.up$gene) {
  endo_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease)) + 
    scale_color_manual(values=c("#4168ff", "#fc7153")) +
    labs(title = i) +
    geom_boxplot(fill = "white") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-0.5, NA)) +
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + theme_classic() + theme(legend.position="none")
  ggsave(filename = paste0("~/vascular/endo_vlnplot/",i,".pdf"),device = "pdf")
}
