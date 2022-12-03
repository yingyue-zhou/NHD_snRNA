### re-cluster of NHD microglia

rm(list = ls())

# load packages
library(tidyverse)
library(Seurat)
library(sctransform)
library(harmony)
library(gplots)
library(ggplot2)
library(dplyr)


# Re-cluster of microglia
load("~/NHD_harmony.Robj")
Idents(merged_seurat) <- "celltype1"
micro <- subset(merged_seurat, idents = "Micro")
micro <- micro %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, assay.use="SCT")
ElbowPlot(micro, ndims = 50)
micro <- micro %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)

micro <- micro %>% FindClusters(resolution = c(0.3,0.4,0.5,0.6,2,3))
Idents(micro) <- "SCT_snn_res.0.4"
DimPlot(object = micro, reduction = "umap", label = T)

# find markers for doublet removal
DefaultAssay(micro) <- "RNA"
micro <- NormalizeData(micro)
micro.markers <- FindAllMarkers(micro, test.use = "MAST", latent.vars = c("age","sex","PMI"))

FeaturePlot(micro, features = c("MOG","VCAN","AQP4","SYT1"))
FeaturePlot(micro, features = "RBFOX1", split.by = "disease")
FeaturePlot(micro, features = "SYT1")

# remove doublets (16 mito,3,5 neuro)
Idents(micro) <- "SCT_snn_res.2"
micro <- subset(micro, idents = c(3,5,16), invert=T)

micro <- micro %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = c(0.3,0.4,0.5,0.6,0.8))

Idents(micro) <- "SCT_snn_res.0.4"
DimPlot(object = micro, reduction = "umap", label = T)

# save microglia UMAP
pal <- c("#a6cee3","#1f78b4","#b2df8a","#e3211c","#fb9a99","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#fddaec")
pdf("~/micro/micro_cluster_umap.pdf", width = 5, height = 5)
DimPlot(object = micro, reduction = "umap", cols = pal, pt.size = 0.01)
dev.off()

# save microglia sub-cluster feature plots
pdf("~/micro/featureplots.pdf", width = 7, height = 5)
FeaturePlot(micro, features = c("TMEM163", "HAMP","P2RY12","CX3CR1","RPL19","RPS8","F13A1","MRC1","CD83","NR4A1", "CSF1R","VCAN"), pt.size = 0.01)
dev.off()

# find markers
DefaultAssay(micro) <- "RNA"
micro <- NormalizeData(micro)
Idents(micro) <- "SCT_snn_res.0.4"
micro.markers <- FindAllMarkers(micro, test.use = "MAST", latent.vars = c("age","sex","PMI"))

# cell type annotation
Idents(micro) <- "SCT_snn_res.0.4"
micro$celltype1 <- NA
micro$celltype1[which(micro$SCT_snn_res.0.4 %in% c(3))] <- "PVM"
micro$celltype1[which(!micro$SCT_snn_res.0.4 %in% c(3))] <- "Micro"
DimPlot(object = micro, reduction = "umap", label = T, group.by = "celltype1")

# save microglia UMAP
pdf("~/micro/micro_cluster_umap_celltype1.pdf", width = 5, height = 5)
DimPlot(object = micro, reduction = "umap", cols = c("#965454","#ead0d1"), pt.size = 0.01, group.by = "celltype1")
dev.off()

save(micro, file = "~/micro/NHD_micro.Robj")



# DE analysis by cell type
DefaultAssay(micro) <- "RNA"
micro <- NormalizeData(micro)
Idents(micro) <- "celltype1"
micro$celltype.disease <- paste(Idents(micro), micro$disease, sep = "_")
Idents(micro) <- "celltype.disease"
levels(micro)
DefaultAssay(micro) <- "RNA"
nhd_c_micro <- FindMarkers(micro, ident.1 = "Micro_NHD", ident.2 = "Micro_Control", test.use = "MAST", latent.vars = c("sex","age","PMI"))
write.table(nhd_c_micro,"~/micro/micro_nhd_c_de_MAST.txt", sep= "\t")

nhd_c_pvm <- FindMarkers(micro, ident.1 = "PVM_NHD", ident.2 = "PVM_Control", test.use = "MAST", latent.vars = c("sex","age","PMI"))
write.table(nhd_c_pvm,"~/micro/pvm_nhd_c_de_MAST.txt", sep= "\t")



# scatter plot of logFC vs avg_exprs
# calculate average expression
DefaultAssay(micro) <- "RNA"
Idents(micro) <- "celltype1"
avg <- AverageExpression(micro, assays = "RNA", slot = "data")
write.table(avg, file="~/micro/micro_avg_RNA.txt", sep = "\t")

avg <- read.table(file="~/micro_avg_RNA.txt", header = T, sep = "\t")
avg <- as.data.frame(avg)
avg <- rownames_to_column(avg, "gene")
micro_avg <- avg %>%
  dplyr::select("gene","RNA.Micro") %>%
  dplyr::rename(Micro = RNA.Micro)

# load DE genes
de <- read.table("~/micro/micro_nhd_c_de_MAST.txt", sep= "\t")
micro.de <- de %>%
  rownames_to_column("gene")

# merge DE and avg_exprs
merge <- merge(micro.de, micro_avg, by.x = "gene")

merge <- merge %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC >0) %>%
  dplyr::select("gene","avg_log2FC","Micro")
colnames(merge) <- c("gene","avg_log2FC","avg_exp")

# plot up-regulated genes
neuron.genes <- c("SPARCL1","SYT1","GRIN2A","SNAP25","VSNL1","CNTNAP2","NEFL","CCK")
up <- merge %>% 
  dplyr::filter(avg_log2FC > 0.5) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) 
sig <- micro.de %>% 
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes) %>%
  top_n(wt = abs(avg_log2FC), n = 60)
anno_list <- c(sig$gene,"SMAD3")

library(ggrepel)
merge  %>% 
  ggplot(aes(x = avg_log2FC, y = log2(avg_exp+1)), label = gene) +
  geom_point(size = 0.5) +
  geom_point(data = up, color = "#e42625") +
  geom_text_repel(data = dplyr::filter(merge, gene %in% anno_list), size = 5,
                  aes(label=gene),box.padding = 0.4, color = "black", max.overlaps = 50) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("NHD Micro vs Control Micro") +
  ylab("log2(average expression+1)") +
  xlab("log2(Fold Change)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=10))

ggsave(filename = "~/micro/micro_up.pdf",device = "pdf", height = 9, width = 9)





### sample-level box plots 
sample.averages <- AverageExpression(micro, add.ident = "sample2", assays = "RNA", slot = "data")
write.table(sample.averages,"~/micro/micro_avg_sample_RNA.txt", sep= "\t")

avg <- read.table(file="~/micro/micro_avg_sample_RNA.txt", header = T)

#extract micro cluster for every sample
avg <- avg %>%
  rownames_to_column("gene")
columns <- c("gene","micro_Control1","micro_Control2","micro_Control3","micro_Control4","micro_Control5",
             "micro_Control6","micro_Control7","micro_Control8","micro_Control9",
             "micro_Control10","micro_Control11","micro_NHD1","micro_NHD2","micro_NHD3")
columns <- gsub("micro_","RNA.Micro_",columns)
avg_micro <- avg[, columns]

# load DE genes
de <- read.table("~/micro/micro_nhd_c_de_MAST.txt", sep= "\t")
micro.de <- de %>%
  rownames_to_column("gene")

# merge DE and average expression
merge <- merge(micro.de, avg_micro, by.x = "gene")

# prepare table for violin plots
micro_vln <- merge %>% 
  pivot_longer(!c(colnames(.)[1:6]), names_to = "sample", values_to = "exprs") %>%
  mutate(disease = ifelse(
    grepl(pattern = "RNA.Micro_NHD1|RNA.Micro_NHD2|RNA.Micro_NHD3", sample),"NHD","Control"
  ))

micro_vln$disease <- as.factor(micro_vln$disease)

# plot box plots
library(ggpubr)
micro.up <- micro.de %>%
  filter(avg_log2FC>0.5 & p_val_adj<0.05) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  top_n(wt = avg_log2FC, n = 10)
for (i in micro.up$gene) {
  micro_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease)) + 
    # geom_violin(trim = F) +
    scale_color_manual(values=c("#4168ff", "#fc7153")) +
    labs(title = i) +
    geom_boxplot(fill = "white") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 5) + theme_classic() + theme(legend.position="none")
  ggsave(filename = paste0("~/micro/boxplot_MAST/",i,".pdf"),device = "pdf", width = 3, height = 5)
}


micro.down <- micro.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < -0.5) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes) %>%
  top_n(wt = avg_log2FC, n = -10)
for (i in micro.down$gene) {
  micro_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=exprs)) + 
    geom_violin(trim = F) +
    labs(title = i) +
    geom_boxplot(width=0.1) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4)
  ggsave(filename = paste0("~/micro/vlnplot/",i,".pdf"),device = "pdf")
}



#### PVM ####
avg <- read.table(file="~/micro/micro_avg_sample_sct.txt", header = T)

#extract cluster micro for every sample
avg <- avg %>%
  rownames_to_column("gene")
celltype = "PVM"
avg_pvm <- avg[, c("gene",paste0(celltype,"_control1"),paste0(celltype,"_control2"),paste0(celltype,"_control3"),paste0(celltype,"_TWCC.JC1.JC1.lib1"),
                   paste0(celltype,"_TWCC.JC2.JC2.lib1"),paste0(celltype,"_TWCC.JC3.JC3.lib1"),
                   paste0(celltype,"_TWCC.JC4.JC4.lib1"),paste0(celltype,"_TWCC.JC5.JC5.lib1"),paste0(celltype,"_TWCC.JC6.JC6.lib1"),paste0(celltype,"_TWCC.JC7.JC7.lib1"),
                   paste0(celltype,"_TWCC.JC8.JC8.lib1"),paste0(celltype,"_TWCC.JC9.JC9.lib1"),paste0(celltype,"_TWCC.JC10.JC10.lib1"),paste0(celltype,"_NHD1"),
                   paste0(celltype,"_NHD2"),paste0(celltype,"_NHD3"))]

de <- read.table("~/micro/pvm_nhd_c_de.txt", sep= "\t")
pvm.de <- de %>%
  rownames_to_column("gene")

merge <- merge(pvm.de, avg_pvm, by.x = "gene")

pvm_vln <- merge %>% 
  pivot_longer(!c(colnames(.)[1:6]), names_to = "sample", values_to = "exprs") %>%
  mutate(disease = ifelse(
    grepl(pattern = "PVM_NHD1|PVM_NHD2|PVM_NHD3", sample),"NHD","Control"
  ))

pvm_vln$disease <- as.factor(pvm_vln$disease)

library(ggpubr)
genes <- c("TGFBR2","IRAK3","LRRK2")
for (i in genes) {
  pvm_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=exprs)) + 
    geom_violin(trim = F) +
    labs(title = i) +
    geom_boxplot(width=0.1) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4)
  ggsave(filename = paste0("~/micro/pvm/",i,".pdf"),device = "pdf")
}
