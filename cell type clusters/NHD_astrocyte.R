### re-cluster of NHD astrocytes

rm(list = ls())

# load packages
library(tidyverse)
library(Seurat)
library(sctransform)
library(harmony)
library(gplots)
library(ggplot2)
library(dplyr)


# Re-cluster of astrocytes
load("~/NHD_harmony.Robj")
Idents(merged_seurat) <- "celltype1"
astro <- subset(merged_seurat, idents = "Astro")
astro <- astro %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, assay.use="SCT")
ElbowPlot(astro, ndims = 50)
astro <- astro %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)
astro <- astro %>% FindClusters(resolution = 3)
DimPlot(object = astro, reduction = "umap", label = T)
DimPlot(object = astro, reduction = "umap", group.by = "orig.ident")
DimPlot(object = astro, reduction = "umap", split.by = "disease")
DimPlot(object = astro, reduction = "umap", split.by = "region")

# find markers for doublet removal
DefaultAssay(astro) <- "RNA"
astro <- NormalizeData(astro)
astro.markers <- FindAllMarkers(astro, test.use = "MAST", latent.vars = c("age","sex","PMI"))
FeaturePlot(astro, features = c("CSF1R","MOG","SYT1"))

# remove contaminated cells (31,33)
Idents(astro) <- "SCT_snn_res.3"
astro <- subset(astro, idents = c(31,33), invert=T)
astro <- astro %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.8)
DimPlot(object = astro, reduction = "umap", label = T)

# remove sample-specific cluster (10)
Idents(astro) <- "SCT_snn_res.0.8"
astro <- subset(astro, idents = c(10), invert=T)
astro <- astro %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.8)
DimPlot(object = astro, reduction = "umap", label = T)

save(astro, file = "~/astro/NHD_astro.Robj")


# DE analysis
DefaultAssay(astro) <- "RNA"
astro <- NormalizeData(astro)
Idents(astro) <- "disease"
astro_de <- FindMarkers(astro, ident.1 = "NHD", ident.2 = "Control", test.use = "MAST", latent.vars = c("age","sex","PMI"), assay = "RNA")
write.table(astro_de,"~/astro/astro_nhd_c_de_RNA.txt", sep= "\t")


# volcano plot
de <- read.table("~/astro/astro_nhd_c_de_RNA.txt", sep= "\t")

neuron.genes <- c("SPARCL1","SYT1","GRIN2A","SNAP25","VSNL1","CNTNAP2","NEFL","CCK")

astro.de <- de %>%
  rownames_to_column("gene")
astro.de <- astro.de %>%
  filter(p_val_adj < 0.05)
up <- astro.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) 
down <- astro.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) 
sig <- astro.de %>% 
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes) %>%
  top_n(wt = abs(avg_log2FC), n = 50)
anno_list <- sig$gene
wound_response <- c("DST","ZFP36L2","CD44","F3","FGF2","FGFR2","FN1","GAP43","GFAP","TNC","CCN1","ITGB1","LYN","SMAD3","MAP3K5","MYH9","NTRK3","PLSCR1",
                    "PROS1","SOD2","TGFB2","ZFP36","YAP1","MACF1","UBASH3B","ANO6","PAM","TRAF3IP2")
not_wound_response <- anno_list[which(!anno_list %in% wound_response)]

astro.de$diffexpressed <- "NO"
astro.de$diffexpressed[astro.de$avg_log2FC > 0.5 & astro.de$p_val_adj < 0.05] <- "UP"
astro.de$diffexpressed[astro.de$avg_log2FC < -0.5 & astro.de$p_val_adj < 0.05] <- "DOWN"

library(ggrepel)
astro.de  %>% 
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(size = 1) +
  geom_text_repel(data = dplyr::filter(astro.de, gene %in% not_wound_response), size = 5,
                  aes(label=gene),box.padding = 0.4, color = "black", max.overlaps = 50) +
  geom_text_repel(data = dplyr::filter(astro.de, gene %in% wound_response), size = 5,
                  aes(label=gene),box.padding = 0.4, color = "#bd1a26", max.overlaps = 50) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("NHD vs. Control") +
  ylab("-log10(adj. p-value)") +
  xlab("log2(Fold Change)") +
  theme_classic()
ggsave(filename = "~/astro/astro_de1_MAST_RNA.pdf",device = "pdf", height = 11, width = 11)




# sample-level box plot
DefaultAssay(astro) <- "RNA"
astro <- NormalizeData(astro)
Idents(astro) <- "sample2"
celltype.averages <- AverageExpression(astro, assays = "RNA", slot = "data")
write.table(celltype.averages,"~/astro/astro_avg_sample_RNA.txt", sep= "\t")

avg <- read.table(file="~/astro/astro_avg_sample_RNA.txt", header = T)
avg <- avg %>%
  rownames_to_column("gene")
micro.name <- c("gene","micro_Control1","micro_Control2","micro_Control3","micro_Control4","micro_Control5",
                "micro_Control6","micro_Control7","micro_Control8","micro_Control9",
                "micro_Control10","micro_Control11","micro_NHD1","micro_NHD2","micro_NHD3")
astro.name <- str_replace_all(micro.name, "micro_", "RNA.")
avg_astro <- avg[, astro.name]

de <- read.table("~/astro/astro_nhd_c_de_RNA.txt", sep= "\t")
astro.de <- de %>%
  rownames_to_column("gene")

merge <- merge(astro.de, avg_astro, by.x = "gene")

astro_vln <- merge %>% 
  pivot_longer(!c(colnames(.)[1:6]), names_to = "sample", values_to = "exprs") %>%
  mutate(disease = ifelse(
    grepl(pattern = "NHD1|NHD2|NHD3", sample),"NHD","Control"
  ))

astro_vln$disease <- as.factor(astro_vln$disease)

library(ggpubr)
genes <- c("ABCC9","GADD45B","BCL6","OSMR","IGFBP7","GFAP","SOD2","NRP1","CCL2","IL1R1","FGF2")
for (i in genes) {
  astro_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease, fill = disease)) + 
    # geom_violin(trim = F) +
    scale_color_manual(values=c("#0a3cfc", "#ff3019")) +
    scale_fill_manual(values=c("#99c6fd", "#fdada2")) +
    labs(title = i) +
    geom_boxplot() + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + 
    theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  + theme(legend.position="none")
  ggsave(filename = paste0("~/astro/vlnplot/",i,".pdf"),device = "pdf", width = 3, height = 5)
}


astro.up <- astro.de %>%
  filter(avg_log2FC>0.5 & p_val_adj<0.05) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  top_n(wt = avg_log2FC, n = 20)
for (i in astro.up$gene) {
  astro_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease, fill = disease)) + 
    scale_color_manual(values=c("#0a3cfc", "#ff3019")) +
    scale_fill_manual(values=c("#99c6fd", "#fdada2")) +
    labs(title = i) +
    geom_boxplot() + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + 
    theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  + theme(legend.position="none")
  ggsave(filename = paste0("~/astro/vlnplot/",i,".pdf"),device = "pdf", width = 3, height = 5)
}


