### re-cluster of NHD oligodendrocytes

rm(list = ls())

# load packages
library(tidyverse)
library(Seurat)
library(sctransform)
library(harmony)
library(gplots)
library(ggplot2)
library(dplyr)


# Re-cluster of oligodendrocytes
load("~/NHD_harmony.Robj")
Idents(merged_seurat) <- "celltype1"
oligo <- subset(merged_seurat, idents = c("Oligo","OPC"))
oligo <- oligo %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, assay.use="SCT")
ElbowPlot(oligo, ndims = 50)
oligo <- oligo %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)
oligo <- oligo %>% FindClusters(resolution = c(3,5))
DimPlot(object = oligo, reduction = "umap", label = T)

Idents(oligo) <- "SCT_snn_res.5"

# find markers
DefaultAssay(oligo) <- "RNA"
oligo <- NormalizeData(oligo)
oligo.markers <- FindAllMarkers(oligo, test.use = "MAST", latent.vars = c("age","sex","PMI"))
FeaturePlot(oligo, features = c("CSF1R","VCAN","SYT1"))
FeaturePlot(oligo, features = c("nFeature_RNA"))

# remove contaminated cells neurons (51)
Idents(oligo) <- "SCT_snn_res.5"
oligo <- subset(oligo, idents = c(51), invert=T)
oligo <- oligo %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.8)

oligo <- oligo %>% FindClusters(resolution = c(1,1.6,2,3,4))
oligo <- oligo %>% FindClusters(resolution = c(6))
Idents(oligo) <- "SCT_snn_res.1"
DimPlot(object = oligo, reduction = "umap", label = T)

Idents(oligo) <- "SCT_snn_res.6"
FeaturePlot(oligo, features = "nFeature_RNA")
FeaturePlot(oligo, features = "FYN", split.by = "disease")
DimPlot(oligo, label=T, cells.highlight=WhichCells(oligo,idents=60), cols.highlight = "darkred", cols= "grey")


# remove contaminated cells astro (60)
Idents(oligo) <- "SCT_snn_res.6"
oligo <- subset(oligo, idents = c(60), invert=T)
oligo <- oligo %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.8)
oligo <- oligo %>% FindClusters(resolution = c(1,1.6,2,3,4))
oligo <- oligo %>% FindClusters(resolution = c(1.2))
Idents(oligo) <- "SCT_snn_res.1.2"
DimPlot(object = oligo, reduction = "umap", label = T)
DimPlot(object = oligo, reduction = "umap", group.by = "orig.ident")


# cell type annotation
Idents(oligo) <- "SCT_snn_res.1.2"
oligo$celltype1 <- NA
oligo$celltype1[which(!oligo$SCT_snn_res.1.2 %in% c(6,13,16))] <- "oligo"
oligo$celltype1[which(oligo$SCT_snn_res.1.2 %in% c(6,13))] <- "OPC"
oligo$celltype1[which(oligo$SCT_snn_res.1.2 == 16)] <- "transition"
DimPlot(object = oligo, reduction = "umap", label = T, group.by = "celltype1")

save(oligo, file = "~/oligo/NHD_oligo.Robj")

# find markers
DefaultAssay(oligo) <- "RNA"
oligo <- NormalizeData(oligo)
Idents(oligo) <- "celltype1"
oligo.markers <- FindAllMarkers(oligo, test.use = "MAST", latent.vars = c("age","sex","PMI"))
write.table(oligo.markers,"~/oligo/oligo.markers.txt", sep= "\t")


# save UMAP
Idents(oligo) <- "celltype1"
levels(oligo)
library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(12,"Paired"))(12)
pal <- c("#a6cee3","#fb9a99","#b2df8a","#e3211c","#b15928","#fdbf6f","#cab2d6","#6a3d9a","#fddaec")

pdf("~/oligo/oligo_umap.pdf", width = 7, height = 5)
DimPlot(object = oligo, reduction = "umap", label = T, cols = pal, pt.size = 0.01)
dev.off()


# marker dotplot
Idents(oligo) <- "celltype1"
levels(oligo) <- c("OPC","transition","oligo") %>% rev()
DefaultAssay(oligo) <- "RNA"
oligo <- NormalizeData(oligo)
pdf("~/oligo/dotplot_oligo_rev.pdf", width = 6, height = 5)
DotPlot(oligo, features = unique(c("PTPRZ1","PDGFRA","GPR17","FYN","MOBP","PLP1"), assay = "RNA") %>% rev()) + RotatedAxis()
dev.off()

# stats
Idents(oligo) <- "celltype1"
nuclei.number <- table(oligo$celltype1, oligo$sample2)
write.table(nuclei.number,"~/oligo/nuclei_number_sample.txt", sep= "\t")


# DE analysis by cell type
Idents(oligo) <- "celltype1"
oligo$celltype.disease <- paste(Idents(oligo), oligo$disease, sep = "_")
Idents(oligo) <- "celltype.disease"
DefaultAssay(oligo) <- "RNA"
oligo <- NormalizeData(oligo)
de <- data.frame()
for (i in unique(oligo@meta.data$celltype1)){
  nhd <- paste(i, "NHD", sep = "_")
  ctrl <- paste(i, "Control", sep = "_")
  de.genes <- FindMarkers(oligo, ident.1 = nhd, ident.2 = ctrl, test.use = "MAST", latent.vars = c("sex","age","PMI"))
  de.genes$cluster <- i
  de.genes$gene <- rownames(de.genes)
  de <- rbind(de, de.genes)
}
write.table(de,"~/oligo/nhd_de_oligo_RNA.txt", sep= "\t")


# volcano plot
de <- read.table("~/oligo/nhd_de_oligo_RNA.txt", sep= "\t")

oligo.de <- de %>%
  filter(cluster=="oligo")
neuron.genes <- c("SPARCL1","SYT1","GRIN2A","SNAP25","VSNL1","CNTNAP2","NEFL","CCK")
up <- oligo.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0)%>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes)
down <- oligo.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes)
sig <- oligo.de %>% 
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes) %>%
  top_n(wt = abs(avg_log2FC), n=50)
anno_list <- sig$gene
gene.label <- c("MBP")

oligo.de$diffexpressed <- "NO"
oligo.de$diffexpressed[oligo.de$avg_log2FC > 0.5 & oligo.de$p_val_adj < 0.05] <- "UP"
oligo.de$diffexpressed[oligo.de$avg_log2FC < (-0.5) & oligo.de$p_val_adj < 0.05] <- "DOWN"

library(ggrepel)
oligo.de  %>% 
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj)), color = diffexpressed) +
  geom_point(size = 0.5) +
  geom_point(data = up, color = "#e42625") +
  geom_point(data = down, color = "#030a8f") +
  geom_text_repel(data = filter(oligo.de, gene %in% anno_list), size = 5,
                  aes(label=gene),box.padding = 0.4, color = "black", max.overlaps = 50) +
  geom_text_repel(data = filter(oligo.de, gene %in% gene.label), size = 5,
                  aes(label=gene),box.padding = 0.4, color = "black", max.overlaps = 50) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("NHD vs. Control") +
  ylab("-log10(adj. p-value)") +
  xlab("log2(Fold Change)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "~/oligo/oligo_de_RNA.pdf",device = "pdf", height = 11, width = 16)




# heatmap
Idents(oligo) <- "celltype1"
celltype.averages <- AverageExpression(oligo, add.ident = "sample2", assays = "RNA", slot = "data")
write.table(celltype.averages,"~/oligo/oligo_avg_sample_RNA.txt", sep= "\t")

avg <- read.table(file="~/oligo/oligo_avg_sample_RNA.txt", header = T)

#extract oligo cluster for every sample
avg <- avg %>%
  rownames_to_column("gene")

celltype = "RNA.oligo"
avg_oligo <- avg[, c("gene",paste0(celltype,"_Control1"),paste0(celltype,"_Control2"),paste0(celltype,"_Control3"),
                     paste0(celltype,"_Control4"),paste0(celltype,"_Control5"),
                     paste0(celltype,"_Control6"),paste0(celltype,"_Control7"),paste0(celltype,"_Control8"),paste0(celltype,"_Control9"),
                     paste0(celltype,"_Control10"),paste0(celltype,"_Control11"),paste0(celltype,"_NHD1"),
                     paste0(celltype,"_NHD2"),paste0(celltype,"_NHD3"))]

de <- read.table("~/oligo/nhd_de_oligo_RNA.txt", sep= "\t")
celltype = "oligo"
oligo.de <- de %>%
  filter(cluster == celltype)

merge <- merge(oligo.de, avg_oligo, by.x = "gene")
neuron.genes <- c("SPARCL1","SYT1","GRIN2A","SNAP25","VSNL1","CNTNAP2","NEFL","CCK")

oligo.up <- oligo.de %>%
  filter(avg_log2FC>0.5 & p_val_adj<0.05) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  top_n(wt = avg_log2FC, n = 40)
oligo.up.gene <- merge %>%
  filter(gene %in% oligo.up$gene) %>%
  arrange(desc(avg_log2FC)) %>%
  select(!c(colnames(.)[2:7])) %>%
  column_to_rownames(var="gene") %>%
  data.matrix()
colpal <- colorRampPalette(c("#070066","white", "#fc8786"))(256)
par(mar=c(0.2,0.2,0.2,0.2))
pdf("~/oligo/oligo_top40_RNA_red.pdf")
heatmap.2(oligo.up.gene, scale = "row", col=colpal, dendrogram = "none", Rowv=FALSE, Colv=FALSE,trace='none',sepwidth=c(0.05,0.05),sepcolor="white",
          colsep=1:ncol(oligo.up.gene), rowsep=1:nrow(oligo.up.gene), density.info = "none", margins = c(5,5))
dev.off()

oligo.down <- oligo.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < (-0.5)) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes) %>%
  top_n(wt = avg_log2FC, n = -40)
oligo.down.gene <- merge %>%
  filter(gene %in% oligo.down$gene) %>%
  arrange(desc(avg_log2FC)) %>%
  select(!c(colnames(.)[2:7])) %>%
  column_to_rownames(var="gene") %>%
  data.matrix()
colpal <- colorRampPalette(c("#070066","white", "#fc8786"))(256)
pdf("~/oligo/oligo_down40_RNA_red.pdf")
heatmap.2(oligo.down.gene, scale = "row", col=colpal, dendrogram = "none", Rowv=FALSE, Colv=FALSE,trace='none',sepwidth=c(0.05,0.05),sepcolor="white",
          colsep=1:ncol(oligo.down.gene), rowsep=1:nrow(oligo.down.gene), density.info = "none", margins = c(5,10))
dev.off()





# box plot
avg <- read.table(file="~/oligo/oligo_avg_sample_RNA.txt", header = T)

#extract oligo cluster for every sample
avg <- avg %>%
  rownames_to_column("gene")

celltype = "RNA.oligo"
avg_oligo <- avg[, c("gene",paste0(celltype,"_Control1"),paste0(celltype,"_Control2"),paste0(celltype,"_Control3"),
                     paste0(celltype,"_Control4"),paste0(celltype,"_Control5"),
                     paste0(celltype,"_Control6"),paste0(celltype,"_Control7"),paste0(celltype,"_Control8"),paste0(celltype,"_Control9"),
                     paste0(celltype,"_Control10"),paste0(celltype,"_Control11"),paste0(celltype,"_NHD1"),
                     paste0(celltype,"_NHD2"),paste0(celltype,"_NHD3"))]

de <- read.table("~/oligo/nhd_de_oligo_RNA.txt", sep= "\t")
celltype = "oligo"
oligo.de <- de %>%
  filter(cluster == celltype)

merge <- merge(oligo.de, avg_oligo, by.x = "gene")
oligo_vln <- merge %>% 
  pivot_longer(!c(colnames(.)[1:7]), names_to = "sample", values_to = "exprs") %>%
  mutate(disease = ifelse(
    grepl(pattern = "oligo_NHD1|oligo_NHD2|oligo_NHD3", sample),"NHD","Control"
  ))

oligo_vln$disease <- as.factor(oligo_vln$disease)

library(ggpubr)
genes <- c("IGF1R","LRP2","MBP","APOD","FAIM2","SELENOP","ADGRB3")
for (i in genes) {
  oligo_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease, fill = disease)) + 
    # geom_violin(trim = F) +
    scale_color_manual(values=c("#0a3cfc", "#ff3019")) +
    scale_fill_manual(values=c("#99c6fd", "#fdada2")) +
    labs(title = i) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") + 
    geom_boxplot() + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + 
    theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  + theme(legend.position="none")
  ggsave(filename = paste0("~/oligo/boxplot/",i,".pdf"),device = "pdf", width = 3, height = 5)
}

oligo.up <- oligo.de %>%
  filter(avg_log2FC>0.5 & p_val_adj<0.05) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("FP|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  top_n(wt = avg_log2FC, n = 20)
for (i in oligo.up$gene) {
  oligo_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease)) + 
    scale_color_manual(values=c("#4168ff", "#fc7153")) +
    labs(title = i) +
    geom_boxplot(fill="white") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + 
    theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  + theme(legend.position="none")
  ggsave(filename = paste0("~/oligo/boxplot/",i,".pdf"),device = "pdf",width = 3,height = 5)
}

oligo.down <- oligo.de %>% 
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < -0.5) %>%
  filter(!grepl("XIST",gene)) %>%
  filter(!grepl("MT-|/.",gene)) %>%
  filter(!grepl("AC|/.",gene)) %>%
  filter(!grepl("AL|/.",gene)) %>%
  filter(!grepl("LINC",gene)) %>%
  filter(!gene %in% neuron.genes) %>%
  top_n(wt = avg_log2FC, n = -10)
for (i in oligo.down$gene) {
  oligo_vln %>%
    filter(gene==i) %>%
    ggplot(aes(x=disease, y=log2(exprs+1), color = disease)) + 
    # geom_violin(trim = F) +
    scale_color_manual(values=c("#4168ff", "#fc7153")) +
    labs(title = i) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") + 
    geom_boxplot(fill="white") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 4) + 
    theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  + theme(legend.position="none")
  ggsave(filename = paste0("~/oligo/boxplot/",i,".pdf"),device = "pdf",width = 3,height = 5)
}



### OPC
load("~/oligo/NHD_oligo.Robj")

DefaultAssay(oligo) <- "RNA"
VlnPlot(oligo, features = "BCAN", split.by = "disease", group.by = "celltype1", cols = c("#b3cde3", "#fbb4ae"), pt.size = 0) +
  geom_boxplot(aes(color = oligo$disease), width =0.1, position = position_dodge(0.9)) + 
  scale_color_manual(values=c("#323695", "#e3211c"))
VlnPlot(oligo, features = "OLIG2", split.by = "disease", group.by = "celltype1")
VlnPlot(oligo, features = "APOD", split.by = "disease", group.by = "celltype1")
VlnPlot(oligo, features = "NLGN4Y", split.by = "disease", group.by = "celltype1")
VlnPlot(oligo, features = "PTGDS", split.by = "disease", group.by = "celltype1")


genes <- c("BCAN","OLIG2","APOD","PTGDS")
for (i in genes) {
  pdf(paste0("~/OPC/vlnplot/",i,".pdf"))
  p=VlnPlot(oligo, features = i, split.by = "disease", group.by = "celltype1", cols = c("#b3cde3", "#fbb4ae"), pt.size = 0) +
    geom_boxplot(aes(color = oligo$disease), width =0.1, position = position_dodge(0.9)) + 
    scale_color_manual(values=c("#323695", "#e3211c"))
  print(p)
  dev.off()
}

for (i in genes) {
  data <- FetchData(oligo, vars = c(i,"disease","celltype1"))
  data %>%
    ggplot(aes_string(x="celltype1", y=i)) +
    geom_violin(
      aes(fill = disease), position = position_dodge(0.6)) +
    scale_fill_manual(values=c("#b3cde3", "#fbb4ae")) +
    stat_compare_means(aes(label = "p.signif"), method = "wilcox", label.x = 1.5, label.y = 5) + 
    theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  ggsave(filename = paste0("~/OPC/vlnplot/",i,"_stat.pdf"),device = "pdf")
  
}

data <- FetchData(oligo, vars = c("BCAN","disease","celltype1"))
head(data)

data$celltype1 <- as.factor(data$celltype1)
library(ggpubr)
data %>%
  ggplot(aes_string(x="celltype1", y=i)) +
  geom_violin(
    aes(fill = disease), position = position_dodge(0.6)) +
  scale_fill_manual(values=c("#b3cde3", "#fbb4ae")) +
  geom_boxplot(
    aes(color = disease), width =0.05, position = position_dodge(0.6)) + 
  scale_color_manual(values=c("#323695", "#e3211c")) +
  stat_compare_means(aes(label = "p.signif"), method = "wilcox", label.x = 1.5, label.y = 5) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


# box plot
avg <- read.table(file="~/oligo/oligo_avg_sample_RNA.txt", header = T)

# extract OPC cluster for every sample
avg <- avg %>%
  rownames_to_column("gene")

celltype = "RNA.OPC"
avg_oligo <- avg[, c("gene",paste0(celltype,"_control2"),paste0(celltype,"_control3"),paste0(celltype,"_TWCC.JC1.JC1.lib1"),
                     paste0(celltype,"_TWCC.JC2.JC2.lib1"),paste0(celltype,"_TWCC.JC3.JC3.lib1"),
                     paste0(celltype,"_TWCC.JC4.JC4.lib1"),paste0(celltype,"_TWCC.JC5.JC5.lib1"),paste0(celltype,"_TWCC.JC6.JC6.lib1"),paste0(celltype,"_TWCC.JC7.JC7.lib1"),
                     paste0(celltype,"_TWCC.JC8.JC8.lib1"),paste0(celltype,"_TWCC.JC9.JC9.lib1"),paste0(celltype,"_NHD1"),
                     paste0(celltype,"_NHD2"),paste0(celltype,"_NHD3"))]

de <- read.table("~/oligo/nhd_de_oligo_RNA.txt", sep= "\t")
celltype = "OPC"
oligo.de <- de %>%
  filter(cluster == celltype)

merge <- merge(oligo.de, avg_oligo, by.x = "gene")
oligo_vln <- merge %>% 
  pivot_longer(!c(colnames(.)[1:7]), names_to = "sample", values_to = "exprs") %>%
  mutate(disease = ifelse(
    grepl(pattern = "NHD1|NHD2|NHD3", sample),"NHD","Control"
  ))

oligo_vln$disease <- as.factor(oligo_vln$disease)

library(ggpubr)
gene.name <- "BCAN"
oligo_vln %>%
  filter(gene==gene.name) %>%
  ggplot(aes(x=disease, y=log2(exprs+1), color = disease)) + 
  scale_color_manual(values=c("#4168ff", "#fc7153")) +
  labs(title = gene.name) +
  geom_boxplot(fill="white") + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  stat_compare_means(aes(label = "p.signif"), label.x = 1.5, label.y = 5) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  + theme(legend.position="none")
ggsave(filename = paste0("~/OPC/",gene.name,".pdf"),device = "pdf",width = 3,height = 5)


