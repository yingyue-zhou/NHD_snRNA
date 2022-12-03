### analysis of NHD, pre-process, metadata

rm(list = ls())

# load packages
library(tidyverse)
library(Seurat)
library(sctransform)
library(harmony)
library(gplots)
library(ggplot2)
library(dplyr)

# load data
# Batch1
dirs1 <- list.dirs(path = "~/batch1",recursive = F) %>%
  list.dirs(recursive = F) %>%
  grep(pattern = "outs", value = TRUE) %>%
  list.dirs(recursive = F)

# Batch2
dirs2 <- list.dirs(path = "~/batch2",recursive = F) %>% 
  list.dirs(recursive = F) %>% 
  grep(pattern = "outs", value = TRUE) %>% 
  list.dirs(recursive = F)

# Batch3
dirs3 <- list.dirs(path = "~/batch3",recursive = F) %>% 
  list.dirs(recursive = F) %>% 
  grep(pattern = "outs", value = TRUE) %>% 
  list.dirs(recursive = F)

# Batch4
dirs4 <- list.dirs(path = "~/batch4",recursive = F) %>% 
  list.dirs(recursive = F) %>% 
  grep(pattern = "outs", value = TRUE) %>% 
  list.dirs(recursive = F)
dirs <- c(dirs1, dirs2, dirs3, dirs4)
dirlist <- grep("outs/filtered", dirs, value = TRUE)
rlist <- list()
slist <- list()
i = 1

# create Seurat object
for(i in 1:length(dirlist)){
  samplename <- dirlist[i] %>%
    dirname() %>%
    dirname() %>%
    basename()
  slist[[i]] <- samplename
  
  
  rlist[[i]] <- Read10X(data.dir = dirlist[i])%>%
    CreateSeuratObject(project = samplename,) %>%
    AddMetaData(metadata = samplename, col.name = "sample")
}

merged_seurat <- merge(x = rlist[[1]], 
                       y = c(rlist[2:32]), add.cell.ids = c(slist[1],slist[2:32]))

# calculate mitochondrial ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# quality control
merged_seurat@meta.data %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

merged_seurat@meta.data %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 1000)

merged_seurat@meta.data %>%   
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 700)

merged_seurat@meta.data %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 0.05)

# filter by UMI counts, gene counts, mitochondrial ratio
merged_seurat <- merged_seurat %>%
  subset(
    subset = (nCount_RNA > 1000 & nCount_RNA < 100000) &
      (nFeature_RNA > 700 & nFeature_RNA < 12000) & 
      (mitoRatio < 0.05)
  )


# add metadata
merged_seurat$disease <- NA
merged_seurat$disease[which(str_detect(merged_seurat$orig.ident, "TWCC-NH1-NH1|TWCC-NH2-NH2|TWCC-N1_97-OC-NHD-1-lib1|MCAM-MMH861_occipital-2-lib1"))] <- "NHD"
merged_seurat$disease[which(!str_detect(merged_seurat$orig.ident, "TWCC-NH1-NH1|TWCC-NH2-NH2|TWCC-N1_97-OC-NHD-1-lib1|MCAM-MMH861_occipital-2-lib1"))] <- "Control"

merged_seurat$sex <- NA
merged_seurat$sex[which(str_detect(merged_seurat$orig.ident, "MCAM-MMH861_occipital-2-lib1|TWCC-C1-C1-lib1|TWCC-JC1-JC1-lib1|TWCC-JC10-JC10-lib1|TWCC-JC3-JC3-lib1|TWCC-JC7-JC7-lib1|TWCC-N1_97-OC-NHD-1-lib1|TWCC-NH1-NH1|TWCC-R8_04-OC-Control-1-lib1"))] <- "female"
merged_seurat$sex[which(!str_detect(merged_seurat$orig.ident, "MCAM-MMH861_occipital-2-lib1|TWCC-C1-C1-lib1|TWCC-JC1-JC1-lib1|TWCC-JC10-JC10-lib1|TWCC-JC3-JC3-lib1|TWCC-JC7-JC7-lib1|TWCC-N1_97-OC-NHD-1-lib1|TWCC-NH1-NH1|TWCC-R8_04-OC-Control-1-lib1"))] <- "male"

merged_seurat$sample2 <- merged_seurat$sample
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-C2-C2|TWCC-R2_07-OC-Control-1-lib1"))] <- "Control1"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-R8_04-OC-Control-1-lib1"))] <- "Control2"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC1-JC1-lib1"))] <- "Control3"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC2-JC2-lib1"))] <- "Control4"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC3-JC3-lib1"))] <- "Control5"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC4-JC4-lib1"))] <- "Control6"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC5-JC5-lib1"))] <- "Control7"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC6-JC6-lib1"))] <- "Control8"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC7-JC7-lib1"))] <- "Control9"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC8-JC8-lib1"))] <- "Control10"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-JC9-JC9-lib1"))] <- "Control11"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-NH1-NH1|TWCC-N1_97-OC-NHD-1-lib1"))] <- "NHD1"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "TWCC-NH2-NH2"))] <- "NHD2"
merged_seurat$sample2[which(str_detect(merged_seurat$orig.ident, "MCAM-MMH861_occipital-2-lib1"))] <- "NHD3"

merged_seurat$age <- NA
merged_seurat$age[which(str_detect(merged_seurat$sample2, "Control1|NHD1|NHD2|NHD3"))] <- "age49-60"
merged_seurat$age[which(str_detect(merged_seurat$sample2, "Control2"))] <- "age61-70"
merged_seurat$age[which(str_detect(merged_seurat$sample2, "Control3|Control4|Control6|Control7|Control9|Control10|Control11"))] <- "age71-80"
merged_seurat$age[which(str_detect(merged_seurat$sample2, "Control5|Control8"))] <- "age81-90"

merged_seurat$PMI <- NA
merged_seurat$PMI[which(str_detect(merged_seurat$sample2, "Control11|NHD3"))] <- "PMI>10"
merged_seurat$PMI[which(!str_detect(merged_seurat$sample2, "Control11|NHD3"))] <- "PMI<10"

save(merged_seurat, file = "~/NHD_preprocess.Robj")
