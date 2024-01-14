
library(Seurat)
library("RColorBrewer")

library(dplyr)
library(patchwork)
library(ggplot2)
library(plyr)

setwd("/work/LAS/geetu-lab/sarwar")


in_dir<- "/work/LAS/geetu-lab/sarwar/data/10x_V2_processed_data"

# ================== All cells=====================

data1 <- read.csv(paste(in_dir,"/MPL39603-matrix.csv", sep=""),
                 header=T,
                row.names = 1,
                check.names = F)

dim(data1) # 21474 X 39603

#View(data1[1:2,1:10])



data1.meta = read.csv(paste(in_dir,"/MPL39603-meta.csv", sep=""),
         check.names = FALSE, header=TRUE,row.names=1)


#View(data1.meta)

levels(as.factor(data1.meta$cluster_14))

# sort


data1_sort <- data1[,match(rownames(data1.meta),colnames(data1))]

#View(data1_sort[1:2,1:5])


# ============== Create Seurat obj

MPL.obj <- CreateSeuratObject(counts = data1_sort,
                              meta.data = data1.meta, 
                              project = "MTR15682",
                              min.cells = 3, min.features = 500)


MPL.obj
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
MPL.obj[["percent.mt"]] <- PercentageFeatureSet(MPL.obj, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(MPL.obj, features = c("nFeature_RNA", 
                              "nCount_RNA", 
                              "percent.mt"), ncol = 3)

# Get number of cells per cluster and per sample of origin
table(MPL.obj@meta.data$orig.ident)

# Get number of cells per cluster and per sample of origin
cells_ids <- data.frame(table(MPL.obj@meta.data$cluster_14, 
                              MPL.obj@meta.data$orig.ident))
#View(cells_ids)
#sum(cells_ids$Freq)

colnames(cells_ids) <- c("Celltypes","Samples","Frequency")


cells_ids$Celltypes  <- revalue(cells_ids$Celltypes, 
                                  c("A" = "A. SpA-TGC", 
                                    "B" = "B. Trophoblast cells", 
                                    "C" = "C. P-TGC",
                                    "D" = "D. Primtive endoderm cell",
                                    "E" = "E. Yolk sac epithelial cell",
                                    "F" = "F. Embryo stem cell",
                                    "G" = "G. Decidual stromal cell",
                                    "H" = "H. Decidula pericyte",
                                    "I" = "I. Endothelial cell",
                                    "J" = "J. Embryo stromal cell",
                                    "K" = "K. Immune cell",
                                    "L" = "L. Megakaryocyte",
                                    "M" = "M. Erythrocyte",
                                    "N" = "N. Hematopoietic cell"))

jpeg(paste(in_dir,"/cells_14.jpeg" ,sep=""),
     width = 1000, height = 480)

ggplot(cells_ids,aes(x=Celltypes,y=Frequency,fill=Samples)) + 
  geom_col(position="dodge") +
  theme_grey(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(name="Number of cells",
                     breaks = seq(0,5500, 1000),
                     limits=c(0,5500)) +
  xlab("Cell-types") 

dev.off()


MPL.obj <- NormalizeData(MPL.obj)

#MPL.obj  <- FindVariableFeatures(MPL.obj, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
MPL.obj <- FindVariableFeatures(MPL.obj, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(MPL.obj)) # 2000

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MPL.obj ), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(MPL.obj )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

MPL.obj <- ScaleData(MPL.obj, features = rownames(MPL.obj))

MPL.obj<- RunPCA(MPL.obj, features = VariableFeatures(object = MPL.obj))

DimPlot(MPL.obj, reduction = "pca")

ElbowPlot(MPL.obj)

#DimHeatmap(MPL.obj, dims = 1:2, cells = 500, balanced = TRUE)

MPL.obj <- RunUMAP(MPL.obj, dims = 1:10)


jpeg(paste(in_dir,"/all_cells_14.jpeg" ,sep=""),
     width = 700, height = 500)

d1<- DimPlot(MPL.obj, reduction = "umap",
        group.by = "cluster_14", label = T,
        pt.size = 1,
        label.size = 4) + ggtitle(NULL)

  d1+ scale_color_manual(labels = c("A. SpA-TGC", 
                                 "B. Trophoblast cells", 
                                 "C. P-TGC",
                                 "D. Primtive endoderm cell",
                                 "E. Yolk sac epithelial cell",
                                 "F. Embryo stem cell",
                                 "G. Decidual stromal cell",
                                 "H. Decidula pericyte",
                                 "I. Endothelial cell",
                                 "J. Embryo stromal cell",
                                 "K. Immune cell",
                                 "L. Megakaryocyte",
                                 "M. Erythrocyte",
                                 "N. Hematopoietic cell"),
                         values = colorRampPalette(brewer.pal(8, "Set1"))(14))


  # https://github.com/satijalab/seurat/issues/1664
  #DimPlot(pbmc_small, label =T, cells.highlight = Cells(pbmc_small)[1:10]) + 
    #scale_color_manual(labels = c("New Value1", "New Value2"), values = c("grey", "red")) +
   # labs(color = "legend title")
dev.off()

# =======================================
rm(list = ls())
gc()


setwd("/work/LAS/geetu-lab/sarwar")


in_dir<- "/work/LAS/geetu-lab/sarwar/data/10x_V2_processed_data"

# ========================== Trophoblast cells ===============

data1 <- read.csv(paste(in_dir,"/MTR15682-matrix.csv", sep=""),
                  header=T,
                  row.names = 1,
                  check.names = F)

dim(data1) # 18909X 15682

#View(data1[1:2,1:10])



data1.meta = read.csv(paste(in_dir,"/MTR15682-meta.csv", sep=""),
                      check.names = FALSE, header=TRUE,row.names=1)


#View(data1.meta)

levels(as.factor(data1.meta$cluster_14))

# sort


data1_sort <- data1[,match(rownames(data1.meta),colnames(data1))]

dim(data1_sort)
#View(data1_sort[1:2,1:5])


# ============== Create Seurat obj

MPL.obj <- CreateSeuratObject(counts = data1_sort,
                              meta.data = data1.meta, 
                              project = "MTR15682",
                              min.cells = 3, min.features = 500)


MPL.obj

# Get number of cells per cluster and per sample of origin
table(MPL.obj@meta.data$orig.ident)

# Get number of cells per cluster and per sample of origin
cells_ids <- data.frame(table(MPL.obj@meta.data$cluster_14, 
                              MPL.obj@meta.data$orig.ident))
#View(cells_ids)
#sum(cells_ids$Freq)

colnames(cells_ids) <- c("Celltypes","Samples","Frequency")


cells_ids$Celltypes  <- revalue(cells_ids$Celltypes, 
                                c("A" = "A. SpA-TGC", 
                                  "B" = "B. Trophoblast cells", 
                                  "C" = "C. P-TGC"))


jpeg(paste(in_dir,"/Trophoblast_filtered_14.jpeg" ,sep=""),
     width = 1000, height = 480)

ggplot(cells_ids,aes(x=Celltypes,y=Frequency,fill=Samples)) + 
  geom_col(position="dodge") +
  theme_grey(base_size = 18) + 
  geom_text(aes(label=Frequency), 
            position=position_dodge(width=0.9),
            vjust=-0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(name="Number of cells",
                     breaks = seq(0,3000, 1000),
                     limits=c(0,3000)) +
  xlab("Cell-types") 

dev.off()


MPL.obj <- NormalizeData(MPL.obj)

MPL.obj  <- FindVariableFeatures(MPL.obj , 
                                 selection.method = "vst", 
                                 nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MPL.obj ), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(MPL.obj )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

MPL.obj <- ScaleData(MPL.obj, features = rownames(MPL.obj))

MPL.obj<- RunPCA(MPL.obj, features = VariableFeatures(object = MPL.obj))

DimPlot(MPL.obj, reduction = "pca")

ElbowPlot(MPL.obj)

#DimHeatmap(MPL.obj, dims = 1:2, cells = 500, balanced = TRUE)

MPL.obj <- RunUMAP(MPL.obj, dims = 1:10)


jpeg(paste(in_dir,"/Trophoblast_filtered_14_umap.jpeg" ,sep=""),
     width = 700, height = 500)

d1<- DimPlot(MPL.obj, reduction = "umap",
             group.by = "cluster_14", label = T,
             pt.size = 1,
             label.size = 7) + ggtitle(NULL)

d1+ scale_color_manual(labels = c("A. SpA-TGC", 
                                  "B. Trophoblast cells", 
                                  "C. P-TGC"),
                       values = c("#E41A1C", "#864F70" ,"#3881AF"))

dev.off()

jpeg(paste(in_dir,"/Trophoblast_filtered_14_umap_split.jpeg" ,sep=""),
     width = 1000, height = 500)

d1<- DimPlot(MPL.obj, reduction = "umap",
             group.by = "cluster_14", label = F,
             split.by = "orig.ident",
             pt.size = 1,
             label.size = 7) + ggtitle(NULL)

d1+ scale_color_manual(labels = c("A. SpA-TGC", 
                                  "B. Trophoblast cells", 
                                  "C. P-TGC"),
values = c("#E41A1C", "#864F70" ,"#3881AF"))

dev.off()

# ============ 17 clusters
levels(as.factor(data1.meta$cluster_21))

# Get number of cells per cluster and per sample of origin
cells_ids <- data.frame(table(MPL.obj@meta.data$cluster_21, 
                              MPL.obj@meta.data$orig.ident))
#View(cells_ids)
#sum(cells_ids$Freq)

colnames(cells_ids) <- c("Celltypes","Samples","Frequency")


cells_ids$Celltypes  <- revalue(cells_ids$Celltypes, 
                                c("1" = "1. TSC and ExE cell",
                                  "10" = "10. Primary P-TGC",
                                  "11" = "11. Secondary P-TGC",
                                  "12" = "12. Secondary P-TGC Precursor",
                                  "14-15" = "14-15. SpT",
                                  "16" = "16. Gly-T",
                                  "17" = "17. SpA-TGC",
                                  "18" = "18. EPC Migratory Cell",
                                  "19" = "19. SynTII Precursor",
                                  "2" = "2. LaTP",
                                  "3" = "3. LaTP 2",
                                  "4-5" = "4-5. SynTI Precursor",
                                  "6-7-8" = "6-7-8. S-TGC Precursor",
                                  "9" = "9. S-TGC",
                                  "E1" = "E1. The progenitor of SpT",
                                  "E2" = "E2. The progenitor of S-TGC Precursor",
                                  "P1-2-3" = "P1-2-3. Biopotentila progenitor"))


jpeg(paste(in_dir,"/Trophoblast_filtered_17.jpeg" ,sep=""),
     width = 1000, height = 600)

ggplot(cells_ids,aes(x=Celltypes,y=Frequency,fill=Samples)) + 
  geom_col(position="dodge") +
  theme_grey(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(name="Number of cells",
                     breaks = seq(0,1500, 300),
                     limits=c(0,1500)) +
  xlab("Cell-types") 

dev.off()


MPL.obj <- NormalizeData(MPL.obj)

MPL.obj  <- FindVariableFeatures(MPL.obj , 
                                 selection.method = "vst", 
                                 nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MPL.obj ), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(MPL.obj )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

MPL.obj <- ScaleData(MPL.obj, features = rownames(MPL.obj))

MPL.obj<- RunPCA(MPL.obj, features = VariableFeatures(object = MPL.obj))

DimPlot(MPL.obj, reduction = "pca")

ElbowPlot(MPL.obj)

#DimHeatmap(MPL.obj, dims = 1:2, cells = 500, balanced = TRUE)

MPL.obj <- RunUMAP(MPL.obj, dims = 1:10)


jpeg(paste(in_dir,"/Trophoblast_filtered_17_umap.jpeg" ,sep=""),
     width = 700, height = 500)

d1<- DimPlot(MPL.obj, reduction = "umap",
             group.by = "cluster_21", label = T,
             pt.size = 1,
             label.size = 6) + ggtitle(NULL)

d1+ scale_color_manual(labels = c("1" = "1. TSC and ExE cell",
                                  "10" = "10. Primary P-TGC",
                                  "11" = "11. Secondary P-TGC",
                                  "12" = "12. Secondary P-TGC Precursor",
                                  "14-15" = "14-15. SpT",
                                  "16" = "16. Gly-T",
                                  "17" = "17. SpA-TGC",
                                  "18" = "18. EPC Migratory Cell",
                                  "19" = "19. SynTII Precursor",
                                  "2" = "2. LaTP",
                                  "3" = "3. LaTP 2",
                                  "4-5" = "4-5. SynTI Precursor",
                                  "6-7-8" = "6-7-8. S-TGC Precursor",
                                  "9" = "9. S-TGC",
                                  "E1" = "E1. The progenitor of SpT",
                                  "E2" = "E2. The progenitor of S-TGC Precursor",
                                  "P1-2-3" = "P1-2-3. Biopotentila progenitor"),
                       values = colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1:9,22:24,18,11:14)])

dev.off()

jpeg(paste(in_dir,"/Trophoblast_filtered_17_umap_split.jpeg" ,sep=""),
     width = 1000, height = 500)

d1<- DimPlot(MPL.obj, reduction = "umap",
             group.by = "cluster_21", label = F,
             split.by = "orig.ident",
             pt.size = 1,
             label.size = 6) + ggtitle(NULL)

d1+ scale_color_manual(labels = c("1" = "1. TSC and ExE cell",
                                  "10" = "10. Primary P-TGC",
                                  "11" = "11. Secondary P-TGC",
                                  "12" = "12. Secondary P-TGC Precursor",
                                  "14-15" = "14-15. SpT",
                                  "16" = "16. Gly-T",
                                  "17" = "17. SpA-TGC",
                                  "18" = "18. EPC Migratory Cell",
                                  "19" = "19. SynTII Precursor",
                                  "2" = "2. LaTP",
                                  "3" = "3. LaTP 2",
                                  "4-5" = "4-5. SynTI Precursor",
                                  "6-7-8" = "6-7-8. S-TGC Precursor",
                                  "9" = "9. S-TGC",
                                  "E1" = "E1. The progenitor of SpT",
                                  "E2" = "E2. The progenitor of S-TGC Precursor",
                                  "P1-2-3" = "P1-2-3. Biopotentila progenitor"),
                       values = colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1:9,22:24,18,11:14)])


dev.off()

# ============ renv
#https://rstudio.github.io/renv/articles/renv.html
#install.packages("renv")
#library(renv)

#renv::init()
#renv::snapshot()
