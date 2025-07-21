library(readr)
library(Signac)
library(Seurat)
library(SCP)
library(GenomicRanges)
library(ggplot2)
library(scales)
library(stringr)
library(patchwork)
library(AnnotationHub)
library(data.table)
library(SingleCellExperiment)
library(scDblFinder)
library(harmony)
library(clustree)
library(cols4all)
library(cowplot)
library(gridExtra)
library(GenomeInfoDb)
library(tidyverse)
library(dplyr)
library(readxl)
library(writexl)
library(BSgenome)
library(Biostrings)
library(BSgenome.Sscrofa.UCSC.susScr11)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(ggseqlogo)
library(ComplexHeatmap)
library(ggplotify)
library(colorRamp2)
library(org.Ss.eg.db)
library(ChIPseeker)
options(ChIPseeker.downstreamDistance = 2000)
library(GenomicFeatures)
library(ggpubr)
library(rtracklayer)
library(coin)
library(biomaRt)
set.seed(1234)

library(future)
nbrOfWorkers()
plan(multisession, workers = 8)
options(future.globals.maxSize = 30 * 1024^3)  # 50G 20*1024^3

####Data loading####
setwd("./atac/cellranger/sample/outs")
counts <- Read10X_h5(filename = "raw_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

peaks.keep <- seqnames(granges(atac)) %in% standardChromosomes(granges(atac))[1:18] #Chr filtered
atac <- atac[as.vector(peaks.keep), ]

ah <- AnnotationHub()
ssc_ensdb_109 <- ah[["AH109771"]]
annotations <- GetGRangesFromEnsDb(ensdb = ssc_ensdb_109)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "susScr11"
Annotation(atac) <- annotations

atac <- NucleosomeSignal(object = atac)
atac <- TSSEnrichment(object = atac, fast = FALSE)
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
atac$high.tss <- ifelse(atac$TSS.enrichment > 3, 'High', 'Low')
atac <- subset(atac, nCount_peaks > 200) # filtered
p1 <- DensityScatter(atac, x = 'nCount_peaks', y = 'TSS.enrichment',
                     log_x = TRUE, quantiles = TRUE) + ggtitle('Density_Scatter')
p2 <- FragmentHistogram(object = atac, group.by = 'nucleosome_group') + ggtitle('Fragment_Histogram')
p3 <- TSSPlot(atac, group.by = 'high.tss') + NoLegend() + ggtitle('TSSPlot')
p4 <- VlnPlot(object = atac,
              features = c('nCount_peaks', 'TSS.enrichment',
                           'nucleosome_signal', 'pct_reads_in_peaks'), # 'blacklist_ratio'
              pt.size = 0, ncol = 4)+ 
  NoLegend() +
  ggtitle(sample)
p <- (p1|p2|p3)/p4 + plot_layout(widths = c(3, 1, 1), heights = unit(c(2.5, 3), c('null', 'null')))
ggsave(paste0(path_tmp,sample,"_QC_metric.pdf"),p, width = 5*3)

cells <- as.data.frame(colnames(atac))
frags <- as.data.frame(rownames(atac))
fwrite(cells, file = paste0(path_tmp,sample,"_raw_cell.txt"), row.names = F, col.names = F)
fwrite(frags, file = paste0(path_tmp,sample,"_raw_frag.txt"), row.names = F, col.names = F)

atac_sc <- as.SingleCellExperiment(atac)
atac_sc <- scDblFinder(atac_sc, aggregateFeatures = TRUE, nfeatures = 25,
                       processing = "normFeatures", 
                       dbr = 0.08, # The expected doublet rate
                       dbr.sd = 0) # Disable the uncertainty around the doublet rate

atac_sc$doublet_logic <- ifelse(atac_sc$scDblFinder.class == "doublet", TRUE, FALSE)
atac$Class <- atac_sc$scDblFinder.class
atac$Score <- atac_sc$scDblFinder.score
atac$Weight <- atac_sc$scDblFinder.weighted
atac <- atac %>% RunTFIDF() %>% FindTopFeatures() %>% RunSVD()
DepthCor(atac)
ggsave(filename = "cor_depth_LSI.pdf",width = 6, height = 5)
atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:30)
total_cell <- length(atac$Class)
removed_cell <- sum(atac$Class == "doublet")
percent  <-  round(removed_cell/total_cell, 3)*100
p5 <- DimPlot(atac, group.by = "Class", cols = c("#FF7676","steelblue")) +
  labs(title = paste0("sample", "\nTotal: ",total_cell, "\nRemoved: ",
                      removed_cell, "\npercent: ", percent,"%"))
pdf("umap_pct.pdf")
p5
dev.off()

write_rds(atac, paste0("scdblfinder/",sample,"_atac_scDblFinder.rds"))

cells <- as.data.frame(colnames(atac))
frags <- as.data.frame(rownames(atac))
fwrite(cells, file = "scDblFinder_cell.txt", row.names = F, col.names = F, quote = F)
fwrite(frags, file = "scDblFinder_frag.txt", row.names = F, col.names = F, quote = F)

####snATAC mergeing####
# base_path
# cellranger_path
# scdblfinder_path
peaks_list <- list()
meta_list <- list()
dbl_list <- list()

for (i in 1:length(samples)){
  # Load peaks
  peaks_list[[samples[i]]] <- read.table(file = file.path(path_cellranger, samples[i],
                                                          "/outs/filtered_peak_bc_matrix/peaks.bed"), col.names = c("chr", "start", "end"))
  # Load metadata
  md <- read.table(
    file = file.path(path_cellranger,samples[i],"/outs/singlecell.csv"),
    stringsAsFactors = FALSE, sep = ",", header = TRUE,
    row.names = 1)[-1, ]
  md <- md[md$passed_filters > 1000, ]
  meta_list[[samples[i]]] <- md
  # Load filtered cell
  dbl_list[[samples[i]]] <- read.table(file = file.path(scdblfinder_path, samples[i], "_scDblFinder_cell.txt"))
}

peaks_list <- lapply(peak_files, load_peaks)
gr_list <- lapply(peaks_list, makeGRangesFromDataFrame)
gr_combined <- do.call(c, unname(gr_list))
combined.peaks <- reduce(gr_combined)
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
write_rds(combined.peaks, file.path(path_combined, "combined.peaks.rds"))

fra_path <- sapply(samples, function(s) paste0(path_cellranger, s, "/outs/fragments.tsv.gz"), USE.NAMES = TRUE)
create_fragments <- function(fra_path, cells) {
  frags <- CreateFragmentObject(path = fra_path, cells = cells)
  counts <- FeatureMatrix(fragments = frags,
                          features = combined.peaks,
                          cells = cells)
  list(frag = frags, count = counts)
}

scdbl_list <- lapply(dbl_list, `[[`, 1)
frag_data <- Map(create_fragments, fra_path, scdbl_list)
create_seurat <- function(count, frag, metadata, sample) {
  assay <- CreateChromatinAssay(count, fragments = frag)
  seurat_obj <- CreateSeuratObject(assay, assay = "peaks", meta.data = metadata)
  seurat_obj$sample <- sample
  seurat_obj
}

frag_arr <- lapply(frag_data, `[[`, 1)
count_arr <- lapply(frag_data, `[[`, 2)
seurat_objects <- Map(create_seurat, count_arr, frag_arr, meta_list, samples)

combined <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[2:length(samples)],
  add.cell.ids = samples)

Annotation(combined) <- annotations
combined <- NucleosomeSignal(object = combined)
combined <- TSSEnrichment(object = combined, fast = FALSE)
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
write_rds(combined,  file.path(path_combined,"combined_unfilter.rds"))

p7_1 <- VlnPlot(object = combined,
                features = c('nCount_peaks', 'TSS.enrichment',
                             'nucleosome_signal', 'pct_reads_in_peaks'), # 'blacklist_ratio'
                pt.size = 0, ncol = 4,
                cols = paletteer::paletteer_d("RColorBrewer::Paired")[1])+ NoLegend() + ggtitle("combined_before_filtered")

combined <- subset(
  x = combined,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2)

p7_2 <- VlnPlot(object = combined,
                features = c('nCount_peaks', 'TSS.enrichment',
                             'nucleosome_signal', 'pct_reads_in_peaks'), # 'blacklist_ratio'
                pt.size = 0, ncol = 4,
                cols = paletteer::paletteer_d("RColorBrewer::Paired")[1])+ NoLegend() + ggtitle("combined_after_filtered")
ggsave(paste0(path_combined,"combined_QC_metric.pdf"),p7_1/p7_2,width = 12, height = 12)
combined <- combined %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q75') %>% RunSVD() # use the top 25% all peaks
write_rds(combined, paste0(path_combined,"combined_filter_normalized.rds"))

# choose dims for lsi
p8_1 <- ElbowPlot(object = combined, ndims = 50, reduction = "lsi")
p8_2 <- DepthCor(combined, n = 30) # assess the correlation between each LSI component and sequencing depth
combined_umap <- RunUMAP(combined, reduction = 'lsi', dims = 2:30)
p8_3 <- DimPlot(combined_umap, group.by = 'sample', pt.size = 0.1)
p8 <- (p8_1|p8_2|p8_3)
ggsave(paste0(path_combined,"lsi_elbow_depthcor_dimplot.pdf"), p8, width = 7*3)

dim_use = 2:30
res = 0.2
combined_integrate <- RunHarmony(combined, group.by.vars = 'sample', reduction.use = 'lsi',
                                 assay.use = 'peaks', project.dim = FALSE, # Project dimension reduction loadings
                                 reduction.save = "harmony_lsi", dims.use = 2:50)

combined_integrate <- FindNeighbors(combined_integrate, reduction = "harmony_lsi", dims = dim_use)

# for(i in  seq(0.1,1,by=0.1)){
#  sce <- FindClusters(object = combined_integrate, verbose = F, resolution = seq(0.1,1,by=0.1))
# }
# clus.tree.out <- clustree(sce)+ theme(legend.position = "bottom")+ 
#   scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "green", high = "red")
# 
# ggsave(paste0(path_combined,"clustetree.pdf"), clus.tree.out, width = 3*5)

combined_integrate <- FindClusters(combined_integrate, resolution = res, cluster.name = "harmony_clusters_lsi")

## UMAP
combined_integrate <- RunUMAP(combined_integrate, reduction = "harmony_lsi", dims = dim_use, reduction.name = "umap_harmony_lsi")
levels(combined_integrate$harmony_clusters_lsi) <- c(0:(length(unique(combined_integrate$harmony_clusters_lsi))-1))

CellDimPlot(combined_integrate, group.by = c("sample"), reduction = "umap_harmony_lsi",
            combine = FALSE, label = F, palcolor = paletteer::paletteer_d("awtools::ppalette"))
ggsave(paste0(path_combined,"combined_peaks_umap_harmony.pdf"), width = 3*3)

CellDimPlot(combined_integrate, group.by = c("harmony_clusters_lsi"), reduction = "umap_harmony_lsi",
            combine = FALSE, label = F)
ggsave(paste0(path_combined,"combined_peaks_umap_cluster_", res, ".pdf"), width = 3*3)

# gene activity
gene.activities <- GeneActivity(combined_integrate)
combined_integrate[['activity']] <- CreateAssayObject(counts = gene.activities)
combined_integrate <- NormalizeData(object = combined_integrate,
                                    assay = 'activity',
                                    normalization.method = 'LogNormalize',
                                    scale.factor = median(combined_integrate$nCount_activity))
DefaultAssay(combined_integrate) <- 'activity'
marker_gene <- fread(paste0(path_base,'marker_gene.csv'), header = T) # Celltype marker
DotPlot(combined_integrate, features = split(marker_gene$gene, marker_gene$celltype)) +
  RotatedAxis()+ scale_color_continuous_c4a_seq('viridis',reverse = T) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(path_combined,"DotPlot_marker_", res, ".pdf"), width = 6*3)
write_rds(combined_integrate, paste0(path_combined,"combined_integrate_celltype.rds"))

combined_integrate <- RenameIdents(combined_integrate,
                                   "0"= "Myofibers",
                                   "1"= "Myofibers",
                                   "2"= "Myofibers",
                                   "3"= "MuSC/Myoblast",
                                   "4"= "FAPs",
                                   "5"= "B Cells",
                                   "6"= "Endothelial",
                                   "7"= "T Cells",
                                   "8" = "Smooth Muscle",
                                   "9" = "MuSC/Myoblast", 
                                   "10" = "Adipocyte") # for example

combined_integrate$celltype <- factor(x =Idents(combined_integrate), levels = sort(levels(combined_integrate@active.ident)))

CellDimPlot(combined_integrate, group.by = c("celltype"), reduction = "umap_harmony_lsi",
            combine = FALSE, label = F)
ggsave(paste0(path_combined,"combined_peaks_umap_celltype.pdf"), width = 3*3)

DotPlot(combined_integrate, features = split(marker_gene$gene, marker_gene$celltype)) + 
  RotatedAxis()+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_continuous_c4a_seq('viridis',reverse = T) + labs(x = "", y = "")
ggsave(paste0(path_combined,"DotPlot_marker.pdf"), width = 6*3)
write_rds(combined_integrate, paste0(path_combined,"combined_integrate_celltype_ann.rds"))

#### corresponding snRNA Data--integreation and imputation ####
snRNA <- read_rds(paste0(path_snRNA,'snRNA_selected.rds'))
snRNA_select <- subset(snRNA, subset = sample %in% smaples)
snRNA_select$sample <- factor(x = snRNA_select$sample,levels = ex_smaple)

# Normalization
snRNA_select <- snRNA_select %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
snRNA_select  <- snRNA_select  %>%
  FindNeighbors(.,reduction = "pca", dims = 1:30, verbose = FALSE) %>% 
  FindClusters(.,resolution= 0.1, verbose = FALSE, cluster.name ="unintegrated_clusters") %>%
  RunUMAP(.,reduction = "pca", dims = 1:30, verbose = FALSE, reduction.name = "unintegrated_umap")
CellDimPlot(snRNA_select, reduction = "unintegrated_umap", group.by = c("sample"))
ggsave(paste0(path_snRNA,"umap_unintegrated.pdf"), width = 3*3)

# Integration
snRNA_int <- RunHarmony(snRNA_select, group.by.vars="sample", ncores=6,plot_convergence = TRUE,
                        reduction.use = 'pca', reduction.save = "harmony_pca", project.dim = F)
snRNA_int  <- snRNA_int  %>%
  FindNeighbors(.,reduction = "harmony_pca", dims = 1:30, verbose = FALSE) %>% 
  FindClusters(.,resolution= 0.1, verbose = FALSE, cluster.name ="harmony_clusters") %>%
  RunUMAP(.,reduction = "harmony_pca", dims = 1:30, verbose = FALSE, reduction.name = "harmony_umap")
CellDimPlot(snRNA_int, reduction = "integrated_umap", ncol=3,
            group.by = c("harmony_clusters","cell_class_toplevel","cell_cluster"))
ggsave(paste0(path_snRNA,"harmony_integrated_umap.pdf"), width = 6*3)
write_rds(snRNA_int, paste0(path_snRNA,"snRNA_int.rds"))

table <- as.data.frame(table(snRNA_int@meta.data$sample,snRNA_int@meta.data$cell_class_toplevel))
names(table) <- c("Samples","Celltype","CellNumber")
ggplot(table, aes(x = Samples, weight = CellNumber, fill = Celltype))+
  geom_bar(position="fill")+ scale_fill_manual(values= paletteer::paletteer_d("RColorBrewer::Paired")) +
  theme_cowplot(font_size = 12) + labs(y="Percentage") + RotatedAxis()
ggsave(paste0(path_snRNA,'celltype_percentage_barplot.pdf'), p8, width = 2*3)

# Loading PCGs snRNA data
path_TSS = "./ref/Sus_scrofa_11.1_109/Sscrofa11.1_109.PCGs.gtf.TSS"
gene.info <- read.table(path_TSS, header = T)
snRNA_int <- read_rds(paste0(path_snRNA,"snRNA_int.rds"))
snRNA_sub <- snRNA_int[gene.info$gene_name[gene.info$gene_name %in% rownames(snRNA_int)], ]
p9_1 <- CellDimPlot(snRNA_sub, reduction = "unintegrated_umap", group.by = c("cell_class_toplevel")) +
  ggtitle("snRNA") 
write_rds(snRNA_sub, paste0(path_intsc,"snRNA_sub_PCG.rds"))

# scATAC data
combined_integrate <- read_rds(paste0(path_combined,"combined_integrate_celltype_ann.rds"))
DefaultAssay(combined_integrate) <- 'activity'
combined_integrate <- combined_integrate %>% FindVariableFeatures() %>% ScaleData()
p9_2 <- CellDimPlot(combined_integrate, reduction = "umap_harmony_lsi", group.by = "celltype") +
  ggtitle("snATAC") 

# Integrating with snRNA-seq data
transfer.anchors <- FindTransferAnchors(
  reference = snRNA_sub,
  query = combined_integrate,
  features = VariableFeatures(snRNA_sub),
  reference.assay = "RNA",
  query.assay = "activity",
  reduction = "cca",
  verbose = T,
  dims = 1:30)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = snRNA_sub$cell_class_toplevel,
  weight.reduction = combined_integrate[["harmony_lsi"]],
  dims = 2:30)

combined_integrate <- AddMetaData(combined_integrate, metadata = predicted.labels)
combined_integrate$predicted.id <- factor(x = combined_integrate$predicted.id,
                                          levels = sort(levels(snRNA_sub$cell_class_toplevel)))
p9_3 <- CellDimPlot(combined_integrate, reduction = "umap_harmony_lsi", group.by = c("predicted.id")) +
  ggtitle("snATAC_predicted_id_raw")

combined_integrate$mapped_celltype <- ifelse(
  combined_integrate$predicted.id %in% c("Endothelial(Blood)", "Endothelial(Lymphatic)"), "Endothelial",
  as.character(combined_integrate$predicted.id))
combined_integrate$mapped_celltype <- factor(combined_integrate$mapped_celltype,
                                             levels = levels(combined_integrate$celltype))
combined_integrate$ann_cor <- combined_integrate$mapped_celltype == combined_integrate$celltype
predictions <- table(combined_integrate$celltype, combined_integrate$mapped_celltype)
predictions <- as.data.frame(predictions/rowSums(predictions))

p9_4 <- ggplot(predictions, aes(x = Var1, y = Var2, fill = Freq)) + geom_tile(color = "white", lwd = 0.5)+ 
  viridis::scale_fill_viridis(option = "magma", name = "Fraction of cells", direction = -1) +
  xlab("Celltype annotation (snRNA)") + ylab("Predicted celltype") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) + coord_fixed()

write_rds(combined_integrate, paste0(path_intsc,"combined_integrate_predict_ann.rds"))

score_ann <- FetchData(combined_integrate, vars = c("prediction.score.max", "ann_cor"))
score_ann$ann_cor <- factor(x = score_ann$ann_cor,
                            levels = sort(unique(score_ann$ann_cor), decreasing = T))
# Density
p9_5 <- ggplot(score_ann, aes(prediction.score.max, fill = ann_cor, colour = ann_cor)) +
  geom_density(alpha = 0.5) +  theme(aspect.ratio = 1.1) + 
  scale_fill_discrete(name = "Annotation Correct",type = c("dodgerblue", "#E31A1CFF"),
                      labels = c(paste0("TRUE (n = ", length(which(combined_integrate$ann_cor == "TRUE")), ")"),
                                 paste0("FALSE (n = ", length(which(combined_integrate$ann_cor == "FALSE")), ")"))) + 
  scale_color_discrete(name = "Annotation Correct", type = c("dodgerblue", "#E31A1CFF"),
                       labels = c(paste0("TRUE (n = ", length(which(combined_integrate$ann_cor == "TRUE")), ")"), 
                                  paste0("FALSE (n = ", length(which(combined_integrate$ann_cor == "FALSE")), ")"))) + 
  xlab("Prediction Score") + ylab("Density") + theme_cowplot() +theme(axis.title.x = element_text(vjust = 0.5))

# score > 0.5, proportion
proportion_0.5 <- round(sum(combined_integrate@meta.data$prediction.score.max > 0.5) / 
                          length(combined_integrate@meta.data$prediction.score.max) * 100, 2)
p9_6 <- ggplot(combined_integrate@meta.data, aes(x=prediction.score.max)) + 
  geom_histogram(position = "identity", alpha = 0.5, fill = "#93deff",colour = "black") +
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = "black") + theme_cowplot() +
  ggtitle(paste0( "Total cell number = ",length(combined_integrate@meta.data$prediction.score.max), "\n",
                  "Prediction score > 0.5 cell number = ", sum(combined_integrate@meta.data$prediction.score.max > 0.5), "\n",
                  "Propotion = ", proportion_0.5, " %")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 12),
        axis.line = element_line(color = "black"), 
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_text(vjust = 0.5)) + xlab("Prediction Score")+ ylab("Count") + coord_flip() 

p9_7 <- CellDimPlot(combined_integrate, group.by = c("ann_cor"), reduction = "umap_harmony_lsi", 
                    palcolor =c("#E31A1CFF","dodgerblue")) + ggtitle("snATAC_predict_matched")

atac_sub <- subset(combined_integrate, prediction.score.max > 0.5)
# table(atac_sub$mapped_celltype)
p9_8 <- CellDimPlot(atac_sub, group.by = c("mapped_celltype"), reduction = "umap_harmony_lsi") +
  ggtitle("snATAC_predict_score > 0.5")

atac_sub = subset(atac_sub, subset = mapped_celltype != c("Adipocyte"))
p9_9 <- CellDimPlot(atac_sub, group.by = c("mapped_celltype"), reduction = "umap_harmony_lsi") +
  ggtitle("snATAC_predict_score > 0.5")

marker_gene <- fread(paste0(path_base,'marker_gene.csv'), header = T) # Celltype marker
marker_gene <- marker_gene[!grepl("Adipocyte", marker_gene$celltype), ]

p9 <- grid.arrange(p9_1, p9_2, p9_3,p9_4, p9_5, p9_6, p9_7, p9_8, p9_9, 
                   heights = c(1, 1.2, 1), ncol = 3)
ggsave(paste0(path_intsc, "atac_integrate_snRNA_predicted_umap.pdf"), p9, width = 20, height = 15)

atac_sub$mapped_celltype <- factor(x = atac_sub$mapped_celltype, levels = sort(unique(atac_sub$mapped_celltype)))
Idents(atac_sub) <- "mapped_celltype"

DotPlot(atac_sub, features = split(marker_gene$gene, marker_gene$celltype)) + 
  RotatedAxis()+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_continuous_c4a_seq('viridis',reverse = T) + labs(x = "", y = "")
ggsave(paste0(path_intsc, "atac_sub_dotplot.pdf"), width = 5*3)
write_rds(atac_sub, paste0(path_intsc,"atac_sub_predict_0.5.rds"))

## Integration Imputation
snRNA_sub <- read_rds( paste0(path_intsc,"snRNA_sub_PCG.rds"))
atac_sub <- read_rds(paste0(path_intsc,"atac_sub_predict_0.5.rds"))

snRNA_sub$tech <- "snRNA"
atac_sub$tech <- "snATAC"
snRNA_sub$mapped_celltype <- ifelse(
  snRNA_sub$cell_class_toplevel %in% c("Endothelial(Blood)", "Endothelial(Lymphatic)"), "Endothelial",
  as.character(snRNA_sub$cell_class_toplevel))
snRNA_sub$mapped_celltype <- factor(snRNA_sub$mapped_celltype,
                                    levels = sort(unique(snRNA_sub$mapped_celltype)))

genes.use <- VariableFeatures(snRNA_sub)
refdata <- GetAssayData(snRNA_sub, assay = "RNA", layer = "data")[genes.use, ]

transfer.anchors <- FindTransferAnchors(
  reference = snRNA_sub,
  query = atac_sub,
  features = VariableFeatures(snRNA_sub),
  reference.assay = "RNA",
  query.assay = "activity",
  reduction = "cca",
  verbose = T,
  dims = 1:30)

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,
                           weight.reduction = atac_sub[["harmony_lsi"]], dims = 2:30)

atac_sub[["RNA"]] <- imputation
multiome <- merge(x =snRNA_sub, y = atac_sub)
multiome <- ScaleData(multiome , features = genes.use, do.scale = FALSE) %>%
  RunPCA(., features = genes.use, verbose = FALSE) %>%
  RunUMAP(., dims = 1:30)
celltype_predict <- c(snRNA_sub$mapped_celltype, atac_sub$mapped_celltype)
celltype_predict <- celltype_predict[colnames(multiome)]
multiome$celltype_predict <- celltype_predict
p9_11 <- CellDimPlot(multiome, group.by = c("tech", "celltype_predict"), reduction = "umap")
ggsave(paste0(path_intsc, "int_umap.pdf"), p9_11, width = 5*3)
write_rds(multiome, paste0(path_intsc,"multiome_coembed.rds"))


######## MACS2-for 501bp peaks #########
# path_mac
# path_macs2_tmp

if (!dir.exists(paste0(path_snRNA,"mac"))) {
  dir.create(paste0(path_snRNA,"mac"))
} else {
  print("Dir already exists!")
}

if (!dir.exists(path_macs2_tmp)) {
  dir.create(path_macs2_tmp)
} else {
  print("Dir already exists!")
}

mac <- read_rds(paste0(path_intsc,"multiome_coembed.rds"))
mac_sub <- subset(mac, tech == "snATAC")
write_rds(mac_sub, paste0(path_mac,"mac_sub.rds"))
DefaultAssay(mac_sub) <- "peaks"
Idents(mac_sub) = "celltype_predict"
tempdir <- function() "./atac/snRNA/mac/macs2_output/"

# Secondly call peaks (merged atac_file), create tmp bed file for greenleaf.sh
peaks <- CallPeaks(object = mac_sub,
                   outdir = path_macs2_tmp, #tempdir()
                   group.by = "celltype_predict",
                   cleanup =FALSE)

# Quantify fragments in each peak
macs2_counts <- FeatureMatrix(fragments = Fragments(mac_sub),
                              features = peaks,
                              cells = colnames(mac_sub))

ah <- AnnotationHub()
ssc_ensdb_113 <- ah[["AH119485"]]
annotations <- GetGRangesFromEnsDb(ensdb = ssc_ensdb_113)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "susScr11"

mac_sub[["macs2_peaks"]] <- CreateChromatinAssay(counts = macs2_counts,
                                                 fragments = Fragments(mac_sub),
                                                 annotation = annotations)

# Get meta.data
if ("celltype_predict" %in% colnames(mac_sub@meta.data)) {
  metadata <- data.frame(
    sample = mac_sub$sample,
    celltype_predict = mac_sub$celltype_predict,
    stringsAsFactors = FALSE)
  metadata$sample = metadata$celltype_predict # only use cell_type
  metadata <- unique(metadata) 
  colnames(metadata) <- c("Sample","Group")
  metadata <- as.data.frame(lapply(metadata, function(x) {
    gsub("[[:space:]/]", "_", x)
  }))
  fwrite(metadata,paste0(path_mac,"metadata.txt"),row.names=F,sep="\t")
} else {
  stop("'celltype_predict' unexist")
}

write_rds(mac_sub, paste0(path_mac,"macs2.rds"))

# https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging

source_file <- "./software/ATAC_IterativeOverlapPeakMerging/greenleaf.sh"
target_file <- paste0(path_mac,"greenleaf.sh")

if (file.copy(source_file, target_file)) {
  print("success")
} else {
  print("failed")
}

lines <- readLines(target_file) %>% 
  gsub("meta_file", paste0(path_mac,"metadata.txt"), .) %>% 
  gsub("output_dir", paste0(path_mac,"macs2_output/"), .)
writeLines(lines, target_file)

# ATAC_IterativeOverlapPeakMerging
system(target_file)
# 501bp peak
atac = mac_sub
granges_501 <- read_rds(paste0(path_macs2_tmp,"All_Samples.fwp.filter.non_overlapping.rds"))

macs2_counts_501bp <- FeatureMatrix(
  fragments = Fragments(atac), # from cellranger fragment result
  features = granges_501,
  cells = colnames(atac))

atac[["macs2_peaks_501bp"]] <- CreateChromatinAssay(
  counts = macs2_counts_501bp,
  fragments = Fragments(atac),
  annotation = annotations)

DefaultAssay(atac) <- 'macs2_peaks_501bp'
atac <- atac %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q75') %>% RunSVD()
write_rds(atac, paste0(path_mac,"macs2_501bp.rds"))

# Visualize the cell-type-specific MACS2 peak calls alongside the 10x Cellranger peak calling
p10 <- CoveragePlot(object = atac, region = "MYH1", ranges = peaks, ranges.title = "MACS2") &
  scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::Paired"))
ggsave(paste0(path_mac,"CoveragePlot_MYH1_macs2.pdf"), p10 ,width=3*3)

# Gene activity for macs2_peaks_501bp
gene.activities.501 <- GeneActivity(atac)
atac[['macs2_501bp_ACTIVITY']] <- CreateAssayObject(counts = gene.activities.501)
atac <- NormalizeData(object = atac,
                      assay = 'macs2_501bp_ACTIVITY',
                      normalization.method = 'LogNormalize',
                      scale.factor = median(atac$nCount_macs2_peaks_501bp))

#### Chromvar ####
# Add motif information
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", 
                                tax_group = 'vertebrates', 
                                all_versions = FALSE))

atac <- AddMotifs(object = atac,
                  genome = BSgenome.Sscrofa.UCSC.susScr11,
                  pfm = pfm)

motif <- LayerData(atac, layer = "motifs")

atac <- RunChromVAR(object = atac, genome = BSgenome.Sscrofa.UCSC.susScr11)
# DefaultAssay(atac) <- 'chromvar'

MotifPlot(object = atac, motifs = "MAF", assay = 'peaks')
FeaturePlot(atac, features= "MAF", order = T) + ggtitle(paste0("Chromvar MAF activity")) +
  scale_color_gradientn(colours = c("#56B4E9", "white", "#e74a32"))
write_rds(atac, paste0(path_mac,"atac_motif_chromvar.rds"))

#### Cicero ####
# path_base
# path_snRNA
# path_mac
atac <- read_rds(paste0(path_mac,"atac_motif_chromvar.rds"))
DefaultAssay(atac) <- "macs2_peaks_501bp"

# convert to CDS format and make the cicero object
atac.cds <- as.cell_data_set(x = atac)
atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)

# genome informations from BSgenome.Sscrofa.UCSC.susScr11 
genome <- seqlengths(BSgenome.Sscrofa.UCSC.susScr11)
genome <- genome[seq(1:18)]
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num = 100)
write_rds(conns, paste0(path_mac,"conns.rds"))

# Find cis-co-accessible networks (CCANs)
ccans <- generate_ccans(conns)

# Add links to a Seurat object
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(atac) <- links
write_rds(atac, paste0(path_mac,"atac_cicero_links.rds"))

# Create a column that identifies which connections belong to a CCAN
ccan1 <- left_join(conns, ccans, by=c("Peak1" = "Peak"))
colnames(ccan1)[4] <- "CCAN1"
ccan2 <- left_join(conns, ccans, by=c("Peak2" = "Peak"))
colnames(ccan2)[4] <- "CCAN2"
df <- cbind(ccan1, CCAN2=ccan2$CCAN2) %>%
  mutate(CCAN = ifelse(CCAN1 == CCAN2, CCAN1, 0)) %>%
  select(-CCAN1, -CCAN2)
fwrite(df, file = paste0(path_mac, "ciceroConns.allcells.csv"), row.names = TRUE)

df_filtered <- df %>% filter(coaccess > 0.2) 
plots <- list()
ct = unique(atac$predicted.id)
for (i in 1:length(unique(atac$predicted.id))) {
  atac_sub <- subset(atac, subset = predicted.id == ct[i])
  peak_counts <- GetAssayData(atac_sub, assay = "macs2_peaks_501bp", layer = "counts")
  accessible_peaks <- rownames(peak_counts)[Matrix::rowSums(peak_counts > 0) > 10]
  # filter accessible_peaks within links
  ct_links <- df_filtered %>%
    filter(Peak1 %in% accessible_peaks & 
             Peak2 %in% accessible_peaks & 
             CCAN != 0)
  Links(atac_sub) <- ConnectionsToLinks(ct_links)
  # CoveragePlot for the subset
  p <- CoveragePlot(
    object = atac_sub,
    region = c("chr12-55185433-55352087"),
    assay = "macs2_peaks_501bp",
    split.by = "predicted.id") & scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::Paired")[i])
  plots[[ct[i]]] <- p
}

combined_plot <- wrap_plots(plots, ncol = 1)
ggsave(paste0(path_mac, "celltype_links_myh3-4_CoveragePlot.pdf"), combined_plot, width = 8, height = 18)

pcover <- CoveragePlot(atac, region = c("chr12-55185433-55352087"), assay = "macs2_peaks_501bp", split.by = "predicted.id") &
  scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::Paired"))
ggsave(paste0(path_mac, "myh3-4_CoveragePlot.pdf"), pcover, width = 8, height = 9)

#### DAR-celltype ####
# loading function
source_all_R <- function(dir) {
  files <- list.files(path = dir, pattern = "\\.R$", full.names = TRUE)
  for (f in files) {
    message(">> Loading ", f)
    tryCatch({
      source(f)
    }, error = function(e) {
      message("!! Error loading ", f, ": ", e$message)
    })
  }
}
source_all_R("./software/chipseeker_anno") #Loading specific promoter annotated function
color <- paletteer::paletteer_d("RColorBrewer::Paired")

# convert gtf/gff3 to TxDb file
txdb <- makeTxDbFromGFF(paste0(path_ref,"Sscrofa11.1_109.gtf"))
head(seqlevels(txdb))
peak_region <- gsub("chr", "", rownames(atac_all)) %>% StringToGRanges(., sep = c("-","-")) 
peakAnno <- annotatePeak_c57(peak_region, tssRegion=c(-2200, 500), TxDb=txdb, annoDb="org.Ss.eg.db")
fwrite(as.data.frame(peakAnno), file = paste0(path_dff,"all_celltype_chipseeker_peakanno.txt"), sep="\t", row.names = F)

pdf(paste0(path_dff,"all_celltype_pie_501bp_peaks_ann.pdf"), width = 10, height = 6)
plotAnnoPie(peakAnno)
dev.off()

pdf(paste0(path_dff,"all_celltype_bar_501bp_peaks_ann.pdf"), width = 8, height =3)
plotAnnoBar(peakAnno)
dev.off()

pdf(paste0(path_dff,"all_celltype_peaks_501bp_peakAnno_DistToTss.pdf"), width = 8, height = 3)
plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci/nrelative to TSS")
dev.off()

create_dirs(path_dff, c("all_celltype"))
path_all <- paste0(path_dff, "all_celltype/")
DefaultAssay(atac_all) <- "macs2_peaks_501bp"
Idents(atac) = "cell_type"
da.peak <- FindAllMarkers(object = atac_all, test.use = 'wilcox', #test.use = 'LR',
                          logfc.threshold = 1, min.pct = 0.1, only.pos = T) %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
fwrite(da.peak, paste0(path_all,"da.peak_wilcox.csv"), quote = F, row.names = F, sep = "\t")
top.da.peak <- da.peak[which(da.peak$p_val_adj <= 0.01), ] # %>% top_n (n = 500,wt = avg_log2FC)
fwrite(top.da.peak, paste0(path_all,"significant0.05_da.peak_wilcox.csv"), quote = F, row.names = F, sep = "\t")
table(top.da.peak$cluster)

peak_region <- gsub("chr", "", top.da.peak$gene) %>% StringToGRanges(., sep = c("-","-")) 
peakAnno <- annotatePeak_c57(peak_region, tssRegion=c(-2200, 500), TxDb=txdb, annoDb="org.Ss.eg.db")
fwrite(as.data.frame(peakAnno), file = paste0(path_all,"all_celltype_da_chipseeker_peakanno.txt"), sep="\t", row.names = F)

pdf(paste0(path_all,"all_celltype_da_pie_501bp_peaks_ann.pdf"), width = 10, height = 6)
plotAnnoPie(peakAnno)
dev.off()

peakAnnoList <- list()
for (cell_type in unique(atac_all$cell_type)){
  cell_type_peaks <- top.da.peak %>% filter(cluster == cell_type)
  peak_region_celltype <- gsub("chr", "", cell_type_peaks$gene) %>% StringToGRanges(., sep = c("-","-")) 
  peakAnno_celltype <- annotatePeak_c57(peak_region_celltype, tssRegion=c(-2200, 500), TxDb=txdb, annoDb="org.Ss.eg.db")
  fwrite(as.data.frame(peakAnno_celltype), file = paste0(path_all,"chipseeker_significant_peakanno_", gsub("[[:space:]/]", "_", cell_type),".txt"),
         sep="\t", row.names = F)
  peakAnnoList[[cell_type]] <- peakAnno_celltype
}

cell_num = as.data.frame(table(top.da.peak$cluster))
names(cell_num) = c("cell_type", "num")
cell_num$cell_type = factor(cell_num$cell_type, levels = cell_num$cell_type[order(cell_num$num, decreasing = TRUE)])
p <- ggbarplot(cell_num, x = "cell_type", y = "num", label = TRUE, lab.vjust = -0.3,fill = "cell_type",
               color = "cell_type", palette = paletteer::paletteer_d("RColorBrewer::Paired")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        axis.title.x = element_blank()) + ylab("Cell-type-specific ACR counts")
ggsave(paste0(path_all,"barplot_celltype_specific_ACR_counts.pdf"),p)


cell_types <- Idents(atac)
peak_matrix <- GetAssayData(atac, assay = "macs2_peaks_501bp", layer = "counts")

peak_counts <- sapply(levels(cell_types), function(ct) {
  cells <- names(cell_types)[cell_types == ct]
  sum(rowSums(peak_matrix[, cells, drop = FALSE] > 0) > 0)
})

result <- data.frame(
  cell_type = names(peak_counts),
  num = as.numeric(peak_counts),
  row.names = NULL)

result$cell_type = factor(result$cell_type, levels = result$cell_type[order(result$num, decreasing = TRUE)])
p <- ggbarplot(result, x = "cell_type", y = "num", label = TRUE, lab.vjust = -0.3,fill = "cell_type",
               color = "cell_type", palette = paletteer::paletteer_d("RColorBrewer::Paired")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        axis.title.x = element_blank()) + ylab("Cell-type-ALL ACR counts")
ggsave(paste0(path_dff,"barplot_celltype_ALL_ACR_counts.pdf"),p)

merged_df <- merge(cell_num, result, by = "cell_type", suffixes = c("_df_peaks", "_peaks"))
merged_df$pct_peaks_celltype <- merged_df$num_df_peaks / merged_df$num_peaks
fwrite(merged_df,paste0(path_dff,"da.peaks_ACR_pct.txt"),sep="\t")

#### myofiber subcluster DAR ####
# path_sub
# path_wilcox
DefaultAssay(atac) <- "macs2_peaks_501bp"
atac_sub <- subset(atac, cell_type == "Myofibers")
atac_sub@reductions <- list()
atac_sub@meta.data <- atac_sub@meta.data[, !grepl("^predict", colnames(atac_sub@meta.data))]
atac_sub <- atac_sub %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = "q75") %>% RunSVD()
atac_sub <- RunHarmony(atac_sub, group.by.vars = 'sample', reduction.use = 'lsi', assay.use = 'macs2_peaks_501bp', project.dim = FALSE, reduction.save = "sub.harmony.lsi")
DefaultAssay(atac_sub) <- "macs2_501bp_ACTIVITY"
atac_sub <- NormalizeData(atac_sub, normalization.method = "LogNormalize", scale.factor = 10000)
atac_sub <- FindVariableFeatures(atac_sub, selection.method = "vst", nfeatures = 1500)
atac_sub <- ScaleData(atac_sub, features = rownames(atac_sub))
atac_sub <- RunPCA(atac_sub, features = VariableFeatures(object = atac_sub))
atac_sub <- RunHarmony(atac_sub, group.by.vars = "sample", reduction.save = "sub.harmony.pca")
atac_sub <- FindNeighbors(atac_sub, reduction = "sub.harmony.pca", dims = 1:40)
atac_sub <- FindClusters(atac_sub,  cluster.name = "sub.harmony.pca_clusters", resolution = 0.1)
atac_sub <- RunUMAP(atac_sub, reduction = "sub.harmony.pca", dims = 1:40, reduction.name = "umap.sub.harmony.pca")
write_rds(atac_sub, paste0(path_snRNA, 'atac_sub.rds'))

rna <- read_rds(paste0('/data/ZhengYundi/result/atac_LW/snRNA/', "snRNA_selected.rds"))
rna_sub <- subset(rna, cell_cluster %in% c("Myofibers_I", "IIA", "IIB"))
rna_sub <- FindVariableFeatures(object = rna_sub, nfeatures = 1500)

# compute anchors between RNA and ATAC
transfer.anchors <- FindTransferAnchors(
  reference = rna_sub,
  query = atac_sub,
  features = VariableFeatures(rna_sub),
  reference.assay = "RNA",
  query.assay = "macs2_501bp_ACTIVITY",
  reduction = "cca",
  verbose = T,
  dims = 1:30) 

predicted.id.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna_sub$cell_cluster,
  weight.reduction = atac_sub[['sub.harmony.lsi']],
  dims = 2:30)

atac_sub <- AddMetaData(object = atac_sub, metadata = predicted.id.labels)

proportion_0.7 <- round(sum(atac_sub@meta.data$prediction.score.max > 0.7) /
                          length(atac_sub@meta.data$prediction.score.max) * 100, 2)
ggplot(atac_sub@meta.data, aes(x=prediction.score.max)) +
  geom_histogram(position = "identity", alpha = 0.5, fill = "#93deff",colour = "black") +
  geom_vline(xintercept = 0.7, linetype = 'dashed', color = "black") + theme_cowplot() +
  ggtitle(paste0( "Total cell number = ",length(atac_sub@meta.data$prediction.score.max), "\n",
                  "Prediction score > 0.7 cell number = ", sum(atac_sub@meta.data$prediction.score.max > 0.7), "\n",
                  "Propotion = ", proportion_0.7, " %")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_text(vjust = 0.5)) + xlab("Prediction Score")+ ylab("Count") + coord_flip()
ggsave(paste0(path_snRNA,"cell_predicted_score_0.7.pdf"), width = 10, height = 10)

atac_sub <- subset(atac_sub, prediction.score.max > 0.7)
Idents(atac_sub) <- "predicted.id"
atac_sub$predicted.id <- factor(x =Idents(atac_sub), levels = c("Myofibers_I", "IIA", "IIB" ))

atac_1 <- subset(atac_sub, tissue =="LDM")
p0 <- CellDimPlot(atac_1, reduction = "umap.sub.harmony.pca", group.by= "predicted.id", theme = ggplot2::theme_classic, theme_args = list(base_size = 16))
ggsave(paste0(path_snRNA,"cell_predicted_umap_LDM.pdf"), p0, width = 8, height = 8)

p1 <- DimPlot(atac_sub, reduction = "umap.sub.harmony.pca", cells.highlight = WhichCells(atac, expression = tissue == "LDM"), cols.highlight = "#62BBC3", cols = "gray", sizes.highlight = F, label.size = 5) + NoLegend() + labs(x = "UMAP 1", y = "UMAP 2", title = "LDM(n=8192)")    
p2 <- DimPlot(atac_sub, reduction = "umap.sub.harmony.pca", cells.highlight = WhichCells(atac, expression = tissue == "PM"), cols.highlight = "#E07F70", cols = "gray", sizes.highlight = F, label.size = 5) + NoLegend() + labs(x = "UMAP 1", y = "UMAP 2", title = "PM(n=6610)")   
ggsave(paste0(path_snRNA, "cell_predicted_umap_tissue.pdf"), p1|p2, width = 5.6*2, height = 6)
write.csv(table(atac_sub$predicted.id, atac_sub$tissue), paste0(path_snRNA,"tissue_cells_num.csv"),quote = F)

marker_genes = c("MYH7", "ATP2A2", "ENSSSCG00000029441", "MYH4")
p3 <- DotPlot(object = atac_sub, features = marker_genes) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(path_snRNA,"cell_predicted_dotplot.pdf"),p3, width= 6,height = 6)
p4 <- VlnPlot(atac_sub, features = marker_genes, pt.size = 0, ncol = 2)
p5 <- FeaturePlot(atac_sub, features= marker_genes, order = T, ncol = 2)
ggsave(paste0(path_snRNA,"cell_predicted_vlnplot_feature.pdf"), p4/p5, width= 10,height = 18)
write_rds(atac_sub, paste0(path_snRNA, 'Myo_sub.rds'))

atac <- read_rds(paste0(path_snRNA,'Myo_sub.rds'))
DefaultAssay(atac) <- "macs2_peaks_501bp"
atac <- RenameIdents(atac,
                     "Myofibers_I"="Myofibers_I",
                     'IIA' = 'Myofibers_II', 
                     'IIB' = 'Myofibers_II') # I vs II, IIA vs IIB
atac$merged_celltype <- factor(x =Idents(atac), levels = unique(Idents(atac)))
Idents(atac) = "merged_celltype"
da.peak <- FindAllMarkers(object = atac, test.use = 'wilcox', #test.use = 'LR',
                          logfc.threshold = 1, min.pct = 0.1, only.pos = T) %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
fwrite(da.peak, paste0(path_wilcox,"da.peak_wilcox.csv"), quote = F, row.names = F, sep = "\t")
top.da.peak <- da.peak[which(da.peak$p_val_adj <= 0.01), ] # %>% top_n (n = 500,wt = avg_log2FC)
fwrite(top.da.peak, paste0(path_wilcox,"significant0.05_da.peak_wilcox.csv"), quote = F, row.names = F, sep = "\t")
table(top.da.peak$cluster)

#significant peaks anno
peak_region <- gsub("chr", "", top.da.peak$gene) %>% StringToGRanges(., sep = c("-","-")) 
peakAnno <- annotatePeak_c57(peak_region, tssRegion=c(-2200, 500), TxDb=txdb, annoDb="org.Ss.eg.db")
fwrite(as.data.frame(peakAnno), file = paste0(path_wilcox,"chipseeker_significant_peakanno.txt"), sep="\t", row.names = F)

pdf(paste0(path_wilcox,"pie_501bp_significant_subcluster_peaks_ann.pdf"), width = 10, height = 6)
plotAnnoPie(peakAnno)
dev.off()

peakAnnoList <- list()
for (cell_type in unique(atac$merged_celltype)){
  cell_type_peaks <- top.da.peak %>% filter(cluster == cell_type)
  peak_region_celltype <- gsub("chr", "", cell_type_peaks$gene) %>% StringToGRanges(., sep = c("-","-")) 
  peakAnno_celltype <- annotatePeak_c57(peak_region_celltype, tssRegion=c(-2200, 500), TxDb=txdb, annoDb="org.Ss.eg.db")
  fwrite(as.data.frame(peakAnno_celltype), file = paste0(path_wilcox,"chipseeker_significant_peakanno_", gsub("[[:space:]/]", "_", cell_type),".txt"),
         sep="\t", row.names = F)
  peakAnnoList[[cell_type]] <- peakAnno_celltype
  pdf(paste0(path_wilcox,"pie_501bp_significant_peaks_ann_",cell_type,".pdf"), width = 10, height = 7)
  plotAnnoPie(peakAnnoList[[cell_type]])
  dev.off()
}

colors <- setNames(c(color[1:length(peakAnno@annoStat$Feature)]), peakAnno@annoStat$Feature)

pdf(paste0(path_wilcox,"celltype_barplot_501bp_significant_subcluster_peaks_ann.pdf"), width = 8, height =3)
plotAnnoBar(peakAnnoList) + scale_fill_manual(values = colors) + guides(fill = guide_legend(reverse = TRUE))
dev.off()

# barplot for cell-specific accessible chromatin regions counts
cell_num = as.data.frame(table(top.da.peak$cluster))
names(cell_num) = c("cell_type", "num")
cell_num$cell_type = factor(cell_num$cell_type, levels = cell_num$cell_type[order(cell_num$num, decreasing = TRUE)])
ggbarplot(cell_num, x = "cell_type", y = "num", label = TRUE, lab.vjust = -0.3,fill = "cell_type",
          color = "cell_type", palette = paletteer::paletteer_d("RColorBrewer::Paired")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        axis.title.x = element_blank()) + ylab("Cell-type-specific ACR counts")
ggsave(paste0(path_wilcox,"barplot_celltype_specific_ACR_counts.pdf"))

# heatmap
# significant peaks z-score
dis_peak <- top.da.peak %>% dplyr::select("gene") %>% distinct()
atac_aver <- AverageExpression(atac_pre, features = dis_peak$gene, assays = "macs2_peaks_501bp")
atac_aver <- as.matrix(atac_aver[["macs2_peaks_501bp"]])
atac_aver_scaled <- t(scale(t(atac_aver)))
# closest peaks
closest_genes <- ClosestFeature(object = atac, regions = rownames(atac_aver_scaled))
marker_gene <- fread(paste0("/data/ZhengYundi/result/atac/",'marker_gene.csv'), header = T)
marker_peaks <- closest_genes[closest_genes$gene_name %in% marker_gene$gene, ] 
marker_peaks <- marker_peaks[!duplicated(marker_peaks$gene_name), ]

row_number <- which(rownames(atac_aver_scaled) %in% marker_peaks$query_region)
labels = paste0(marker_peaks$gene_name,", ","peak region = ", marker_peaks$query_region)
row_anno <- rowAnnotation(foo = anno_mark(row_number, labels = labels,
                                          link_width = unit(5, "mm"), 
                                          padding = unit(2, "mm"),
                                          extend = unit(0, "mm")))
col_bar <- colorRamp2(c(-2, 0, 2), c("#56B4E9", "white", "#e74a32"))
p11_2 <- as.grob(Heatmap(atac_aver_scaled, name = "ACRs z-score", col = col_bar,
                         cluster_rows = F, cluster_columns = F,show_row_names = F, right_annotation = row_anno))
ggsave(paste0(path_wilcox,"Heatmap_maker_peaks.pdf"), p11_2, width= 10, height = 12)

# find peaks open in cells
open.peaks <- AccessiblePeaks(atac)
meta.feature <- GetAssayData(atac, assay = "macs2_peaks_501bp", layer = "meta.features")
# match the overall GC content in the peak set
peaks.matched <- MatchRegionStats( meta.feature = meta.feature[open.peaks, ],
                                   query.feature = na.omit(meta.feature[top.da.peak$gene, ]),
                                   n = nrow(meta.feature))

for (i in 1:length(unique(top.da.peak$cluster))) {
  celltype <- unique(top.da.peak$cluster)[i]
  safe_celltype <- gsub("[[:space:]/]", "_", celltype)
  DAP_sub <- subset(top.da.peak, cluster == celltype)
  fwrite(DAP_sub, file = file.path(paste0(path_wilcox,"DAR/", safe_celltype, "_DAR.csv")), row.names = FALSE, quote = FALSE, sep = "\t")
  # test enrichment in specific cells
  enriched.motifs <- FindMotifs(object = atac, background = peaks.matched, features = DAP_sub$gene)
  # sort enriched.motifs with fold.enrichment
  enriched.motifs <- enriched.motifs[order(enriched.motifs$fold.enrichment, decreasing = TRUE), ]
  fwrite(enriched.motifs, file = file.path(path_wilcox,"motif/", paste0(safe_celltype, "_motif.csv")), sep = "\t", quote = F)
  top.enriched.motifs <- enriched.motifs[which(enriched.motifs$pvalue < 0.05), ]
  fwrite(top.enriched.motifs, file = file.path(path_wilcox,"motif/", paste0(safe_celltype, "_p0.05_motif.csv")), sep = "\t", quote = F)
  MotifPlot(object = atac, motifs = head(rownames(top.enriched.motifs), n = 10), ncol = 2)
  ggsave(filename = file.path(path_wilcox,"motif/", paste0(safe_celltype, "_motif.pdf")), height = 8)
  # find closest gene
  closest_genes <- ClosestFeature(atac, regions = DAP_sub$gene)
  fwrite(closest_genes, file = paste0(path_wilcox,"closest_gene/", safe_celltype, "_closest_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  only.gene.names <- na.omit(closest_genes$gene_name) %>% str_subset(.,pattern="^.+$") %>% sort() %>% unique()
  fwrite(as.data.frame(only.gene.names), file = paste0(path_wilcox,"closest_gene/", safe_celltype, "_closest_genename.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
  # keep closest gene distance <= 2000 bp
  closest_genes_2kb <- subset(closest_genes,distance <= 2000)
  fwrite(closest_genes_2kb, file = paste0(path_wilcox,"closest_gene/", safe_celltype, "_closest_genes_2kb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  only.gene.names_2kb <- na.omit(closest_genes_2kb$gene_name) %>% str_subset(.,pattern="^.+$") %>% sort() %>% unique()
  fwrite(as.data.frame(only.gene.names_2kb), file = paste0(path_wilcox,"closest_gene/", safe_celltype, "_closest_genename_2kb.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
  # Transfer to human gene set 
  human <- read.table("./database/TF_for_SCENIC/human.pig.ortholog.109.one2one.txt", header = TRUE, sep = "\t", fill = TRUE, comment.char = "")
  pig <- fread(paste0(path_wilcox,"closest_gene/", safe_celltype, "_closest_genename.txt"))
  colnames(pig) = "Gene.name"
  pig2human = merge(pig, human) 
  pig2human_genename <- pig2human$Human.gene.name %>% str_subset(.,pattern="^.+$") %>% sort() %>% unique()
  fwrite(as.data.frame(pig2human_genename), file = paste0(path_wilcox, "metascape/", safe_celltype, "_METASCAPE.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
  # Transfer to 2kb human gene set
  pig <- fread(paste0(path_wilcox,"closest_gene/", safe_celltype, "_closest_genename_2kb.txt"))
  colnames(pig) = "Gene.name"
  pig2human = merge(pig, human)
  pig2human_genename_2k  <- pig2human$Human.gene.name %>% str_subset(.,pattern="^.+$") %>% sort() %>% unique()
  fwrite(as.data.frame(pig2human_genename_2k), file = paste0(path_wilcox, "metascape/", safe_celltype, "_METASCAPE_2kb.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
  # # merge DAR file and chipseeker
  chip <- fread(paste0(path_wilcox,"chipseeker_significant_peakanno_", safe_celltype, ".txt"))
  chip$gene <- paste0("chr", paste(chip$seqnames, chip$start, chip$end, sep = "-"))
  mer <- merge(chip, DAP_sub)
  fwrite(mer, file = paste0(path_wilcox, "DAR/", safe_celltype, "_merge_chipseek_DAR.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

#### DAR between tissue ####
atac <- read_rds(paste0(path_mac,"atac_motif_chromvar.rds"))
DefaultAssay(atac) = "macs2_peaks_501bp"
Idents(atac) = "cell_type"
atac$breed <- factor(atac$breed, levels = c("LD", "RC"))

DE = run_de(atac, de_family = "pseudobulk", replicate_col = "sample", de_method = "edgeR", de_type = "QLF",
            cell_type_col = "cell_type", label_col = "breed", n_threads = 16)
fwrite(DE, file = paste0(path_edger,"edgeR_QLF_DAR.csv"), row.names = F, quote = F, sep = "\t")

DE_sig = subset(DE, p_val_adj < 0.01 & abs(avg_logFC) > 0.25)
fwrite(DE_sig, file = paste0(path_edger, "edgeR_DAR_sig.csv"), row.names = F, quote = F, sep = "\t")

DE_diff <- as.data.frame(DE_sig[order(DE_sig$p_val_adj, DE_sig$avg_logFC, decreasing = c(FALSE, TRUE)), ])
DE_diff[which(DE_diff$avg_logFC >= 0.25 & DE_diff$p_val_adj < 0.01),'sig'] <- 'up'
DE_diff[which(DE_diff$avg_logFC <= 0.25 & DE_diff$p_val_adj < 0.01),'sig'] <- 'down'
DE_diff[which(abs(DE_diff$avg_logFC) <= 0.25 | DE_diff$p_val_adj => 0.01),'sig'] <- 'none'
fwrite(DE_diff, file = paste0(path_edger, "edgeR_DAR_sig_ann.csv"), row.names = F, quote = F, sep = "\t")

top.da.peak = DE_diff
summarize <- DE_diff %>% group_by(cell_type, sig) %>% summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = sig, values_from = count, values_fill = 0)
fwrite(summarize, file = paste0(path_edger, "counts_edgeR_DAR_summarize.csv"), row.names = F, quote = F, sep = "\t")

open.peaks <- AccessiblePeaks(atac)
meta.feature <- GetAssayData(atac, assay = "macs2_peaks_501bp", layer = "meta.features")
peaks.matched <- MatchRegionStats( meta.feature = meta.feature[open.peaks, ],
                                   query.feature = na.omit(meta.feature[top.da.peak$gene, ]),
                                   n = nrow(meta.feature))

for (i in 1:length(unique(top.da.peak$cell_type))) {
  celltype <- unique(top.da.peak$cell_type)[i]
  safe_celltype <- gsub("[[:space:]/]", "_", celltype)
  DAP_sub <- subset(top.da.peak, cell_type == celltype)
  DAP_sub_up <- subset(DAP_sub, sig == "up")
  DAP_sub_down <- subset(DAP_sub, sig == "down")
  fwrite(DAP_sub, file = paste0(path_edger,"DAR/", safe_celltype, "_DAR.csv"), row.names = FALSE, quote = FALSE, sep = "\t")
  fwrite(DAP_sub_up, file = paste0(path_edger,"DAR/", safe_celltype, "_up", "_DAR.csv"), row.names = FALSE, quote = FALSE, sep = "\t")
  fwrite(DAP_sub_down, file = paste0(path_edger,"DAR/", safe_celltype,"_down", "_DAR.csv"), row.names = FALSE, quote = FALSE, sep = "\t")
  for (dap in list(DAP_sub_up,DAP_sub_down)){
    if (length(dap$gene) > 0) {
      enriched.motifs <- FindMotifs(object = atac, background = peaks.matched, features = dap$gene)
      enriched.motifs <- enriched.motifs[order(enriched.motifs$fold.enrichment, decreasing = TRUE), ]
      fwrite(enriched.motifs, file = file.path(path_edger,"motif/", paste0(safe_celltype, "_", unique(dap$sig), "_motif.csv")), sep = "\t", quote = F)
      top.enriched.motifs <- enriched.motifs[which(enriched.motifs$pvalue < 0.05), ]
      fwrite(top.enriched.motifs, file = file.path(path_edger,"motif/", paste0(safe_celltype, "_", unique(dap$sig), "_p0.05_motif.csv")), sep = "\t", quote = F)
      MotifPlot(object = atac, motifs = head(rownames(top.enriched.motifs), n = 10), ncol = 2)
      ggsave(filename = file.path(path_edger,"motif/", paste0(safe_celltype, "_", unique(dap$sig), "_motif.pdf")), height = 8)
      closest_genes <- ClosestFeature(atac, regions = dap$gene)
      fwrite(closest_genes, file = file.path(path_edger,"closest_gene/", paste0(safe_celltype, "_", unique(dap$sig), "_closest_genes.txt")), sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
      closest_genes_2kb <- subset(closest_genes, distance <= 2000)
      fwrite(closest_genes_2kb, file = file.path(path_edger,"closest_gene/", paste0(safe_celltype, "_", unique(dap$sig), "_closest_genes_2kb.txt")), sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
      only.gene.names <- na.omit(closest_genes$gene_name) %>% str_subset(.,pattern="^.+$") %>% sort()
      fwrite(as.data.frame(only.gene.names), file = paste0(path_edger,"closest_gene/", safe_celltype, "_", unique(dap$sig), "_closest_genename.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
      only.gene.names_2kb <- na.omit(closest_genes_2kb$gene_name) %>% str_subset(.,pattern="^.+$") %>% sort()
      fwrite(as.data.frame(only.gene.names_2kb), file = paste0(path_edger,"closest_gene/", safe_celltype,  "_", unique(dap$sig), "_closest_genename_2kb.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
      human <- read.table("./database/TF_for_SCENIC/human.pig.ortholog.109.one2one.txt", header = TRUE, sep = "\t", fill = TRUE, comment.char = "")
      pig <- fread(paste0(path_edger,"closest_gene/", safe_celltype, "_", unique(dap$sig), "_closest_genename.txt"))
      colnames(pig) = "Gene.name"
      pig2human = merge(pig, human) 
      pig2human_genename <- pig2human$Human.gene.name %>% str_subset(.,pattern="^.+$") %>% sort()
      fwrite(as.data.frame(pig2human_genename), file = paste0(path_edger, "metascape/", safe_celltype, "_", unique(dap$sig), "_METASCAPE.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
      pig <- fread(paste0(path_edger,"closest_gene/", safe_celltype, "_", unique(dap$sig), "_closest_genename_2kb.txt"))
      colnames(pig) = "Gene.name"
      pig2human = merge(pig, human)
      pig2human_genename_2k  <- pig2human$Human.gene.name %>% str_subset(.,pattern="^.+$") %>% sort()
      fwrite(as.data.frame(pig2human_genename_2k), file = paste0(path_edger, "metascape/", safe_celltype, "_", unique(dap$sig), "_METASCAPE_2kb.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
    } else {
      message("No regions in ", celltype, " DAP ", unique(dap$sig), " process")}
  }
}

txdb <- makeTxDbFromGFF(paste0(path_ref,"Sscrofa11.1_109.gtf"))
DAP_gr <- DE_sig$gene
DAP_gr_rmchr <- gsub('chr', '', DAP_gr)
DAP_gr_rmchr <- StringToGRanges(DAP_gr_rmchr, sep = c("-","-"))
peakAnno <- annotatePeak_c57(DAP_gr_rmchr, tssRegion=c(-2200, 500), TxDb=txdb, annoDb="org.Ss.eg.db")
pdf("diff_breed_peakAnno_Pie.pdf", width = 6, height = 6)
plotAnnoPie(peakAnno)
dev.off()

peakAnno_info <- as.data.frame(peakAnno)
write.table(peakAnno_info, file="DAP_chipseeker_info_all.txt", sep="\t", row.names=FALSE, quote=FALSE)

CellTypes <- unique(edgeR_DAR_sig$cell_type)
peakAnnoList <- lapply(CellTypes, function(CellType) {
  cell_type_peaks <- edgeR_DAR_sig %>% filter(CellType == cell_type)
  DAP_gr <- cell_type_peaks$gene
  DAP_gr_rmchr <- gsub('chr', '', DAP_gr)
  DAP_gr_rmchr <- StringToGRanges(DAP_gr_rmchr, sep = c("-", "-"))
  peakAnno <- annotatePeak_c57(DAP_gr_rmchr, tssRegion = c(-2200, 500), TxDb = txdb, annoDb="org.Ss.eg.db")
  peakAnno_info <- as.data.frame(peakAnno)
  safe_celltype <- gsub("[^a-zA-Z0-9]", "_", CellType)
  write.table(peakAnno_info, file=paste0("DAP_chipseeker_info_", safe_celltype, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
  return(peakAnno)
})
names(peakAnnoList) <- CellTypes
color <- paletteer::paletteer_d("RColorBrewer::Paired")
colors <- setNames(c(color[1:length(peakAnno@annoStat$Feature)]), peakAnno@annoStat$Feature)

pdf('diff_breed_peakAnno_Bar_cluster.pdf',width = 3*3)
plotAnnoBar(peakAnnoList)+ scale_fill_manual(values = colors) + guides(fill = guide_legend(reverse = TRUE))
dev.off()

#### Footprint ####
atac <- read_rds(paste0(path_mac,"atac_motif_chromvar.rds"))
DefaultAssay(atac) <- "macs2_peaks_501bp"
Idents(atac) = "celltype_predict"

for (i in 1:length(unique(atac$celltype_predict))) {
  celltypes <- unique(atac$celltype_predict)[i]
  safe_celltype <- gsub("[[:space:]/]", "_", celltypes) 
  feature <-  fread(paste0(path_dff,"motif/",safe_celltype,"_p0.05_motif.csv"))
  feature_motif <- head(feature$motif.name, n = 2)
  color_mapping <- setNames(rep("gray", length(unique(atac$celltype_predict))), unique(atac$celltype_predict))
  color_mapping[as.character(celltypes)] <- "#2F318D"
  for (motif in feature_motif){
    atac_fp = atac
    atac_fp <- Footprint(object = atac_fp, motif.name = motif, in.peaks = TRUE, genome = BSgenome.Sscrofa.UCSC.susScr11)
    p <- PlotFootprint(atac_fp, features = motif, label = FALSE, show.expected = FALSE, label.idents = celltypes) + 
      scale_fill_manual(values =  color_mapping) +
      scale_color_manual(values =  color_mapping) +
      theme(legend.position = "none")
    ggsave(paste0(path_foot,safe_celltype,"_", motif, "_footprint.pdf"), p, width = 4, height=4)
  }
}