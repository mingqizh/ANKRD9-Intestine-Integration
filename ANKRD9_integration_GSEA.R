library(SeuratDisk)
library(Seurat)
library(WGCNA)

# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
# This .d5seurat object can then be read in manually
seuratObject <- LoadH5Seurat(".h5seurat", assays = "RNA")
gene1 = 'ANKRD9'
pbmc = seuratObject
pdf(file = 'all cells UMAP GTEx Small Intestine.pdf')
p2 <- DimPlot(pbmc, reduction = "umap",  group.by = 'free_annotation',label = TRUE,
              repel = TRUE) +  NoLegend()
print(p2)
dev.off()
pdf(file = paste0(gene1, ' expression UMAP GTEx Small Intestine.pdf'))
FeaturePlot(pbmc, features = c(gene1),  reduction = "umap",) +  NoLegend()
dev.off()
Idents(pbmc) = pbmc@meta.data$free_annotation
table(Idents(pbmc))
pdf(file = paste0(gene1, ' expression Violin GTEx Small Intestine.pdf'))
VlnPlot(object = pbmc, features = c(gene1), group.by = 'free_annotation' ) +  NoLegend()
dev.off()
run_paths_incell = function(cell_subset){
  pbmc_small.c1 <- subset(pbmc, idents = cell_subset)
  cnts_mat = as.data.frame(pbmc_small.c1[["RNA"]]@counts)
  colnames(cnts_mat) = pbmc_small.c1@meta.data$free_annotation
  cnts_mat[1:5,1:5]
  cnts_mat1 = cnts_mat[!rowSums(cnts_mat == 0) >= (length(colnames(cnts_mat))*.8), , drop = FALSE]
  colnames(cnts_mat)[1:10]
  gene_subset = gene1
  cc1 = as.data.frame(t(cnts_mat1))
  cc2 = cc1[,colnames(cc1)==gene_subset,]
  cc3 = cc1[,!colnames(cc1)==gene_subset,]
  bics_table = bicorAndPvalue(cc2, cc3, use = 'p')
  cor_table = reshape2::melt(bics_table$bicor)
  cor_table$Var1=NULL
  colnames(cor_table) = c('gene_symbol', 'bicor')
  cor_table$pvalue = reshape2::melt(bics_table$p)$value
  head(cor_table)
  cor_table = cor_table[order(cor_table$pvalue, decreasing = F),]
  deg_set1 = cor_table[cor_table$pvalue<0.2,]
  library(clusterProfiler)
  library(enrichR)
  library(cowplot)
  library(enrichplot)
  library(ggplot2)
  #BiocManager::install('org.Hs.eg.db', force = T)
  library(org.Hs.eg.db)
  fc_dko = deg_set1$bicor
  names(fc_dko) <- deg_set1$gene_symbol
  fc_dko <- sort(fc_dko, decreasing = TRUE)
  fc_dko = fc_dko[!duplicated(names(fc_dko))]
  organism = "org.Hs.eg.db"
  #Note. exponent, minGSsize and maxGSsize tailored to each cell type given the variability in cell expression
  gse <-gseGO(
    geneList=fc_dko,
    ont = "ALL",
    OrgDb= organism,
    keyType = "SYMBOL",
    exponent = 1,
    minGSSize = 2,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH")
  gse_results = as.data.frame(gse@result)
  write.csv(gse_results, file = paste0('full pathway enrichment file for ', gene_subset, ' in ', cell_subset, '.csv'), row.names = F)
  ## gsego https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/gseGO
  pdf(file =paste0('GSEA output ', gene_subset, ' in ', cell_subset, '.pdf'))
  g3 = dotplot(gse, showCategory=15, split=".sign", color="pvalue") + facet_grid(.~.sign) + ggtitle(paste0( 'GSEA terms')) + theme(axis.text.y=element_text(size=4))
  print(g3)
  dev.off()
  pdf(file = paste0('GSEA network graph ', gene_subset, ' in ', cell_subset, '.pdf'))
  x2<- pairwise_termsim(gse)
  g4 = emapplot(x2, showCategory = 15, color = "pvalue", cex_label_category=1)+ ggtitle(paste0('GSEA Network graph'))
  print(g4)
  dev.off()
  edox2 <- pairwise_termsim(gse)
  p1 <- treeplot(edox2)
  p2 <- treeplot(edox2, hclust_method = "average")
  gg42 = aplot::plot_list(p1, p2, tag_levels='A')
  pdf(file = paste0('pathway tree graph ', gene_subset, ' in ', cell_subset, '.pdf'))
  print(gg42)
  dev.off()
  pdf(file = paste0('Ridgeplot ',gene_subset,' in ', cell_subset,  '.pdf'))
  g1 = ridgeplot(gse, showCategory = 20) + labs(x = "logFC") + ggtitle('Distribution of fold changes')+ theme(axis.text.y=element_text(size=5))
  print(g1)
  dev.off()
  #
  gse@result$Description[1:20]
  ##cnetplot gene concept networks
  pdf(file = paste0('cnet gene concept network ', gene_subset,' in ', cell_subset,  '.pdf'), width = 11, height = 8)
  g22 = cnetplot(gse, foldChange=fc_dko,
                 showCategory = gse@result$Description[1:3],
                 circular = TRUE, colorEdge = TRUE) +
    ggtitle('Gene-Concept Network from GSEA-GO')
  print(g22)
  dev.off()
}
table(Idents(pbmc))
run_paths_incell('immature enterocyte')
run_paths_incell('mature enterocyte')
########get pc paths. Change cell_subset for cell of interest
cell_subset = 'mature enterocyte'
pbmc_small.c1 <- subset(pbmc, idents = cell_subset)
cnts_mat = as.data.frame(pbmc_small.c1[["RNA"]]@counts)
colnames(cnts_mat) = pbmc_small.c1@meta.data$free_annotation
cnts_mat[1:5,1:5]
#adjust criteria based on cell type to remove missing data
cnts_mat1 = cnts_mat[!rowSums(cnts_mat == 0) >= (length(colnames(cnts_mat))*.8), , drop = FALSE]
colnames(cnts_mat)[1:10]
gene_subset = gene1
cc1 = as.data.frame(t(cnts_mat1))
cc2 = cc1[,colnames(cc1)==gene_subset,]
cc3 = cc1[,!colnames(cc1)==gene_subset,]
bics_table = bicorAndPvalue(cc2, cc3, use = 'p')
cor_table = reshape2::melt(bics_table$bicor)
cor_table$Var1=NULL
colnames(cor_table) = c('gene_symbol', 'bicor')
cor_table$pvalue = reshape2::melt(bics_table$p)$value
head(cor_table)
cor_table = cor_table[order(cor_table$pvalue, decreasing = F),]
cc1 = cor_table[cor_table$pvalue<1e-50,]
cc1 = na.omit(cc1)
write.csv(cc1, file = 'gene set for ankrd9.csv', row.names = F)