## make cleaned up dataset for Anya

## from mouse_cortex_analysis
annot <- read.csv('data/S1R1_with_structure_id_name.csv.gz', row.names=2)
annot <- annot[,-1]
dim(annot)
head(annot)

celltype <- read.csv('cell_type_clustering_11-18-2021.csv.gz', row.names = 1)
## python writing out names got a little messed up
rownames(celltype) <- sapply(rownames(celltype), function(x) as.numeric(strsplit(x, '-')[[1]][1]))
## double check that names match up
table(rownames(celltype) %in% rownames(pos))
head(celltype)
ct <- celltype$leiden
names(ct) <- rownames(celltype)
ct <- as.factor(ct)
table(ct)

## plot
par(las=2, mar=rep(1,4))
pos <- annot[, c('center_x', 'center_y')]
MERINGUE::plotEmbedding(pos, groups=ct, cex=0.1, shuffle.colors = FALSE)

struc <- annot$acronym
names(struc) <- rownames(annot)
MERINGUE::plotEmbedding(pos, groups=struc, cex=0.1, shuffle.colors = FALSE)

## get gene expression
gexp <- read.csv('data/datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S1R1.csv.gz', row.names=1)
gexp[1:5,1:5]
dim(gexp)

## make one simple table
df <- data.frame(annot, celltype[rownames(annot),], gexp[rownames(annot),])

write.csv(df, file='~/Desktop/S2R2_with_structure_id_name_and_cell_type_clustering_11-18-2021.csv')

## double check we can get info back out
df <- read.csv('~/Desktop/S2R2_with_structure_id_name_and_cell_type_clustering_11-18-2021.csv.gz', row.names=1)
head(df)
pos <- df[, c('center_x', 'center_y')]
vi <- df$acronym == 'CA3' & df$leiden %in% c('Oligodendrocyte Progenitor Cells', 'Oligodendrocyte Progenitor Cells(1)')
table(vi)
plot(pos, pch=".")
points(pos[vi,], col='red', pch=16)

