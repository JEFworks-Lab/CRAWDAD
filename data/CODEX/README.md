Notes about data:

See `metadata.xlsx` for additional information about the CODEX datasets

For spleen, only used datasets with 28 marker panels (so 6 total, 2 paired datasets from 3 different donors)

For annotations, used Maigan's (University of Florida group) for `HBM389.PKHL.939`. Then transferred the annotations to the other datasets using Linear Discriminant Analysis (LDA).

Harmonized paired datasets before transfer: HBM556.KSFB.592 and HBM568.NGPL.345, and HBM825.PBVN.284 and HBM342.FSLD.938

LDA transfers were done between the following reference and query datasets:

HBM389.PKHL.936 -\> HBM772.XXCD.697

HBM389.PKHL.936 -\> HBM556.KSFB.592

HBM389.PKHL.936 -\> HBM825.PBVN.284

HBM556.KSFB.592 -\> HBM568.NGPL.345

HBM825.PBVN.284 -\> HBM342.FSLD.938

The column names in each `*.meta.csv.gz` file indicate these label transfers.

As it happens, when comparing the protein expressions of each annotateed cluster after label transfer, there was high correlation between the same annotations of the reference and query dataset, which means that cell types in the query dataset assigned a given cell type annotation had a similar protein expression profile as the cells with that annotation in the reference dataset. Harmonization actually didn't have a significant effect in this particular case, so probably wasn't entirely necessary but doesn't hurt.

Maigan's original cell type annotations based on normalized cluster protein expression, and her confidence of each label assignment, can be found in `ProposedAnnotations.csv`

The protein expression values in each `*.exp.csv.gz` file are the raw protein expression values from the Akoya output for each dataset. I generated the original clusters by normalizing the protein expression by `cellarea` and applying a `log10 + 1` transform. Then I assigned cells to communities using `igraph::cluster_louvain` with `k=50`. Protein expression values were summed together for cells in a given community, and then the values were scaled. The scaled expression values were visually inspected on a heatmap to assign annotations. The 28 protein markers were chosen to select for specific cell types and for many of the clusters, just a few identifiable markers were enriched providing some confidence in the cell type label assignment.

Similar procedure was done with the thymus samples (also 6 samples, 2 paired samples from 3 donors)

HBM288.XSQZ.639 was the reference dataset annotated by Maigan for the thymus samples
