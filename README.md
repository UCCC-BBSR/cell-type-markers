# cell-type-markers

 Gene lists and gene expression matrix references for annotating cell types. Their primary use is for annotating cell types in single-cell RNA-seq data.

## How to use this repository

This repository contains gene lists and references for annotating cell types in single-cell RNA-seq data. The gene lists are stored in CSV files and can be used with the `clustifyr` package to annotate cell types in Seurat objects. The gene expression matrix references are stored in RDS files and can be used directly with `clustifyr`.

**Please see the example in the `example-usage-clustifyr` folder.**

### How to contribute

The structure of this repo is not solidified. Let's keep it simple and start a discussion about how to organize and optimize the resources. If you have suggestions or contributions, please feel free to open an issue or pull request. The most basic way to contribute is to add your gene lists or gene expression matrix references in the appropriate folders. Ultimately we should come up with a standardized format for the gene lists and references.

If you have other sources of gene lists or references that you think would be useful, please share them! We can create a folder for curated resources or links to external databases.

### Gene Lists

`markers-gene-lists`

Curated gene lists for various cell types are stored here in CSV format. Suggested columns for the CSV files include:

| Gene | Cell Type | Species | Description | Source | Pubmed ID / DOI | Notes |
|------|-----------|---------|-------------|--------|------------------|-------|
| CD3E | T Cell | Human | CD3 epsilon subunit of T-cell receptor complex | PanglaoDB | 12345678 | Marker for T cells |
| CD19 | B Cell | Human | B-lymphocyte antigen CD19 | CellMarker | 23456789 | Marker for B cells |

### Gene Expression Matrix References

`markers-expression-objects`

RDS files can be stored here if they are small (< 100 MB). Larger files will need a different storage solution like Google Drive or Dropbox.
These could be in the form of Seurat objects, SingleCell Experiment objects, or any other format that can be easily read into R.

### Scripts

`helper-scripts`

**WORK IN PROGRESS**
These are helper scripts for generating, curating, formatting, and using the gene lists and gene expression matrix references.

Example scripts could include:

- `generate_gene_list.R`: A script to generate a gene list from a publication or database.
- `format_gene_list.R`: A script to format gene lists into the desired CSV structure.
- `annotate_cell_types.R`: A script to annotate cell types using the `clustifyr` package and the gene lists or expression matrices.
- `visualize_annotations.R`: A script to visualize the results of cell type annotations on a Seurat object.
- `download_refs.R`: A script to download and prepare reference datasets for use with `clustifyr`.
- `check_gene_lists.R`: A script to check the integrity and format of gene lists before use.
- `merge_gene_lists.R`: A script to merge multiple gene lists into a single comprehensive list.
- `update_gene_lists.R`: A script to update existing gene lists with new information or genes.

---

## Annotating Cell Types with clustifyr

**Please see the example in the `example-usage-clustifyr` folder.**

The `clustifyr` package can be used to annotate cell types using reference gene lists or a gene expression matrix from a Seurat object. Below are two common approaches:

### Using a gene list stored in a CSV file

To annotate cell types using a pre-defined gene list stored in a CSV file:

```r
library(clustifyr)
library(readr)

# Load gene list
gene_markers <- read_csv("path/to/gene_list.csv")

# Available metrics include: "hyper", "jaccard", "spearman", "gsea"
list_res <- clustify_lists(
    input = so, # matrix of normalized single-cell RNA-seq counts
    cluster_col = "RNA_snn_res.0.5", # name of column in meta.data containing cell clusters
    marker = cell_type_markers, # list of known marker genes
    marker_inmatrix = FALSE,
    metric = "jaccard", # test to use for assigning cell types
    obj_out = FALSE # return Seurat object
)

plot_cor_heatmap(
    cor_mat = list_res, # matrix of correlation coefficients from clustify_lists()
    cluster_rows = TRUE, # cluster by row
    cluster_columns = TRUE, # cluster by column
    legend_title = "jaccard" # title of heatmap legend
)

so_res <- clustify_lists(
    input = so, # matrix of normalized single-cell RNA-seq counts
    cluster_col = "RNA_snn_res.0.5", # name of column in meta.data containing cell clusters
    marker = cell_type_markers, # list of known marker genes
    marker_inmatrix = FALSE,
    metric = "jaccard", # test to use for assigning cell types
    obj_out = TRUE # return Seurat object
)

# clustifyr stores the cell type assignments in the metadata column "type"
so_res <- SetIdent(so_res, value = "type")
```

### Using a Seurat Object Gene Expression Matrix

To automatically annotate Seurat clusters based on gene expression and annotated cell types from a reference Seurat object.

```r
library(Seurat)
library(clustifyr)

# Load the Seurat object
so_ref <- readRDS("so_ref.rds")

ref_mat <- seurat_ref(
  seurat_object = so_ref,
  cluster_col = "cell_types"
)

so_res2 <- clustify(
  input = so,                      # unannotated Seurat object
  ref_mat = ref_mat,               # matrix of average expression per reference cell type
  cluster_col = "seurat_clusters", # or whatever your cluster column is
  obj_out = FALSE                   # return Seurat object with metadata updated
)

plot_cor_heatmap(
  cor_mat = so_res2,
  cluster_rows = TRUE,
  cluster_columns = TRUE)

so_res2 <- clustify(
  input = so,                      # unannotated Seurat object
  ref_mat = ref_mat,               # matrix of average expression per reference cell type
  cluster_col = "seurat_clusters", # or whatever your cluster column is
  obj_out = TRUE                   # return Seurat object with metadata updated
)
```

---

## Support and Resources

https://panglaodb.se/

http://biocc.hrbmu.edu.cn/CellMarker/
http://xteam.xbio.top/CellMarker/
http://xteam.xbio.top/ACT/

https://www.immgen.org/Databrowser19/DatabrowserPage.html

http://cloud.capitalbiotech.com/SingleCellBase/

https://ngdc.cncb.ac.cn/databasecommons/database/id/6110

https://functionome.geneontology.org/
https://www.nature.com/articles/s41586-025-08592-0

https://github.com/rnabioco/clustifyr
https://rnabioco.github.io/clustifyrdata/articles/download_refs.html
