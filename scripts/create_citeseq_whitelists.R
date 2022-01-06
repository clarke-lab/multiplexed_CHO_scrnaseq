#!/usr/bin/env Rscript --vanilla
#### Description: imports kallisto bustools counts for both
####              multiplexed libraries and outputs whitelist           
####              for cells for cite-seq counts
####
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

package_list <- c("Seurat", "Matrix","DropletUtils")

invisible(lapply(package_list, require, character.only = TRUE, quietly=TRUE))

root_dir <-c("./")

read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cho.combined.cds in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

# first library
cho_tagged_lib1_mat <- read_count_output("kallisto_counts/multiplexed_scRNAseq_data_1/counts_unfiltered/", 
                                          name = "cells_x_genes")

tot_counts <- Matrix::colSums(cho_tagged_lib1_mat)
summary(tot_counts)
bc_rank <- barcodeRanks(cho_tagged_lib1_mat, lower = 10)

options(repr.plot.width=9, repr.plot.height=6)

cho_tagged_lib1_mat <- cho_tagged_lib1_mat[, tot_counts > metadata(bc_rank)$inflection]

print("SBO library 1")
print(dim(cho_tagged_lib1_mat))

write.table(colnames(cho_tagged_lib1_mat), 
            file="CITE-seq_counts/SBO_library_1/whitelist.tsv", 
            append = FALSE, 
            sep = "\t",
            row.names = F, 
            col.names = F, 
            quote = F)

# second library
cho_tagged_lib2_mat <- read_count_output("kallisto_counts/multiplexed_scRNAseq_data_2/counts_unfiltered/", 
                                          name = "cells_x_genes")

tot_counts <- Matrix::colSums(cho_tagged_lib2_mat)
summary(tot_counts)
bc_rank <- barcodeRanks(cho_tagged_lib1_mat, lower = 10)

options(repr.plot.width=9, repr.plot.height=6)

cho_tagged_lib2_mat <- cho_tagged_lib2_mat[, tot_counts > metadata(bc_rank)$inflection]

print("SBO library 2")
print(dim(cho_tagged_lib2_mat))

write.table(colnames(cho_tagged_lib2_mat), 
            file="CITE-seq_counts/SBO_library_2/whitelist.tsv", 
            append = FALSE, 
            sep = "\t",
            row.names = F, 
            col.names = F, 
            quote = F)