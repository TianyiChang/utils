library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# import CDS annot and MGE coord tables
import_combine_dts <- function(relative_path, pattern) {

    path2files <- list.files(
        path = str_c(args[1], {{relative_path}}),
        pattern = {{pattern}},
        full.names = TRUE,
        recursive = TRUE)

    # read all files into a single data.table
    combined_dt <- rbindlist(lapply(path2files, fread))

}

cds_annot_dt <- import_combine_dts(
    "/candid_mge_aa_chunks",
    "split_batch_\\d+.txt")

# Check if FSlink directory exists and has files
fslink_path <- file.path(args[1], "FSlink")
if (dir.exists(fslink_path) && length(list.files(fslink_path, pattern = "_fslink.tsv$")) > 0) {
    fslinks <- import_combine_dts(
        "/FSlink",
        "_fslink.tsv")
} else {
    # Create an empty data.table with the required column
    fslinks <- data.table(cds_id = character(), filename = character())
}

#=======================#
# Process CDS annot table
#=======================#

# sort the data.table by gene, score (descending), and evalue
cds_annot_dt <- cds_annot_dt[order(gene, -score, evalue)]

# select the top row for each gene
cds_annot_dt <- cds_annot_dt[, .SD[1], by = gene]

cds_annot_dt <- cds_annot_dt %>% 
    select(cds_id = gene, pfam_id)

# append pfam metadata and file names
pfam_metad <- read_tsv(args[2])

cds_annot_dt <- cds_annot_dt %>% 
    mutate(pfam_id = str_replace(pfam_id, '\\.\\d+$', '')) %>% 
    left_join(pfam_metad, by = 'pfam_id') %>% 
    left_join(fslinks, by = 'cds_id') %>% 
    distinct()

write_tsv(cds_annot_dt, args[3])
