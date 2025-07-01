library(dplyr)

# Set the paths
input_dir <- "/path/to/alignment_files" # DIR_ALIGN
output_dir <- "/path/to/orthogonal_datasets_for_sqanti3" # DIR_ORTHOGONAL
output_file <- paste0(output_dir,"/GUSTO_filtered_SJ.out.tab")

# List all SJ.out.tab files
sj_files <- list.files(path = input_dir, pattern = "SJ.out.tab$", full.names = TRUE)

# Function to read and assign column names
read_sj_file <- function(file) {
  read.table(file, header = FALSE, sep = "\t",
             col.names = c("chr", "start", "end", "strand", "motif", "annotation", 
                           "uniq_reads", "multi_reads", "max_overhang"))
}

# Read all files and bind them into one dataframe
all_sj <- do.call(rbind, lapply(sj_files, read_sj_file))

# Merge by columns 1–6 and summarize
merged_sj <- all_sj %>%
  group_by(chr, start, end, strand, motif, annotation) %>%
  summarise(
    uniq_reads = sum(uniq_reads),
    multi_reads = sum(multi_reads),
    max_overhang = max(max_overhang),
    .groups = "drop"
  )

merged_sj_out <- merged_sj[merged_sj$chr%in%c(paste0("chr",seq.int(1,22,1)),"chrX","chrY","chrM"),]

filtered_sj <- merged_sj_out %>%
  mutate(total_support = uniq_reads + multi_reads) %>%
  filter(total_support >= 2, max_overhang >= 10)
filtered_sj$total_support <- NULL

# Write to output file
write.table(filtered_sj, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)