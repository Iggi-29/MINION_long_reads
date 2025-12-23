#######################
### GFF3 filtration ###
#######################

library(readr)
library(dplyr)

### Load original GFF3
raw_gff3 <- snakemake@input[["gff3_filtered"]]
final_gff3 <- snakemake@output[["gff3_filtered_MANE"]]

gff <- readr::read_tsv(raw_gff3, skip = 7, col_names = FALSE)
gff_header <- readLines(con = raw_gff3, n = 7)

colnames(gff) <- c("chr","source_of_ann","type","start","end","dot1","strand","dunno","tag")
valid_chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                "chr21","chr22","chrX","chrY")
### filter the GFF3 file
gff_filtered <- gff %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(keep = ifelse(type == "gene","YES, gene",
                              ifelse(grepl(x = tag, pattern = "Mane|MANE"),"YES, MANE","NO!"))) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(chr %in% valid_chrs)

gff_filtered <- gff_filtered %>% 
  dplyr::filter(grepl(x = keep, pattern = "YES")) %>% 
  dplyr::select(-c(keep))

writeLines(gff_header, con = final_gff3)
readr::write_tsv(gff_filtered, final_gff3, append = TRUE, col_names = FALSE)

readr::write_tsv(file = final_gff3,
                 x = gff_filtered, col_names = FALSE)
