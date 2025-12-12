#######################
### GFF3 filtration ###
#######################

library(readr)
library(dplyr)

### Load original GFF3
raw_gff3 <- snakemake@input[["gff3_filtered"]]
final_gff3 <- snakemake@output[["gff3_filtered_MANE"]]

gff <- readr::read_tsv(raw_gff, skip = 7, col_names = FALSE)
gff_header <- readLines(con = raw_gff, n = 7)

colnames(gff) <- c("chr","source_of_ann","type","start","end","dot1","strand","dunno","tag")

### filter the GFF3 file
gff_filtered <- gff %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(keep = ifelse(type == "gene","YES, gene",
                              ifelse(grepl(x = tag, pattern = "Mane|MANE"),"YES, MANE","NO!"))) %>% 
  dplyr::ungroup()

gff_filtered <- gff_filtered %>% 
  dplyr::filter(grepl(x = keep, pattern = "YES")) %>% 
  dplyr::select(-c(keep))

writeLines(gff_header, con = final_gff3)
readr::write_tsv(gff_filtered, final_gff3, append = TRUE, col_names = FALSE)

readr::write_tsv(file = final_gff3,
                 x = gff_filtered, col_names = FALSE)
