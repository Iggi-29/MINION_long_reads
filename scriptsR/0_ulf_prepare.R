########################################
### Prepare the data for ulf scripts ###
########################################

#### libraries
library(regioneR)

#### snakemake parameters
output_file_name = snakemake@output[["done_step"]]
flair_corrected_bed = snakemake@params[["flair_corrected_bed"]]
exons_folder = snakemake@parms[["exons_folder"]]

#### Do the work
### list all exons df in exons_folder

### Corrected bed file into GRanged object
auxDF <- read.table(flair_corrected_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
correctedRanges = regioneR::toGRanges(auxDF[,c(1:3)])

#### write the output
cmd <- paste0("touch ",output_file_name)
system(cmd)