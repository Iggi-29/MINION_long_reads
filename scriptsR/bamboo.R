######################################
### Bambu long reads data analysis ###
######################################

#### libraries
library(bambu)
library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

#### pick the raw data
# bam_file <- c("/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS4_LRS4_barcode01.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS4_LRS4_barcode02.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS4_LRS4_barcode03.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS4_LRS4_barcode04.flair.aligned.bam")
# bam_file <- c("/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS5_LRS5_barcode05.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS5_LRS5_barcode06.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS5_LRS5_barcode07.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS5_LRS5_barcode08.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS5_LRS5_barcode09.flair.aligned.bam")
bam_file <- c("/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS6_LRS6_barcode12.flair.aligned.bam",
              "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS6_LRS6_barcode13.flair.aligned.bam",
              "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS6_LRS6_barcode14.flair.aligned.bam",
              "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS6_LRS6_barcode15.flair.aligned.bam")
# bam_file <- c("/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS7_LRS7_barcode16.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS7_LRS7_barcode17.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS7_LRS7_barcode18.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS7_LRS7_barcode19.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS7_LRS7_barcode20.flair.aligned.bam")
# bam_file <- c("/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS8_LRS8_barcode21.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS8_LRS8_barcode22.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS8_LRS8_barcode23.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS8_LRS8_barcode24.flair.aligned.bam",
#               "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/alignments/flair/LRS8_LRS8_barcode01.flair.aligned.bam")


fa_file <- "/imppc/labs/eclab/Resources/MINION_ddbb/Ref_files/GRCh38.primary_assembly.genome.fa"
gtf_file_name <- "/imppc/labs/eclab/Resources/MINION_ddbb/Ref_files/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf"
# gtf_file_name <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/ref_files/filtered.gff3"
genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this_ENIGMA.txt"

#### Biomart data importation
genes_to_work_on <- readLines(con = genes_to_work_on)
### connect to the service
mart <- biomaRt::useEnsembl(biomart = "ensembl", version = 115)
datasets <- biomaRt::listDatasets(mart = mart)
ensembl <- biomaRt::useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
### get the filters
filts <- biomaRt::listFilters(mart = ensembl)
attrs <- biomaRt::listAttributes(mart = ensembl)
### get the annotations
annotation_of_genes <- biomaRt::getBM(mart = ensembl,
                                      attributes = c("ensembl_gene_id","ensembl_gene_id_version",
                                                     "hgnc_symbol","ensembl_transcript_id_version"),
                                      filters = c("hgnc_symbol"), 
                                      values = genes_to_work_on)
annotation_of_genes_select <- biomaRt::getBM(mart = ensembl,
                                             attributes = c("hgnc_symbol","ensembl_transcript_id_version","transcript_mane_select"), 
                                             filters = c("hgnc_symbol"), values = genes_to_work_on)
annotation_of_genes_select <- annotation_of_genes_select %>% 
  dplyr::filter(transcript_mane_select != "")

annotation_of_genes <- annotation_of_genes %>% 
  dplyr::mutate(is_MANE_select = ifelse(ensembl_transcript_id_version %in% annotation_of_genes_select$ensembl_transcript_id_version,
                                        "yes","no"))


#### pick the gtf and filter it 
gtf_file_raw <- rtracklayer::import(con = gtf_file_name)
gtf_file_raw <- gtf_file_raw[gtf_file_raw$type == "exon"]
gtf_file_raw <- gtf_file_raw[
  mcols(gtf_file_raw)$gene_name %in%
    annotation_of_genes$hgnc_symbol
]

annotation_of_genes <- annotation_of_genes %>% 
  dplyr::filter(ensembl_transcript_id_version %in% gtf_file_raw$transcript_id)

#### Generate the Bambu data
tictoc::tic()
gtf_file <- bambu::prepareAnnotations(x = gtf_file_name)
tictoc::toc()

tictoc::tic()
se <- bambu::bambu(reads = bam_file, annotations = gtf_file, genome = fa_file, NDR = 1)
tictoc::toc()

#### Check the Bambu data
# check the quantification levels
assayNames(se)
# check the metadata of our exeriment
head(colData(se))
# The annotation of each found transcript 
head(mcols(rowRanges(se)))
# Check the number of transcripts
table(mcols(rowRanges(se))$GENEID) |> summary()

#### Filter junk isoforms
## Extract assays
counts <- assay(se, "counts")
flc <- assay(se, "fullLengthCounts")

## Filtering rule
keep_tx <- (rowSums(flc > 0) >= 1) &      # full-length support in ≥1 samples
  (rowSums(flc) >= 3) &          # ≥3 full-length reads in total
  (rowSums(counts) >= 25)        # sufficient overall expression

## Apply filter
se_robust <- se
# se_robust <- se[keep_tx, ]

#### Compare isoforms
rr <- rowRanges(se_robust)

#### Genes with canonical isoforms
annotated_isoforms <- as.data.frame(rr@elementMetadata)

#### Filter the isoforms data
annotated_isoforms <- annotated_isoforms %>% 
  ## filter for non-novel genes
  dplyr::filter(novelGene == FALSE) %>% 
  ## add check if isoforms are in the genes of interest and filter
  dplyr::mutate(normal_gene_id = gsub(pattern = "\\..*$", replacement = "", x = GENEID)) %>% 
  dplyr::mutate(Gene_in_interest = normal_gene_id %in% annotation_of_genes$ensembl_gene_id) %>% 
  dplyr::relocate(Gene_in_interest, normal_gene_id, .before = GENEID) %>% 
  dplyr::rename(GENEID_ensembl = normal_gene_id) %>% 
  dplyr::filter(Gene_in_interest == TRUE) %>% 
  ## add the annotations of the genes 
  dplyr::mutate(Gene_symbol = annotation_of_genes$hgnc_symbol[match(GENEID_ensembl,
                                                                    annotation_of_genes$ensembl_gene_id)]) %>% 
  dplyr::relocate(Gene_symbol)

final_data_list <- list()
gtf_by_tx <- split(gtf_file_raw, mcols(gtf_file_raw)$transcript_id)
#### For each gene, do the check
for (i in 1:length(unique(annotated_isoforms$Gene_symbol))) {
  ### data_now
  gene_now <- unique(annotated_isoforms$Gene_symbol)[i]
  data_now <- annotated_isoforms %>% 
    dplyr::filter(Gene_symbol == gene_now)
  if (all(data_now$TXNAME %in% gtf_file_raw$transcript_id)) {
    cat(paste0("All isoforms of ",gene_now," are canonical No comparison will be performed","\n"))
  } else if (!all(data_now$TXNAME %in% gtf_file_raw$transcript_id)) {
    cat(paste0("Not all isoforms of ",gene_now," are canonical, comparsion will be performed","\n"))
    # pick reference GRAnges to compare
    transcripts_to_compare_with <- annotation_of_genes %>%
      dplyr::filter(hgnc_symbol == gene_now) %>%
      dplyr::pull(ensembl_transcript_id_version)
    
    ranges_of_subject <- rr[mcols(rr)$TXNAME %in% transcripts_to_compare_with]
    ranges_of_query <- rr[mcols(rr)$TXNAME %in% data_now$TXNAME[!data_now$TXNAME %in% transcripts_to_compare_with]]
    
    for (j in 1:length(names(ranges_of_query))) {
      ## data now
      query_now <- names(ranges_of_query)[j]
      ranges_of_query_now <- ranges_of_query[mcols(ranges_of_query)$TXNAME %in% query_now]
      for (k in 1:length(names(ranges_of_subject))) {
        ## data now
        subject_now <- names(ranges_of_subject)[k]
        ranges_of_subject_now <- ranges_of_subject[mcols(ranges_of_subject)$TXNAME %in% subject_now]
        ## final data
        name_final <- paste0(query_now,"_",subject_now)
        comparison_now <- bambu::compareTranscripts(query = ranges_of_query_now,
                                                    subject = ranges_of_subject_now)
        final_data_list[[name_final]] <- comparison_now
      }
    }
  }
  }

final_data <- do.call(rbind, final_data_list)
final_data <- merge(x = final_data, y = annotation_of_genes, 
                    by.x = "subjectId", by.y = "ensembl_transcript_id_version", 
                    all.x = TRUE)
final_data <- final_data %>% 
  dplyr::relocate(ensembl_gene_id_version, hgnc_symbol) %>% 
  dplyr::select(-c(ensembl_gene_id)) %>% 
  dplyr::relocate(is_MANE_select, .after = subjectId)

##### get the expression
final_data_ids <- unique(x = c(final_data$queryId, final_data$subjectId))
filtered_counts <- counts[grep(pattern = paste0(x = final_data_ids, collapse = "|"), x = rownames(counts)),,
                          drop = FALSE]

##### get the data to visualize
se_export <- se_robust[mcols(rowRanges(se_robust))$TXNAME %in% final_data_ids,]
rr_export <- rowRanges(se_export)

se_export_canonical <- se_export[
  mcols(rowRanges(se_export))$TXNAME %in% c(
    final_data_ids[final_data_ids %in% annotation_of_genes_select$ensembl_transcript_id_version],
    final_data_ids[!final_data_ids %in% annotation_of_genes$ensembl_transcript_id_version])]
rr_export_canonical <- rowRanges(se_export_canonical)

bambu::writeToGTF(annotation = rr_export,
  file = "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/bambu/gtf/gtf_for_lrs6.gtf"
)

bambu::writeToGTF(annotation = rr_export_canonical,
                  file = "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/bambu/gtf_canonical/gtf_for_lrs6_canonical.gtf"
)

openxlsx::write.xlsx("/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/bambu/lrs6_isoform_data.xlsx", x = final_data)
