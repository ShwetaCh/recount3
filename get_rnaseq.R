## References
#https://www.bioconductor.org/packages/devel/bioc/vignettes/recount3/inst/doc/recount3-quickstart.html#41_Available_data
#https://www.raynamharris.com/blog/recount3/
#https://rdrr.io/bioc/recount3/f/vignettes/recount3-quickstart.Rmd
#https://rdrr.io/bioc/recount/man/getTPM.html

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("recount3")

## Check that you have a valid Bioconductor installation
BiocManager::valid()

library(recount3)
## Find all available human projects
human_projects <- available_projects()
human_projects

proj_info <- subset(human_projects, project == "LIVER" & project_type == "data_sources")
proj_info

rse_gene_liver <- create_rse(proj_info)

## Once you have your RSE object, you can transform the raw coverage
## base-pair coverage counts using transform_counts().
## For RPKM, TPM or read outputs, check the details in transform_counts().

## Scale the counts using the AUC
assay(rse_gene_liver, "counts") <- transform_counts(rse_gene_liver)

## If needed, install recount, the R/Bioconductor package for recount2:
# BiocManager::install("recount")
library(recount)
## Compute TPMs
assays(rse_gene_liver)$TPM <- recount::getTPM(rse_gene_liver, length_var = "bp_length")
colSums(assay(rse_gene_liver, "TPM")) / 1e6 ## Should all be equal to 1
assayNames(rse_gene_liver)

## Information about how this RSE object was made
metadata(rse_gene_liver)
assays(rse_gene_liver)$counts
assays(rse_gene_liver)$TPM

## Number of genes by number of samples
dim(rse_gene_liver)
#> [1] 63856    251

## Information about the genes/Annotaion
gene_annot = rowRanges(rse_gene_liver) %>% as.data.frame()
## samples and quality metrics
colData(rse_gene_liver)

## Make it a data frame and remove miRNA files & any files that don’t have an accession number linking the source files
colData0 <- colData(rse_gene_liver) %>% as.data.frame()
colData0 <- filter(colData0, gtex.run_acc != "NA", gtex.smnabtcht != "RNA isolation_PAXgene Tissue miRNA") 
countData0 <- assays(rse_gene_liver)$counts %>% as.data.frame() %>% dplyr::select(rownames(colData0))
TPMData0 <- assays(rse_gene_liver)$TPM %>% as.data.frame() %>% dplyr::select(rownames(colData0))
## Check if this matches
length(rownames(colData0) == length(colnames(countData0))) == length(colnames(TPMData0))
head(rownames(colData0) == colnames(countData0))
head(rownames(colData0) == colnames(TPMData0))

## Look at data by age and biological condition if any
library(tidyverse)
colData0 %>%
  dplyr::group_by(gtex.sex, gtex.smtsd) %>%
  dplyr::summarise(cohort_size = length(gtex.sex)) %>%
  ggplot(aes(x = gtex.sex,  y = cohort_size, fill = gtex.smtsd)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Sex", y = "Cohort size", fill = "Tissue",
       subtitle = "GTEx data obtained using recount3 ") +
  theme_linedraw(base_size = 15) +
  geom_text(aes(label = cohort_size),
            position = position_dodge(width = .9),
            vjust = -0.25)

## Create a mapping from Ensembl gene name to hgnc_symbol
library(biomaRt)
library(tidyr)
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_info <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),  
                   mart = ensembl) %>% 
  mutate_all(na_if, "") %>%
  drop_na(.) %>%
  unique(.) %>%
  mutate(ensembl_gene_id = paste0(ensembl_gene_id, ".1", sep = ""))
head(gene_info)

## create long data for counts
countData_long <- countData0 %>%
  mutate(ensembl_gene_id = rownames(.)) %>%
  pivot_longer(-ensembl_gene_id, 
               names_to = "external_id", values_to = "counts") %>%
  dplyr::inner_join(gene_info, .,by = "ensembl_gene_id") %>% 
  dplyr::full_join(colData0, ., by = "external_id") %>%
  arrange(desc(counts))
head(countData_long)

countData_long %>%
  filter( hgnc_symbol == "UGT1A9") %>%
  ggplot(aes(x = as.character(gtex.sex), y = counts, fill = gtex.smtsd)) +
  geom_boxplot() +
  #scale_y_log10() +
  labs(y = 'counts', x = "Sex", fill = "Tissue") +
  theme_linedraw(base_size = 15) +
  scale_y_continuous(labels = scales::label_number_si()) 

## create long data for TPM
tmpData_long <- TPMData0 %>%
  mutate(ensembl_gene_id = rownames(.)) %>%
  pivot_longer(-ensembl_gene_id, 
               names_to = "external_id", values_to = "TPM") %>%
  dplyr::inner_join(gene_info, .,by = "ensembl_gene_id") %>% 
  dplyr::full_join(colData0, ., by = "external_id") %>%
  arrange(desc(TPM))
head(countData_long)

tmpData_long %>%
  filter( .,hgnc_symbol == "UGT1A9") %>%
  ggplot(aes(x = as.character(gtex.sex), y = TPM, fill = gtex.smtsd)) +
  geom_boxplot() +
  #scale_y_log10() +
  labs(y = 'TPM', x = "Sex", fill = "Tissue") +
  theme_linedraw(base_size = 15) +
  scale_y_continuous(labels = scales::label_number_si()) 

## Remove mitochondrial genes
tmpData_long_fil = filter(tmpData_long, !hgnc_symbol %like% "MT-", ) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, TPM, gtex.subjid, gtex.sex, gtex.smtsd); 
dim(tmpData_long_fil)
names(tmpData_long_fil)
head(tmpData_long_fil)

tmpData_long_fil = left_join(tmpData_long_fil, gene_annot, by = c('ensembl_gene_id' = 'gene_id')) %>%
  filter(gene_type == "protein_coding")
dim(tmpData_long_fil)
names(tmpData_long_fil)
head(tmpData_long_fil)

tmpData_long_fil.aggr = tmpData_long_fil %>% 
  group_by(hgnc_symbol,gtex.sex) %>% 
  dplyr::summarise(.,avgTPM = mean(TPM)) %>%
  arrange(desc(avgTPM))
head(tmpData_long_fil.aggr)

tmpData_long_fil.aggr = tmpData_long_fil %>% 
  group_by(hgnc_symbol) %>% 
  dplyr::summarise(.,avgTPM = mean(TPM)) %>%
  arrange(desc(avgTPM)) %>%
  filter(avgTPM >=0.1)
head(tmpData_long_fil.aggr)
dim(tmpData_long_fil.aggr)

