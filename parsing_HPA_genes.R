###########################################################################################
############################ HPA Genes ######################################
###########################################################################################
library(readxl)
library(dplyr)
library(data.table)

mydir = "~/Box/schavan/Bill/"
myfiles = list.files(mydir, full.names = TRUE)
myfiles

ligep_paper = read_excel("/Users/schavan/Box/schavan/Bill//LiGEP_paper.xlsx", col_names = FALSE); 
dim(ligep_paper)
colnames(ligep_paper)
ligep_paper = dplyr::rename(ligep_paper, 'gene_name' = colnames(ligep_paper)[1])
length(ligep_paper$gene_name) #93

gtex_19k = read_excel("/Users/schavan/Box/schavan/Bill//GTEX-19K-mean.xlsx", col_names = TRUE);
dim(gtex_19k)
length(gtex_19k$gene_name) #19582

known_iedb = read_excel("/Users/schavan/Box/schavan/Bill//Known Autoantigens IEDB truncated 20221213.xlsx")
known_iedb = dplyr::rename(known_iedb, 'gene_name' = 'Gene')
length(known_iedb$gene_name) #19

liver_exclu = read_excel("/Users/schavan/Box/schavan/Bill//Protein atlas liver only genes 20221213.xlsx")
liver_exclu = dplyr::rename(liver_exclu, 'gene_name' = 'Gene_name')
length(liver_exclu$gene_name) #32

mydir = "~/Box/schavan/Bill/HumanProteinAtlasGenes"
myfiles = list.files(mydir, full.names = TRUE)
myfiles

## Specialized epithelial cells broken down by “Cell Type Enriched”, “Group Enriched”, and “Cell Type Enhanced”
## The specialized epithelial cell transcriptome

hpa_sc_cell_type_cholangio = fread("/Users/schavan/Box/schavan/Bill/HumanProteinAtlasGenes/cell_type_category_rna_Cholangiocytes_Cell_208.tsv")
hpa_sc_cell_type_cholangio = dplyr::rename(hpa_sc_cell_type_cholangio, 'gene_name' = 'Gene')
length(hpa_sc_cell_type_cholangio$gene_name) #208

hpa_sc_cell_type_hepato = fread("/Users/schavan/Box/schavan/Bill/HumanProteinAtlasGenes/cell_type_category_rna_Hepatocytes_Cell_843.tsv")
hpa_sc_cell_type_hepato = dplyr::rename(hpa_sc_cell_type_hepato, 'gene_name' = 'Gene')
length(hpa_sc_cell_type_hepato$gene_name) #843

## Tissue Cell types
## Genes with predicted cell type specificity within liver

hpa_en_liver_any = fread("/Users/schavan/Box/schavan/Bill/HumanProteinAtlasGenes/ce_enriched_liver_any_Very_2396.tsv")
hpa_en_liver_any = dplyr::rename(hpa_en_liver_any, 'gene_name' = 'Gene')
length(hpa_en_liver_any$gene_name) #2396

hpa_lcten_cholangio = fread("/Users/schavan/Box/schavan/Bill/HumanProteinAtlasGenes/ce_enriched_liver_Cholangiocyte_Very_57.tsv")              
hpa_lcten_cholangio = dplyr::rename(hpa_lcten_cholangio, 'gene_name' = 'Gene')
length(hpa_lcten_cholangio$gene_name) #57

hpa_lcten_hepato = fread("/Users/schavan/Box/schavan/Bill/HumanProteinAtlasGenes/ce_enriched_liver_Hepatocyte_1_Hepatocyte_2_Very_1264.tsv")
hpa_lcten_hepato = dplyr::rename(hpa_lcten_hepato, 'gene_name' = 'Gene')
length(hpa_lcten_hepato$gene_name) #1264

hpa_lcten_hepatic_stellate = fread("/Users/schavan/Box/schavan/Bill/HumanProteinAtlasGenes/ce_enriched_liver_Hepatic_316.tsv")
hpa_lcten_hepatic_stellate = dplyr::rename(hpa_lcten_hepatic_stellate, 'gene_name' = 'Gene')
length(hpa_lcten_hepatic_stellate$gene_name) #316

hpa_lcten_endothelial = fread("/Users/schavan/Box/schavan/Bill/HumanProteinAtlasGenes/ce_enriched_liver_Sinusoid_330.tsv")
hpa_lcten_endothelial = dplyr::rename(hpa_lcten_endothelial, 'gene_name' = 'Gene')
length(hpa_lcten_endothelial$gene_name) #330

library("ggvenn")
library("ggVennDiagram")

pdf("~/Box/schavan/Bill/HumanProteinAtlasGenes/HPA_liver_genes_overlaps.pdf")
A <-list(SE_Cholangio = hpa_sc_cell_type_cholangio$gene_name,
         SE_Hepato = hpa_sc_cell_type_hepato$gene_name)
ggvenn(A, set_name_size = 4, text_size = 3)
ggVennDiagram(A, set_name_size = 4, text_size = 3)

#Nothing in common
A <-list(
         WL_Cholangio = hpa_lcten_cholangio$gene_name,
         WL_Hepato = hpa_lcten_hepato$gene_name,
         WL_Stellate = hpa_lcten_hepatic_stellate$gene_name,
         WL_Endothelial = hpa_lcten_endothelial$gene_name)
ggvenn(A, set_name_size = 4, text_size = 3)
ggVennDiagram(A, set_name_size = 4, text_size = 3)

#############

WL_subtypes <-c(
 hpa_lcten_cholangio$gene_name,
 hpa_lcten_hepato$gene_name,
 hpa_lcten_hepatic_stellate$gene_name,
 hpa_lcten_endothelial$gene_name)

setdiff(hpa_en_liver_any$gene_name, WL_subtypes)
A <-list(
  WL_Any = hpa_en_liver_any$gene_name,
  WL_4subtypes = WL_subtypes)
ggvenn(A, set_name_size = 4, text_size = 3)
ggVennDiagram(A, set_name_size = 4, text_size = 3)

#############

A <-list(
  WL_Any = hpa_en_liver_any$gene_name,
  WL_Cholangio = hpa_lcten_cholangio$gene_name,
  WL_Hepato = hpa_lcten_hepato$gene_name)
ggvenn(A, set_name_size = 4, text_size = 3)
ggVennDiagram(A, set_name_size = 4, text_size = 3)

A <-list(
  WL_Any = hpa_en_liver_any$gene_name,
  WL_Stellate = hpa_lcten_hepatic_stellate$gene_name,
  WL_Endothelial = hpa_lcten_endothelial$gene_name)
ggvenn(A, set_name_size = 4, text_size = 3)
ggVennDiagram(A, set_name_size = 4, text_size = 3)

#############

A <-list(WL_Any = hpa_en_liver_any$gene_name,
         SE_Cholangio = hpa_sc_cell_type_cholangio$gene_name,
         SE_Hepato = hpa_sc_cell_type_hepato$gene_name)
ggvenn(A, set_name_size = 4, text_size = 3)
ggVennDiagram(A, set_name_size = 4, text_size = 3)

A <-list(
         WL_4subtypes = WL_subtypes,
         SE_Cholangio = hpa_sc_cell_type_cholangio$gene_name,
         SE_Hepato = hpa_sc_cell_type_hepato$gene_name)
ggvenn(A, set_name_size = 4, text_size = 3)
ggVennDiagram(A, set_name_size = 4, text_size = 3)

A <-list(WL_Cholangio = hpa_lcten_cholangio$gene_name,
         WL_Hepato = hpa_lcten_hepato$gene_name,
         SE_Cholangio = hpa_sc_cell_type_cholangio$gene_name,
         SE_Hepato = hpa_sc_cell_type_hepato$gene_name)
ggvenn(A, set_name_size = 4, text_size = 3)
ggVennDiagram(A, set_name_size = 4, text_size = 3)

dev.off()

#############

A <-list(ligep_paper = ligep_paper$gene_name, 
         gtex_19k = gtex_19k$gene_name, 
         known_iedb = toupper(known_iedb$gene_name))
ggvenn(A, set_name_size = 4, text_size = 3)

setdiff(known_iedb$gene_name, gtex_19k$gene_name)
#[1] "POLG_HCV77" "H2AZ1"      "UL98"       
setdiff(ligep_paper$gene_name, gtex_19k$gene_name)
#[1] "C5orf27"
intersect("C5ORF27",gtex_19k$gene_name)
A <-list(liver_exclu = liver_exclu$gene_name,
         ligep_paper = ligep_paper$gene_name, 
         gtex_19k = gtex_19k$gene_name, 
         known_iedb = toupper(known_iedb$gene_name))
ggvenn(A, set_name_size = 4, text_size = 3)
setdiff(liver_exclu$gene_name, gtex_19k$gene_name) #1 "RTL4"

#############

A <-list(WL_Any = hpa_en_liver_any$gene_name,
         liver_exclu = liver_exclu$gene_name,
         gtex_19k = gtex_19k$gene_name, 
         known_iedb = toupper(known_iedb$gene_name))
ggvenn(A, set_name_size = 4, text_size = 3)
setdiff(hpa_en_liver_any$gene_name, gtex_19k$gene_name) #199

A <-list(WL_4subtypes = WL_subtypes,
         liver_exclu = liver_exclu$gene_name, 
         gtex_19k = gtex_19k$gene_name, 
         known_iedb = toupper(known_iedb$gene_name))
ggvenn(A, set_name_size = 4, text_size = 3)
setdiff(WL_subtypes, gtex_19k$gene_name) #50

#############

A <-list(SE_Any = c(SE_Cholangio = hpa_sc_cell_type_cholangio$gene_name,
                    SE_Hepato = hpa_sc_cell_type_hepato$gene_name),
         liver_exclu = liver_exclu$gene_name,
         gtex_19k = gtex_19k$gene_name, 
         known_iedb = toupper(known_iedb$gene_name))
ggvenn(A, set_name_size = 4, text_size = 3)
setdiff(SE_Any$gene_name, gtex_19k$gene_name)

#############         
         
