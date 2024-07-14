# Count of gene copies of functions of interest
setwd("Manuscript_1/Figure 2/")
library(dplyr)
library(readxl)
library(tidyverse)
library("writexl")

a01_vec_genomeID<-c(
  "g__Alistipes FM_02",
  "g__Alistipes FM_01",
  "g__Alistipes FM_04",
  "g__Alistipes FM_05",
  "g__Alistipes FM_03",
  "g__Alistipes FM_13",
  "g__Alistipes FM_12",
  "g__Alistipes FM_11",
  "g__Alistipes FM_10",
  "g__Alistipes FM_09",
  "g__Alistipes FM_08",
  "g__Alistipes FM_07",
  "g__Alistipes FM_06",
  "g__Alistipes FM_20",
  "g__Alistipes FM_18",
  "g__Alistipes FM_19",
  "g__Alistipes FM_16",
  "g__Alistipes FM_15",
  "g__Alistipes FM_17",
  "g__Alistipes FM_14",
  "f__Rikenellaceae FM_21",
  "f__Rikenellaceae FM_29",
  "f__Rikenellaceae FM_28",
  "f__Rikenellaceae FM_27",
  "f__Rikenellaceae FM_26",
  "f__Rikenellaceae FM_25",
  "g__Mucinivorans FM_24",
  "g__Mucinivorans FM_23",
  "g__Mucinivorans FM_22",
  "g__Aphodosoma FM_31",
  "g__HGM05232 FM_30",
  "g__Akkermansia FM_36",
  "g__CAKUIA01 FM_35",
  "o__Christensenellales FM_59",
  "o__Christensenellales FM_60",
  "o__Christensenellales FM_58",
  "o__Christensenellales FM_57",
  "g__WQYD01 FM_56",
  "f__Lachnospiraceae FM_68",
  "f__Lachnospiraceae FM_67",
  "f__Lachnospiraceae FM_66",
  "g__RGIG6307 FM_65",
  "f__Lachnospiraceae FM_64",
  "f__Lachnospiraceae FM_63",
  "g__Anaerotignum FM_62",
  "f__CAG-274 FM_61",
  "g__JAGHEK01 FM_55",
  "g__Faecalibacterium FM_53",
  "g__UBA6857 FM_54",
  "g__SIG603 FM_50",
  "g__SIG603 FM_49",
  "f__Butyricicoccaceae FM_52",
  "f__Butyricicoccaceae FM_51",
  "g__Lawsonibacter FM_48",
  "g__Lawsonibacter FM_47",
  "f__Oscillospiraceae FM_46",
  "g__DUWA01 FM_45",
  "f__UBA660 FM_44",
  "g__RQZE01 FM_43",
  "f__Anaeroplasmataceae FM_42",
  "f__Anaeroplasmataceae FM_41",
  "g__CALURL01 FM_40",
  "f__Treponemataceae FM_38",
  "g__Nitratidesulfovibrio FM_37",
  "o__RF32 FM_39",
  "g__Ferrimonas FM_34",
  "s__Vibrio pomeroyi FM_33",
  "f__Psittacicellaceae FM_32")

a02_tax_color<-c("Bacteroidia"="#B2DF8A",
             "Bacilli"="#FDBF6F",
             "Clostridia"="#FF7F00",
             "Spirochaetia"="#FB9A99",
             "Desulfovibrionia"="#FFFF99",
             "Verrucomicrobiae"="#B15928",
             "Alphaproteobacteria"="#CAB2D6",
             "Gammaproteobacteria"="#6A3D9A",
             "Vampirovibrionia"="#1F78B4")

# Aggregated annotation table
a03_annotation<-read.csv("../Annotation_table.csv", header = TRUE)%>%
  select(MAG,query_name,description_uniprot,description_uniref,description_kegg,sseqid_pfam,description_pfam,sseqid_tigrfam,description_tigrfam)

##############################################################
##############################################################
##############################################################
##############################################################
# CAZy families and CAZy families containing EC similarity

# CAZy families and EC functions of interest
b01_subst_meta<-read_excel("Polysaccharide_meta_interest.xlsx")%>%
  select(Fam_CAZy,EC_cazyEC,name,substrate_broad,substrate_simple,substrate_order)%>%
  rename(metabolism=substrate_broad)

# matching CAZy families containing unique activities
b02_family_match<-b01_subst_meta%>%
  filter(!grepl("3.|4.",EC_cazyEC))%>%
  select(-c(EC_cazyEC))

b03_gene_count_fam<-read.csv("../Annotation_table.csv", header = TRUE)%>%
  select(genomeID,query_name,Fam_sub_CAZy,EC_cazyEC)%>%
  mutate(Fam_sub_CAZy = strsplit(Fam_sub_CAZy, "\\|"))%>%
  unnest(Fam_sub_CAZy)%>%
  filter(grepl("GH|PL",Fam_sub_CAZy))%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_sub_CAZy))%>%
  inner_join(b02_family_match,by = "Fam_CAZy",relationship="many-to-many")
colnames(b03_gene_count_fam)
# matching CAZy families containing EC
b04_family_EC_match<-b01_subst_meta%>%
  filter(grepl("3.|4.",EC_cazyEC))

b05_gene_count_fam_EC<-read.csv("../Annotation_table.csv", header = TRUE)%>%
  select(genomeID,query_name,Fam_sub_CAZy,EC_cazyEC)%>%
  mutate(Fam_sub_CAZy = strsplit(Fam_sub_CAZy, "\\|"))%>%
  unnest(Fam_sub_CAZy)%>%
  filter(grepl("GH|PL",Fam_sub_CAZy))%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_sub_CAZy))%>%
  inner_join(b04_family_EC_match,by = c("Fam_CAZy","EC_cazyEC"),relationship="many-to-many")
colnames(b05_gene_count_fam_EC)
# bind fam based and EC base
b06_matches<-rbind(b03_gene_count_fam,b05_gene_count_fam_EC)

# counts
b07_count_fam_ec<-read.csv("../Annotation_table.csv", header = TRUE)%>%
  select(genomeID)%>%
  distinct()%>%
  left_join(b06_matches,by = "genomeID")%>%
  select(genomeID,query_name,Fam_CAZy,EC_cazyEC,name,metabolism,substrate_simple,substrate_order)%>%
  group_by(genomeID,Fam_CAZy,EC_cazyEC)%>%
  mutate(count=n())%>%
  ungroup()%>%
  select(genomeID,Fam_CAZy,name,EC_cazyEC,count,metabolism,substrate_order)%>%
  distinct()%>%
  arrange(match(genomeID, a01_vec_genomeID))%>%
  pivot_wider(id_cols = c("metabolism","substrate_order","Fam_CAZy","EC_cazyEC","name"), names_from = genomeID, values_from=count, values_fill = 0)%>%
  arrange(substrate_order)%>%
  filter(grepl("alginate|laminarin|FCSP|carrageenan|agarose|galactan|starch",metabolism))%>%
  pivot_longer(cols = -c("metabolism","substrate_order","Fam_CAZy","EC_cazyEC","name"), names_to = "genomeID", values_to = "count")%>%
  group_by(genomeID,name)%>%
  mutate(enzyme_count=sum(count))%>%
  select(genomeID,metabolism,substrate_order,name,enzyme_count)%>%
  distinct()
write_xlsx(b07_count_fam_ec,"01a_enzyme_interest.xlsx")
