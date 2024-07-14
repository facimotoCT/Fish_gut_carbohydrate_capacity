# Matrix of EC similar CAZymes, mannitol genes and sulfatases of interest

setwd("Manuscript_1/Figure 2/")

library(dplyr)
library(readxl)
library(tidyverse)
vec_genomeID<-c(
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

tax_color<-c("Bacteroidia"="#B2DF8A",
             "Bacilli"="#FDBF6F",
             "Clostridia"="#FF7F00",
             "Spirochaetia"="#FB9A99",
             "Desulfovibrionia"="#FFFF99",
             "Verrucomicrobiae"="#B15928",
             "Alphaproteobacteria"="#CAB2D6",
             "Gammaproteobacteria"="#6A3D9A",
             "Vampirovibrionia"="#1F78B4")


####################################
####################################
####################################
####################################
# ECs of interest
b01_CAZyme_copies<-read_xlsx("01a_enzyme_interest.xlsx")

# count based on substrate
b02_meta_enzyme<-read.csv("../Annotation_table.csv",header=TRUE)%>%
  select(Class,genomeID,MAG)%>%
  distinct()%>%
  left_join(b01_CAZyme_copies,by = "genomeID")%>%
  mutate(Class_plot = gsub(".*__", "", Class))%>%
  arrange(substrate_order)
b03_vec_enzyme<-unique(b02_meta_enzyme$name)

####################################
####################################
####################################
####################################
# SULFATASES of interest
c01_sulf_count<-read.csv("../Annotation_table.csv",header=TRUE)%>%
  select(Class,genomeID,MAG,query_name,Fam_SULFATLAS)%>%
  na.omit()%>%
  group_by(MAG,Fam_SULFATLAS)%>%
  mutate(count_sulf=n())%>%
  ungroup()%>%
  select(Class,genomeID,MAG,Fam_SULFATLAS,count_sulf)%>%
  distinct()

c02_sulf_no_meta_count<-read.csv("../Annotation_table.csv",header=TRUE)%>%
  select(MAG)%>%
  distinct()%>%
  left_join(c01_sulf_count,by = "MAG")%>%
  ungroup()%>%
  mutate(
    metabolism = case_when( 
      (str_detect(Fam_SULFATLAS, str_c("S")) )~ 'sulfur')) %>%
  select(MAG,Fam_SULFATLAS,metabolism,count_sulf)%>%
  distinct()%>%
  pivot_wider(id_cols = c("metabolism","Fam_SULFATLAS"), names_from = MAG, values_from=count_sulf, values_fill = 0)%>%
  na.omit()%>%
  pivot_longer(cols = -c("metabolism","Fam_SULFATLAS"), names_to = "MAG", values_to = "count_sulf")%>%
  arrange(desc(count_sulf))

# count based on substrate
c03_meta_sulf<-read.csv("../Annotation_table.csv",header=TRUE)%>%
  select(Class,genomeID,MAG)%>%
  distinct()%>%
  left_join(c02_sulf_no_meta_count,by = "MAG")%>%
  mutate(Class_plot = gsub(".*__", "", Class))%>%
  arrange(desc(count_sulf))%>%
  filter(grepl("S1_7\\b|S1_8\\b|S1_15\\b|S1_16\\b|S1_17\\b|S1_19\\b|S1_25\\b",Fam_SULFATLAS))
c05_vec_sul_fam<-unique(c03_meta_sulf$Fam_SULFATLAS)

####################################
####################################
####################################
####################################
# Mannitol genes in operons
d01_man_count<-read.csv("02d_mannitol_genes.csv",header=TRUE)%>%
  select(MAG,query_name,mannitol_genes)%>%
  na.omit()%>%
  group_by(MAG,mannitol_genes)%>%
  mutate(count_man=n())%>%
  ungroup()%>%
  select(MAG,mannitol_genes,count_man)%>%
  distinct()

d02_man_no_meta_count<-read.csv("../Annotation_table.csv",header=TRUE)%>%
  select(MAG)%>%
  distinct()%>%
  left_join(d01_man_count,by = "MAG")%>%
  ungroup()%>%
  mutate(
    metabolism = case_when( 
      (str_detect(mannitol_genes, str_c("fructo|mann")) )~ 'mannitol')) %>%
  select(MAG,mannitol_genes,metabolism,count_man)%>%
  distinct()%>%
  pivot_wider(id_cols = c("metabolism","mannitol_genes"), names_from = MAG, values_from=count_man, values_fill = 0)%>%
  na.omit()%>%
  pivot_longer(cols = -c("metabolism","mannitol_genes"), names_to = "MAG", values_to = "count_man")%>%
  arrange(desc(count_man))

# count based on substrate
d03_meta_man<-read.csv("../Annotation_table.csv",header=TRUE)%>%
  select(Class,genomeID,MAG)%>%
  distinct()%>%
  left_join(d02_man_no_meta_count,by = "MAG")%>%
  mutate(Class_plot = gsub(".*__", "", Class))%>%
  arrange(desc(count_man))

d04_vec_man_fam<-c("mannitol PTS system","mannitol-1-phosphate 5-dehydrogenase","mannitol 2-dehydrogenase","fructokinase")  

####################################
####################################
####################################
####################################
# plot all together
e01_all_genes_screening_ec<-b02_meta_enzyme%>%
  select(Class_plot,genomeID,name,enzyme_count,metabolism)%>%
  rename(Class=Class_plot,gene=name,count=enzyme_count)
e01_all_genes_screening_sulf<-c03_meta_sulf%>%
  select(Class_plot,genomeID,Fam_SULFATLAS,count_sulf,metabolism)%>%
  rename(Class=Class_plot,gene=Fam_SULFATLAS,count=count_sulf)
e01_all_genes_screening_man<-d03_meta_man%>%
  select(Class_plot,genomeID,mannitol_genes,count_man,metabolism)%>%
  rename(Class=Class_plot,gene=mannitol_genes,count=count_man)

e01_all<-rbind(e01_all_genes_screening_ec,e01_all_genes_screening_sulf,e01_all_genes_screening_man)

Summary_table<-e01_all%>%
  select(genomeID,metabolism,gene,count)%>%
  arrange(match(genomeID, vec_genomeID))%>%
  pivot_wider(id_cols = c("metabolism","gene"), names_from = genomeID, values_from=count, values_fill = 0)
writexl::write_xlsx(Summary_table,"03a_genes_interest_EC_man_sulf.xlsx")

e02_vec_all_genes <- c(b03_vec_enzyme,c05_vec_sul_fam,d04_vec_man_fam)
e03_vec_metabolism<-c("alginate","laminarin","FCSP","carrageenan","agarose","porphyran","galactan","starch","sulfur","mannitol")
# Define the desired order of the y-axis labels
e06_all_genes <- ggplot(e01_all, aes(x = factor(genomeID, levels = as.character(vec_genomeID)),
                                     y = factor(gene, levels = as.character(rev(e02_vec_all_genes))),
                                     fill = ifelse(count > 0, Class, "black"))) + 
    facet_grid(vars(factor(metabolism, levels = e03_vec_metabolism)),
               vars(factor(Class,levels=c("Bacteroidia","Verrucomicrobiae","Clostridia",
                               "Bacilli","Vampirovibrionia","Spirochaetia",
                               "Desulfovibrionia","Alphaproteobacteria","Gammaproteobacteria"))),
               scales = "free",
               space = "free",
               labeller = label_wrap_gen(multi_line = TRUE)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(tax_color, "black"), breaks = c("Bacteroidia","Verrucomicrobiae","Clostridia",
                                                             "Bacilli","Vampirovibrionia","Spirochaetia",
                                                             "Desulfovibrionia","Alphaproteobacteria","Gammaproteobacteria")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 7, face = "bold", hjust = 0, vjust = 0.5, margin = margin(1, 1, 1, 1, "mm")),
    legend.position = "bottom",
    legend.title = element_text(size = 6, face = "bold"), #change legend title font size
    legend.text = element_text(size = 5), #change legend text font size
    legend.justification = "centre",
    legend.title.align = 1,
    legend.direction = "horizontal",
    legend.key.size = unit(2, "mm"), # Adjust the size of the legend keys
    legend.key.width = unit(2, "mm"), # Adjust the width of the legend keys
    legend.key.height = unit(2, "mm"), # Adjust the height of the legend keys
    
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = 6, angle = 0, face="bold"),

    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    axis.ticks.y = element_line(size = 0.2),
    axis.ticks.x = element_line(size = 0.2),

    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0,face="bold", size = 6, 
                                     hjust = 0.5, vjust = 0.5, 
                                     margin = margin(1, 1, 1, 1, "mm")),
    strip.placement = "inside",
    plot.margin = margin(1,1,1,1, unit = "mm"),
    plot.background = element_blank(),
    panel.spacing = unit(1, "mm"),
    panel.grid = element_line(size = 0.1),
    panel.border = element_rect(size = 0.2)
    ) +
  labs(
     fill = expression(paste("Class"))  # Example of Greek character (beta) in legend title
     ) 
e06_all_genes  
ggsave(filename = "03b_screening_EC_sulf_man.pdf", e06_all_genes,device = cairo_pdf,  width = 200, height = 160 , units = "mm")

