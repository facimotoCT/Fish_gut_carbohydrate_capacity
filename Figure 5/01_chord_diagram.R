# CGC expression
setwd("Manuscript_1/Figure 5/")
library("dplyr")
library("tidyr")
library("stringr")
library(patchwork)
library(ggplot2)
library(ggrepel)
library(scales)
library(rstatix)
library(ggnewscale)
library(ggpubr)
library(readxl)
library(tidyverse)
library(gtools)

# Expression (Data from recently generated large dataset from 2017 and 2020)
a01_TPM<-read.delim("../local_files/3F0502_normalised_summary_count_table.tsv")%>%
  select(geneID,ends_with("_III_TPM"),ends_with("_IV_TPM"),ends_with("_V_TPM"))%>%
  rename_with(~ gsub("rep_community_concatenated_vs_", "", .), contains("rep_community_concatenated_vs_"))%>%
  rename(query_name=geneID)%>%
  pivot_longer(-c("query_name"), names_to = "Sample", values_to = "TPM")%>%
  mutate(Sample = gsub("_TPM", "", Sample))%>%
  mutate(Section = gsub(".*_", "", Sample))%>%
  mutate(Fish = gsub("_.*", "", Sample))%>%
  select(query_name,Sample,Fish,Section,TPM)%>%
  mutate(Total_TPM=sum(TPM))%>%
  group_by(Section)%>%
  mutate(Section_TPM=sum(TPM))%>%
  group_by(Sample)%>%
  mutate(Sample_TPM=sum(TPM))%>%
  filter(!grepl("III",Section))
  
# Annotation (Data from recently generated large dataset from 2017 and 2020)
a02_annotation<-read.csv("../local_files/annotation.csv")%>%
  filter(!grepl("III", genome))

# genome relative expression
a03_genome_expression<-a02_annotation%>%
  select(Class,genomeID,query_name)%>%
  inner_join(a01_TPM,by="query_name")%>%
  mutate(Total_TPM=sum(TPM))%>%
  group_by(Section)%>%
  mutate(Section_TPM=sum(TPM))%>%
  group_by(Section)%>%
  mutate(Sample_TPM=sum(TPM))%>%
  group_by(genomeID,Sample)%>%
  mutate(genomeID_Sample_TPM=sum(TPM))%>%
  select(Total_TPM,Section,Section_TPM,Sample,Sample_TPM,Class,genomeID,genomeID_Sample_TPM)%>%
  distinct()
  
a04_class_col<-c("Bacteroidia"="#B2DF8A",
  "Bacilli"="#FDBF6F",
  "Clostridia"="#FF7F00",
  "Spirochaetia"="#FB9A99",
  "Desulfovibrionia"="#FFFF99",
  "Verrucomicrobiae"="#B15928",
  "Alphaproteobacteria"="#CAB2D6",
  "Gammaproteobacteria"="#6A3D9A",
  "Vampirovibrionia"="#1F78B4",
  "Brachyspirae"="#E31A1C",
  "Fusobacteriia"="#33A02C",
  "Dehalobacteriia"="#A6CEE3")

##############################################################################
##############################################################################
##############################################################################
##############################################################################
# CAZyme expression section, class
b01_CAZy_Sulfatase_CGC <- a02_annotation %>%
  select(Class, query_name, Fam_CAZy, CGC_ID, gene_type_DBCAN, Fam_SULFATLAS) %>%
  mutate(CAZyme = if_else(!is.na(Fam_CAZy), "CAZyme", NA_character_)) %>%
  mutate(Sulfatase = if_else(!is.na(Fam_SULFATLAS), "Sulfatase", NA_character_)) %>%
  mutate(CGC = if_else(!is.na(CGC_ID), "CGC", NA_character_)) %>%
  mutate("CGC CAZyme" = if_else(gene_type_DBCAN == "CAZyme", "CGC CAZyme", NA_character_))%>%
  mutate("CGC TP" = if_else(gene_type_DBCAN == "TC", "CGC TP", NA_character_))%>%
  mutate("CGC TF" = if_else(gene_type_DBCAN == "TF", "CGC TF", NA_character_))%>%
  mutate("CGC STP" = if_else(gene_type_DBCAN == "STP", "CGC STP", NA_character_))%>%
  mutate("CGC other genes" = if_else(gene_type_DBCAN == "null", "CGC other genes", NA_character_))%>%
  select(Class, query_name, CAZyme, Sulfatase, CGC, "CGC CAZyme","CGC TP","CGC TF","CGC STP","CGC other genes") %>%
  distinct()%>%
  pivot_longer(cols=-c("Class", "query_name"),names_to = "types",values_to = "Element")%>%
  na.omit()%>%
  select(Class, query_name,Element)%>%
  distinct()
  
b02_element_expression<-a01_TPM%>%
  group_by(Section)%>%
  mutate(Section_TPM=sum(TPM))%>%
  inner_join(b01_CAZy_Sulfatase_CGC,by="query_name",relationship="many-to-many")%>%
  group_by(Section,Class,Element)%>%
  mutate(Section_Class_Element_TPM=sum(TPM))%>%
  mutate(Section_Class_Element_relative_expression=(Section_Class_Element_TPM/Section_TPM)*100)%>%
  select(Class, Element,Section,Section_Class_Element_relative_expression)%>%
  distinct()%>%
  group_by(Class)%>%
  mutate(Class_order=sum(Section_Class_Element_relative_expression))%>%
  arrange(desc(Class_order))%>%
  mutate(Class = gsub("c__", "", Class))%>%
  filter(grepl("CAZyme|Sulfatase|CGC CAZymes",Element ))

b03_Carbo_elements<-ggplot(b02_element_expression, aes(x = factor(Element,
                                 levels = 
                                   c("CAZyme","Sulfatase","CGC",
                                     "CGC CAZyme", "CGC TP",
                                     "CGC TF","CGC STP","CGC other genes")),   
                      y = Section_Class_Element_relative_expression,
                      fill = factor(Class,levels = unique(b02_element_expression$Class)))) +
  facet_grid(cols = vars(factor(Section, levels = c("III","IV","V"))),
             scales = "free",
             space = "free",
             labeller = label_wrap_gen(multi_line = TRUE)) + 
  geom_bar(stat = "identity", 
           show.legend = TRUE,size=0.1,
           color="black") +
  scale_fill_manual(values=a04_class_col) +
  scale_y_continuous(
    name = "Relative expression",
    position = "left",
    limits = c(0,2.5),
    expand = c(0,0),
    breaks=seq(0,4,1)) +
  scale_x_discrete(
    name = "",
    position = "bottom") +
  theme_bw() +
  labs(title = "",
       fill = "Class")+
  theme(
    title = element_blank(),
    legend.position = "right",
    legend.spacing.x = unit(1, "mm"),
    legend.title = element_text(size = 6, face = "bold"), #change legend title font size
    legend.text = element_text(size = 5), #change legend text font size
    legend.justification = "left",
    legend.title.align = 0,
    legend.direction = "vertical",
    legend.key.size = unit(2, "mm"), # Adjust the size of the legend keys
    legend.key.width = unit(2, "mm"), # Adjust the width of the legend keys
    legend.key.height = unit(2, "mm"), # Adjust the height of the legend keys
    
    axis.text.x = element_text(size = 5, angle= 90,hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 5, angle= 0,hjust = 1, vjust = 0.5),
    axis.ticks.y = element_line(size = 0.2),
    axis.ticks.x = element_line(size = 0.2),
    axis.title.x = element_text(size = 7, face="bold",angle= 0,hjust = 0.5, vjust = 0.5),
    axis.title.y = element_text(size = 7, face="bold",angle= 90,hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 7,face="bold", angle= 0,hjust = 0.5, vjust = 0.5),
    plot.margin = margin(1,1,1,1, unit = "mm"),
    plot.background = element_blank(),
    panel.spacing = unit(1, "mm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.2)
  ) 
b03_Carbo_elements
ggsave(filename = "01a_carbohydrate_elements.pdf", b03_Carbo_elements, width = 75, height = 50 , units = "mm")

##############################################################################
##############################################################################
##############################################################################
##############################################################################
# Highly expressed CAZymes
c01_Top_expressing_CAZymes<-a02_annotation%>%
  select(Fam_CAZy,query_name)%>%
  na.omit()%>%
  inner_join(a01_TPM,by="query_name")%>%
  group_by(Sample)%>%
  mutate(CAZy_expression_Sample=sum(TPM))%>%
  mutate(Fam_CAZy = strsplit(Fam_CAZy, "\\|"))%>%
  unnest(Fam_CAZy)%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_CAZy))%>%
  filter(grepl("GH|PL|CE",Fam_CAZy))%>%
  ungroup()%>%
  group_by(Fam_CAZy,Sample)%>%
  mutate(Fam_CAZy_expression_Sample=sum(TPM))%>%
  select(Fam_CAZy,Sample,Fam_CAZy_expression_Sample,CAZy_expression_Sample)%>%
  distinct()%>%
  mutate(Fam_CAZy_impact_Sample=(Fam_CAZy_expression_Sample/CAZy_expression_Sample)*100)%>%
  filter(Fam_CAZy_impact_Sample>0)%>%
  group_by(Sample)%>%
  arrange(desc(CAZy_expression_Sample),desc(Fam_CAZy_impact_Sample))%>%
  mutate(Perc_cumulative = cumsum(Fam_CAZy_impact_Sample))%>%
  filter(Perc_cumulative < 50 )
length(unique(c01_Top_expressing_CAZymes$Fam_CAZy))

# differentially expressed
c02_families_KW<-a02_annotation%>%
  select(Fam_CAZy,query_name)%>%
  na.omit()%>%
  mutate(Fam_CAZy = strsplit(Fam_CAZy, "\\|"))%>%
  unnest(Fam_CAZy)%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_CAZy))%>%
  filter(grepl("GH|PL|CE",Fam_CAZy))%>%
  distinct()%>%
  inner_join(a01_TPM,by="query_name",relationship="many-to-many")%>%
  select(Fam_CAZy,query_name,Section,Sample,TPM)%>%
  distinct()%>%
  group_by(Fam_CAZy,Sample)%>%
  mutate(Fam_CAZy_Sample=sum(TPM))%>%
  select(Sample,Section,Fam_CAZy,Fam_CAZy_Sample)%>%
  distinct()%>%
  ungroup()%>%
  group_by(Fam_CAZy)%>%
  mutate(KW_Fam_Section=kruskal.test(Fam_CAZy_Sample~Section)$stat, p_KW_Fam_Section = kruskal.test(Fam_CAZy_Sample~Section)$p.value)%>%
  select(Sample,Section,Fam_CAZy,Fam_CAZy_Sample,KW_Fam_Section,p_KW_Fam_Section)%>%
  distinct()

# Fam Dunn 005
c03_dunn_CAZy_fam_TPM_005<-c02_families_KW%>%
  filter(p_KW_Fam_Section < 0.05)%>%
  select(Sample,Section,Fam_CAZy,Fam_CAZy_Sample)%>%
  group_by(Fam_CAZy)%>%
  dunn_test(Fam_CAZy_Sample ~ Section)%>%
  filter(p.adj < 0.05)%>%
  arrange(p.adj)
length(unique(c03_dunn_CAZy_fam_TPM_005$Fam_CAZy))

# Fam Dunn 001
c04_dunn_CAZy_fam_TPM_001<-c02_families_KW%>%
  filter(p_KW_Fam_Section < 0.01)%>%
  select(Sample,Section,Fam_CAZy,Fam_CAZy_Sample)%>%
  group_by(Fam_CAZy)%>%
  dunn_test(Fam_CAZy_Sample ~ Section)%>%
  filter(p.adj < 0.01)%>%
  arrange(p.adj)
length(unique(c04_dunn_CAZy_fam_TPM_001$Fam_CAZy))

c05_matrix_top_expressing<-a02_annotation%>%
  select(Fam_CAZy,query_name)%>%
  na.omit()%>%
  inner_join(a01_TPM,by="query_name")%>%
  mutate(Fam_CAZy = strsplit(Fam_CAZy, "\\|"))%>%
  unnest(Fam_CAZy)%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_CAZy))%>%
  group_by(Fam_CAZy,Sample)%>%
  mutate(Fam_expression_Sample=sum(TPM))%>%
  ungroup()%>%
  select(Sample,Fam_CAZy,Fam_expression_Sample)%>%
  distinct()%>%
  pivot_wider(id_cols="Fam_CAZy",names_from = Sample, values_from = Fam_expression_Sample,values_fill = 0)%>%
  #filter(Fam_CAZy %in% unique(c03_dunn_CAZy_fam_TPM_005$Fam_CAZy))%>%
  filter(Fam_CAZy %in% unique(c01_Top_expressing_CAZymes$Fam_CAZy))%>%
  column_to_rownames("Fam_CAZy")
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# CGCs containing these highly expressing degradative CAZy families
d01_CGCs_containing_fams<-a02_annotation%>%
  select(CGC_ID,Fam_CAZy,query_name)%>%
  na.omit()%>%
  mutate(Fam_CAZy = strsplit(Fam_CAZy, "\\|"))%>%
  unnest(Fam_CAZy)%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_CAZy))%>%
  filter(Fam_CAZy %in% unique(c01_Top_expressing_CAZymes$Fam_CAZy))%>%
  select(CGC_ID)%>%
  distinct()

# CGCs expressing at least 2 degradative CAZymes in the same sample
d02_concert_expression<-a02_annotation%>%
  select(CGC_ID,query_name,gene_type_DBCAN,Fam_CAZy)%>%
  filter(grepl("CGC",CGC_ID))%>%
  # Total genes in CGC
  group_by(CGC_ID)%>%
  mutate(CGC_genes=n())%>%
  # Total gene type
  group_by(CGC_ID,gene_type_DBCAN)%>%
  mutate(CGC_type_genes=n())%>%
  # Degradative per CGC
  mutate(CAZyme_activity = if_else(grepl("GH|PL|CE", Fam_CAZy), "degradative", "other"))%>%
  group_by(CGC_ID,CAZyme_activity)%>%
  mutate(CGC_deg=n())%>%
  # Expressed genes per CGC in each sample
  inner_join(a01_TPM,by="query_name")%>%
  mutate(Expression = if_else(TPM > 0, "expressed", "no"))%>%
  group_by(CGC_ID,Sample,Expression)%>%
  mutate(CGC_genes_exp_sample=n())%>%
  # Expressed degradative per CGC in each sample
  group_by(CGC_ID,Sample,CAZyme_activity,Expression)%>%
  mutate(Deg_exp_sample=n())%>%
  select(CGC_ID,CGC_genes,gene_type_DBCAN,CGC_type_genes,CAZyme_activity,CGC_deg,Sample,Expression,CGC_genes_exp_sample,Deg_exp_sample)%>%
  distinct()%>%
  mutate(Percentage_CGC_expressed=(CGC_genes_exp_sample/CGC_genes)*100,
         Percentage_deg_CGC_expressed=(Deg_exp_sample/CGC_deg)*100)%>%
  filter(grepl("CAZyme",gene_type_DBCAN))%>%
  filter(grepl("degradative",CAZyme_activity))%>%
  filter(grepl("expressed",Expression))%>%
  filter(Deg_exp_sample > 1)%>%
  filter(Percentage_deg_CGC_expressed > 99)%>%
  filter(Percentage_CGC_expressed > 99)%>%
  ungroup()%>%
  select(CGC_ID)%>%
  distinct()

# CGCs highly expressed per sample
d03_highly_expressing_CGCs<-a02_annotation%>%
  select(Class,CGC_ID,query_name)%>%
  na.omit()%>%
  inner_join(a01_TPM,by="query_name")%>%
  group_by(CGC_ID,Sample)%>%
  mutate(CGC_Sample_expression=sum(TPM))%>%
  ungroup()%>%
  select(Class,CGC_ID,Sample,CGC_Sample_expression)%>%
  distinct()%>%
  group_by(Sample)%>%
  mutate(Total_CGC_expression=sum(CGC_Sample_expression))%>%
  mutate(CGC_impact_Sample=(CGC_Sample_expression/Total_CGC_expression)*100)%>%
  filter(CGC_impact_Sample>0)%>%
  group_by(Sample)%>%
  arrange(desc(Total_CGC_expression),desc(CGC_impact_Sample))%>%
  mutate(Perc_cumulative = cumsum(CGC_impact_Sample))%>%
  filter(Perc_cumulative < 25 )
length(unique(d03_highly_expressing_CGCs$CGC_ID))

d04_CGC_content<-a02_annotation%>%
  filter(grepl("CGC",CGC_ID))%>%
  select(Class,CGC_ID,query_name,Fam_CAZy,Fam_SULFATLAS,gene_type_DBCAN)%>%
  mutate(CAZyme = if_else(!is.na(Fam_CAZy), Fam_CAZy, NA_character_))%>%
  mutate(Sulfatase = if_else(!is.na(Fam_SULFATLAS), Fam_SULFATLAS, NA_character_)) %>%
  mutate(Transporter = if_else(gene_type_DBCAN == "TC", "TC", NA_character_))%>%
  mutate(Transcription_factor = if_else(gene_type_DBCAN == "TF", "TF", NA_character_))%>%
  mutate(Signal_transduction_protein = if_else(gene_type_DBCAN == "STP", "STP", NA_character_))%>%
  unite(CGC_annotation,CAZyme,Sulfatase,Transporter,Transcription_factor,Signal_transduction_protein, na.rm = TRUE, sep = ' ')%>%
  mutate(CGC_annotation = if_else(nzchar(CGC_annotation), CGC_annotation, "null"))%>%
  group_by(CGC_ID)%>%
  arrange(CGC_ID,gtools::mixedorder(query_name))%>%
  mutate(CGC_content = paste0(CGC_annotation, collapse = ":::"))%>%
  filter(grepl("GH|PL|CE",Fam_CAZy))%>%
  ungroup()%>%
  select(CGC_ID,CGC_content,Fam_CAZy)%>%
  distinct()%>%
  group_by(CGC_ID)%>%
  arrange(CGC_ID,Fam_CAZy)%>%
  mutate(Degradative_CAZyme_content = paste0(Fam_CAZy, collapse = ":::"))%>%
  ungroup()%>%
  select(CGC_ID,CGC_content,Degradative_CAZyme_content)%>%
  distinct()
  
d05_classify<-d04_CGC_content%>%
  filter(CGC_ID %in% unique(d03_highly_expressing_CGCs$CGC_ID))%>%
  filter(CGC_ID %in% unique(d02_concert_expression$CGC_ID))%>%
  #################
  # Alginate
  mutate(
  Substrate = case_when(
    ##### ALGINATE
    (str_detect(CGC_content, str_c("PL6|PL7|PL34|PL17|PL15|PL38")))~ 'alginate',
    
    ##### FCSP
    # pairs
    (str_detect(CGC_content, str_c("GH29\\b")))~ 'FCSP',
    (str_detect(CGC_content, str_c("GH30_4\\b")))~ 'FCSP',
    (str_detect(CGC_content, str_c("GH151\\b")) & str_detect(CGC_content, str_c("GH116\\b")))~ 'FCSP',
    (str_detect(CGC_content, str_c("GH141\\b")) & str_detect(CGC_content, str_c("GH117\\b|GH16_12\\b|GH81\\b")))~ 'FCSP',
    
    ##### LAMINARIN
    # pairs
    (str_detect(CGC_content, str_c("GH16_3")) & str_detect(CGC_content, str_c("GH31\\b|GH3\\b")))~ 'laminarin',
    (str_detect(CGC_content, str_c("GH149")) & str_detect(CGC_content, str_c("GH1\\b|GH3\\b")))~ 'laminarin',
    (str_detect(CGC_content, str_c("GH30_3")))~ 'laminarin',
    (str_detect(CGC_content, str_c("GH161")) & str_detect(CGC_content, str_c("GH94\\b")))~ 'laminarin',
    (str_detect(CGC_content, str_c("GH128")) & str_detect(CGC_content, str_c("GH3\\b")))~ 'laminarin',
    
    ##### CARRAGEENAN
    # pairs
    (str_detect(CGC_content, str_c("GH150\\b")) & str_detect(CGC_content, str_c("GH150|GH43_29|GH2\\b")))~ 'carrageenan',
    (str_detect(CGC_content, str_c("GH167\\b")) & str_detect(CGC_content, str_c("GH1\\b|GH2\\b|GH3\\b")))~ 'carrageenan',
    (str_detect(CGC_content, str_c("GH16_13\\b")) & str_detect(CGC_content, str_c("GH2\\b")))~ 'carrageenan',
    (str_detect(CGC_content, str_c("GH16_3")) & str_detect(CGC_content, str_c("GH2\\b")))~ 'carrageenan',
    
    ##### STARCH
    # pairs
    (str_detect(CGC_content, str_c("GH57\\b")) & str_detect(CGC_content, str_c("GH133\\b")))~ 'starch',
    (str_detect(CGC_content, str_c("GH13\\b|GH13_")) & str_detect(CGC_content, str_c("GH77\\b|GH31\\b|GH97\\b|GH133\\b|GH57\\b|GH2\\b|GH3\\b|GH13_36\\b|GH13_9\\b|GH13_28\\b|GH13_31\\b")))~ 'starch',
    
    ##### AGAROSE/PORPHYRAN
    # pairs
    (str_detect(CGC_content, str_c("GH86\\b|GH16_11\\b|GH16_12\\b|GH16_14\\b|GH16_11\\b|GH16_16\\b")))~ 'agarose',
  ))%>%
  select(CGC_ID,CGC_content,Substrate)%>%
  distinct()%>%
  na.omit()

##############################################################################
##############################################################################
##############################################################################
##############################################################################
# Chord diagram
library("circlize")
e01_Section_CGC<-a02_annotation%>%
  select(Class,Genus,genomeID,query_name,CGC_ID,gene_type_DBCAN,Fam_CAZy,Fam_SULFATLAS)%>%
  filter(grepl("CGC",CGC_ID))%>%
  inner_join(a01_TPM, by ="query_name")%>%
  group_by(Section)%>%
  mutate(CGC_Section_TPM=sum(TPM))%>%
  filter(grepl("GH|PL|CE|AA|GT|CBM",Fam_CAZy))%>%
  group_by(Section)%>%
  mutate(CAZyme_Section_TPM=sum(TPM))%>%
  mutate(Fam_CAZy = strsplit(Fam_CAZy, "\\|"))%>%
  unnest(Fam_CAZy)%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_CAZy))%>%
  filter(grepl("GH|PL|CE",Fam_CAZy))%>%
  filter(CGC_ID %in% unique(d05_classify$CGC_ID))%>%
  left_join(d05_classify,by="CGC_ID")%>%
  group_by(Section,Class,Substrate,Fam_CAZy)%>%
  mutate(Substrate_Class_CGC_Fam_Section=sum(TPM))%>%
  mutate(Substrate_Class_CGC_Fam_Section_relative_expression=(Substrate_Class_CGC_Fam_Section/CGC_Section_TPM)*100
         )%>%
  ungroup()%>%
  select(Substrate,Class,Fam_CAZy,Section,Substrate_Class_CGC_Fam_Section_relative_expression)%>%
  filter(Substrate_Class_CGC_Fam_Section_relative_expression > 0)%>%
  distinct()%>%
  mutate(Class = gsub("c__", "", Class))%>%
  mutate(`CAZy_Substrate` = paste0(as.character(Fam_CAZy), ":::", as.character(Substrate)))%>%
  mutate(`Section_Class` = paste0(as.character(Section), ":::", as.character(Class)))
length(unique(e01_Section_CGC$Fam_CAZy))

length(unique(e01_Section_CGC$Substrate))

# Sectors
e02_class_section_vec<-e01_Section_CGC%>%
  group_by(Section)%>%
  mutate(Section_order=sum(Substrate_Class_CGC_Fam_Section_relative_expression))%>%
  group_by(Class)%>%
  mutate(Class_order=sum(Substrate_Class_CGC_Fam_Section_relative_expression))%>%
  ungroup()%>%
  select(Class,Class_order,Section,Section_order,Section_Class)%>%
  distinct()%>%
  arrange(desc(Section_order),desc(Class_order))

e03_Substrate_Fam_vec<-e01_Section_CGC%>%
  group_by(Substrate)%>%
  mutate(Substrate_order=sum(Substrate_Class_CGC_Fam_Section_relative_expression))%>%
  group_by(Fam_CAZy)%>%
  mutate(Fam_order=sum(Substrate_Class_CGC_Fam_Section_relative_expression))%>%
  arrange(match(Substrate, c("alginate","laminarin","FCSP","carrageenan","agarose","starch")),desc(Fam_order))

e04_vector_sector <- unique(c(e02_class_section_vec$Section_Class, e03_Substrate_Fam_vec$CAZy_Substrate))
e04_vector_sector
e05_label_sector <- structure(gsub(".*:::", "", e04_vector_sector), names = e04_vector_sector)
e05_label_sector

e06_color_sector_meta<-read_xlsx("sectors_colors.xlsx")
e07_link_cols <- setNames(e06_color_sector_meta$color, e06_color_sector_meta$sector)
e07_link_cols

##############################################################################
f02_IV_chord<-e01_Section_CGC%>%
  filter(grepl("IV",Section))
library(circlize)
circos.clear()
cairo_pdf("01b_chord_diagram_IV.pdf")  # Adjust dimensions as needed
chordDiagramFromDataFrame(f02_IV_chord %>% 
                            select(Section_Class,CAZy_Substrate),
                          order = e04_vector_sector,
                          grid.col = structure(e07_link_cols[match(gsub(".*:::", "", e04_vector_sector),
                                                                   labels(e07_link_cols))],  
                                               names = e04_vector_sector),
                          link.lty = e01_Section_CGC$Substrate_Class_CGC_Fam_Section_relative_expression, 
                          target.prop.height = 1,
                    
                          directional = 1,
                          annotationTrack = c("grid", "axis"),
                          preAllocateTracks = 1)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(
    CELL_META$xcenter,
    CELL_META$ylim[1]+cm_h(1.1),
    sub(":::.*", "", CELL_META$sector.index),
    facing = "clockwise",
    niceFacing = FALSE,cex = 0.8,
    adj = c(0, 0.5),
  )
}, bg.border = NA)

dev.off()
while (!is.null(dev.list()))  dev.off()
circos.clear()

##############################################################################
f03_V_chord<-e01_Section_CGC%>%
  filter(grepl("\\bV",Section))

circos.clear()
cairo_pdf("01c_chord_diagram_V.pdf")  # Adjust dimensions as needed
chordDiagramFromDataFrame(f03_V_chord %>% 
                            select(Section_Class,CAZy_Substrate),
                          order = e04_vector_sector,
                          grid.col = structure(e07_link_cols[match(gsub(".*:::", "", e04_vector_sector),
                                                                   labels(e07_link_cols))],  
                                               names = e04_vector_sector),
                          link.lty = e01_Section_CGC$Substrate_Class_CGC_Fam_Section_relative_expression, 
                          target.prop.height = 1,
                          
                          directional = 1,
                          annotationTrack = c("grid", "axis"),
                          preAllocateTracks = 1)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(
    CELL_META$xcenter,
    CELL_META$ylim[1]+cm_h(1.1),
    sub(":::.*", "", CELL_META$sector.index),
    facing = "clockwise",
    niceFacing = FALSE,cex = 0.8,
    adj = c(0, 0.5),
  )
}, bg.border = NA)

dev.off()
while (!is.null(dev.list()))  dev.off()
circos.clear()

##############################################################################
##############################################################################
##############################################################################
##############################################################################
# Summary table
f04_table_carbo_elements<-b02_element_expression%>%
  pivot_wider(id_cols = c("Class","Class_order"),names_from = c("Element","Section"), values_from = Section_Class_Element_relative_expression, values_fill = 0)%>%
  select(Class,Class_order,starts_with("CAZyme"),starts_with("Sulfatase"),starts_with("CGC"))

f05_CGC_tops<-a02_annotation%>%
  select(Class,Genus,genomeID,query_name,CGC_ID,gene_type_DBCAN,Fam_CAZy,Fam_SULFATLAS)%>%
  filter(grepl("CGC",CGC_ID))%>%
  inner_join(a01_TPM, by ="query_name")%>%
  group_by(Section)%>%
  mutate(CGC_Section_TPM=sum(TPM))%>%
  mutate(Fam_CAZy = strsplit(Fam_CAZy, "\\|"))%>%
  unnest(Fam_CAZy)%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_CAZy))%>%
  filter(grepl("GH|PL|CE",Fam_CAZy))%>%
  mutate(Fam_CAZy = if_else(!CGC_ID %in% unique(d03_highly_expressing_CGCs$CGC_ID), "degradative CAZyme within lower quarters CGCs", 
                            if_else(!CGC_ID %in% unique(d05_classify$CGC_ID), "degradative CAZyme within unclassified CGCs",
                                    Fam_CAZy)))%>%
  left_join(d05_classify,by="CGC_ID")%>%
  mutate(Substrate = if_else(!CGC_ID %in% unique(d03_highly_expressing_CGCs$CGC_ID), "unknown", 
                             if_else(!CGC_ID %in% unique(d05_classify$CGC_ID), "unknown",
                                     Substrate)))%>%
  select(Sample,Section,Substrate,query_name,Fam_CAZy,Class,CGC_Section_TPM,TPM)%>%
  distinct()%>%
  group_by(Section,Class,Substrate,Fam_CAZy)%>%
  mutate(Substrate_Class_CGC_Fam_Section=sum(TPM))%>%
  mutate(Substrate_Class_CGC_Fam_Section_relative_expression=(Substrate_Class_CGC_Fam_Section/CGC_Section_TPM)*100)%>%
  ungroup()%>%
  select(Substrate,Class,Fam_CAZy,Section,Substrate_Class_CGC_Fam_Section_relative_expression)%>%
  filter(Substrate_Class_CGC_Fam_Section_relative_expression > 0)%>%
  distinct()%>%
  mutate(Class = gsub("c__", "", Class))%>%
  mutate(`CAZy_Substrate` = paste0(as.character(Fam_CAZy), ":::", as.character(Substrate)))%>%
  mutate(`Section_Class` = paste0(as.character(Section), ":::", as.character(Class)))%>%
  select(Class,Section,Fam_CAZy,Substrate_Class_CGC_Fam_Section_relative_expression,Substrate)%>%
  pivot_wider(id_cols = c("Class","Fam_CAZy","Substrate"),names_from = Section, values_from = Substrate_Class_CGC_Fam_Section_relative_expression, values_fill = 0)%>%
  right_join(f04_table_carbo_elements,by="Class")%>%
  arrange(desc(Class_order),
          match(Substrate, c("alginate","laminarin","FCSP","carrageenan","agarose","starch","unknown")),
          match(Fam_CAZy, c((e03_Substrate_Fam_vec$Fam_CAZy),"degradative CAZyme within lower quarters CGCs","degradative CAZyme within unclassified CGCs")))%>%
  select(Class,Class_order,starts_with("CAZyme"),starts_with("Sulfatase"),starts_with("CGC"),Fam_CAZy,Substrate,IV,V)

library(readxl)
writexl::write_xlsx(f05_CGC_tops,"01d_CGCs_expression.xlsx")
