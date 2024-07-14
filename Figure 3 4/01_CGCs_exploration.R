# Explore CGC counts and classify
setwd("Manuscript_1/Figure 3 4/")
library(dplyr)
library(readxl)
library(tidyverse)

##############################################################
##############################################################
##############################################################
# Tidy metadata from CAZy
a01_fam_act<-read.delim("../local_files/CAZyDB.08062022.fam-activities.txt",header =FALSE)%>%
  dplyr::rename(Fam_CAZy=V1,Fam_activity=V2)%>%
  select(Fam_CAZy,Fam_activity)

# filter all the CAZy annotation within the metagenome
a02_screening_fam_and_sub<-read.csv("../Annotation_table.csv",header=TRUE)%>%
  select(Fam_sub_CAZy)%>%
  na.omit()%>%
  mutate(Fam_sub_CAZy = strsplit(Fam_sub_CAZy, "\\|"))%>%
  unnest(Fam_sub_CAZy)%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_sub_CAZy))%>%
  distinct()%>%
  select(Fam_CAZy,Fam_sub_CAZy)

# Families of interest (based on thei presence within the metagenome)
a03_Fam_tidy <- a02_screening_fam_and_sub %>%
  select(Fam_CAZy) %>%
  left_join(a01_fam_act, by = "Fam_CAZy")%>%
  mutate(Fam_activity = gsub("Iota;-carrageenase \\(EC 3.2.1.162\\)", "Lambda-carrageenase \\(EC 3.2.1.162\\)", Fam_activity))%>% # fix typos of greek namings
  mutate(Fam_activity = gsub(";-", "-", Fam_activity))%>% # fix typos of greek namings
  mutate(Fam_activity = gsub("; ", "::", Fam_activity))%>%
  mutate(Name = strsplit(Fam_activity, "::")) %>%
  unnest(Name)%>%
  mutate_all(str_trim)%>%
  mutate_all(na_if, "")%>%
  na.omit()%>%
  distinct()%>%
  group_by(Fam_CAZy)%>%
  mutate(n_activities_per_family=n())%>%
  mutate(EC_1= sub(".*?([0-9]\\..*)", "\\1", Name))%>%
  mutate(EC_2= gsub("\\).*| .*", "", EC_1))%>%
  mutate(EC= gsub(".*[A-z][A-z].*", NA, EC_2))%>%
  select(Fam_CAZy,Fam_activity,Name,n_activities_per_family,EC)%>%
  distinct()%>%
  group_by(Fam_CAZy,EC)%>%
  mutate(unespecific_act_per_fam=n())%>%
  filter(unespecific_act_per_fam < 2)%>%
  filter(!grepl("amino acid transporter",Name))

##############################################################
##############################################################
##############################################################
# CGCs
# appending to genomes
b01_annot<-read.csv("../Annotation_table.csv")%>%
  select(Class,genomeID,MAG,query_name,Fam_sub_CAZy,EC_dbcan,Fam_SULFATLAS,pident_cazyEC,EC_cazyEC,sseqid_pfam,sseqid_tigrfam)%>%
  filter(grepl("FM",MAG))

b02_cgc_numbering<-read.delim("../local_files/020305_CGC_annotation.txt")%>%
  dplyr::rename(contig=contigID_DBCAN)%>%
  select(query_name,contig,CGC_number_DBCAN)%>%
  na.omit()%>%
  left_join(b01_annot,by="query_name")%>%
  select(MAG,contig,CGC_number_DBCAN)%>%
  arrange(MAG,contig,CGC_number_DBCAN)%>%
  distinct()%>%
  group_by(MAG)%>% 
  mutate(id = row_number())%>%
  mutate(id = str_replace(id, "(^[0-9]$)", "0\\1"))%>%
  mutate(genome_number = str_replace(MAG, ".*_", "CGC"))%>%
  mutate(`CGC_ID` = paste0(as.character(genome_number), "_", as.character(id)))%>%
  mutate(`genome_contig_cgc` = paste0(as.character(MAG), "_", as.character(contig), "_", as.character(CGC_number_DBCAN)))%>%
  ungroup()%>%
  select(genome_contig_cgc,CGC_ID)%>%
  distinct()

b03_CGCIDs<-read.delim("../local_files/020305_CGC_annotation.txt")%>%
  dplyr::rename(contig=contigID_DBCAN)%>%
  select(query_name,contig,CGC_number_DBCAN,gene_type_DBCAN)%>%
  right_join(b01_annot,by = "query_name")%>%
  select(MAG,contig,query_name,CGC_number_DBCAN,gene_type_DBCAN)%>%
  distinct()%>%
  na.omit()%>%
  mutate(`genome_contig_cgc` = paste0(as.character(MAG), "_", as.character(contig), "_", as.character(CGC_number_DBCAN)))%>%
  select(MAG,query_name,genome_contig_cgc,gene_type_DBCAN)%>%
  left_join(b02_cgc_numbering,by = "genome_contig_cgc")%>%
  mutate(CGC_gene = str_replace(query_name, ".*_", ""))%>%
  mutate(CGC_gene = str_replace(CGC_gene, "(^[0-9][0-9]$)", "0\\1"))%>%
  mutate(CGC_gene = str_replace(CGC_gene, "(^[0-9]$)", "00\\1"))%>%
  mutate(`CGC_gene` = paste0(as.character(CGC_ID), "_", as.character(CGC_gene), ""))%>%
  left_join(b01_annot,by = c("MAG","query_name"))%>%
  select(Class,genomeID,query_name,CGC_ID,CGC_gene,gene_type_DBCAN)

b04_cgc_content <- b01_annot %>%
  select(Class,genomeID,query_name,Fam_sub_CAZy, Fam_SULFATLAS)%>%
  right_join(b03_CGCIDs,by=c("Class","genomeID","query_name"))%>%
  unite(Full_CAZy_sulf, Fam_sub_CAZy, Fam_SULFATLAS, na.rm = TRUE, sep = '+')%>%
  mutate(gene_type_DBCAN_modified = gsub("CAZyme|null", NA, gene_type_DBCAN)) %>%
  mutate_all(na_if, "") %>%
  unite(Full_annot, gene_type_DBCAN_modified, Full_CAZy_sulf, na.rm = TRUE, sep = '+') %>%
  mutate_all(na_if, "") %>%
  mutate(Full_annot = replace(Full_annot, is.na(Full_annot), "unk"))%>%
  group_by(CGC_ID)%>%
  mutate(CGC_content = paste0(Full_annot, collapse = "::"))%>%
  distinct()%>%
  left_join(b01_annot,by=c("Class","genomeID","query_name"))%>%
  select(Class,genomeID,query_name,Full_annot,Fam_sub_CAZy,EC_dbcan,pident_cazyEC,EC_cazyEC,sseqid_pfam,sseqid_tigrfam,CGC_ID,CGC_gene,CGC_content)
write.csv(b04_cgc_content, '01a_CGC_content_long.csv', row.names = FALSE)

##############################################################
##############################################################
##############################################################
# total and deg
# CGCs counts exploration
c01_CGCs<-b04_cgc_content%>%
  mutate(Class = gsub(".*__", "", Class))
  
# CGCs Total per Phylum and average per phylum
c02_CGC_total<-c01_CGCs%>%
  select(Class,genomeID,CGC_ID)%>%
  distinct()%>%
  group_by(genomeID)%>%
  mutate(count_CGC_per_genome=n())%>%
  ungroup()%>%
  select(Class,genomeID,count_CGC_per_genome)%>%
  distinct()%>%
  # Total CGC count
  ungroup()%>%
  mutate(count_CGC=sum(count_CGC_per_genome))%>%
  # Total CGC per Class
  group_by(Class)%>%
  mutate(count_CGC_per_Class=sum(count_CGC_per_genome))%>%
  ungroup()%>%
  select(Class,genomeID,count_CGC_per_genome,count_CGC_per_Class,count_CGC)%>%
  distinct()%>%
  # Mean SD CGC per Class
  group_by(Class)%>%
  mutate(mean_CGC_per_Class=mean(count_CGC_per_genome),
         sd_CGC_per_Class=sd(count_CGC_per_genome))%>%
  mutate(CGC_type = "Total")

# CGCs degradative action (containing more than 2 degradative CAZymes GH, PL or CE) and average per phylum
c03_CGCs_deg2<-c01_CGCs%>%
  select(Class,genomeID,query_name,CGC_ID,Fam_sub_CAZy)%>%
  filter(grepl("GH|PL|CE",Fam_sub_CAZy))%>%
  group_by(CGC_ID)%>%
  mutate(Deg_count=n())%>%
  filter(Deg_count > 1)%>%
  ungroup()%>%
  select(Class,genomeID,CGC_ID)%>%
  distinct()%>%
  group_by(genomeID)%>%
  mutate(count_CGC_per_genome=n())%>%
  ungroup()%>%
  select(Class,genomeID,count_CGC_per_genome)%>%
  distinct()%>%
  # Deg 2 CGC count
  ungroup()%>%
  mutate(count_CGC=sum(count_CGC_per_genome))%>%
  # Deg 2 CGC per Class
  group_by(Class)%>%
  mutate(count_CGC_per_Class=sum(count_CGC_per_genome))%>%
  ungroup()%>%
  select(Class,genomeID,count_CGC_per_genome,count_CGC_per_Class,count_CGC)%>%
  distinct()%>%
  # Mean SD per Class
  group_by(Class)%>%
  mutate(mean_CGC_per_Class=mean(count_CGC_per_genome),
         sd_CGC_per_Class=sd(count_CGC_per_genome))%>%
  mutate(CGC_type = "Degradative")

# Join to dataframes
c04_CGC_summary<-rbind(c02_CGC_total,c03_CGCs_deg2)

# Plot only Bacteroidia, Clostridia and Gammaproteobacteria
c05_CGC_summary_abundant<-c04_CGC_summary%>%
  filter(grepl("Bacteroidia|Clostridia", Class))%>%
  filter(grepl("Total|Degradative\\b",CGC_type))%>%
  mutate(Class = factor(Class, levels = c("Bacteroidia","Clostridia")))

c06_class_color<-c("Bacteroidia"="#B2DF8A",
            "Clostridia"="#FF7F00")

c07_class_order<-c("Bacteroidia", "Clostridia")

c08_CGC_counts<-ggplot(c05_CGC_summary_abundant, aes(x = Class, y = mean_CGC_per_Class, fill = Class)) +
  geom_errorbar(aes(ymin = mean_CGC_per_Class - sd_CGC_per_Class, ymax = mean_CGC_per_Class + sd_CGC_per_Class), width = 0.3, size = 0.2, position = position_dodge(width = 0.9)) +
  #geom_line(aes(x = Class, y = mean_CGC_per_Class, group = CGC_type), color = "gray", size = 0.2) +
  geom_point(position = position_jitterdodge(jitter.width = 1, dodge.width = 1),
             #aes(y = count_CGC_per_genome, fill = factor(Class, levels = b08_phylum_order)), pch = 21,
             aes(y = count_CGC_per_genome), pch = 21,
             size = 1, stroke = 0.1,
             show.legend = FALSE, alpha = 0.8) +
  geom_text(aes(label = count_CGC_per_Class, y = 48), vjust = -0.5, size = 2, position = position_dodge(width = 1)) +
  geom_point(shape = 17, color = "black", size = 1, position = position_dodge(width = 0.9), show.legend = FALSE) +
  scale_fill_manual(values = c06_class_color) +
  scale_y_continuous(name = "Average CGC per Class", limits = c(0, 50), breaks = seq(0, 50, 10)) +
  facet_grid(cols = vars(factor(CGC_type, levels = c("Total","Degradative"))), scales = "free", space = "free", labeller = label_wrap_gen(multi_line = TRUE)) +
  theme_bw() +
  theme(
    axis.title.y = element_text(angle = 90, size = 6, face = "bold"),
    axis.text.y = element_text(angle = 0, size = 5),
    strip.text.x = element_text(angle = 0, size = 6, face = "bold"),
    axis.title.x = element_text(angle = 0, size = 6, face = "bold"),
    axis.text.x = element_text(angle = 0, size = 5),
    panel.grid.major = element_line(size = 0.1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.5),
    axis.ticks = element_line(size = 0.2))
c08_CGC_counts
ggsave(filename = "01b_CGC_summary.pdf", c08_CGC_counts,width = 70, height = 70 , units = "mm")

# CGCs outside Bacteroidia and Clostridia
c09_count_other_class<-c04_CGC_summary%>%
  filter(!grepl("Bacteroidia|Clostridia", Class))%>%
  filter(grepl("Total|Degradative\\b",CGC_type))%>%
  select(Class,genomeID,count_CGC_per_genome,count_CGC_per_Class,CGC_type)%>%
  distinct()%>%
  arrange(desc(CGC_type),desc(count_CGC_per_Class))

c10_other_phyla<-ggplot(c09_count_other_class, aes(x=factor(CGC_type,levels = c("Total","Degradative")),y = count_CGC_per_Class, fill = CGC_type)) +
  facet_grid(cols = vars(factor(Class,
                                levels=c("Gammaproteobacteria","Verrucomicrobiae",
                                         "Desulfovibrionia","Bacilli",
                                         "Spirochaetia","Alphaproteobacteria",
                                         "Vampirovibrionia")))
             , scales = "free_y", space = "free", labeller = label_wrap_gen(multi_line = TRUE)) +
  geom_bar(stat = "identity",position = "dodge", size = 0.1) +
  scale_fill_manual(values = c("Total"="#303A52",
                               "Degradative"="#94ACBF"), name = "CGC type") +
  scale_y_continuous(name = "CGC count", limits = c(0, 70), breaks = seq(0, 65, 5), expand = c(0, 1)) +
  scale_x_discrete(name = "Class") +
  theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 5),
        axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 7, face="bold"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 7, face="bold"),
        axis.ticks.y = element_line(size = 0.25),
        axis.ticks.x = element_blank(),
        
        strip.background = element_blank(),
        strip.text = element_blank(),
        rect = element_rect(size = 0.1),
        line = element_line(size = 0.1),
        plot.margin = margin(1,1,1,1, unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 0.5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6, face = "bold"),
        legend.key.size = unit(4, 'mm'),
        legend.position = "right") +
  geom_text(aes(label = Class, vjust = 1,hjust=0),x="Degradative",y=30,angle=90, size = 2, show.legend = F) +
  geom_text(aes(label = count_CGC_per_Class, y = count_CGC_per_Class + 2), vjust = 0, size = 2, position = position_dodge(width = 1))
c10_other_phyla
ggsave(filename = "01c_CGC_other_summary.pdf", c10_other_phyla,width = 120, height = 70 , units = "mm")

library(patchwork)
c11_CGC_mean_distr <-  c08_CGC_counts + c10_other_phyla + plot_layout(width = c(70,120))
c11_CGC_mean_distr
ggsave("01d_CGC_mean_and_distr.pdf", c11_CGC_mean_distr, width = 180, height = 70,units = "mm")

##############################################################
##############################################################
##############################################################
##############################################################
# Classify CGCs based on their CAZyme content
# inspected cazymes with similarity and what cazymes pair with them
# from the classified cgcs, picked other families with similar function
# inspect if the pfam and tigrfam annotation are consistent
d01_CGCs_deg2_classify<-c01_CGCs%>%
  select(Class,genomeID,query_name,CGC_ID,Fam_sub_CAZy)%>%
  filter(grepl("GH|PL|CE",Fam_sub_CAZy))%>%
  group_by(CGC_ID)%>%
  mutate(Deg_count=n())%>%
  filter(Deg_count > 1)%>%
  ungroup()%>%
  select(CGC_ID)%>%
  distinct()%>%
  left_join(c01_CGCs,by = "CGC_ID")%>%
  select(Class,genomeID,CGC_ID,CGC_content)%>%
  distinct()

# classify based on pairs indicate pathway
d02_CGCs_specificity<-d01_CGCs_deg2_classify%>%
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
write.csv(d02_CGCs_specificity, '01e_CGC_ID_substrate_specificity.csv', row.names = FALSE)

##############################################################
##############################################################
##############################################################
##############################################################
# CAZy families percentage composition within CGC
e01_fam_CAZy_top<-d02_CGCs_specificity%>%
  select(Substrate,CGC_ID)%>%
  left_join(c01_CGCs,by = "CGC_ID")%>%
  select(Substrate,CGC_ID,query_name,Fam_sub_CAZy)%>%
  na.omit()%>%
  mutate(Fam_CAZy = strsplit(Fam_sub_CAZy, "\\|"))%>%
  unnest(Fam_CAZy)%>%
  mutate(Fam_CAZy = gsub("_.*", "", Fam_CAZy))%>%
  filter(grepl("GH|PL|CE",Fam_CAZy))%>%
  select(Substrate,CGC_ID,query_name,Fam_CAZy)%>%
  distinct()%>%
  # Total genes annotated as CAZymes per substrate
  group_by(Substrate)%>%
  mutate(Substrate_gene_total=n())%>%
  # Total CAZy fams annotated per substrate
  group_by(Substrate,Fam_CAZy)%>%
  mutate(Substrate_Fam_count=n())%>%
  select(Substrate,CGC_ID,Fam_CAZy,Substrate_gene_total,Substrate_Fam_count)%>%
  distinct()%>%
  group_by(Substrate,Fam_CAZy)%>%
  mutate(CGCs_with_fam=n())%>%
  select(Substrate,Fam_CAZy,Substrate_gene_total,Substrate_Fam_count,CGCs_with_fam)%>%
  arrange(Substrate,desc(Substrate_gene_total),desc(CGCs_with_fam),desc(Substrate_Fam_count))%>%
  distinct()%>%
  mutate(Percentage=((Substrate_Fam_count/Substrate_gene_total)*100))%>%
  select(Substrate,Fam_CAZy,Substrate_gene_total,Substrate_Fam_count,Percentage)%>%
  distinct()%>%
  arrange(Substrate,desc(Percentage))%>%
  group_by(Substrate)%>%
  filter(row_number() <= n() / 3)%>%
  mutate(Algae = ifelse(Substrate %in% c("alginate", "laminarin", "FCSP"), "Brown", "Red"))%>%
  arrange(Substrate,desc(Percentage))%>%
  group_by(Substrate)%>%
  mutate(Perc_cumulative = cumsum(Percentage))%>%
  mutate(Perc_cumulative_mid = cumsum(Percentage)-(Percentage/2))

Supplementary_table<-e01_fam_CAZy_top%>%
  select(Algae,Substrate,Fam_CAZy,Substrate_gene_total,Substrate_Fam_count,Percentage)%>%
  arrange(Algae,Substrate,desc(Percentage))%>%
  rename(Total_degradative_CAZymes_within_CGC=Substrate_gene_total,Degradative_CAZy_families_within_CGC=Substrate_Fam_count,CAZy_family_percentage=Percentage,CAZy_family=Fam_CAZy)
write.csv(Supplementary_table, '01f_Suppplementary_cazyme_frequency.csv', row.names = FALSE)

  