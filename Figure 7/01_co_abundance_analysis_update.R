# Abundance based on groups
setwd("Manuscript_1/Figure 7")
library(broom)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)

##################################################################################
# load data
a01_annotation<-read.csv("../Annotation_table.csv")

a02_taxonomy<-a01_annotation%>%
  select(Phylum,Class,Order,Family,Genus,genomeID,MAG)%>%
  distinct()%>%
  mutate(Family = coalesce(Family, Order))%>%
  mutate(Genus = coalesce(Genus, Family))%>%
  mutate(Phylum = gsub("p__", "", Phylum))%>%
  mutate(Class = gsub("c__", "", Class))%>%
  mutate(Order = gsub("o__", "", Order))%>%
  mutate(Family = gsub("f__", "", Family))%>%
  mutate(Genus = gsub("g__", "", Genus))

a03_coabund_groups<-read_excel("Co_abundant_group_assignment.xlsx")%>%
  select(genomeID,Group)

tax_color<-c("Bacteroidia"="#B2DF8A",
             "Bacilli"="#FDBF6F",
             "Clostridia"="#FF7F00",
             "Spirochaetia"="#FB9A99",
             "Desulfovibrionia"="#FFFF99",
             "Verrucomicrobiae"="#B15928",
             "Alphaproteobacteria"="#CAB2D6",
             "Gammaproteobacteria"="#6A3D9A",
             "Vampirovibrionia"="#1F78B4")

a04_refined_cazymes<-read_xlsx("../Figure 3 4/03a_CAZymes_path_refined.xlsx")

a05_CGC_gene_class<-read.csv("../Figure 3 4/01e_CGC_ID_substrate_specificity.csv")%>%
  select(CGC_ID,Substrate)

a06_CGC_all<-read.csv("../Figure 3 4/01a_CGC_content_long.csv")%>%
  select(genomeID,query_name,CGC_ID,CGC_content)

a07_order_groups<-c(
  "f__Anaeroplasmataceae FM_41",
  "g__Alistipes FM_07",
  "f__Psittacicellaceae FM_32",
  "g__JAGHEK01 FM_55",
  "g__Ferrimonas FM_34",
  "f__Anaeroplasmataceae FM_42",
  "g__DUWA01 FM_45",
  "f__CAG-274 FM_61",
  "g__Anaerotignum FM_62",
  "g__Lawsonibacter FM_47",
  "g__Lawsonibacter FM_48",
  "f__UBA660 FM_44",
  "f__Lachnospiraceae FM_68",
  "f__Butyricicoccaceae FM_51",
  "g__Faecalibacterium FM_53",
  "f__Lachnospiraceae FM_63",
  "g__RQZE01 FM_43",
  "f__Lachnospiraceae FM_64",
  "g__Akkermansia FM_36",
  "f__Oscillospiraceae FM_46",
  "g__UBA6857 FM_54",
  "g__RGIG6307 FM_65",
  "g__SIG603 FM_50",
  "g__Alistipes FM_02",
  "g__Alistipes FM_05",
  "f__Treponemataceae FM_38",
  "g__Alistipes FM_09",
  "g__Alistipes FM_08",
  "g__Alistipes FM_16",
  "g__SIG603 FM_49",
  "f__Lachnospiraceae FM_66",
  "g__Nitratidesulfovibrio FM_37",
  "g__Alistipes FM_03",
  "g__HGM05232 FM_30",
  "f__Rikenellaceae FM_28",
  "g__Mucinivorans FM_24",
  "f__Rikenellaceae FM_29",
  "g__Alistipes FM_18",
  "g__Mucinivorans FM_22",
  "g__Alistipes FM_14",
  "f__Lachnospiraceae FM_67",
  "g__Alistipes FM_01",
  "g__Alistipes FM_17",
  "g__Alistipes FM_06",
  "g__Alistipes FM_12",
  "f__Rikenellaceae FM_27",
  "f__Butyricicoccaceae FM_52",
  "g__Alistipes FM_13",
  "o__Christensenellales FM_57",
  "g__Alistipes FM_10",
  "g__Alistipes FM_15",
  "g__Aphodosoma FM_31",
  "f__Rikenellaceae FM_21",
  "s__Vibrio pomeroyi FM_33",
  "g__WQYD01 FM_56",
  "g__CALURL01 FM_40",
  "o__Christensenellales FM_59",
  "g__Alistipes FM_11",
  "g__Mucinivorans FM_23",
  "g__Alistipes FM_19",
  "g__Alistipes FM_04",
  "g__CAKUIA01 FM_35",
  "o__Christensenellales FM_60",
  "o__RF32 FM_39",
  "f__Rikenellaceae FM_25",
  "f__Rikenellaceae FM_26",
  "g__Alistipes FM_20",
  "o__Christensenellales FM_58")

#################################################################################
#################################################################################
#################################################################################
# relative abundance of groups
##################################################################################
# ABUNDANCE AND METADATA
# Relative Abundance
b01_relative_abundance <- read.delim("../local_files/020205_normalised_summary_count_table_genomeSummary.tsv")%>%
  select(genomeID, ends_with("_pileup_counts_lengthNorm_libNorm"))%>%
  dplyr::rename(
    MAG=genomeID,
    G117_V=S4_pileup_counts_lengthNorm_libNorm,
    G121_V=S2_pileup_counts_lengthNorm_libNorm,
    G124_V=S6_pileup_counts_lengthNorm_libNorm,
    G125_V=S3_pileup_counts_lengthNorm_libNorm,
    G117_IV=S8_pileup_counts_lengthNorm_libNorm,
    G121_IV=S7_pileup_counts_lengthNorm_libNorm,
    G124_IV=S1_pileup_counts_lengthNorm_libNorm,
    G125_IV=S5_pileup_counts_lengthNorm_libNorm)%>%
  left_join(a02_taxonomy,by = "MAG")%>%
  mutate(across(starts_with("G1"),~ ((.) / sum(.)) * 100)) %>%
  left_join(a03_coabund_groups,by = "genomeID")%>%
  pivot_longer(cols=-c("Group","Phylum","Class","Order","Family","Genus","genomeID","MAG"),names_to = "Sample",values_to = "Relative_abundance")%>%
  separate(Sample, into = c("Fish", "Section"), sep = "_",remove = FALSE)%>%
  select(Group,Phylum,Class,Genus,genomeID,Sample,Fish,Section,Relative_abundance)

b02_sum_relative_abundance<-read.delim("../local_files/020205_normalised_summary_count_table_genomeSummary.tsv")%>%
  select(genomeID, ends_with("_pileup_counts_lengthNorm_libNorm"))%>%
  dplyr::rename(
    MAG=genomeID,
    G117_V=S4_pileup_counts_lengthNorm_libNorm,
    G121_V=S2_pileup_counts_lengthNorm_libNorm,
    G124_V=S6_pileup_counts_lengthNorm_libNorm,
    G125_V=S3_pileup_counts_lengthNorm_libNorm,
    G117_IV=S8_pileup_counts_lengthNorm_libNorm,
    G121_IV=S7_pileup_counts_lengthNorm_libNorm,
    G124_IV=S1_pileup_counts_lengthNorm_libNorm,
    G125_IV=S5_pileup_counts_lengthNorm_libNorm)%>%
  left_join(a02_taxonomy,by = "MAG")%>%
  left_join(a03_coabund_groups,by = "genomeID")%>%
  pivot_longer(cols=-c("Group","Phylum","Class","Order","Family","Genus","genomeID","MAG"),names_to = "Sample",values_to = "Counts")%>%
  separate(Sample, into = c("Fish", "Section"), sep = "_",remove = FALSE)%>%
  select(Group,Phylum,Class,Genus,genomeID,Fish,Section,Counts)%>%
  group_by(Section)%>%
  mutate(Total_section=sum(Counts))%>%
  group_by(genomeID,Section)%>%
  mutate(Sum_genomeID_Section=sum(Counts))%>%
  select(Group,Phylum,Class,Genus,genomeID,Section,Sum_genomeID_Section,Total_section)%>%
  distinct()%>%
  mutate(sum_Relative_abundance=(Sum_genomeID_Section/Total_section)*100)%>%
  select(Group,Phylum,Class,Genus,genomeID,Section,sum_Relative_abundance)%>%
  distinct()
  
#################################################################################
#################################################################################
#################################################################################
# Relative abundance of groups and their coexistence
c01_group_abundance<-b01_relative_abundance%>%
  group_by(Group,Sample)%>%
  mutate(Group_Sample_abundance=sum(Relative_abundance))%>%
  select(Group,Sample,Fish,Section,Group_Sample_abundance)%>%
  distinct()
library(ggplot2)
ggplot(c01_group_abundance, aes(x = Fish, y = Group_Sample_abundance, fill = Group)) +
  geom_bar(show.legend = FALSE, position = position_stack(reverse = TRUE), stat = "identity") +
  facet_grid(cols = vars(factor(Section)), scales = "free", space = "free", labeller = label_wrap_gen(multi_line = TRUE)) +
  #scale_fill_manual(values = phylum_color) +
  scale_y_continuous(name = "", limits = c(0, 101), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, size = 5, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(angle = 90, size = 5, hjust = 0.5, vjust = 1),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", size = 0.1, linetype = "solid"),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.background = element_rect(colour = "black", size = 0.1),
    axis.ticks = element_line(size = 0.1),
    plot.margin = margin(5, 5, 5, 5),  # Adjust the plot margin here
    panel.spacing = unit(0, "lines")
  )

#################################################################################
#################################################################################
#################################################################################
# Barplots highlighting Phylum relative abundance within Groups across samples
d01_Group_Class<-b01_relative_abundance%>%
  group_by(Class,Group,Sample)%>%
  mutate(Class_Group_Sample_abundance=sum(Relative_abundance))%>%
  select(Class,Group,Sample,Fish,Section,Class_Group_Sample_abundance)%>%
  distinct()

library(ggplot2)
d02_isolated_group_abundance<-ggplot(d01_Group_Class, aes(x = Fish, y = Class_Group_Sample_abundance, fill = Class)) +
  geom_bar(show.legend = FALSE, position = position_stack(reverse = TRUE), stat = "identity", color = "black", size = 0.1) +
  facet_grid(cols = vars(factor(Section)),rows = vars(factor(Group)), scales = "free_x", space = "free", labeller = label_wrap_gen(multi_line = TRUE)) +
  scale_fill_manual(values = tax_color) +
  scale_y_continuous(name = "", limits = c(0, 98), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_discrete(name = "") +
  theme_bw() +
  labs(title = "b)") +
  theme(
    plot.title = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(angle = 90, size = 5, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5),
    axis.ticks = element_line(size = 0.1,),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_rect(size = 0.5),
    strip.text.x = element_text(angle=0,size = 5, hjust = 0.5, vjust = 0.5,margin = margin(1,1,1,1, "mm")),
    strip.text.y = element_text(angle=0,size = 5, hjust = 0.5, vjust = 0.5,margin = margin(1,1,1,1, "mm")),
    plot.margin = margin(1,1,1,1,unit = "mm"),  # Adjust the plot margin here
    panel.spacing = unit(1, "mm"),
    panel.grid = element_blank(),
    #panel.grid = element_line(color = "gray", size = 0.1),  # Adjust the gridline size and color here
    panel.border = element_rect(color = "black", size = 0.5)
  )
d02_isolated_group_abundance  
ggsave(filename = "01a_Groups_highlight_vert.pdf", d02_isolated_group_abundance, width = 30, height = 280 , units = "mm")

#################################################################################
#################################################################################
#################################################################################
# Mean relative abundance
e01a_group_mean_relative_abundance<-b01_relative_abundance%>%
  group_by(Group,Section)%>%
  mutate(mean_Group_Sec=mean(Relative_abundance),
         sd_Group_Sec=sd(Relative_abundance))%>%
  ungroup()%>%
  select(Group,Section,mean_Group_Sec,sd_Group_Sec)%>%
  distinct()

e01b_group_cla_mean_relative_abundance<-b01_relative_abundance%>%
  group_by(Group,Phylum,Class,Section)%>%
  mutate(mean_Group_cla_Sec=mean(Relative_abundance),
         sd_Group_cla_Sec=sd(Relative_abundance))%>%
  ungroup()%>%
  select(Group,Phylum,Class,Section,mean_Group_cla_Sec,sd_Group_cla_Sec)%>%
  distinct()

e01c_group_gen_mean_relative_abundance<-b01_relative_abundance%>%
  group_by(Group,Phylum,Class,Genus,Section)%>%
  mutate(mean_Group_Gen_Sec=mean(Relative_abundance),
         sd_Group_Gen_Sec=sd(Relative_abundance))%>%
  ungroup()%>%
  select(Group,Phylum,Class,Genus,Section,mean_Group_Gen_Sec,sd_Group_Gen_Sec)%>%
  distinct()

#################################################################################
# Sum relative abundance
e02a_group_sum_relative_abundance<-b02_sum_relative_abundance%>%
  group_by(Group,Section)%>%
  mutate(sum_Group_Sec=sum(sum_Relative_abundance))%>%
  ungroup()%>%
  select(Group,Section,sum_Group_Sec)%>%
  distinct()

e02b_group_cla_sum_relative_abundance<-b02_sum_relative_abundance%>%
  group_by(Group,Phylum,Class,Section)%>%
  mutate(sum_Group_Class_Sec=sum(sum_Relative_abundance))%>%
  ungroup()%>%
  select(Group,Phylum,Class,Section,sum_Group_Class_Sec)%>%
  distinct()

e02c_group_gen_sum_relative_abundance<-b02_sum_relative_abundance%>%
  group_by(Group,Phylum,Class,Genus,Section)%>%
  mutate(sum_Group_Gen_Sec=sum(sum_Relative_abundance))%>%
  ungroup()%>%
  select(Group,Phylum,Class,Genus,Section,sum_Group_Gen_Sec)%>%
  distinct()

#################################################################################
# Log fold change
# using Mean
# group LogFC
e03a_LOGFC_mean_group<-e01a_group_mean_relative_abundance%>%
  select(Group,Section,mean_Group_Sec)%>%
  pivot_wider(id_cols="Group", names_from = Section,values_from = mean_Group_Sec)%>%
  mutate(log10fold_change_mean=log10((V+1e-6)/(IV+1e-6)))%>%
  pivot_longer(cols = -c("Group","log10fold_change_mean"), names_to = "Section", values_to = "mean_Group_Sec")%>%
  arrange(log10fold_change_mean)

# group Class LogFC
e03b_LOGFC_mean_cla_group<-e01b_group_cla_mean_relative_abundance%>%
  select(Group,Phylum,Class,Section,mean_Group_cla_Sec)%>%
  pivot_wider(id_cols=c("Group","Phylum","Class"), names_from = Section,values_from = mean_Group_cla_Sec)%>%
  mutate(log10fold_change_mean=log10((V+1e-6)/(IV+1e-6)))%>%
  pivot_longer(cols = -c("Group","Phylum","Class","log10fold_change_mean"), names_to = "Section", values_to = "mean_Group_cla_Sec")%>%
  arrange(log10fold_change_mean)

# group Genus LogFC
e03c_LOGFC_mean_Gen_group<-e01c_group_gen_mean_relative_abundance%>%
  select(Group,Phylum,Class,Genus,Section,mean_Group_Gen_Sec)%>%
  pivot_wider(id_cols=c("Group","Phylum","Class","Genus"), names_from = Section,values_from = mean_Group_Gen_Sec)%>%
  mutate(log10fold_change_mean=log10((V+1e-6)/(IV+1e-6)))%>%
  pivot_longer(cols = -c("Group","Phylum","Class","Genus","log10fold_change_mean"), names_to = "Section", values_to = "mean_Group_Gen_Sec")%>%
  arrange(log10fold_change_mean)

# using Sum
# group LogFC
e04a_LOGFC_sum_group<-e02a_group_sum_relative_abundance%>%
  select(Group,Section,sum_Group_Sec)%>%
  pivot_wider(id_cols="Group", names_from = Section,values_from = sum_Group_Sec)%>%
  mutate(log10fold_change_sum=log10((V+1e-6)/(IV+1e-6)))%>%
  pivot_longer(cols = -c("Group","log10fold_change_sum"), names_to = "Section", values_to = "sum_Group_Sec")%>%
  arrange(log10fold_change_sum)

# group Phylum LogFC
e04b_LOGFC_sum_Phy_group<-e02b_group_cla_sum_relative_abundance%>%
  select(Group,Phylum,Class,Section,sum_Group_Class_Sec)%>%
  pivot_wider(id_cols=c("Group","Phylum","Class"), names_from = Section,values_from = sum_Group_Class_Sec)%>%
  mutate(log10fold_change_sum=log10((V+1e-6)/(IV+1e-6)))%>%
  pivot_longer(cols = -c("Group","Phylum","Class","log10fold_change_sum"), names_to = "Section", values_to = "sum_Group_Class_Sec")%>%
  arrange(log10fold_change_sum)

# group Genus LogFC
e04c_LOGFC_sum_Gen_group<-e02c_group_gen_sum_relative_abundance%>%
  select(Group,Phylum,Class,Genus,Section,sum_Group_Gen_Sec)%>%
  pivot_wider(id_cols=c("Group","Phylum","Class","Genus"), names_from = Section,values_from = sum_Group_Gen_Sec)%>%
  mutate(log10fold_change_sum=log10((V+1e-6)/(IV+1e-6)))%>%
  pivot_longer(cols = -c("Group","Phylum","Class","Genus","log10fold_change_sum"), names_to = "Section", values_to = "sum_Group_Gen_Sec")%>%
  arrange(log10fold_change_sum)

#################################################################################
# join sum and mean
# group
e05a_sum_mean_group<-e01a_group_mean_relative_abundance%>%
  left_join(e03a_LOGFC_mean_group,by = c("Group","Section","mean_Group_Sec"))%>%
  left_join(e04a_LOGFC_sum_group,by = c("Group","Section"))%>%
  select(Group,log10fold_change_sum,Section,sum_Group_Sec)%>%
  arrange(log10fold_change_sum)

e06_area_log_fold<-ggplot(e05a_sum_mean_group, aes(x = log10fold_change_sum,  y = sum_Group_Sec, group = Section)) +
  geom_area(aes(fill = Section), colour = "black", size = 0.2, alpha = 0.6, position = "identity", show.legend = FALSE) +
  geom_point(aes(fill = Section), shape = 16, size = 1, show.legend = FALSE) +
  scale_x_reverse(limits = c(0.62, -0.99),name = "Group log 10 fold change",breaks = seq(-0.95, 0.6, 0.1), expand = c(0, 0)) +
  coord_flip() +
  scale_y_continuous(name = "Sum group relative abundance (% Section)", limits = c(0, 80), breaks = seq(0, 80, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c("IV" = "#f061eb", "V" = "#205233"), name = "Section") +
  geom_label(aes(label = paste0(Group,"_",Section, " = ", sprintf("%.1f", sum_Group_Sec), " %"), 
                 group = Section, 
                 x = ifelse(Section == "IV", log10fold_change_sum - 0.012, log10fold_change_sum + 0.012), 
                 y = 62),
             fill = "white", colour = "black", size = 2, angle = 0) +
  geom_vline(aes(xintercept=log10fold_change_sum),
             color="black", linetype="dashed", size=0.2) +
  theme_bw() +
  labs(title = "a)") +
  theme(
    plot.title = element_text(size = 7, face = "bold"),
    axis.ticks = element_line(size = 0.1),
    axis.text.x = element_text(angle = 0, size = 5, hjust = 0.5, vjust =0.5),
    axis.text.y = element_text(angle = 0, size = 5, hjust = 0.5, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(1,1,1,1,unit = "mm"),  # Adjust the plot margin here
    panel.grid = element_blank(),
#    panel.grid = element_line(color = "gray", size = 0.1),  # Adjust the gridline size and color here
    panel.border = element_rect(color = "black", size = 0.5)
  ) 
e06_area_log_fold
ggsave(filename = "01b_Groups_logfold_vert.pdf", e06_area_log_fold, width = 65, height = 280 , units = "mm")

library(patchwork)
e07_log_relative_abundance <- e06_area_log_fold + d02_isolated_group_abundance + plot_layout(width = c(65, 35))
e07_log_relative_abundance
ggsave(filename = "01c_Groups_logfold_vert_abundance.pdf", e07_log_relative_abundance, width = 100, height = 280 , units = "mm")

#################################################################################
#################################################################################
#################################################################################
# Enzymatic activity distribution
# total number of genomes within group
# genomes containing activity within group
f01_Genomes_containing_in_Group<-a03_coabund_groups%>%
  group_by(Group)%>%
  mutate(Total_genomes_group=n())%>%
  left_join(a04_refined_cazymes,by = "genomeID")%>%
  select(Group,genomeID,EC_refined,Total_genomes_group)%>%
  distinct()%>%
  na.omit()%>%
  group_by(Group,EC_refined)%>%
  mutate(genomes_containing_EC_per_group=n())%>%
  select(Group,EC_refined,genomes_containing_EC_per_group,Total_genomes_group)%>%
  distinct()%>%
  pivot_wider(id_cols = c("Group","Total_genomes_group"),names_from =EC_refined,values_from =genomes_containing_EC_per_group,values_fill = 0)%>%
  pivot_longer(cols = -c("Group","Total_genomes_group"),names_to = "EC_refined",values_to = "genomes_containing_EC_per_group")%>%
  select(Group,genomes_containing_EC_per_group,Total_genomes_group,EC_refined)%>%
  distinct()

f02_Function<-a04_refined_cazymes%>%
  select(Fam_refined,EC_refined,Name)%>%
  distinct()%>%
  na.omit()%>%
  group_by(EC_refined)%>%
  mutate(Fams = paste0(Fam_refined, collapse = ", "))%>%
  mutate(`Function` = paste0(as.character(Fams), "  ", as.character(Name),"  ", as.character(EC_refined)))%>%
  select(Function,Name,EC_refined)%>%
  distinct()%>%
  left_join(a04_refined_cazymes, by =c("Name","EC_refined"))%>%
  select(Function,Name,EC_refined,Activity,Substrate)%>%
  distinct()

EC_vec<-c("4.2.2.11","4.2.2.3","4.2.2.26",
          "3.2.1.75","3.2.1.39",
          "3.2.1.21","3.2.1.51","3.2.1.38","3.2.1.159",
          "3.2.1.83","3.2.1.162","3.2.1.81","3.2.1.23","3.2.1.49","3.2.1.80")

f04_plot_bar<-a04_refined_cazymes%>%
  select(genomeID,query_name,EC_refined)%>%
  distinct()%>%
  group_by(genomeID,EC_refined)%>%
  mutate(EC_per_genome=n())%>%
  select(genomeID,EC_refined,EC_per_genome)%>%
  distinct()%>%
  right_join(a02_taxonomy,by = "genomeID")%>%
  right_join(a03_coabund_groups,by = "genomeID")%>%
  select(Group,Class,genomeID,EC_refined,EC_per_genome)%>%
  na.omit()%>%
  # pivot_wider(id_cols = c("Group","Phylum","Genus","genomeID"),names_from =EC_refined,values_from =EC_per_genome,values_fill = 0)%>%
  # select(-c("NA"))%>%
  # pivot_longer(cols = -c("Group","Phylum","Genus","genomeID"),names_to = "EC_refined",values_to = "EC_per_genome")%>%
  # Count EC per group
  group_by(Group,EC_refined)%>%
  mutate(Group_EC_sum=sum(EC_per_genome))%>%
  # Count EC per group per Phylum
  group_by(Class,Group,EC_refined)%>%
  mutate(Cla_Group_EC_sum=sum(EC_per_genome))%>%
  left_join(f01_Genomes_containing_in_Group,by = c("Group","EC_refined"))%>%
  left_join(f02_Function,by = "EC_refined", relationship = "many-to-many")%>%
  mutate(Activity = gsub("MW", "", Activity))%>%
  arrange(factor(EC_refined, levels= EC_vec))%>%
  group_by(Group)%>%
  mutate(CAZyme_per_group=sum(EC_per_genome))

f05_substrate_vec<-c("alginate","laminarin","FCSP","carrageenan")
f07_substrate_color<-c("laminarin"="#A16928",
                       "alginate"="#C29B64",
                       "FCSP"="#E0CFA2",
                       "carrageenan"="#CBD5BC")
f09_order_vec_name<-unique(f04_plot_bar$Name)
f09_order_vec_name
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Plot group capacity separately
library(ggnewscale)
library(ggbreak)
library(stringr)

f10_barplot_group_enzyme<-ggplot(f04_plot_bar, aes(y = factor(Name, levels = rev(f09_order_vec_name)), x = EC_per_genome, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", show.legend = FALSE) +
  facet_grid(rows = vars(factor(Group)), scales = "free_x", space = "free", labeller = label_wrap_gen(multi_line = TRUE)) + 
  scale_fill_manual(values = tax_color) +
  scale_x_continuous(name = "", limits = c(-2, 175),breaks = seq(0, 160, 20), expand = c(0, 0), position = "bottom") +
  scale_x_break(breaks = c(75, 125),space = 0.05) +
  scale_y_discrete(position = "right", labels = function(x) str_wrap(x, width = 15)) +
  theme_bw() +
  labs(title = "c)") +
  theme(
    plot.title = element_text(size = 7, face = "bold"),
    strip.text = element_blank(),
    axis.text.y = element_text(
      angle = 0, size = 5, hjust = 1, vjust = 0.5),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x.bottom = element_text(angle = 0, size = 6, hjust = 0.5, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.1),
    plot.margin = margin(1,1,1,1,unit = "mm"),  # Adjust the plot margin here
    panel.spacing = unit(1, "mm"),
    panel.border = element_rect(color = "black", size = 0.5),
    panel.grid = element_blank()
    ) +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(y = factor(Name, levels = as.character(rev(f09_order_vec_name))), x = -1, fill = Substrate),
            show.legend = FALSE,
            width = 2,
            height = 1, color = "black") +
  scale_fill_manual(values = f07_substrate_color) +
  geom_text(aes(label = ifelse(after_stat(x) != 0, after_stat(x), ""), 
                group = Group, 
                y = Name,fontface = "bold"), 
            stat = 'summary', 
            fun = sum, 
            hjust = -0.4, 
            vjust = 0.5, 
            size = 1.8) +
  geom_label(aes(label= genomes_containing_EC_per_group, 
                 y = Name, 
                 x = 166),
             hjust = 1, 
             vjust = 0.5, 
             size = 1.8) +
  geom_label(aes(label = paste0(Activity), 
                 x = 172,
                 group = Group),
           fill = "white", 
           hjust=0.5, 
           vjust=0.5, 
           colour = "black", 
           size = 1.8, 
           angle = 0) +
  geom_label(aes(label = paste0("MAGs in ",Group, " = ",Total_genomes_group), 
                 x = 75, 
                 y = "guluronate-specific alginate lyase"), 
             hjust = 1, 
             vjust = 1, 
             size = 1.8, 
             show.legend = FALSE) +
  geom_label(aes(label = paste0("CAZymes in ",Group, " = ",CAZyme_per_group), 
                 x = 75, 
                 y = "mannuronate-specific alginate lyase"), 
             hjust = 1, 
             vjust = 1, 
             size = 1.8, 
             show.legend = FALSE)
  
f10_barplot_group_enzyme
ggsave(filename = "01d_capacity_groups.pdf", f10_barplot_group_enzyme, width = 100, height = 400 , units = "mm",device = cairo_pdf)

#################################################################################
# gather plots
f11_log_relative_abundance_pool <- e07_log_relative_abundance + f10_barplot_group_enzyme + plot_layout(width = c(100,100))
f11_log_relative_abundance_pool
ggsave(filename = "01e_group_logfc_relabund_pool_genes.pdf", f11_log_relative_abundance_pool, width = 200, height = 280 , units = "mm",device = cairo_pdf)

#################################################################################
# Table summary
coabund_vec<-a03_coabund_groups$genomeID

g01_summary_table_abund<-b01_relative_abundance%>%
  arrange(factor(genomeID, levels= coabund_vec),Fish)%>%
  pivot_wider(id_cols = c("Group","Phylum","Class","Genus","genomeID"),names_from =Sample,values_from =Relative_abundance,values_fill = 0)%>%
  select(Group,Phylum,Class,Genus,genomeID,ends_with("_IV"),ends_with("_V"))

g02_summary_table_LFC<-e05a_sum_mean_group%>%
  rename(Sum_Group_Section=sum_Group_Sec,Log10_FC_Sum_Group_Section=log10fold_change_sum)%>%
  pivot_wider(id_cols = c("Group","Log10_FC_Sum_Group_Section"),names_from =Section,values_from =Sum_Group_Section,values_fill = 0)%>%
  left_join(g01_summary_table_abund,by = "Group")%>%
  rename(Sum_Group_IV=IV,Sum_Group_V=V)%>%
  arrange(factor(genomeID, levels= coabund_vec))

g03_CAZymes<-a04_refined_cazymes%>%
  select(genomeID,query_name,EC_refined)%>%
  left_join(a06_CGC_all,by = c("genomeID","query_name"))%>%
  select(-c("CGC_content"))%>%
  mutate(across(starts_with("CGC"), ~ifelse(is.na(.), "out", "in")))%>%
  rename(CGC_localization=CGC_ID)%>%
  distinct()%>%
  group_by(genomeID,CGC_localization,EC_refined)%>%
  mutate(EC_per_genome=n())%>%
  select(genomeID,EC_refined,CGC_localization,EC_per_genome)%>%
  distinct()%>%
  right_join(a02_taxonomy,by = "genomeID")%>%
  select(genomeID,EC_refined,CGC_localization,EC_per_genome)%>%
  pivot_wider(id_cols = c("EC_refined","CGC_localization"),names_from =genomeID,values_from =EC_per_genome,values_fill = 0)%>%
  na.omit()%>%
  pivot_longer(cols = -c("EC_refined","CGC_localization"),names_to = "genomeID",values_to = "EC_per_genome")%>%
  pivot_wider(id_cols = c("genomeID","CGC_localization"),names_from =EC_refined,values_from =EC_per_genome,values_fill = 0)%>%
  na.omit()%>%
  pivot_longer(cols = -c("genomeID","CGC_localization"),names_to = "EC_refined",values_to = "EC_per_genome")%>%
  left_join(a03_coabund_groups,by="genomeID")%>%
  select(Group,genomeID,EC_refined,CGC_localization,EC_per_genome)%>%
  left_join(f01_Genomes_containing_in_Group,by = c("Group","EC_refined"))%>%
  left_join(f02_Function,by = "EC_refined", relationship = "many-to-many")%>%
  mutate(`Function` = paste0(as.character(Substrate), " ", as.character(Activity), ": ", as.character(Function)))%>%
  arrange(CGC_localization,factor(EC_refined, levels= EC_vec))%>%
  ungroup()%>%
  select(Group,Total_genomes_group,genomeID,CGC_localization,EC_per_genome,Function)%>%
  rename(number_of_MAGs_in_group=Total_genomes_group)%>%
  pivot_wider(id_cols = c("Group","number_of_MAGs_in_group","genomeID"),names_from =c("Function","CGC_localization"),values_from = EC_per_genome,values_fill = 0)%>%
  arrange(factor(genomeID, levels= coabund_vec))%>%
  inner_join(g02_summary_table_LFC, by = c("Group","genomeID"))%>%
  select(Group,number_of_MAGs_in_group,Sum_Group_IV,Sum_Group_V,Log10_FC_Sum_Group_Section,Phylum,Class,Genus,genomeID,starts_with("G"),
         starts_with("alginate HMW"),starts_with("laminarin HMW"),starts_with("FCSP HMW"),starts_with("carrageenan HMW"),
         starts_with("alginate LMW"),starts_with("laminarin LMW"),starts_with("FCSP LMW"),starts_with("carrageenan LMW"))
         
library(readr)
write_excel_csv(g03_CAZymes, '01f_Summary_function.csv')
