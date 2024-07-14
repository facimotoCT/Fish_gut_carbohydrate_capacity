# Correlation and functional capacity
setwd("Manuscript_1/Figure 6/")
library(broom)
library(dplyr)
library(tidyr)

##################################################################################
# ABUNDANCE AND METADATA
# COUNTS PER GENOMES
a01_counts_norm <- read.delim("../local_files/020205_normalised_summary_count_table_genomeSummary.tsv")%>%
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
    G125_IV=S5_pileup_counts_lengthNorm_libNorm)
  
# TAXONOMY
a02_taxonomy<-read.csv("../Annotation_table.csv")%>%
  select(Phylum,Class,Order,Family,Genus,genomeID,MAG)%>%
  distinct()%>%
  mutate(Family = coalesce(Family, Order))%>%
  mutate(Genus = coalesce(Genus, Family))%>%
  mutate(Phylum = gsub("p__", "", Phylum))%>%
  mutate(Class = gsub("c__", "", Class))%>%
  mutate(Order = gsub("o__", "", Order))%>%
  mutate(Family = gsub("f__", "", Family))%>%
  mutate(Genus = gsub("g__", "", Genus))

tax_color<-c("Bacteroidia"="#B2DF8A",
             "Bacilli"="#FDBF6F",
             "Clostridia"="#FF7F00",
             "Spirochaetia"="#FB9A99",
             "Desulfovibrionia"="#FFFF99",
             "Verrucomicrobiae"="#B15928",
             "Alphaproteobacteria"="#CAB2D6",
             "Gammaproteobacteria"="#6A3D9A",
             "Vampirovibrionia"="#1F78B4")

#################################################################################
#################################################################################
#################################################################################
##################################################################################
# append taxonomy to counts
# Relative abundance
b01_relative_abundance<-a02_taxonomy%>%
  left_join(a01_counts_norm,by = "MAG")%>%
  mutate(across(starts_with("G1"),~ ((.) / sum(.)) * 100))%>%
  pivot_longer(cols=-c("Phylum","Class","Order","Family","Genus","genomeID","MAG"),names_to = "Sample",values_to = "Relative_abundance")%>%
  separate(Sample, into = c("Fish", "Section"), sep = "_",remove = FALSE)

# Normalized counts abundance
b02_norm_counts<-a02_taxonomy%>%
  left_join(a01_counts_norm,by = "MAG")%>%
  pivot_longer(cols=-c("Phylum","Class","Order","Family","Genus","genomeID","MAG"),names_to = "Sample",values_to = "Normalized_counts")%>%
  separate(Sample, into = c("Fish", "Section"), sep = "_",remove = FALSE)

#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Correlation of abundances
library(tidyverse)
library(Hmisc)
library(corrplot)
# genomeID
c01_genomeID_coabundance<-b01_relative_abundance%>%
  select(genomeID,Sample,Relative_abundance)%>%
  pivot_wider(id_cols = Sample,names_from = genomeID,values_from = Relative_abundance, values_fill = 0)%>%
  column_to_rownames("Sample")

c02_genomeID_res<-rcorr(as.matrix(c01_genomeID_coabundance[,1:68]), type=c("pearson"))

# Create the correlation plot
library(RColorBrewer)
pdf(file="01a_corr_hclust_complete.pdf",width = 8,height = 8)
corrplot(c02_genomeID_res$r,
         method = "color",
         type = "full",
         order = "hclust",
         sig.level = 0.01,
         insig = "blank",
         tl.cex = 0.4, tl.col = "black",
         cl.cex = 0.6,
         tl.pos = "l",
         cl.pos = "b",
         cl.ratio = 0.1,
         col=brewer.pal(n=10, name="RdYlBu"),
         addgrid.col = "gray")
dev.off()

# reorder highly correlated groups and their log fold changes
c03_order_groups<-c(
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

c05_matrix_R<-as.matrix(c02_genomeID_res$r)
c06_matrix_ordered_R<-c05_matrix_R[c03_order_groups,c03_order_groups]
pdf(file="01b_corr_hclust_complete_order_groups.pdf",width = 8,height = 8)
corrplot(c06_matrix_ordered_R,
         method = "color",
         type = "lower",
         order = "original",
         insig = "blank",
         tl.cex = 0.4, tl.col = "black",
         cl.cex = 0.6,
         tl.pos = "l",
         cl.pos = "b",
         cl.ratio = 0.1,
         col=brewer.pal(n=10, name="RdYlBu"),
         addgrid.col = "gray")
dev.off()

c05_matrix_P<-as.matrix(c02_genomeID_res$P)
c06_matrix_ordered_P<-c05_matrix_P[c03_order_groups,c03_order_groups]
c06_matrix_ordered_P[is.na(c06_matrix_ordered_P)] <- 1
pdf(file="01c_corr_hclust_complete_order_groups_P_values.pdf",width = 8,height = 8)
corrplot(c06_matrix_ordered_R,
         method = "color",
         type = "upper",
         order = "original",
         p.mat=c06_matrix_ordered_P,
         sig.level = 0.05,
         insig = "blank",
         tl.cex = 0.4, tl.col = "black",
         cl.cex = 0.6,
         tl.pos = "l",
         cl.pos = "b",
         cl.ratio = 0.1,
         col=brewer.pal(n=10, name="RdYlBu"),
         addgrid.col = "gray")
dev.off()

#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Log fold change
d01_sum_genomeID<-b02_norm_counts%>%
  select(Phylum,Class,genomeID,Fish,Section,Normalized_counts)%>%
  group_by(genomeID,Section)%>%
  mutate(sum_genome_counts_sec=sum(Normalized_counts))%>%
  ungroup()%>%
  group_by(Section)%>%
  mutate(sum_counts_sec=sum(Normalized_counts))%>%
  mutate(Relative_abundance_sum_sec=((sum_genome_counts_sec/sum_counts_sec)*100))%>%
  select(Phylum,Class,genomeID,Section,Relative_abundance_sum_sec)%>%
  distinct()%>%
  pivot_wider(id_cols = c("Phylum","Class","genomeID"), names_from = Section,values_from = Relative_abundance_sum_sec)%>%
  mutate(log2fold_change=log2((V+1e-6)/(IV+1e-6)),
         log10fold_change=log10((V+1e-6)/(IV+1e-6)))

# Plot log10 fold change in geom_points solid and open circles
library(ggplot2)
library(ggnewscale)

d02_L10_FC <- ggplot(d01_sum_genomeID, aes(x = log10fold_change, y = factor(genomeID, levels = rev(c03_order_groups)), fill = Class)) +
  geom_point(alpha = 0.8, aes(color = ifelse(log10fold_change < 0, "IV", "V"),
                              size = ifelse(log10fold_change < 0, IV, V)), show.legend = TRUE) +
  scale_color_manual(values = c("IV" = "#f061eb", "V" = "#205233")) +
  scale_x_continuous(name = "log10FC", limits = c(-4, 6), breaks = seq(-3, 6, 1), expand = c(0, 0)) +
  scale_size(range = c(0, 5), breaks = c(0, 2.5, 5, 10, 20)) +
  guides(size = guide_legend(title = "Point Size")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 6),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.justification = "left"
  ) +
  labs(size = "Count sum", color = "Section") +
  geom_tile(aes(y = factor(genomeID, levels = rev(c03_order_groups)), x = -3.5), color = "black") +
  scale_fill_manual(values = tax_color) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") 
d02_L10_FC
ggsave(filename = "02d_fold_change.pdf", d02_L10_FC, width = 50, height = 200 , units = "mm")

#################################################################################
# enrichment of refined CAZymes
library(readxl)
d03_CAZymes_path<-read_xlsx("../Figure 3 4/03a_CAZymes_path_refined.xlsx") # CAZymes manually curated

d04_count_CAZyme<-d03_CAZymes_path%>%
  select(genomeID,Fam_refined,EC_refined,Substrate,Activity)%>%
  group_by(genomeID,Substrate,Activity)%>%
  mutate(count_refined_CAZyme=n())%>%
  select(genomeID,Substrate,Activity,count_refined_CAZyme)%>%
  distinct()%>%
  right_join(a02_taxonomy,by = "genomeID")%>%
  pivot_wider(id_cols=c("Substrate","Activity"),names_from = genomeID,values_from = count_refined_CAZyme,values_fill = 0)%>%
  na.omit()%>%
  pivot_longer(cols=-c("Substrate","Activity"),names_to = "genomeID",values_to = "count_refined_CAZyme")

d05_vec_Substrates<-c("alginate", "laminarin","FCSP","carrageenan")
d04_substrate_color<-c("laminarin"="#A16928",
                       "alginate"="#C29B64",
                       "FCSP"="#E0CFA2",
                       "carrageenan"="#CBD5BC")


d06_CAZyme_path<-ggplot(d04_count_CAZyme, aes(x = factor(Activity), y = factor(genomeID, levels = rev(c03_order_groups)), fill = Substrate)) +
  geom_point(alpha = 0.8, aes(color = Substrate, size = count_refined_CAZyme)) +
  scale_color_manual(values = d04_substrate_color) +
  scale_size(range = c(0, 10),breaks = c(0,1,2,3,4,5)) +
  facet_grid(cols = vars(factor(Substrate, levels = d05_vec_Substrates)), 
             scales = "free_x", space = "free", labeller = label_wrap_gen(multi_line = TRUE)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    legend.position="bottom", legend.box = "horizontal", legend.justification = "right",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6)
  )  
d06_CAZyme_path
library(patchwork)
d07_logfc_CAZy<-d02_L10_FC + d06_CAZyme_path + plot_layout(widths = c(1,2 ))
d07_logfc_CAZy
ggsave(filename = "01e_fold_change_cazyme_path.pdf", d07_logfc_CAZy, width = 100, height = 200,units = "mm")
