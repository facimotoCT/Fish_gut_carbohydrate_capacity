# Explore CGC classified and relative abundance
setwd("Manuscript_1/Figure 3 4/")
library(dplyr)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
##############################################################
##############################################################
##############################################################
# CGCs counts exploration
a01_CGCs<-read.csv("01a_CGC_content_long.csv")%>%
  select(-c("Class"))

a02_CGC_classified<-read.csv("01e_CGC_ID_substrate_specificity.csv")%>%
  select(CGC_ID,Substrate)%>%
  distinct()%>%
  left_join(a01_CGCs,by = "CGC_ID")%>%
  select(genomeID,CGC_ID,Substrate)%>%
  distinct()

# taxonomy
a03_taxonomy<-read.csv("../Annotation_table.csv")%>%
  select(Phylum,Class,Order,Family,Genus,genomeID,MAG)%>%
  distinct()%>%
  mutate(Family = coalesce(Family, Order))%>%
  mutate(Genus = coalesce(Genus, Family))%>%
  mutate(Phylum = gsub("p__", "", Phylum))%>%
  mutate(Class = gsub("c__", "", Class))%>%
  mutate(Order = gsub("o__", "", Order))%>%
  mutate(Family = gsub("f__", "", Family))%>%
  mutate(Genus = gsub("g__", "", Genus))

# counts abundance
a04_relative_abundance <- read.delim("../local_files/020205_normalised_summary_count_table_genomeSummary.tsv")%>%
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
  pivot_longer(cols =-c("MAG"), names_to = "Sample",values_to = "counts")%>%
  group_by(Sample)%>%
  mutate(Total_sample=sum(counts))%>%
  mutate(Relative_abundance=(counts/Total_sample)*100)%>%
  separate(Sample, into = c("Fish", "Section"), sep = "_",remove = FALSE)%>%
  ungroup()%>%
  select(MAG,Fish,Section,Relative_abundance)

##############################################################
##############################################################
##############################################################
##############################################################
# counts abundance only the ones containing CGCs of interest
b01_relative_abundance_CGC<-a03_taxonomy%>%
  select(Phylum,Class,Genus,genomeID,MAG)%>%
  distinct()%>%
  left_join(a04_relative_abundance,by="MAG")%>%
  filter(genomeID %in% unique(a02_CGC_classified$genomeID))

##############################################################
##############################################################
##############################################################
##############################################################
# Mean Genus
c01_genus_stats<-a03_taxonomy%>%
  select(Phylum,Class,Genus,genomeID,MAG)%>%
  distinct()%>%
  left_join(a04_relative_abundance,by="MAG")%>%
  select(Phylum,Class,Genus,genomeID,Fish,Section,Relative_abundance)%>%
  group_by(Genus,Fish,Section)%>%
  mutate(Genus_relative_abundance=sum(Relative_abundance))%>%
  select(Genus,Section,Fish,Genus_relative_abundance)%>%
    distinct()%>%
  group_by(Genus,Section)%>%
  mutate(Mean_Genus=mean(Genus_relative_abundance),
         SD_Genus=sd(Genus_relative_abundance))%>%
  group_by(Genus)%>%
  mutate(KW_Genus=kruskal.test(Genus_relative_abundance~Section)$stat, p_KW_Genus = kruskal.test(Genus_relative_abundance~Section)$p.value)

# CGCs distribution total, degradative and predicted
c02_CGC_total_tax<-a03_taxonomy%>%
  select(Class,Genus,genomeID)%>%
  distinct()%>%
  left_join(a01_CGCs,by = "genomeID")%>%
  select(Class,Genus,genomeID,CGC_ID)%>%
  distinct()%>%
  group_by(genomeID)%>%
  mutate(CGC_per_genome=n())%>%
  ungroup()%>%
  select(Class,Genus,genomeID,CGC_per_genome)%>%
  distinct()%>%
  # Total CGC count
  ungroup()%>%
  mutate(Total_CGC=sum(CGC_per_genome))%>%
  select(Class,Genus,genomeID,Total_CGC,CGC_per_genome)%>%
  # Total CGC per Class
  group_by(Class)%>%
  mutate(CGC_per_Class=sum(CGC_per_genome))%>%
  ungroup()%>%
  select(Class,Genus,genomeID,Total_CGC,CGC_per_Class,CGC_per_genome)%>%
  distinct()%>%
  # Mean SD CGC per Class
  group_by(Class)%>%
  mutate(Mean_CGC_per_Class=mean(CGC_per_genome),
         SD_CGC_per_Class=sd(CGC_per_genome))%>%
  # Total CGC per Genus
  group_by(Genus)%>%
  mutate(CGC_per_Genus=sum(CGC_per_genome))%>%
  ungroup()%>%
  select(Class,Genus,genomeID,Total_CGC,Mean_CGC_per_Class,SD_CGC_per_Class,CGC_per_Class,CGC_per_Genus,CGC_per_genome)%>%
  distinct()%>%
  # Mean SD CGC per Genus
  group_by(Genus)%>%
  mutate(Mean_CGC_per_Genus=mean(CGC_per_genome),
         SD_CGC_per_Genus=sd(CGC_per_genome))

# CGCs degradative action (containing more than 2 degradative CAZymes GH, PL or CE) and average per phylum
c03_CGC_deg_tax<-a03_taxonomy%>%
  select(Class,Genus,genomeID)%>%
  distinct()%>%
  left_join(a01_CGCs,by = "genomeID")%>%
  
  # Select only degradative and containing 2 CAZymes
  select(Class,Genus,genomeID,query_name,CGC_ID,Fam_sub_CAZy)%>%
  filter(grepl("GH|PL|CE",Fam_sub_CAZy))%>%
  group_by(CGC_ID)%>%
  mutate(Deg_count=n())%>%
  filter(Deg_count > 1)%>%
  ungroup()%>%
  
  select(Class,Genus,genomeID,CGC_ID)%>%
  distinct()%>%
  group_by(genomeID)%>%
  mutate(deg_CGC_per_genome=n())%>%
  ungroup()%>%
  select(Class,Genus,genomeID,deg_CGC_per_genome)%>%
  distinct()%>%
  # Total deg CGC count
  ungroup()%>%
  mutate(Total_deg_CGC=sum(deg_CGC_per_genome))%>%
  select(Class,Genus,genomeID,Total_deg_CGC,deg_CGC_per_genome)%>%
  # Total CGC per Class
  group_by(Class)%>%
  mutate(deg_CGC_per_Class=sum(deg_CGC_per_genome))%>%
  ungroup()%>%
  select(Class,Genus,genomeID,Total_deg_CGC,deg_CGC_per_Class,deg_CGC_per_genome)%>%
  distinct()%>%
  # Mean SD CGC per Class
  group_by(Class)%>%
  mutate(Mean_deg_CGC_per_Class=mean(deg_CGC_per_genome),
         SD_deg_CGC_per_Class=sd(deg_CGC_per_genome))%>%
  # Total CGC per Genus
  group_by(Genus)%>%
  mutate(deg_CGC_per_Genus=sum(deg_CGC_per_genome))%>%
  ungroup()%>%
  select(Class,Genus,genomeID,Total_deg_CGC,deg_CGC_per_Class,Mean_deg_CGC_per_Class,SD_deg_CGC_per_Class,deg_CGC_per_Genus,deg_CGC_per_genome)%>%
  distinct()%>%
  # Mean SD CGC per Genus
  group_by(Genus)%>%
  mutate(Mean_deg_CGC_per_Genus=mean(deg_CGC_per_genome),
         SD_deg_CGC_per_Genus=sd(deg_CGC_per_genome))

# CGCs degradative predicted
c04_CGC_pred_tax<-a03_taxonomy%>%
  select(Class,Genus,genomeID)%>%
  distinct()%>%
  left_join(a02_CGC_classified,by = "genomeID")%>%

  # count substrate per genome
  group_by(genomeID,Substrate)%>%
  mutate(Substrate_CGC_genome=n())%>%
  select(Class,Genus,genomeID,Substrate,Substrate_CGC_genome)%>%
  distinct()%>%
  pivot_wider(id_cols=c("Class","Genus","genomeID"),names_from = Substrate,values_from = Substrate_CGC_genome,values_fill = 0)%>%
  select(-c("NA"))%>%
  
  # predicted CGC per genome
  mutate(pred_CGC_per_genome = (laminarin+alginate+FCSP+carrageenan+starch+agarose))%>%
  
  # Total deg pred CGC count
  ungroup()%>%
  mutate(Total_pred_CGC=sum(pred_CGC_per_genome))%>%
  select(Class,Genus,genomeID,laminarin,alginate,FCSP,carrageenan,starch,agarose,Total_pred_CGC,pred_CGC_per_genome)%>%
  # Total CGC per Class
  group_by(Class)%>%
  mutate(pred_CGC_per_Class=sum(pred_CGC_per_genome))%>%
  ungroup()%>%
  select(Class,Genus,genomeID,laminarin,alginate,FCSP,carrageenan,starch,agarose,Total_pred_CGC,pred_CGC_per_Class,pred_CGC_per_genome)%>%
  distinct()%>%
  # Mean SD CGC per Class
  group_by(Class)%>%
  mutate(Mean_pred_CGC_per_Class=mean(pred_CGC_per_genome),
         SD_pred_CGC_per_Class=sd(pred_CGC_per_genome))%>%
  # Total CGC per Genus
  group_by(Genus)%>%
  mutate(pred_CGC_per_Genus=sum(pred_CGC_per_genome))%>%
  ungroup()%>%
  select(Class,Genus,genomeID,laminarin,alginate,FCSP,carrageenan,starch,agarose,Total_pred_CGC,pred_CGC_per_Class,Mean_pred_CGC_per_Class,SD_pred_CGC_per_Class,pred_CGC_per_Genus,pred_CGC_per_genome)%>%
  distinct()%>%
  # Mean SD CGC per Genus
  group_by(Genus)%>%
  mutate(Mean_pred_CGC_per_Genus=mean(pred_CGC_per_genome),
         SD_pred_CGC_per_Genus=sd(pred_CGC_per_genome))

# Collapse tables
c05_summary_Genus<-c01_genus_stats%>%
  select(Genus,Section,Mean_Genus,SD_Genus,KW_Genus,p_KW_Genus)%>%
  distinct()%>%
  pivot_wider(id_cols=c("Genus","KW_Genus","p_KW_Genus"),names_from = Section,values_from = c("Mean_Genus","SD_Genus"))%>%
  select(Genus,Mean_Genus_IV,SD_Genus_IV,Mean_Genus_V,SD_Genus_V,KW_Genus,p_KW_Genus)%>%
  left_join(c01_genus_stats,by=c("Genus","KW_Genus","p_KW_Genus"))%>%
  select(-c("Mean_Genus","SD_Genus"))%>%
  arrange(Fish,Section)%>%
  mutate(Fish = gsub("(G[0-9])", "Genus_\\1", Fish))%>%
  pivot_wider(id_cols=c("Genus","Mean_Genus_IV","SD_Genus_IV","Mean_Genus_V","SD_Genus_V","KW_Genus","p_KW_Genus"),names_from = c("Fish","Section") ,values_from = c("Genus_relative_abundance"))%>%
  select(Genus,starts_with("Genus_") & ends_with("IV"),starts_with("Genus_") & ends_with("_V"),Mean_Genus_IV,SD_Genus_IV,Mean_Genus_V,SD_Genus_V,KW_Genus,p_KW_Genus)

c06_Genus_genome_sum<-a03_taxonomy%>%
  left_join(a04_relative_abundance,by="MAG")%>%
  select(Phylum,Class,Genus,genomeID,Fish,Section,Relative_abundance)%>%
  arrange(Fish,Section)%>%
  pivot_wider(id_cols=c("Phylum","Class","Genus","genomeID"),names_from = c("Fish","Section") ,values_from = c("Relative_abundance"))%>%
  left_join(c05_summary_Genus,by="Genus")%>%
  left_join(c02_CGC_total_tax,by=c("Class","Genus","genomeID"))%>%
  left_join(c03_CGC_deg_tax,by=c("Class","Genus","genomeID"))%>%
  left_join(c04_CGC_pred_tax,by=c("Class","Genus","genomeID"))%>%
  select(Phylum,Class,Genus,
       Genus_G117_IV,Genus_G121_IV,Genus_G124_IV,Genus_G125_IV,
       Genus_G117_V,Genus_G121_V,Genus_G124_V,Genus_G125_V,
       Mean_Genus_IV,SD_Genus_IV,Mean_Genus_V,SD_Genus_V,KW_Genus,p_KW_Genus,
       genomeID,
       G117_IV,G121_IV,G124_IV,G125_IV,
       G117_V,G121_V,G124_V,G125_V,
       Total_CGC,CGC_per_Class,Mean_CGC_per_Class,SD_CGC_per_Class,CGC_per_Genus,Mean_CGC_per_Genus,SD_CGC_per_Genus,CGC_per_genome,
       Total_deg_CGC,deg_CGC_per_Class,Mean_deg_CGC_per_Class,SD_deg_CGC_per_Class,deg_CGC_per_Genus,Mean_deg_CGC_per_Genus,SD_deg_CGC_per_Genus,deg_CGC_per_genome,
       Total_pred_CGC,pred_CGC_per_Class,Mean_pred_CGC_per_Class,SD_pred_CGC_per_Class,pred_CGC_per_Genus,Mean_pred_CGC_per_Genus,SD_pred_CGC_per_Genus,pred_CGC_per_genome,
       laminarin,alginate,FCSP,carrageenan,starch,agarose)%>%
  mutate(KW_pred_CGC_Genus=kruskal.test(pred_CGC_per_genome~Genus)$stat, p_KW_pred_CGC_Genus = kruskal.test(pred_CGC_per_genome~Genus)$p.value)

##############################################################
##############################################################
##############################################################
##############################################################
# Plots
d01_genus_CGCs_substrate<-c06_Genus_genome_sum%>%
  select(Class,Genus,genomeID,alginate,laminarin,FCSP,carrageenan,starch,agarose)%>%
  pivot_longer(cols = -c("Class","Genus","genomeID"),names_to = "Substrate",values_to="Substrate_genome")%>%
  group_by(Genus,Substrate)%>%
  mutate(Sum_CGC_substrate_genus=sum(Substrate_genome))%>%
  select(Class,Genus,Substrate,Sum_CGC_substrate_genus)%>%
  distinct()%>%
  group_by(Genus)%>%
  mutate(Genus_CGCs_total=sum(Sum_CGC_substrate_genus))%>%
  arrange(desc(Genus_CGCs_total))

substrate_vec<-c("laminarin", "alginate","FCSP", "carrageenan", "starch", "agarose")

# substrate colors
library(rcartocolor)
# Generate color palette
rcartocolor::display_carto_all(colorblind_friendly = TRUE)
substrate_color_palette <- rcartocolor::carto_pal(n=6,"Earth")
names(substrate_color_palette) <- substrate_vec
substrate_color_palette

class_colors<-c("Bacteroidia"="#B2DF8A",
                "Bacilli"="#FDBF6F",
                "Clostridia"="#FF7F00",
                "Spirochaetia"="#FB9A99",
                "Desulfovibrionia"="#FFFF99",
                "Verrucomicrobiae"="#B15928",
                "Alphaproteobacteria"="#CAB2D6",
                "Gammaproteobacteria"="#6A3D9A",
                "Vampirovibrionia"="#1F78B4")

# only totals
d02_totals<-d01_genus_CGCs_substrate%>%
  select(Class,Genus,Substrate,Sum_CGC_substrate_genus,Genus_CGCs_total)%>%
  filter(Genus_CGCs_total > 0)

### order x axis
d03_rich_genus<-unique(d02_totals$Genus)

# plot distribution of classified CGCs
library(ggplot2)
library(ggnewscale)
d04_CGCs_classified<-ggplot(d02_totals, aes(x = factor(Genus, levels = as.character(d03_rich_genus)), y = Sum_CGC_substrate_genus, fill = Substrate)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.1) +
  scale_fill_manual(values = substrate_color_palette, limits = substrate_vec) +
  scale_y_continuous(name = "Predicted CGCs", limits = c(-2, 80), breaks = seq(0, 76, 10), expand = c(0, 1)) +
  scale_x_discrete(name = "Genus") +
  theme_bw() +
  labs(title = "a) Distribution of predicted CGCs") +
  theme(
    plot.title = element_text(size = 7, face = "bold"),
    axis.title.y = element_text(angle = 90, size = 6, face = "bold"),
    axis.text.y = element_text(angle = 0, size = 5),
    axis.title.x = element_text(angle = 0, size = 6, face = "bold"),
    axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.5),
    
    axis.ticks = element_line(size = 0.25),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6, face = "bold"),
    legend.position="bottom", 
    legend.box = "horizontal", 
    legend.justification = "left",
    legend.key.size = unit(6, 'mm')
  )  +
  guides(fill = guide_legend(nrow = 2)) +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = factor(Genus, levels = as.character(d03_rich_genus)), y = -1.5, fill = Class), color = "black") +
  scale_fill_manual(values = class_colors) +
  geom_text(aes(label = after_stat(y), group = Genus), stat = 'summary', fun = sum, vjust = -1, size = 1.8) +
  guides(fill = guide_legend(nrow = 2)) +
  theme(
    plot.title = element_text(size = 7, face = "bold"),
    axis.title.y = element_text(angle = 90, size = 6, face = "bold"),
    axis.text.y = element_text(angle = 0, size = 5),
    axis.title.x = element_text(angle = 0, size = 6, face = "bold"),
    axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.5),
    
    axis.ticks = element_line(size = 0.25),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6, face = "bold"),
    legend.position="bottom", 
    legend.box = "horizontal", 
    legend.justification = "left",
    legend.key.size = unit(3, 'mm')
  ) 

d04_CGCs_classified
ggsave(filename = "02a_CGC_substrate.pdf", d04_CGCs_classified, width = 150, height = 100 , units = "mm")

##############################################################
# Plot sum abundance, mean and sd lines
library(scales)
library(ggpubr)
e05_subset_genus_stats<-c01_genus_stats%>%
  filter(Genus %in% d03_rich_genus)

e06_mean_abund<-ggplot(e05_subset_genus_stats, aes(x = factor(Genus, levels = as.character(d03_rich_genus)), y = Mean_Genus, fill = Section)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             aes(y = Genus_relative_abundance), pch = 21,
             size = 0.5, stroke = 0.1,
             show.legend = FALSE, alpha = 0.5) + 
  geom_point(shape = 21, color = "black", size = 1,stroke=0.1, position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_errorbar(aes(ymin = Mean_Genus - SD_Genus, ymax = Mean_Genus + SD_Genus), width = 0.2, size = 0.1, position = position_dodge(width = 0.8,preserve = "total"),show.legend = FALSE) +
  scale_y_continuous(name = "Relative abundance", limits = c(-3,90), breaks = seq(0, 90, 20)) +
  scale_x_discrete(name="Genus") +
  scale_fill_manual(values = c("IV" = "#f061eb", "V" = "#205233"), name = "Section") +
  theme_bw() +
  # stat_compare_means(method="wilcox.test", p.adjust = "BH",
  #                    comparisons = c("IV","V"), group = factor(Genus),
  #                    label.y = 80, 
  #                    #                     hide.ns = TRUE, 
  #                    label = "p.signif",
  #                    bracket.size = 0.2,size=1.8) 
  labs(title = "b)") +
    theme(
      plot.title = element_text(size = 7, face = "bold"),
    axis.ticks = element_line(size = 0.25),
    axis.title.y = element_text(angle = 90, size = 6, face = "bold"),
    axis.text.y = element_text(angle = 0, size = 5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5),
    strip.text.y.right = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5),
    panel.grid.major = element_line(size = 0.1),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(size = 0.5),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6, face = "bold"),
    legend.position="bottom", legend.box = "horizontal", legend.justification = "right"
    ) 
e06_mean_abund
ggsave(filename = "02b_abundance.pdf", e06_mean_abund, width = 100, height = 80 , units = "mm")

# Mean Alistipes vs other Bacteroidia vs Clostridia vs other
e07_Alis_CGC_mean_other_genus<-c06_Genus_genome_sum%>%
  select(Class,Genus,genomeID,pred_CGC_per_genome)%>%
  distinct()%>%
  group_by(Class)%>%
  mutate(Sum_CGC_Class=sum(pred_CGC_per_genome))%>%
  filter(Sum_CGC_Class > 0)%>%
  mutate(ngenomes_family=n())%>%
  group_by(Genus)%>%
  mutate(ngenomes_Genus=n())%>%
  mutate(Taxa = case_when(
    ngenomes_Genus > 2 ~ Genus,
    TRUE ~ NA_character_
  )) %>%
  na.omit()%>%
  filter(!grepl("o__Chris",Taxa))%>%
  group_by(Taxa)%>%
  mutate(mean_Taxa=mean(pred_CGC_per_genome),
         SD_Taxa=mean(pred_CGC_per_genome)
         
         )

e08_mean_pred<-ggplot(e07_Alis_CGC_mean_other_genus, aes(fill = Class, x = factor(Taxa,levels = c("Alistipes","f__Rikenellaceae","f__Lachnospiraceae","Mucinivorans")), y = mean_Taxa)) +
  geom_bar(aes(fill = Class, y = mean_Taxa),
           stat = "identity", width = 1, size = 0.2, color = "black", show.legend = F, position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin = mean_Taxa - 0 , ymax = mean_Taxa + SD_Taxa), width = 0.3, size = 0.2, color = "black", position = position_dodge(width = 1)) +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(name="Mean of Predicted CGCs",position = "left",)+
  theme_bw() +
  labs(title = "c)") +
  theme(
    plot.title = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(angle= 90, size = 5,  hjust = 1, vjust = 0.5),
    axis.ticks.x = element_line(size = 0.25),
    axis.title.x = element_blank(),
    
    axis.text.y = element_text(size = 5,  hjust = 0.5, vjust = 0.5),
    axis.ticks.y = element_line(size = 0.25),
    axis.title.y = element_text(angle = 90, size = 6, face = "bold"),
    
#    plot.margin = margin(1,1,1,1, unit = "mm"),
    panel.spacing = unit(1, "mm"),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.5),
    plot.background = element_blank()
  )
e08_mean_pred
ggsave(filename = "02c_mean_CGC_predicted.pdf", e08_mean_pred, width = 30, height = 80 , units = "mm")

library(patchwork)
e09_inset_abundance_CGCs_origin<-e06_mean_abund + e08_mean_pred + plot_layout(widths = c(90,10))
e09_inset_abundance_CGCs_origin
e10_inset_abundance_CGCs_origin<-d04_CGCs_classified + inset_element(e09_inset_abundance_CGCs_origin, 0.07, 0.1, 0.99, 0.93)
e10_inset_abundance_CGCs_origin

ggsave(filename = "02d_CGC_substrate_genus_abundance.pdf", e10_inset_abundance_CGCs_origin, width = 150, height = 120 , units = "mm")

# Supplementary table
write.csv(c06_Genus_genome_sum, '02e_table_summary_CGC_abundance.csv', row.names = FALSE)

