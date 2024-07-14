# CAZyme, sulfatase density
setwd("C:/Users/cfac668/OneDrive - The University of Auckland/Chapter_1/Revisiting/Review_analysis_V2_CF13_KH9_KDC3_ISME_submission/Figure 1/")

library(dplyr)
library(readxl)
library(tidyr)

##################################################################################
# Taxonomy
tax<-read.csv("01_taxonomy_genomeID_refs_tree.csv")

# CAZy density
a01_CAZy_density<-read.csv("../Annotation_table.csv")%>%
  select(Class,MAG,query_name,Fam_sub_CAZy,Genome_length_Mbp)%>%
  na.omit()%>%
  group_by(MAG)%>%
  mutate(CAZy_count=n())%>%
  ungroup()%>%
  select(Class,MAG,CAZy_count,Genome_length_Mbp)%>%
  distinct()%>%
  mutate(CAZy_density = CAZy_count/Genome_length_Mbp)%>%
  mutate(Class = sub(".*__", "", Class))%>%
  filter(grepl("Bacteroidia|Clostridia|Bacilli|Gammaproteobacteria",Class))%>%
  arrange(desc(CAZy_density))

# Plotting the boxplot with significances
library(ggplot2)
library(ggpubr)
# a02_cols<-c("Bacteroidota"="#1AFF1A",
#             "Proteobacteria"="#EEDD88",
#             "Firmicutes"="#FFC20A",
#             "Firmicutes_A"="#EE8866")

a02_cols<-c("Bacteroidia"="#B2DF8A",
            "Bacilli"="#FDBF6F",
            "Clostridia"="#FF7F00",
            "Spirochaetia"="#FB9A99",
            "Desulfovibrionia"="#FFFF99",
            "Verrucomicrobiae"="#B15928",
            "Alphaproteobacteria"="#CAB2D6",
            "Gammaproteobacteria"="#6A3D9A",
            "Vampirovibrionia"="#1F78B4")

a03_comparison<-list(
  c("Bacteroidia", "Clostridia"),
  c("Bacteroidia", "Bacilli"),
  c("Bacteroidia", "Gammaproteobacteria"),
  c("Clostridia", "Bacilli"),
  c("Clostridia", "Gammaproteobacteria"),
  c("Bacilli", "Gammaproteobacteria"))

a04_order_Class<-unique(a01_CAZy_density$Class)

a05_cazyme_mbp<-ggplot(a01_CAZy_density, aes(fill = Class, y = CAZy_density, x = factor(Class,a04_order_Class))) +
  geom_boxplot(position = position_dodge(width = 1), show.legend = FALSE, size = 0.3,
               aes(fill = Class, x = factor(Class,a04_order_Class))) +
  geom_point(position = position_jitterdodge(jitter.width = 1, dodge.width = 1),
             aes(fill = factor(Class,a04_order_Class)), pch = 21,
             size = 1, stroke = 0.2,
             show.legend = FALSE, alpha = 0.5) +
  scale_fill_manual(values = a02_cols) +
  theme_bw() +
  scale_y_continuous(name = "CAZyme / Mbp", limits = c(0, 85), breaks = seq(0, 85, 20)) +
  labs(title = "b) CAZyme density", x = "Class", y = "CAZyme / Mbp", face = "bold") +
  theme(panel.grid.major = element_line(size = 0.25),
        plot.title = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 7, face="bold"),
        axis.ticks.y = element_line(size = 0.25),
        axis.ticks.x = element_line(size = 0.25),
        rect = element_rect(size = 0.1),
        line = element_line(size = 0.1),
        plot.margin = margin(1,1,1,1, unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.grid = element_line(size = 0.1),
        panel.border = element_rect(size = 0.5)) +
  stat_compare_means(method="wilcox.test", comparisons = a03_comparison, label.y = c(65, 70, 75), hide.ns = TRUE,bracket.size = 0.2,size=2)
a05_cazyme_mbp
ggsave(filename = "02a_CAZyme_density.pdf", a05_cazyme_mbp, width = 70, height = 70 , units = "mm")

#################################################################################
# NMDS considering CAZy family density
library(tidyverse)
b01_CAZy_dens<-read.csv("../Annotation_table.csv")%>%
  select(query_name,Class,MAG,Fam_sub_CAZy,Genome_length_Mbp)%>%
  na.omit()%>%
  mutate(Class = sub(".*__", "", Class))%>%
  filter(grepl("Bacteroidia|Clostridia|Bacilli",Class))%>%
  group_by(MAG,Fam_sub_CAZy)%>%
  mutate(count=n())%>%
  mutate(density=count/Genome_length_Mbp)%>%
  select(Class,MAG,Fam_sub_CAZy,density)%>%
  distinct()%>%
  pivot_wider(id_cols=c("Class","MAG"),names_from = Fam_sub_CAZy, values_from = density,values_fill = 0)

b02_CAZy_dens_NMDS<-b01_CAZy_dens%>%
  select(-c(Class))%>%
  column_to_rownames("MAG")

################################################################################################
# NMDS BRAY CURTIS binary
# Calculate dissimilarities ----
## Weighted dissimilarities
#c01_bray_KO_weighted <- b01_KO_density %>% vegdist(x = ., method = "bray") 
#c01_bray_CAZY_weighted <- b01_CAZy_density %>% vegdist(x = ., method = "bray") 
#c01_bray_SULF_weighted <- b01_SULF_density %>% vegdist(x = ., method = "bray") 

library(vegan)
## Weighted dissimilarities
b02_bray_CAZY_weighted <- b02_CAZy_dens_NMDS %>% vegdist(x = ., method = "bray")

# Collect dissimilarities into a list for collective processing ----
b03_bray_CAZy_list <- list("b02_bray_CAZY_weighted" = b02_bray_CAZY_weighted)

# Perform non-metric multidimensional scaling (nMDS) ----
b04_CAZy_nmds_list <- map(b03_bray_CAZy_list, function(bray) metaMDS(bray, trymax = 999))

# Extract data from nMDS ----
## Obtain coordinates
b05_CAZy_scrs_list <- map(b04_CAZy_nmds_list, function(nmds) {scores(nmds, display = "sites") %>%as.data.frame() %>%rownames_to_column("sample")})

## Collect nMDS scores in a single data frame
b06_CAZy_scrs_all <- bind_rows(b05_CAZy_scrs_list, .id = "data_type")

## Collect nMDS statistics (stress values)
b07_CAZy_stress_values <- map(b05_CAZy_scrs_list, function(nmds) {data.frame("label" = paste("Stress =", nmds$stress))}) %>%bind_rows(.id = "data_type")

# Plot nMDS using ggplot2 ----
## Set up colours
b08_metadata<-b01_CAZy_dens%>%select(MAG,Class)%>%distinct()%>%mutate(Group = as.factor(Class))%>%select(-c(Class))

## Append sample grouping to scores
b09_CAZy_scrs_all <- left_join(b06_CAZy_scrs_all, b08_metadata, by = c("sample" = "MAG"))

## Create panel headers
b10_CAZy_panel_labels <- c("b02_bray_CAZY_weighted" = "CAZyme density\n(weighted Bray-Curtis)")

# centroids
b11_CAZy_centroid<-b09_CAZy_scrs_all%>%group_by(Group)%>%summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))

## Call ggplot
library(ggplot2)
b12_CAZy_density_bray_nmds<-ggplot(data = b09_CAZy_scrs_all, aes(x = NMDS1, y = NMDS2, color=Group)) +
  geom_point(data=b09_CAZy_scrs_all,size = 1) +
  scale_color_manual(values =c("Bacteroidia"="#B2DF8A",
                               "Bacilli"="#FDBF6F",
                               "Clostridia"="#FF7F00")) +
  theme_bw() +
  labs(title = "c) CAZy family density") +
  theme(plot.title = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 5),
        axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 7, face="bold"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 7, face="bold"),
        axis.ticks.y = element_line(size = 0.25),
        axis.ticks.x = element_line(size = 0.25),
        rect = element_rect(size = 0.1),
        line = element_line(size = 0.1),
        plot.margin = margin(1,1,1,1, unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.5)) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")
b12_CAZy_density_bray_nmds
ggsave("02b_CAZyme_density_NMDS.pdf", b12_CAZy_density_bray_nmds, width = 70, height = 70 , units = "mm")

library(patchwork)
dens_bray_nmds <-  a05_cazyme_mbp / b12_CAZy_density_bray_nmds + plot_layout(heights = c(70,70))
dens_bray_nmds
ggsave("02c_density_nmds.pdf", dens_bray_nmds, width = 70, height = 140,units = "mm")

