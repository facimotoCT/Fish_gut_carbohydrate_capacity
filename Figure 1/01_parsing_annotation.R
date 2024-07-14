# make annotation for unrooted host and species tree with all the alistipes to highlight the diversification of the Alistipes in fish
setwd("C:/Users/cfac668/OneDrive - The University of Auckland/Chapter_1/Revisiting/Review_analysis_V2_CF13_KH9_KDC3_ISME_submission/Figure 1/")

################################################################################
# gtdbtk
library(dplyr)
library(tidyverse)
a01_tax <- read.delim("../local_files/02010401_gtdbtk.bac120.summary.tsv") %>%
  select(user_genome, classification) %>%
  rename(genome = user_genome) %>%
  separate(classification, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(Species = gsub("s__\\b", "", Species)) %>%
  mutate(Genus = gsub("g__\\b", "", Genus)) %>%
  mutate(Family = gsub("f__\\b", "", Family)) %>%
  mutate(Order = gsub("o__\\b", "", Order)) %>%
  mutate_all(na_if, "") %>%
  mutate(genomeID = ifelse(is.na(Species),
                           ifelse(is.na(Genus),
                                  ifelse(is.na(Family), Order, Family), Genus), Species))%>%
  mutate(genomeID = ifelse(substr(genome, 1, 2) == "FM", paste(genomeID, genome, sep = " "), genomeID))%>%
  rename(MAG=genome)
write.csv(a01_tax, '01_taxonomy_genomeID_refs_tree.csv', row.names = FALSE)

################################################################################
# checkm
a02_checkm_refined<-read.delim("../local_files/0109_refined-after-duplicate-removal.bins_for_drep.checkm")%>%
  select(Bin.Id,Completeness,Contamination,Strain.heterogeneity) %>%
  rename(MAG=Bin.Id)

################################################################################
# refined bins Bin.Id
library(readxl)
a03_refined<-read_excel("../local_files/0109_refined_to_genome.xlsx")%>%
  left_join(a02_checkm_refined,by = "MAG")%>%
  left_join(a01_tax,by = c("MAG"))
write.csv(a03_refined, '02_taxonomy_checkm_refined.csv', row.names = FALSE)

# NCBI match
a04_NCBI<-read.delim("../local_files/taxonomy_map_summary_GTDB_NCBI.tsv")

################################################################################
# Calculate count per taxonomic level
b01_tax <- a03_refined %>%
  select(MAG, Phylum, Class, Order, Family, Genus) %>%
  group_by(Phylum) %>%
  mutate(Phylum_count = n()) %>%
  group_by(Phylum, Class) %>%
  mutate(Class_count = n()) %>%
  group_by(Phylum, Class, Order) %>%
  mutate(Order_count = n()) %>%
  group_by(Phylum, Class, Order, Family) %>%
  mutate(Family_count = n()) %>%
  group_by(Phylum, Class, Order, Family, Genus) %>%
  mutate(Genus_count = n()) %>%
  mutate(Phylum = ifelse(is.na(Phylum), "unknown", Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), "unknown", Class)) %>%
  mutate(Order = ifelse(is.na(Order), "unknown", Order)) %>%
  mutate(Family = ifelse(is.na(Family), "unknown", Family)) %>%
  mutate(Genus = ifelse(is.na(Genus), "unknown", Genus)) %>%
  arrange(desc(Phylum_count),desc(Class_count))%>%
  mutate(Phylum_name = gsub("p__", "", Phylum))%>%
  mutate(Class_name = gsub("c__", "", Class))
  


b02_Phylum_order<-(unique(b01_tax$Phylum_name))
b02_Class_order<-(unique(b01_tax$Class_name))
b02_Class_order
# b03_cols<-c("Bacteroidia"="#1AFF1A",
#             "Pseudomonadota"="#EEDD88",
#             "Verrucomicrobiae"="#AAAA00",
#             "Desulfovibrionia"="#99DDFF",
#             "Spirochaetia"="#FFAABB",
#             "Vampirovibrionia"="#44BB99",
#             "Bacilli"="#FFC20A",
#             "Clostridia"="#EE8866")

b03_cols<-c("Bacteroidia"="#B2DF8A",
            "Bacilli"="#FDBF6F",
            "Clostridia"="#FF7F00",
            "Spirochaetia"="#FB9A99",
            "Desulfovibrionia"="#FFFF99",
            "Verrucomicrobiae"="#B15928",
            "Alphaproteobacteria"="#CAB2D6",
            "Gammaproteobacteria"="#6A3D9A",
            "Vampirovibrionia"="#1F78B4")

b04_genome_count<-ggplot(b01_tax, aes(fill = Class_name, x = Class_count, y = Class_name, color = Class_name)) +
  geom_bar(aes(y = factor(Class_name, rev(b02_Class_order))),
           position = "dodge",
           stat = "identity",
           color = "black", size = 0.3, show.legend = FALSE) +
  geom_text(aes(label = Class_count), hjust = -0.5, size = 3, color = "black") +  # Add text labels with black color
  scale_fill_manual(values = b03_cols) +
#  scale_x_continuous(name = "MAGs", limits = c(0,40), breaks = seq(0, 40, 5)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 10, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"),
    axis.title.y = element_blank(),  # Remove y-axis title
    rect = element_rect(size = 0.2),
    line = element_line(size = 0.2)
  )
b04_genome_count
ggsave("06_genome_count.pdf", b04_genome_count, width = 3, height = 3)


