# Mannitol cluster search based on the most frequent of pfams surrounding mannitol associated genes
# clusters will be identified if any of the annotations defined prior are found less than 5 genes apart
# these will be called mannitol operons

setwd("C:/Users/cfac668/OneDrive - The University of Auckland/Chapter_1/Revisiting/Review_analysis_V2_CF17_KH12_KDC4_ISME_communications/R_codes_submitted/Figure 2/")
library(dplyr)
library(tidyr)
library(stringr)

a01_annotation<-read.csv("../Annotation_table.csv", header = TRUE)%>%
  select(MAG,query_name,description_uniprot,description_uniref,description_kegg,sseqid_pfam,description_pfam, sseqid_tigrfam,description_tigrfam)

# Search for genes associated to mannitol annotations
a02_mannitol_screenning <- a01_annotation %>%
  mutate(
    Fructokinase_suggestive = case_when(
      str_detect(description_uniprot, "Fructokinase|fructokinase|Hexokinase|hexokinase|PFK|PfkB") ~ 'Fructokinase',
      str_detect(description_uniref, "Fructokinase|fructokinase|Hexokinase|hexokinase|PFK|PfkB") ~ 'Fructokinase',
      str_detect(description_kegg, "Fructokinase|fructokinase|Hexokinase|hexokinase|PFK|PfkB") ~ 'Fructokinase',
      str_detect(description_pfam, "Fructokinase|fructokinase|Hexokinase|hexokinase|PFK|PfkB") ~ 'Fructokinase',
    ))%>%
  mutate(
    Transporter_suggestive = case_when(
      str_detect(description_uniprot, "PTS|Phosphotransferase|Phosphoenolpyruvate|transporter|Transporter") ~ 'Transporter',
      str_detect(description_uniref, "PTS|Phosphotransferase|Phosphoenolpyruvate|transporter|Transporter") ~ 'Transporter',
      str_detect(description_kegg, "PTS|Phosphotransferase|Phosphoenolpyruvate|transporter|Transporter") ~ 'Transporter',
      str_detect(description_pfam, "PTS|Phosphotransferase|Phosphoenolpyruvate|transporter|Transporter") ~ 'Transporter',
    ))%>%
  mutate(
    Mannitol_suggestive = case_when(
      str_detect(description_uniprot, "Mannitol dehydrogenase|Mannitol-1-phosphate 5-dehydrogenase|Mannitol 2-dehydrogenase") ~ 'Mannitol associated',
      str_detect(description_uniref, "Mannitol dehydrogenase|Mannitol-1-phosphate 5-dehydrogenase|Mannitol 2-dehydrogenase") ~ 'Mannitol associated',
      str_detect(description_kegg, "mannitol dehydrogenase|Mannitol dehydrogenase|Mannitol-1-phosphate/altronate dehydrogenases|mannitol-1-phosphate 5-dehydrogenase") ~ 'Mannitol associated',
      str_detect(description_pfam, "mannitol|Mannitol") ~ 'Mannitol associated',
    ))%>%
  select(MAG,query_name,Mannitol_suggestive,Fructokinase_suggestive,Transporter_suggestive)%>%
  rename(bin=MAG,contig_gene=query_name)
write.csv(a02_mannitol_screenning,"02a_mannitol_associated_to_screen.csv",row.names = FALSE)

########################################################################
# Search for genes surrounding these annotations (adapted from Mike Hoggard script for operon search)
# use a sliding window of 5 containing at least 2 genes associated/annotated as mannitol

# Function to identify gene clusters and return new table
identify_gene_clusters <- function(input_table, sliding_window, min_cluster_size) {
  # Load (and/or install and load) required libraries
  dynamic_require <- function(x){
    for( i in x ){
      if( ! require( i , character.only = TRUE ) ){
        install.packages( i , dependencies = TRUE , repos = "http://cran.us.r-project.org", quiet = TRUE)
        require( i , character.only = TRUE )
      }
    }
  }
  dynamic_require( c("dplyr" , "tibble" , "readr", "tidyr", "stringr") )
  # Read in table; convert to long format (databases in one column; annotations in another column), filter out rows with no annotations; split contig_gene into contig and gene_number columns
  annotations_df <- read_csv(input_table, col_types = cols(.default = "c")) %>%
    pivot_longer(cols = -c(bin, contig_gene), names_to = "database", values_to = "annotation") %>% 
    filter(!is.na(annotation)) %>%  
    mutate(contigID = gsub("_[0-9]+$", "", contig_gene)) %>%
    mutate(gene_number = as.numeric(gsub(".*_([0-9]+$)", "\\1", contig_gene))) %>%
    arrange(., bin, contigID, gene_number)
  # initiate cluster_number
  cluster_number = 1
  # Loop row-by-row
  ## For each: set cluster_number; check if contig id same as previous row; if so, check if gene_number within x genes of previous row; if so, keep cluster number and move to next row; if not, update to new cluster number. 
  for (row in 1:nrow(annotations_df)) {
    if (row == 1) {
      annotations_df[row, "cluster_number"] <- cluster_number
    } else {
      if (annotations_df[row, "contigID"] == annotations_df[(row-1), "contigID"]) {
        if ((annotations_df[row, "gene_number"] - annotations_df[(row-1), "gene_number"]) < sliding_window) {
          annotations_df[row, "cluster_number"] <- cluster_number
        } else {
          cluster_number = cluster_number + 1
          annotations_df[row, "cluster_number"] <- cluster_number
        }
      } else {
        cluster_number = cluster_number + 1
        annotations_df[row, "cluster_number"] <- cluster_number
      }
    }
  }
  # Filter to retain only gene clusters with more entries than min_cluster_size, and pivot_wider back into original format (database and annotation split over separate columns)
  annotations_df <- annotations_df %>% 
    group_by(cluster_number) %>% 
    pivot_wider(names_from = database, values_from = annotation) %>%
    filter(n() >= min_cluster_size) 
  # Create subset table of clusters that contain any of the refined annotated pfam and keyword from m2dh m2dh or PTS
  subset_df <- annotations_df %>% 
    pivot_longer(cols = c("Mannitol_suggestive","Fructokinase_suggestive","Transporter_suggestive"), names_to = "gene", values_to = "annotation") %>%
    group_by(cluster_number) %>% 
    filter(any(!is.na(annotation)) & any(str_detect(annotation, pattern = "Mannitol"))) %>%
    pivot_wider(names_from = gene, values_from = annotation)
  # # # Create list of output tables
  tables_list <- list(annotations_df, subset_df)
  # Return new tables
  return(tables_list)
}

# Run function
clusters_tables_list <- identify_gene_clusters('02a_mannitol_associated_to_screen.csv', 5, 1)
write_csv(as.data.frame(clusters_tables_list[1]), '02b_mannitol_clusters_result.csv')
write_csv(as.data.frame(clusters_tables_list[2]), '02c_mannitol_subset.csv')

########################################################################
########################################################################
########################################################################
########################################################################
# filter clusters containing mannitol genes
b01_clusters<-read.csv("02c_mannitol_subset.csv")%>%
  rename(query_name=contig_gene)%>%
  left_join(a01_annotation,by = "query_name")%>%
  mutate(
    mannitol_genes = case_when(
      # PTS system
      (str_detect(Transporter_suggestive, str_c("Transporter")) & str_detect(sseqid_pfam, str_c("PF00359.22|PF02378.18|PF02302.17|PF03480.13|PF04290.12|PF06808.12")))~ 'mannitol PTS system',
      (str_detect(Transporter_suggestive, str_c("Transporter")) & str_detect(sseqid_tigrfam, str_c("TIGR00848|TIGR01419|TIGR00851|TIGR01995|TIGR00787|TIGR00786|TIGR02123")))~ 'mannitol PTS system',
      
      # mannitol 1 phosphate 5-dehydrogenase
      (str_detect(Mannitol_suggestive, str_c("Mannitol associated")) & str_detect(description_uniprot, str_c("Mannitol-1-phosphate 5-dehydrogenase")))~ 'mannitol-1-phosphate 5-dehydrogenase',
      (str_detect(Mannitol_suggestive, str_c("Mannitol associated")) & str_detect(description_uniref, str_c("Mannitol-1-phosphate 5-dehydrogenase")))~ 'mannitol-1-phosphate 5-dehydrogenase',
      (str_detect(Mannitol_suggestive, str_c("Mannitol associated")) & str_detect(description_kegg, str_c("mannitol-1-phosphate 5-dehydrogenase")))~ 'mannitol-1-phosphate 5-dehydrogenase',
      (str_detect(Mannitol_suggestive, str_c("Mannitol associated")) & str_detect(sseqid_pfam, str_c("PF00107.26|PF08240.12")))~ 'mannitol-1-phosphate 5-dehydrogenase',
      (str_detect(Mannitol_suggestive, str_c("Mannitol associated")) & str_detect(sseqid_tigrfam, str_c("TIGR00692|TIGR01202|TIGR01751|TIGR02818|TIGR02824|TIGR03451|TIGR03989|TIGR02819")))~ 'mannitol-1-phosphate 5-dehydrogenase',
      
      # mannitol 2 dehydrogenase
      (str_detect(Mannitol_suggestive, str_c("Mannitol associated")) & str_detect(sseqid_pfam, str_c("PF01232.23|PF08125.13")))~ 'mannitol 2-dehydrogenase',
      
      # fructokinase
      (str_detect(Fructokinase_suggestive, str_c("Fructokinase")) & str_detect(sseqid_pfam, str_c("PF00294.24|PF08543.12|PF00365.20|PF01256.17|PF01467.26|PF00294.24|PF08543.12|PF01513.21|PF00365.20")))~ 'fructokinase',
      (str_detect(Fructokinase_suggestive, str_c("Fructokinase")) & str_detect(sseqid_tigrfam, str_c("TIGR00097|TIGR00125|TIGR00196|TIGR01231|TIGR01518|TIGR02152|TIGR02198|TIGR02199|TIGR03168|TIGR03828|TIGR04382|TIGR00687|TIGR02477|TIGR02478|TIGR02482|TIGR02483")))~ 'fructokinase',
    ))%>%
  select(bin,query_name,mannitol_genes)%>%
  na.omit()%>%
  rename(MAG=bin)
write.csv(b01_clusters,"02d_mannitol_genes.csv",row.names = FALSE)
