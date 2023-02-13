
library(msigdbr)
library(tidyverse)
library(data.table)

setwd("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Proteomic analysis/")
df <- read_csv("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Proteomic analysis/20220822_Proteomic_10k.csv")

#annotation####
ensembl <- read_csv("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Analysis/ensembl_annot_human.csv")

entrez <- read_csv("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Analysis/gene_result.csv")

df_ID <- left_join(df, entrez %>% dplyr::select(Symbol, GeneID), by = c("Gene" = "Symbol"))
df_ID$GeneID[df_ID$Gene == "AKAP2"] <- "445815"
df_ID$GeneID[df_ID$Gene == "ANKRD20A4P"] <- "728747"
df_ID$GeneID[df_ID$Gene == "ARNTL"] <- "406"
df_ID$GeneID[df_ID$Gene == "ARNTL2"] <- "56938"
df_ID$GeneID[df_ID$Gene == "C12orf45"] <- "121053"
df_ID$GeneID[df_ID$Gene == "C16orf70"] <- "80262"
df_ID$GeneID[df_ID$Gene == "C1orf147"] <- "574431"
df_ID$GeneID[df_ID$Gene == "C5orf51"] <- "285636"
df_ID$GeneID[df_ID$Gene == "CBWD6"] <- "644019"
df_ID$GeneID[df_ID$Gene == "DDX58"] <- "23586"
df_ID$GeneID[df_ID$Gene == "DIEXF"] <- "27042"
df_ID$GeneID[df_ID$Gene == "EEF1AKNMT"] <- "51603"
df_ID$GeneID[df_ID$Gene == "ERVK-24"] <- "100862684"
df_ID$GeneID[df_ID$Gene == "EURL"] <- "54149"
df_ID$GeneID[df_ID$Gene == "FAM126A"] <- "84668"
df_ID$GeneID[df_ID$Gene == "FAM160B1"] <- "57700"
df_ID$GeneID[df_ID$Gene == "FAM207A"] <- "85395"
df_ID$GeneID[df_ID$Gene == "GATD3B"] <- "102724023"
df_ID$GeneID[df_ID$Gene == "GBA"] <- "2629"
df_ID$GeneID[df_ID$Gene == "H3-2"] <- "440686"
df_ID$GeneID[df_ID$Gene == "H4-16"] <- "121504"
#df_ID$GeneID[df_ID$Gene == "HERVK_113"] <- ""
df_ID$GeneID[df_ID$Gene == "HNRNPA1P48"] <- "642659"
df_ID$GeneID[df_ID$Gene == "KCT2"] <- "56951"
#df_ID$GeneID[df_ID$Gene == "L1RE1"] <- ""
df_ID$GeneID[df_ID$Gene == "MPP5"] <- "64398"
df_ID$GeneID[df_ID$Gene == "MPP6"] <- "51678"
df_ID$GeneID[df_ID$Gene == "OCC1"] <- "387882"
df_ID$GeneID[df_ID$Gene == "PALM2"] <- "445815"
df_ID$GeneID[df_ID$Gene == "PHB"] <- "5245"
df_ID$GeneID[df_ID$Gene == "RPL9P9"] <- "254948"
df_ID$GeneID[df_ID$Gene == "SKIV2L"] <- "6499"
df_ID$GeneID[df_ID$Gene == "SMAP"] <- "10944"
df_ID$GeneID[df_ID$Gene == "TTC37"] <- "9652"
df_ID$GeneID[df_ID$Gene == "U2AF1L5"] <- "102724594"
df_ID$GeneID[df_ID$Gene == "WDR61"] <- "80349"
df_ID$GeneID[df_ID$Gene == "WDR92"] <- "116143"
df_ID$GeneID[df_ID$Gene == "ZADH2"] <- "284273"

sum(is.na(df_ID))

df_ID <- left_join(df_ID, ensembl %>% dplyr::select(ensembl_gene_id, external_gene_name), by = c("Gene" = "external_gene_name"))

df_ID$ensembl_gene_id[df_ID$Gene == "RPL9P9"] <- "ENSG00000237550"
df_ID$ensembl_gene_id[df_ID$Gene == "GATD3B"] <- "ENSG00000280071"
df_ID$ensembl_gene_id[df_ID$Gene == "LOC102724159"] <- "ENSG00000275464"
df_ID$ensembl_gene_id[df_ID$Gene == "U2AF1L5"] <- "ENSG00000275895"
df_ID$ensembl_gene_id[df_ID$Gene == "EEF1AKNMT"] <- "ENSG00000010165"
df_ID$ensembl_gene_id[df_ID$Gene == "SMAP"] <- "ENSG00000110696"
#df_ID$ensembl_gene_id[df_ID$Gene == "ERVK-24"] <- ""
df_ID$ensembl_gene_id[df_ID$Gene == "WDR92"] <- "ENSG00000243667"
df_ID$ensembl_gene_id[df_ID$Gene == "DIEXF"] <- "ENSG00000117597"
df_ID$ensembl_gene_id[df_ID$Gene == "FAM160B1"] <- "ENSG00000151553"
df_ID$ensembl_gene_id[df_ID$Gene == "MPP6"] <- "ENSG00000105926"
df_ID$ensembl_gene_id[df_ID$Gene == "C16orf70"] <- "ENSG00000125149"
df_ID$ensembl_gene_id[df_ID$Gene == "MPP5"] <- "ENSG00000072415"
df_ID$ensembl_gene_id[df_ID$Gene == "FAM207A"] <- "ENSG00000160256"
df_ID$ensembl_gene_id[df_ID$Gene == "AKAP2"] <- "ENSG00000157654"
df_ID$ensembl_gene_id[df_ID$Gene == "KCT2"] <- "ENSG00000113583"
df_ID$ensembl_gene_id[df_ID$Gene == "PALM2"] <- "ENSG00000157654"
df_ID$ensembl_gene_id[df_ID$Gene == "LOC400499"] <- "ENSG00000188897"
df_ID$ensembl_gene_id[df_ID$Gene == "OCC1"] <- "ENSG00000235162"
df_ID$ensembl_gene_id[df_ID$Gene == "C12orf45"] <- "ENSG00000151131"
df_ID$ensembl_gene_id[df_ID$Gene == "LOC101059915"] <- "ENSG00000283599"
df_ID$ensembl_gene_id[df_ID$Gene == "EURL"] <- "ENSG00000154642"

write.csv(df_ID, "L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Proteomic analysis/DF_ID.csv")

#divide into groups####
grouping <- function(data, cell_ln, term_in){
  out <- data %>%
    filter(cell_line == cell_ln) %>%
    filter(term %like% term_in)
  return(out)
}
MM649_shBRN2 <- grouping(df_ID, 'MM649', 'IDB')
MM649_shMITF <- grouping(df_ID, 'MM649', 'IDM')
MM649_9B1M <- grouping(df_ID, 'MM649', 'ID9B1M')
MM649_5B5M <- grouping(df_ID, 'MM649', 'ID5B5M')
MM649_1B9M <- grouping(df_ID, 'MM649', 'ID1B9M')
MM96L_shBRN2 <- grouping(df_ID, 'MM96L', 'IDB')
MM96L_shMITF <- grouping(df_ID, 'MM96L', 'IDM')
MM96L_9B1M <- grouping(df_ID, 'MM96L', 'ID9B1M')
MM96L_5B5M <- grouping(df_ID, 'MM96L', 'ID5B5M')
MM96L_1B9M <- grouping(df_ID, 'MM96L', 'ID1B9M')

keytypes(org.Hs.eg.db)


#gsea
library(fgsea)

all_gene_sets = msigdbr(species = "Homo sapiens")
all.ensembl_1 <- dplyr::select(all_gene_sets, gs_name, gs_cat, gs_subcat) %>%
  distinct()
all.ensembl <- all_gene_sets %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(ensembl_gene))) %>%
  deframe()

set.seed(6)
gsea_2 <- function (data){
  list <- data %>%
    dplyr::select(ensembl_gene_id, statistic) %>%
    arrange(desc(statistic))
  
  FC.vec <- list$statistic
  names(FC.vec) <- list$ensembl_gene_id
  
  gsea.allgeneset <- fgseaMultilevel(pathways = all.ensembl,
                                 stats = FC.vec,
                                 scoreType = "std",
                                 nPermSimple = 1000) %>%
    left_join(., 
              all.ensembl_1 %>% dplyr::select(gs_name, gs_cat, gs_subcat), 
              by = c("pathway" = "gs_name")) 
  
}

MM649_shBRN2_gsea_2 <- gsea_2(MM649_shBRN2) 
MM649_shMITF_gsea_2 <- gsea_2(MM649_shMITF) 
MM649_9B1M_gsea_2 <- gsea_2(MM649_9B1M) 
MM649_5B5M_gsea_2 <- gsea_2(MM649_5B5M) 
MM649_1B9M_gsea_2 <- gsea_2(MM649_1B9M) 
MM96L_shBRN2_gsea_2 <- gsea_2(MM96L_shBRN2) 
MM96L_shMITF_gsea_2 <- gsea_2(MM96L_shMITF) 
MM96L_9B1M_gsea_2 <- gsea_2(MM96L_9B1M)
MM96L_5B5M_gsea_2 <- gsea_2(MM96L_5B5M) 
MM96L_1B9M_gsea_2 <- gsea_2(MM96L_1B9M) 

gsea_2 <- Map(cbind, gsea, group = names(gsea)) %>%
  bind_rows() %>%
  separate(group, c("cell_line", "group"), sep = "_") %>%
  dplyr::select(cell_line, group, pathway, gs_cat, gs_subcat, pval, padj, ES, NES) %>%
  filter(pval <= 0.05) %>%
  filter(gs_cat %in% c("C2", "H")) %>%
  filter(gs_subcat %in% c("", "CP:KEGG", "CP:PID", "CP")) 

setwd("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Proteomic analysis/Pathway analysis/Finalised")

filtering_2 <- function (data) {
  output <- data %>%
    filter(padj <= 0.05) %>%
    filter(gs_cat %in% c("C2", "H")) %>%
    filter(gs_subcat %in% c("", "CP:KEGG", "CP:PID", "CP")) %>%
    return(output)
}

MM649_shBRN2_filtered_2 <- filtering_2(MM649_shBRN2_gsea_2) 
MM649_shMITF_filtered_2 <- filtering_2(MM649_shMITF_gsea_2) 
MM649_9B1M_filtered_2 <- filtering_2(MM649_9B1M_gsea_2) 
MM649_5B5M_filtered_2 <- filtering_2(MM649_5B5M_gsea_2) 
MM649_1B9M_filtered_2 <- filtering_2(MM649_1B9M_gsea_2) 
MM96L_shBRN2_filtered_2 <- filtering_2(MM96L_shBRN2_gsea_2) 
MM96L_shMITF_filtered_2 <- filtering_2(MM96L_shMITF_gsea_2) 
MM96L_9B1M_filtered_2 <- filtering_2(MM96L_9B1M_gsea_2)
MM96L_5B5M_filtered_2 <- filtering_2(MM96L_5B5M_gsea_2) 
MM96L_1B9M_filtered_2 <- filtering_2(MM96L_1B9M_gsea_2)

MM649 <- list(MM649_shBRN2 = MM649_shBRN2_filtered_2, MM649_shMITF = MM649_shMITF_filtered_2, 
              MM649_9B1M = MM649_9B1M_filtered_2, MM649_5B5M = MM649_5B5M_filtered_2, 
              MM649_1B9M = MM649_1B9M_filtered_2) %>%
  bind_rows(.id = "Group") %>%
  select(-leadingEdge)

write.csv(MM649, file= "MM649_GSEA.csv", row.names =F)

MM96L <- list(MM96L_shBRN2 = MM96L_shBRN2_filtered_2, MM96L_shMITF = MM96L_shMITF_filtered_2, 
              MM96L_9B1M = MM96L_9B1M_filtered_2, MM96L_5B5M = MM96L_5B5M_filtered_2, 
              MM96L_1B9M = MM96L_1B9M_filtered_2) %>%
  bind_rows(.id = "Group") %>%
  select(-leadingEdge)

write.csv(MM96L, file= "MM96L_GSEA.csv", row.names =F)

windowsFonts(Arial = windowsFont("Arial"))
library(ggplot2)
library(ggforce)

#MM649
plot_MM649 <- ggplot(MM649, aes(x = reorder(pathway, -padj), y = NES, size = size), binwidth = 10) +
  geom_point(aes(colour = padj), position = position_dodge2(w = 0.75)) +
  scale_color_gradient (low = "red", high = "blue") +
  scale_radius(limits = c(0, NA), range = c(0, 10)) +
  labs(y = "NES", x = "Pathway", fill = "p.adjust") +
  facet_grid(Group ~ ., scales = "free", space = "free") +
  ggtitle("MM649") +
  theme(text = element_text(family = "Arial"), plot.margin = ) +
  scale_y_continuous(limits=c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  coord_flip() 

ggsave("MM649.tiff", 
       plot = plot_MM649 , device = "tiff", width = 8, height = 18, units = "in", dpi = 600)

plot_MM649

plot_MM96L <- ggplot(MM96L, aes(x = reorder(pathway, -padj), y = NES, size = size), binwidth = 10) +
  geom_point(aes(colour = padj), position = position_dodge2(w = 0.75)) +
  scale_color_gradient (low = "red", high = "blue") +
  scale_radius(limits = c(0, NA), range = c(0, 10)) +
  labs(y = "NES", x = "Pathway", fill = "p.adjust") +
  facet_grid(Group ~ ., scales = "free", space = "free") +
  ggtitle("MM96L") +
  theme(text = element_text(family = "Arial"), plot.margin = ) +
  scale_y_continuous(limits=c(-3, 3.5), breaks = seq(-3, 3.5, by = 1)) +
  coord_flip() 

ggsave("MM96L.tiff", 
       plot = plot_MM96L , device = "tiff", width = 10, height = 40, units = "in", dpi = 600)

plot_MM96L

