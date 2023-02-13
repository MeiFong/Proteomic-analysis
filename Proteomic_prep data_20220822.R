
library(readxl)
library(data.table)
library(tidyverse)

#load data sheet####
path <- "L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Proteomic analysis/20220701_Spheroids_AllProteins_MF.xlsx"

dataframe <- path %>%
  excel_sheets() %>%
  set_names() %>%
  map_dfr(read_excel, path = path, .id = "cell_line") %>%
  dplyr::select(cell_line, Gene, 'Exp. q-value', "Normalized Abundance 1B9M 1", "Normalized Abundance 1B9M 2", 
         "Normalized Abundance 1B9M 3", "Normalized Abundance 1B9M T 1", "Normalized Abundance 1B9M T 2",
         "Normalized Abundance 1B9M T 3", "Normalized Abundance 5B5M 1", "Normalized Abundance 5B5M 2", 
         "Normalized Abundance 5B5M 3", "Normalized Abundance 5B5M T 1", "Normalized Abundance 5B5M T 2",
         "Normalized Abundance 5B5M T 3", "Normalized Abundance 9B1M 1", "Normalized Abundance 9B1M 2", 
         "Normalized Abundance 9B1M 3", "Normalized Abundance 9B1M T 1", "Normalized Abundance 9B1M T 2", 
         "Normalized Abundance 9B1M T 3", "Normalized Abundance B 1", "Normalized Abundance B 2","Normalized Abundance B 3",
         "Normalized Abundance B T 1", "Normalized Abundance B T 2", "Normalized Abundance B T 3", "Normalized Abundance M 1",
         "Normalized Abundance M 2", "Normalized Abundance M 3", "Normalized Abundance M T 1","Normalized Abundance M T 2",
         "Normalized Abundance M T 3", "Normalized Abundance Eml 1", "Normalized Abundance Eml 2", 
         "Normalized Abundance Eml 3", "Normalized Abundance Eml T 1", "Normalized Abundance Eml T 2", 
         "Normalized Abundance Eml T 3") %>%
  filter(!if_all(4:39, is.na)) %>%
  pivot_longer(c("Normalized Abundance 1B9M 1", "Normalized Abundance 1B9M 2", "Normalized Abundance 1B9M 3", "Normalized Abundance 1B9M T 1",
                 "Normalized Abundance 1B9M T 2", "Normalized Abundance 1B9M T 3", "Normalized Abundance 5B5M 1", "Normalized Abundance 5B5M 2",
                 "Normalized Abundance 5B5M 3", "Normalized Abundance 5B5M T 1", "Normalized Abundance 5B5M T 2", "Normalized Abundance 5B5M T 3", 
                 "Normalized Abundance 9B1M 1", "Normalized Abundance 9B1M 2", "Normalized Abundance 9B1M 3", "Normalized Abundance 9B1M T 1", 
                 "Normalized Abundance 9B1M T 2", "Normalized Abundance 9B1M T 3", "Normalized Abundance B 1", "Normalized Abundance B 2", 
                 "Normalized Abundance B 3", "Normalized Abundance B T 1", "Normalized Abundance B T 2", "Normalized Abundance B T 3", 
                 "Normalized Abundance M 1", "Normalized Abundance M 2", "Normalized Abundance M 3", "Normalized Abundance M T 1", 
                 "Normalized Abundance M T 2", "Normalized Abundance M T 3", "Normalized Abundance Eml 1", "Normalized Abundance Eml 2", 
                 "Normalized Abundance Eml 3", "Normalized Abundance Eml T 1", "Normalized Abundance Eml T 2", "Normalized Abundance Eml T 3"), 
               names_to = "group", values_to = "normalised_abundance") %>% 
  separate(group, c("1", "2", "ID", "group", "3"), sep = " ") %>%
  dplyr::select(-'1', -'2', -'3') %>%
  mutate(group = str_replace(group, "[1-3]", "C"))

dataframe_gl <- dataframe %>%
  dplyr::select(Gene) %>%
  distinct()

#data: replace NA with 10k####
dataframe_2 <- dataframe %>%
  group_by(cell_line, Gene) %>%
  mutate(Missing = sum(is.na(normalised_abundance)), 
         Valid = sum(!is.na(normalised_abundance)), 
         percentMissing = Missing/(Missing+Valid)) %>% 
  #filter(percentMissing  <= 0.8) %>%
  replace(is.na(.), 10000) %>%
  mutate(logAbundance = log2(normalised_abundance)) %>%
  group_by (cell_line, Gene) %>%
  mutate(group = factor(group, levels = c("C", "T")),
         ID = factor(ID, levels = c("Eml", "B", "M", "9B1M", "5B5M", "1B9M"))) %>%
  nest() %>%
  mutate(model = map(data, ~glm(logAbundance ~ ID * group, family = gaussian, data = .)),
         model_tidy = map(model, broom::tidy)) %>%
  unnest(model_tidy) %>%
  ungroup %>%
  dplyr::select(-data, -model) %>%
  filter(term %like% ":groupT")

write_csv(dataframe_2, "L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Proteomic analysis/20220822_Proteomic_10k.csv")

#Clear environment
rm(list = ls(all.names = TRUE))