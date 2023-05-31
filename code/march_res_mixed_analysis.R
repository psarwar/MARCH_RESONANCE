## MARCH and RESONANCE joint ordination
##look at files inthe working directory
list.files()
# occasionally might need to clear environment for sanity
rm(list =ls())

# load required packages
library(tidyverse)
library(vegan)
library(glue)

## Merge metaphlan outputs
# terminal; env metaphlan; merge_metaphlan_tables.py *_profile.tsv > MARCH_merged_ab_table.tsv
# couldnt join the RES and MARCH profiles directly because ph metaphlan version issues
# merged_metaphlan_tables2.py is a separate exec script associated with a version of metaphlan
# it is in the folder with all the resonance metaphlan profiles
./merge_metaphlan_tablesv2.py *_profile.tsv > mergedab.tsv

# Call in the two abundance tables
res_ab <- read_tsv("./raw_data_res/seq_res_mergedab.tsv", skip=1)
#the following dont have samples 71; 144; 166, 193
march_ab <- read_tsv("./raw_data_march/march_mergedab.tsv", skip = 1)

# have samples as rows and taxa as column
res_trans <- res_ab %>% select(-NCBI_tax_id) %>% pivot_longer(cols = -1) %>%
  pivot_wider(names_from = "clade_name", values_from = "value") %>%
  select(name, contains("s__")) 

march_trans <- march_ab %>% select(-NCBI_tax_id) %>% pivot_longer(cols = -1) %>%
  pivot_wider(names_from = "clade_name", values_from = "value") %>% 
  select(name, contains("s__"))

# join and replace all NAs with 0
marchres_ab <- full_join(res_trans, march_trans) %>%
  mutate(seqname = str_replace_all(name, "_profile", "")) %>% select(-name) %>%
  select(seqname, everything())

####Mixed Ordination of MARCH and RESONANCE ####
## Calculate Bray-Curtis
bc_input <- marchres_ab %>%
  mutate_if(is.numeric, coalesce, 0) %>%
  column_to_rownames("seqname")
  
set.seed(1994) # for reproducibility 

# Bray-Curtis Distance calculation
marchres_bc <- vegdist(bc_input, method = "bray")

# PCoA calculation
marchres_pcoa <- cmdscale(marchres_bc, k=2, eig = T, add = T)
positions <- marchres_pcoa$points
colnames(positions) <- c("PCoA1", "PCoA2")

## Calculate Axis variance
percent_explained <- 100 * marchres_pcoa$eig / sum(marchres_pcoa$eig)
percent_explained[1:2]
percent_explained_round <- format(round(percent_explained[1:2], digits = 2), nsmall = 1, trim = T)

##Plot PCoA
labs <- c(glue("PCoA 1 ({percent_explained_round[1]}%)"),
          glue("PCoA 2 ({percent_explained_round[2]}%)"))

plot_pcoa <- positions %>%
  as_tibble(rownames = "seqname") %>%
  left_join(., joined_md, by = "seqname") %>%
  left_join(., marchres_ab, by = "seqname") %>% filter(!is.na(child_age_month)) %>%
  ggplot(aes(x=PCoA1, y=PCoA2, shape = cohort)) +
  geom_point(aes(color = child_age_month), 
             position = position_jitterdodge(dodge.width=0.9), size = 3) +
  scale_color_gradient(low = "#FDE725FF", high = "#440154FF") +
  labs(x = labs[1], y = labs[2], shape = "Cohort", color = "Age (months)") +
  stat_ellipse() + 
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))


plot_pcoa
ggsave(plot_pcoa, file = "./pictures/bc_pcoa_age_cohort.jpeg", device = "jpeg")

## PERMANOVA
##distance matrix for adonis input
set.seed(1994)
adonis_input <- as.matrix(marchres_bc)
adonis_input <- as_tibble(adonis_input, rownames = "seqname") %>%
  inner_join(., joined_md, by = "seqname")
adonis_input$cohor

dist <- adonis_input %>%
  select(all_of(.[["seqname"]])) %>%
  as.dist()

perm1 <- adonis2(dist~cohort, data = adonis_input, permutations = 9999, na.action = na.omit)
perm2 <- adonis2(dist~child_age_month, data = adonis_input, permutations = 9999, na.action = na.omit)
perm3 <- adonis2(dist~cohort*child_age_month, data = adonis_input, permutations = 9999, na.action = na.omit)
perm4 <- adonis2(dist~child_age_month+cohort, data = adonis_input, permutations = 9999, na.action = na.omit)
perm3
#### Stacked Barplot of HMO metabolizer ####
## metadata

# Call in filenames of all metaphlan outputs
filenames <- read_csv("./raw_data_res/filenames.csv", col_names = FALSE) %>%
  separate(X1, into = c("sample", NA), sep = "_" , remove = FALSE) %>%
  mutate(seqname = str_replace_all(X1, "_profile.tsv", "")) %>% select(-X1)

res_md <- read_csv("./processed/metadata_res_seq.csv") %>% 
  left_join(., filenames, by = "sample") %>%
  rename(eczema = eczema_studyID, child_age_month = childAgeMonths) %>%
  mutate(cohort = "RES")

march_md <- read_csv("./processed/march_seq_metadata_allvar.csv") %>%
  rename(child_id = SAMPLEID, sample = specimen_ID) %>%
  mutate(cohort = "MARCH")

joined_md <- full_join(march_md, res_md, by = c("sample", "feedtype", "eczema", 
                                                "child_age_month", "seqname", "cohort"))
## Select the columns with only HMO metabolizers
## list of HMO metabolizers picked from eczema review paper Table 2

hmo_met <- joined_md %>%
  left_join(., marchres_ab, by = "seqname") %>% 
  select(seqname, sample, cohort, eczema, feedtype, child_age_month, 
         contains(c("s__Bifidobacterium_bifidum", "s__Bifidobacterium_breve", 
                    "s__Bifidobacterium_longum", "s__Bifidobacterium_catenulatum",
                    "s__Bifidobacterium_animalis", "s__Bacteroides_fragilis", 
                    "s__Lactobacilus_casei", "s__Lactobacilus_acidophilus", 
                    "s__Escherichia_coli", "s__Klebsiella_pneumoniae", 
                    "s__Faecalibacterium_prausnitzii", "s__Ruminococcus_gnavus", 
                    "s__Akkermansia_muciniphila"))) %>%
  pivot_longer(cols = contains("s__"), names_to = "HMO_met", 
               values_to = "relab") %>%
  separate(HMO_met, into = c(NA, "HMO_met"), sep = "s__" , remove = FALSE) %>%
  mutate(feedtype = factor(feedtype, levels = c("BreastFed", "Mixed", "FormulaFed")))

taxon_eczema <- hmo_met %>% group_by(cohort, eczema, HMO_met) %>% 
  reframe(mean_relab = mean(relab)) %>% filter(!is.na(eczema)) %>% 
  ggplot(aes(x=eczema, y= mean_relab, fill = HMO_met)) +
  geom_col() + facet_wrap(~cohort) + theme_classic()

hmo_met %>% group_by(cohort, feedtype, HMO_met) %>% 
  reframe(mean_relab = mean(relab)) %>% filter(!is.na(feedtype)) %>% 
  ggplot(aes(x=feedtype, y= mean_relab, fill = HMO_met)) +
  geom_col() + facet_wrap(~cohort) + theme_classic()

taxon_cohort <- hmo_met %>% group_by(cohort, HMO_met) %>% 
  mutate(cohort = case_when(cohort == "RES" ~ "RESONANCE", TRUE ~ as.character(cohort))) %>%
  mutate(HMO_met = str_replace(HMO_met, "_", " ")) %>%
  reframe(mean_relab = mean(relab)) %>% 
  ggplot(aes(x=cohort, y= mean_relab, fill = HMO_met)) +
  geom_col() + 
  labs(x = NULL, y= "Mean Relative Abundance (%)", fill = "HMO metabolizing bacteria") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))


taxon_cohort
