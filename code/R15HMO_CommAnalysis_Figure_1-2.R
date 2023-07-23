## Community Diversity Analysis of MARCH and RESONANCE Figures 1 and 2
# occasionally might need to clear environment for sanity
rm(list =ls())

# load required packages
library(tidyverse)
library(viridis)
library(vegan)
library(glue)

####INPUT FILES####
## Call in the taxon relative abundances from metaphlan
# merge the resonance and march metaphlan outputs from metaphlan v3.1
# put all profiles in one folder and use metphlans merge tables function/script
# ./merge_metaphlan_tablesv2.py *_profile.tsv > marchres_merged_relab.tsv
relab <- read_tsv("./raw_data_mixed/marchres_merged_relab.tsv", skip = 1)
## Call in the RESONANCE metadata
marchres_md <- read_csv("./metadata_exports/marchres_combinedmd_ids.csv")
# Call in the metadata for the whole community samples including those not sequenced
res_wholemd <- read_csv("./metadata_exports/commdem_resmd.csv") %>% 
  mutate(cohort = "RESONANCE") %>%
  rename(eczema_timepoint = eczema) %>%
  rename(eczema = eczema_studyID)
march_wholemd <- read_csv("./metadata_exports/commdem_march3momd.csv")%>% 
  mutate(cohort = "MARCH")

####CHI^2 BREASTFEEDING AND ECZEMA####
## RESONANCE + MARCH
# Make new breastfeeding vs. non breastfeeding column
marchres_wholemd <- march_wholemd %>% rename (sample = SAMPLEID) %>%
  full_join(res_wholemd, .) %>%
  mutate(feedtype_bfvnbf = case_match(feedtype, c("Mixed", "FormulaFed") ~ "Non-BreastFed",
                                      .default = feedtype)) %>%
  mutate(feedtype = ordered(feedtype, levels = c("BreastFed", "Mixed", "FormulaFed")))

#Pearson's Chi^2
chi2<- table(marchres_wholemd$feedtype_bfvnbf, marchres_wholemd$eczema)
chisq.test(chi2)

p1 <- marchres_wholemd %>% filter(!is.na(eczema)) %>%
  filter(!is.na(feedtype)) %>%
  mutate(eczema = factor(eczema, levels = c("Yes", "No"))) %>%
  ggplot(aes(x=eczema, fill = feedtype)) +
  geom_bar(position = "stack") +
  scale_fill_viridis(discrete = T) +
  labs(x="Eczema", y= "No. of Samples", fill = "Feeding Status") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))

ggsave(filename="f1_stacked_bfec.pdf",
       plot = p1,
       device="pdf",path="./images",
       width= 5,
       height=5,
       units="in",
       dpi=500)

####ALPHA DIVERSITY####

# make the input matrix for vegdist
bc_input <- relab %>%
  select(-NCBI_tax_id) %>%
  pivot_longer(cols = -1) %>%
  pivot_wider(names_from = "clade_name", values_from = "value") %>%
  separate(name, into = c("uid", NA), sep = "_") %>%
  column_to_rownames("uid")

# Bray-Curtis Distance calculation
# for reproducibility 
set.seed(1994)
marchres_bc <- vegdist(bc_input, method = "bray")

marchres_ad <- relab %>% 
  select(-NCBI_tax_id) %>%
  pivot_longer(cols = -1) %>% 
  pivot_wider(names_from = "clade_name", values_from = "value") %>%
  select(name, contains("s__")) %>%
  separate(name, into = c("uid", NA), sep = "_" )%>%
  pivot_longer(cols = -1) %>%
  group_by(uid) %>%
  summarize(shannon = diversity(value,index = "shannon"), 
            simpson = diversity(value, index = "simpson"),
            inverse_simpson = diversity(value, index = "invsimpson"))

marchres_admd <- full_join(marchres_md, marchres_ad) %>%
  mutate(feedtype = factor(feedtype, 
                           levels = c("BreastFed", "Mixed", "FormulaFed")),
         eczema = factor(eczema, levels = c("Yes", "No")))
## Eczema
p2 <- marchres_admd %>% 
  drop_na(eczema) %>% 
  ggplot(aes(x = eczema, y = shannon, fill = eczema)) + 
  geom_boxplot(width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
  scale_fill_manual(values=c("#e00201","#8faadc")) +
  labs(x= "Eczema", y= "Shannon Diversity", fill = "Eczema Status") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))

ggsave(filename="f1_ad_eczema.pdf",
       plot = p2,
       device="pdf",path="./images",
       width= 5,
       height=5,
       units="in",
       dpi=500)

anova_ad_eczema <- aov(shannon ~ eczema, data = marchres_admd)
summary(anova_ad_eczema)

## FeedType 

p3 <- marchres_admd %>% 
  drop_na(feedtype) %>% 
  ggplot(aes(x = feedtype, y = shannon, fill = feedtype)) + 
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
  scale_fill_viridis(discrete = T) +
  labs(x= "Feeding Status", y= "Shannon Diversity", fill = "Feeding Status") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))

ggsave(filename="f1_ad_feeding.pdf",
       plot = p3,
       device="pdf",path="./images",
       width= 5,
       height=5,
       units="in",
       dpi=500)

anova_ad_feeding <- aov(shannon ~ feedtype, data = marchres_admd)
summary(anova_ad_feeding)

####PCOA####
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

relab_uid <- bc_input %>% as_tibble(rownames = "uid")

pcoa_table <- positions %>%
  as_tibble(rownames = "uid") %>%
  full_join(., marchres_md) %>%
  full_join(., relab_uid)
### Cohort By Age
p4 <- pcoa_table %>%
  drop_na(child_age_month) %>%
  ggplot(aes(x=PCoA1, y=PCoA2, shape = cohort)) +
  geom_point(aes(color = child_age_month), size = 4) +
  scale_color_gradient(low = "#FDE725FF", high = "#440154FF") +
  labs(x = labs[1], y = labs[2], shape = "Cohort", color = "Age (months)") +
  stat_ellipse() + 
  theme_classic() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

ggsave(filename="f1_pcoa_cohortbyage.pdf",
       plot = p4,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)

##Cohort By B.longum Relab
p5 <- pcoa_table %>%
  rename(b_longum = `k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|s__Bifidobacterium_longum`) %>%
  drop_na(b_longum) %>%
  ggplot(aes(x=PCoA1, y=PCoA2, shape = cohort)) +
  geom_point(aes(color = b_longum), size = 4) +
  scale_color_gradient(low = "#FDE725FF", high = "#800000") +
  labs(x = labs[1], y = labs[2], shape = "Cohort", color = "Relative Abundance\nB. longum") +
  stat_ellipse() + 
  theme_classic() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

ggsave(filename="f1_pcoa_cohortbyblongumrelab.pdf",
       plot = p5,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)
