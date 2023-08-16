## Community Diversity Analysis of MARCH and RESONANCE Figures 1 and 2
# occasionally might need to clear environment for sanity
rm(list =ls())

# load required packages
library(tidyverse)
library(ggtext)
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

brks <- c(0, 0.25, 0.5, 0.75, 1)

p1_percentages <- marchres_wholemd %>% filter(!is.na(eczema)) %>%
  filter(!is.na(feedtype)) %>%
  mutate(eczema = factor(eczema, levels = c("Yes", "No"))) %>%
  ggplot(aes(x=eczema, fill = factor(feedtype))) +
  geom_bar(position = "fill") +
  scale_fill_viridis(discrete = T) +
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  labs(x="Eczema", y= "Percentage", fill = "Feeding Status") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))

ggsave(filename="f1_stacked_bfec_percentages.pdf",
       plot = p1_percentages,
       device="pdf",path="./images",
       width= 5,
       height=5,
       units="in",
       dpi=500)

####ALPHA DIVERSITY####

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
  ggplot(aes(x=PCoA1, y=PCoA2, linetype = cohort)) +
  stat_ellipse() +
  geom_point(aes(shape = cohort, color = child_age_month), size = 4) +
  scale_color_gradient(low = "#440154FF", high = "#FDE725FF") +
  labs(x = labs[1], y = labs[2], shape = "Cohort", color = "Age (months)") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), 
        axis.title = element_text(size = 23))

ggsave(filename="f1_pcoa_cohortbyage.pdf",
       plot = p4,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)

### Cohort By Alpha Diversity

p5 <- pcoa_table %>% left_join(.,marchres_ad) %>%
  ggplot(aes(x=PCoA1, y=PCoA2)) +
  geom_point(aes(color = shannon, shape = cohort), size = 4) +
  scale_color_viridis(option = "magma") +
  labs(x = labs[1], y = labs[2], shape = "Cohort", color = "Shannon\nDiversity") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), 
        axis.title = element_text(size = 23))

ggsave(filename="sf1_pcoa_cohortbyshannon.pdf",
       plot = p5,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)

##Cohort By B.longum relab
p6 <- pcoa_table %>%
  rename(b_longum = `k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|s__Bifidobacterium_longum`) %>%
  ggplot(aes(x=PCoA1, y=PCoA2)) +
  geom_point(aes(color = b_longum, shape = cohort), size = 4) +
  scale_fill_viridis() +
  labs(x = labs[1], y = labs[2], shape = "Cohort", color = "Relative Abundance\nB. longum") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), 
        axis.title = element_text(size = 23))

ggsave(filename="sf1_pcoa_cohortbyblongumrelab.pdf",
       plot = p6,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)

##Cohort By B.breve relab
p7 <- pcoa_table %>%
  rename(b_breve = `k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|s__Bifidobacterium_breve`) %>%
  ggplot(aes(x=PCoA1, y=PCoA2, linetype = feedtype)) +
  geom_point(aes(color = b_breve, shape = cohort), size = 4) +
  stat_ellipse() +
  scale_fill_viridis() +
  labs(x = labs[1], y = labs[2], shape = "Cohort", color = "Relative Abundance\nB. breve") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), 
        axis.title = element_text(size = 23))

ggsave(filename="sf1_pcoa_cohortbybbreverelab_feedtypeellipse.pdf",
       plot = p7,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)

## Cohort by feedtype

p8 <- pcoa_table %>%
  drop_na(feedtype) %>%
  mutate(feedtype = factor(feedtype, 
                           levels = c("BreastFed", "Mixed", "FormulaFed"))) %>%
  ggplot(aes(x=PCoA1, y=PCoA2)) +
  geom_point(aes(color = feedtype, shape = cohort), size = 4) +
  #stat_ellipse() +
  scale_color_viridis(discrete = T) +
  labs(x = labs[1], y = labs[2], shape = "Cohort", color = "Feeding Status") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), 
        axis.title = element_text(size = 23))

p8
ggsave(filename="sf1_pcoa_cohortbyfeeding.pdf",
       plot = p8,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)

####PERMANOVA####
##Cohorts Combined
##distance matrix for adonis input
set.seed(1994)
adonis_input <- marchres_bc %>% as.matrix %>%
  as_tibble(rownames = "uid") %>%
  full_join(marchres_md,.) %>%
  mutate(feedtype_bfvnbf = case_match(feedtype, c("Mixed", "FormulaFed") ~ "Non-BreastFed",
                                      .default = feedtype))

dist <- adonis_input %>%
  select(all_of(.[["uid"]])) %>%
  as.dist()

perm1 <- adonis2(dist~cohort, data = adonis_input, 
                 permutations = 9999, na.action = na.omit)
perm2 <- adonis2(dist~child_age_month, data = adonis_input, 
                 permutations = 9999, na.action = na.omit)
perm3 <- adonis2(dist~eczema, data = adonis_input, 
                       permutations = 9999, na.action = na.omit)
perm4 <- adonis2(dist~feedtype_bfvnbf, data = adonis_input, 
                       permutations = 9999, na.action = na.omit)
perm5 <- adonis2(dist~pet_owned, data = adonis_input, 
                       permutations = 9999, na.action = na.omit)
perm6 <- adonis2(dist~household_children, data = adonis_input, 
                       permutations = 9999, na.action = na.omit)
perm7 <- adonis2(dist~maternal_education, data = adonis_input, 
                       permutations = 9999, na.action = na.omit)
perm8 <- adonis2(dist~delivery, data = adonis_input, 
                 permutations = 9999, na.action = na.omit)
perm9 <- adonis2(dist~childcare, data = adonis_input, 
                 permutations = 9999, na.action = na.omit)
perm10 <- adonis2(dist~eczema*child_age_month, data = adonis_input, 
                 permutations = 9999, na.action = na.omit)
perm11 <- adonis2(dist~feedtype_bfvnbf*child_age_month, data = adonis_input, 
                 permutations = 9999, na.action = na.omit)

str(perm1)
str(perm2)
str(perm3)
str(perm4)
str(perm5)
str(perm6)
str(perm7)
str(perm8)
str(perm9)
str(perm10)
str(perm11)

#when update whole permanova need to re-run all the variables since they are called the same things for each sample subset run for ease of join later
variable <- c("Cohort","Age", "Eczema", "BreastFed", "Pets", "Siblings", "Maternal Education", "Delivery Mode", "Childcare", "Eczema*Age", "BreastFed*Age")
Df <- c(perm1$Df[1], perm2$Df[1], perm3$Df[1], perm4$Df[1], perm5$Df[1], 
        perm6$Df[1], perm7$Df[1], perm8$Df[1], perm9$Df[1], perm10$Df[1], 
        perm11$Df[1])
SumofSqs <- c(perm1$SumOfSqs[1], perm2$SumOfSqs[1], perm3$SumOfSqs[1], 
              perm4$SumOfSqs[1], perm5$SumOfSqs[1], perm6$SumOfSqs[1], 
              perm7$SumOfSqs[1], perm8$SumOfSqs[1], perm9$SumOfSqs[1], 
              perm10$SumOfSqs[1], perm11$SumOfSqs[1])
R2 <- c(perm1$R2[1], perm2$R2[1], perm3$R2[1], perm4$R2[1], perm5$R2[1], 
        perm6$R2[1], perm7$R2[1], perm8$R2[1], perm9$R2[1], 
        perm10$R2[1], perm11$R2[1])
`F` <- c(perm1$F[1], perm2$F[1], perm3$F[1], perm4$F[1], perm5$F[1], 
         perm6$F[1], perm7$F[1], perm8$F[1], perm9$F[1], 
         perm10$F[1], perm11$F[1])
`Pr(>F)` <- c(perm1$`Pr(>F)`[1], perm2$`Pr(>F)`[1], perm3$`Pr(>F)`[1], 
              perm4$`Pr(>F)`[1], perm5$`Pr(>F)`[1], perm6$`Pr(>F)`[1], 
              perm7$`Pr(>F)`[1], perm8$`Pr(>F)`[1], perm9$`Pr(>F)`[1], 
              perm10$`Pr(>F)`[1], perm11$`Pr(>F)`[1])

permanova_df <- tibble(variable, Df, SumofSqs, R2, `F`, `Pr(>F)`)
write_csv(permanova_df, file = "permamova_summary.csv")

## MARCH
## distance matrix for adonis input
march_uids <- marchres_md %>%
  filter(cohort == "MARCH") %>%
  .$uid

adonis_input_march <- marchres_bc %>% as.matrix %>%
  as_tibble(rownames = "uid") %>%
  select(uid, march_uids) %>%
  filter(uid %in% march_uids) %>%
  left_join(., marchres_md) %>%
  mutate(feedtype_bfvnbf = case_match(feedtype, c("Mixed", "FormulaFed") ~ "Non-BreastFed",
                                      .default = feedtype))
dist_march <- adonis_input_march %>%
  select(all_of(.[["uid"]])) %>%
  as.dist()

perm2_mar <- adonis2(dist_march~child_age_month, data = adonis_input_march, 
                 permutations = 9999, na.action = na.omit)
perm3_mar <- adonis2(dist_march~eczema, data = adonis_input_march, 
                 permutations = 9999, na.action = na.omit)
perm4_mar <- adonis2(dist_march~feedtype_bfvnbf, data = adonis_input_march, 
                 permutations = 9999, na.action = na.omit)
perm5_mar <- adonis2(dist_march~pet_owned, data = adonis_input_march, 
                 permutations = 9999, na.action = na.omit)
perm6_mar <- adonis2(dist_march~household_children, data = adonis_input_march, 
                 permutations = 9999, na.action = na.omit)
perm7_mar <- adonis2(dist_march~maternal_education, data = adonis_input_march, 
                 permutations = 9999, na.action = na.omit)
perm8_mar <- adonis2(dist_march~delivery, data = adonis_input_march, 
                 permutations = 9999, na.action = na.omit)
perm9_mar <- adonis2(dist_march~childcare, data = adonis_input_march, 
                 permutations = 9999, na.action = na.omit)
perm10_mar <- adonis2(dist_march~eczema*child_age_month, data = adonis_input_march, 
                  permutations = 9999, na.action = na.omit)
perm11_mar <- adonis2(dist_march~feedtype_bfvnbf*child_age_month, data = adonis_input_march, 
                  permutations = 9999, na.action = na.omit)

str(perm2_mar)
str(perm3_mar)
str(perm4_mar)
str(perm5_mar)
str(perm6_mar)
str(perm7_mar)
str(perm8_mar)
str(perm9_mar)
str(perm10_mar)
str(perm11_mar)

#when update whole permanova need to re-run all the variables since they are called the same things for each sample subset run for ease of join later
variable <- c("Age", "Eczema", "BreastFed", "Pets", "Siblings", "Maternal Education", "Delivery Mode", "Childcare", "Eczema*Age", "BreastFed*Age")
Df <- c(perm2_mar$Df[1], perm3_mar$Df[1], perm4_mar$Df[1], perm5_mar$Df[1], 
        perm6_mar$Df[1], perm7_mar$Df[1], perm8_mar$Df[1], perm9_mar$Df[1], perm10_mar$Df[1], 
        perm11_mar$Df[1])
SumofSqs <- c(perm2_mar$SumOfSqs[1], perm3_mar$SumOfSqs[1], 
              perm4_mar$SumOfSqs[1], perm5_mar$SumOfSqs[1], perm6_mar$SumOfSqs[1], 
              perm7_mar$SumOfSqs[1], perm8_mar$SumOfSqs[1], perm9_mar$SumOfSqs[1], 
              perm10_mar$SumOfSqs[1], perm11_mar$SumOfSqs[1])
R2 <- c(perm2_mar$R2[1], perm3_mar$R2[1], perm4_mar$R2[1], perm5_mar$R2[1], 
        perm6_mar$R2[1], perm7_mar$R2[1], perm8_mar$R2[1], perm9_mar$R2[1], 
        perm10_mar$R2[1], perm11_mar$R2[1])
`F` <- c(perm2_mar$F[1], perm3_mar$F[1], perm4_mar$F[1], perm5_mar$F[1], 
         perm6_mar$F[1], perm7_mar$F[1], perm8_mar$F[1], perm9_mar$F[1], 
         perm10_mar$F[1], perm11_mar$F[1])
`Pr(>F)` <- c(perm2_mar$`Pr(>F)`[1], perm3_mar$`Pr(>F)`[1], 
              perm4_mar$`Pr(>F)`[1], perm5_mar$`Pr(>F)`[1], perm6_mar$`Pr(>F)`[1], 
              perm7_mar$`Pr(>F)`[1], perm8_mar$`Pr(>F)`[1], perm9_mar$`Pr(>F)`[1], 
              perm10_mar$`Pr(>F)`[1], perm11_mar$`Pr(>F)`[1])

permanova_df_march <- tibble(variable, Df, SumofSqs, R2, `F`, `Pr(>F)`)

## RESONANCE
## distance matrix for adonis input
res_uids <- marchres_md %>%
    filter(cohort == "RES") %>%
  .$uid

adonis_input_res <- marchres_bc %>% as.matrix %>%
  as_tibble(rownames = "uid") %>%
  select(uid, res_uids) %>%
  filter(uid %in% res_uids) %>%
  left_join(., marchres_md) %>%
  mutate(feedtype_bfvnbf = case_match(feedtype, c("Mixed", "FormulaFed") ~ "Non-BreastFed",
                                      .default = feedtype))
dist_res <- adonis_input_res %>%
  select(all_of(.[["uid"]])) %>%
  as.dist()

perm2_res <- adonis2(dist_res~child_age_month, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm3_res <- adonis2(dist_res~eczema, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm4_res <- adonis2(dist_res~feedtype_bfvnbf, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm5_res <- adonis2(dist_res~pet_owned, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm6_res <- adonis2(dist_res~household_children, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm7_res <- adonis2(dist_res~maternal_education, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm8_res <- adonis2(dist_res~delivery, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm9_res <- adonis2(dist_res~childcare, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm10_res <- adonis2(dist_res~eczema*child_age_month, data = adonis_input_res, 
                      permutations = 9999, na.action = na.omit)
perm11_res <- adonis2(dist_res~feedtype_bfvnbf*child_age_month, data = adonis_input_res, 
                      permutations = 9999, na.action = na.omit)

str(perm2_res)
str(perm3_res)
str(perm4_res)
str(perm5_res)
str(perm6_res)
str(perm7_res)
str(perm8_res)
str(perm9_res)
str(perm10_res)
str(perm11_res)

#when update whole permanova need to re-run all the variables since they are called the same things for each sample subset run for ease of join later
variable <- c("Age", "Eczema", "BreastFed", "Pets", "Siblings", "Maternal Education", "Delivery Mode", "Childcare", "Eczema*Age", "BreastFed*Age")
Df <- c(perm2_res$Df[1], perm3_res$Df[1], perm4_res$Df[1], NA, 
        perm6_res$Df[1], perm7_res$Df[1], perm8_res$Df[1], perm9_res$Df[1], perm10_res$Df[1], 
        perm11_res$Df[1])
SumofSqs <- c(perm2_res$SumOfSqs[1], perm3_res$SumOfSqs[1], 
              perm4_res$SumOfSqs[1], NA, perm6_res$SumOfSqs[1], 
              perm7_res$SumOfSqs[1], perm8_res$SumOfSqs[1], perm9_res$SumOfSqs[1], 
              perm10_res$SumOfSqs[1], perm11_res$SumOfSqs[1])
R2 <- c(perm2_res$R2[1], perm3_res$R2[1], perm4_res$R2[1], NA, 
        perm6_res$R2[1], perm7_res$R2[1], perm8_res$R2[1], perm9_res$R2[1], 
        perm10_res$R2[1], perm11_res$R2[1])
`F` <- c(perm2_res$F[1], perm3_res$F[1], perm4_res$F[1], NA, 
         perm6_res$F[1], perm7_res$F[1], perm8_res$F[1], perm9_res$F[1], 
         perm10_res$F[1], perm11_res$F[1])
`Pr(>F)` <- c(perm2_res$`Pr(>F)`[1], perm3_res$`Pr(>F)`[1], 
              perm4_res$`Pr(>F)`[1], NA, perm6_res$`Pr(>F)`[1], 
              perm7_res$`Pr(>F)`[1], perm8_res$`Pr(>F)`[1], perm9_res$`Pr(>F)`[1], 
              perm10_res$`Pr(>F)`[1], perm11_res$`Pr(>F)`[1])

permanova_df_res <- tibble(variable, Df, SumofSqs, R2, `F`, `Pr(>F)`)

## Combine all the dataframes

permanova_df <- permanova_df %>%
  mutate(dataset = "Combined")

permanova_df_march <- permanova_df_march %>%
  mutate(dataset = "MARCH")

permanova_df_res <- permanova_df_res %>%
  mutate(dataset = "RESONANCE")

permanova_combined <- permanova_df %>% 
  full_join(., permanova_df_march) %>%
  full_join(., permanova_df_res) %>%
  mutate(significant = case_when(`Pr(>F)`> 0.05 ~ "No",
                                 `Pr(>F)`< 0.05 ~ "Yes"))

write_csv(permanova_combined, file = "full_permanova.csv")

##Bubble Plot

bubble_p <- permanova_combined %>%
  ggplot(aes(x=dataset, y=variable, size = R2, color = significant))+ 
  geom_point()+
  labs(y="",x="",colour="",
       size="R^2 Value")+
  scale_y_discrete(limits = rev(c("Cohort", "Age", "BreastFed", "BreastFed*Age", 
                              "Eczema", "Eczema*Age", "Childcare", "Siblings",
                              "Maternal Education", "Delivery Mode", "Pets"))) +
  scale_color_manual(values = c("Yes" = "#229954", "No" = "grey")) +
  theme_classic()

ggsave(filename="f1_permanova_bubble.pdf",
       plot = bubble_p,
       device="pdf",path="./images",
       width= 5,
       height=5,
       units="in",
       dpi=500)

#+
  scale_y_discrete(limits = rev)
theme(panel.background = element_rect(colour = "black", size=1),
      legend.position="right",
      legend.title=element_text(size=10, color="black",face="bold"),
      legend.text=element_text(size=10),
      axis.text.y=element_text(size=10),
      axis.text.x=element_text(size=8, angle = 90, vjust = 0.66),
      axis.ticks.x=element_blank());


####HMO METABOLIZERS RELATIVE ABUDNANCE####

## Select the columns with only HMO metabolizers
## list of HMO metabolizers picked from eczema review paper Table 2
hmo_met <- relab %>%
  select(-NCBI_tax_id) %>%
  pivot_longer(cols = -1) %>%
  pivot_wider(names_from = "clade_name", values_from = "value") %>%
  separate(name, into = c("uid", NA), sep = "_") %>%
  select(uid, contains(c("s__Bifidobacterium_bifidum", "s__Bifidobacterium_breve", 
                         "s__Bifidobacterium_longum", "s__Bifidobacterium_catenulatum",
                         "s__Bifidobacterium_animalis", "s__Bacteroides_fragilis", 
                         "s__Lactobacilus_casei", "s__Lactobacilus_acidophilus", 
                         "s__Escherichia_coli", "s__Klebsiella_pneumoniae", 
                         "s__Faecalibacterium_prausnitzii", "s__Ruminococcus_gnavus", 
                         "s__Akkermansia_muciniphila"))) %>%
  pivot_longer(cols = contains("s__"), names_to = "HMO_metabolizing_bacteria", 
               values_to = "relab") %>%
  separate(HMO_metabolizing_bacteria, into = c(NA, "HMO_metabolizing_bacteria"), sep = "s__") %>%
  left_join(., marchres_md) %>%
  mutate(feedtype = factor(feedtype, levels = c("BreastFed", "Mixed", "FormulaFed")))

## Relative Abundance Plot 
relab_eczema_p <- hmo_met %>%
  drop_na(eczema) %>%
  filter(relab>0) %>%
  ggplot(aes(x = HMO_metabolizing_bacteria, y=relab, fill = eczema)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.5, size = 0.5, 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             dodge.width = 0.6, seed = 511)) +
  labs(x=NULL, y = "Relative Abundance") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, face = "italic"),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 15))

ggsave(filename="f2_relativeabundance_HMOmet.pdf",
       plot = relab_eczema_p,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)

## Relative Abundance log transformed plot
hmo_met_eczema_log <- hmomet_relab_eczema <- hmo_met %>%
  filter(relab>0, !is.na(eczema)) %>%
  ggplot(aes(x = HMO_metabolizing_bacteria, y=log(relab), fill = eczema)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.5, size = 0.5, 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             dodge.width = 0.6, seed = 511)) +
  labs(x=NULL, y = "Relative Abundance (log)") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, face="italic"),
                                   axis.text.y = element_text(size = 10),
                                   axis.title = element_text(size = 15))

ggsave(filename="f2_relativeabundance_HMOmet_logtransformed.pdf",
       plot = hmo_met_eczema_log,
       device="pdf",path="./images",
       width= 12,
       height=10,
       units="in",
       dpi=500)
## Mean Relative Abundance Plot
hmomet_relab_p <- hmo_met %>% 
  group_by(cohort, HMO_metabolizing_bacteria) %>% 
  mutate(HMO_metabolizing_bacteria = str_replace(HMO_metabolizing_bacteria, "_", " ")) %>%
  reframe(mean_relab = mean(relab)) %>% 
  ggplot(aes(x=cohort, y= mean_relab, fill = HMO_metabolizing_bacteria)) +
  geom_col() + 
  labs(x = NULL, y= "Mean Relative Abundance (%)", fill = "HMO metabolizing bacteria") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))
