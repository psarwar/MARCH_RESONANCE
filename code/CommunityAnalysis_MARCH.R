# Community Level Analysis for MARCH dataset
## contains beta and alpha diversity indices and maaslin run correlations

## look at files in the working directory
list.files()
## clear the R environment
rm(list = ls())

library(vegan)
library(tidyverse)
library(broom)
library(ggrepel)
library(vcdExtra)

#### Resolve Tengfei's QC issue ####
# Tengfei reported some samples as having low QC
# but some of those passed our metaphlan pipeline
# find which ones of Tengfei's reported low QC samples passed our pipeline
# for those samples double check the kneaddata reads output

## Call in all the file with all the sequences

seqnames <- read_csv("./raw_data_march/march_sequence_records.csv")
seqnames <- seqnames %>%
  mutate(tengfei_passQC = str_remove(tengfei_passQC, "_L001_R1_001.fastq.gz"),
         all_sequenced = str_remove(all_sequenced, "_L001_R1_001.fastq.gz"),
         metaphlan_output_vkc = str_remove(metaphlan_output_vkc, "_profile.tsv"))

vec1 <- seqnames$all_sequenced

vec2 <- seqnames$metaphlan_output_vkc

vec3 <- seqnames$tengfei_passQC

setdiff(vec1, vec2)

setdiff(vec1, vec3)

setdiff(vec2, vec3)

setdiff(vec3, vec2)

## our metaphlan outputs are missing sample no. 16; 197; 23
## tengfei's low qc samples are: 71; 133; 144; 166; 193
## according to tengfei's missing qc file duplicate pairs: 71 - 14; 133 - 166; 144 - 139
## 193 is an unknown sample

## 16; 197; 23 do not have kneaddata logs and kevin confirmed that this is expected

# rsync -avP ada:/lovelace/sequencing/processed/mgx/kneaddata/71-3mos_S81_kneaddata.log ./
# rsync -avP ada:/lovelace/sequencing/processed/mgx/kneaddata/133_3mos_S37_kneaddata.log ./
# rsync -avP ada:/lovelace/sequencing/processed/mgx/kneaddata/144_3mos_S48_kneaddata.log ./
# rsync -avP ada:/lovelace/sequencing/processed/mgx/kneaddata/166_3mos_S70_kneaddata.log ./
# rsync -avP ada:/lovelace/sequencing/processed/mgx/kneaddata/193-3mos_S1_kneaddata.log ./
# rsync -avP ada:/lovelace/sequencing/processed/mgx/kneaddata/14-3mos_S62_kneaddata.log ./
# rsync -avP ada:/lovelace/sequencing/processed/mgx/kneaddata/139_3mos_S43_kneaddata.log ./


## Compared the read counts of all seqeunces and 14 has more read counts than 71 - keep 14
## 133 has higher final read count than 166 - keep 133
## 139 has higher final read count than 144 - keep 139
## 139 and 14 are also the duplicates that did not have a QC issue in Tengfei's analysis
## For the final analysis get rid of samples 71; 166; 144; 193 from the vkc metaphlan outputsß


#### Import in raw data/metadata ####
# Merge all outputs from metaphlan using a metaphlan subscript
# merge_metaphlan_tables.py *_profile.tsv > MARCH_merged_ab_table.tsv

## Call in the metadata
# Call in the two merged abundance tables
res_relab <- read_tsv("./raw_data_res/res_merged_relab.tsv", skip=1)
#the following dont have samples 71; 144; 166, 193
march_relab <- read_tsv("./raw_data_march/march_merged_relab.tsv", skip = 1)

## Call in the crosswalk file between the old and new filenames
# biospecimen: old seauence name
# uid : new sequence name
# filename: new sequence name with well information
# S_well: well number
filename_code <- read_csv("./old_new_seqID_mgx.csv")

## Call in the metadata for the sequenced samples
# Use only the sequences identified in the metadata files
# they are the quality controlled cross-sectional samples (see: Demographics_RES/MARCH code files for details)
res_md <- read_csv("./metadata_exports/seqdem_resmd.csv") 
march_md <- read_csv("./metadata_exports/seqdem_marchmd.csv") %>%
  mutate(biospecimen = str_replace_all(seqname, "_S.*", ""))

# Call in the metadata for the whole community samples including those not sequenced
res_wholemd <- read_csv("./metadata_exports/commdem_resmd.csv")
march_wholemd <- read_csv("./metadata_exports/commdem_march3momd.csv")

# Find the new sequence IDs for our MARCH and RESONANCE dataset
res_md <- filename_code %>% 
  select(uid, filename, biospecimen, S_well) %>% 
  rename(sample = biospecimen) %>%
  left_join(res_md, ., by = "sample") %>%
  mutate(cohort = "RESONANCE")

march_md <- filename_code %>% 
  select(uid, filename, biospecimen, S_well) %>%
  distinct(biospecimen, .keep_all = TRUE) %>%
  left_join(march_md, ., by = "biospecimen") %>%
  mutate(cohort = "MARCH")

## list of uid
sample_uid <- c(res_md$filename,march_md$filename)
sample_uid <- paste(sample_uid, ".sam.bz2", sep = "")
writeLines(sample_uid, con = "deniz_samfiles.txt")

## list of MARCH uids
march_uid <- c(march_md$filename)
march_uid <- paste(march_uid, "_profile.tsv ./march_metaphlan_profiles", sep = "")
march_uid <- paste("mv ./", march_uid, sep = "")
writeLines(march_uid, con = "move_marchfiles.txt")

####Chi^2 Test of Breastfeeding by Eczema at the Whole Community Level####
# Make new breastfeeding vs. non breastfeeding column
march_wholemd <- march_wholemd %>% 
  mutate(cohort = "MARCH") %>%
  mutate(feedtype = ordered(feedtype, levels = c("BreastFed", "Mixed", "FormulaFed"))) %>%
  mutate(breastfeed_type = 
           case_when(feedtype == "Mixed" | feedtype == "FormulaFed" ~ "Non-BreastFed", 
                     feedtype == "BreastFed" ~ "BreastFed"))

res_wholemd <- res_wholemd %>% 
  mutate(cohort = "RESONANCE") %>%
  mutate(feedtype = ordered(feedtype, levels = c("BreastFed", "Mixed", "FormulaFed"))) %>%
  mutate(breastfeed_type = 
           case_when(feedtype == "Mixed" | feedtype == "FormulaFed" ~ "Non-BreastFed", 
                     feedtype == "BreastFed" ~ "BreastFed"))

## RESONANCE + MARCH
marchres_wholemd <- march_wholemd %>% rename (sample = SAMPLEID) %>%
  full_join(res_wholemd, .)

plot1 <- marchres_wholemd %>% filter(!is.na(eczema)) %>%
  filter(!is.na(feedtype)) %>%
  ggplot(aes(x=eczema,fill=feedtype)) +
  geom_bar(position = "dodge") +
  labs(x="Eczema", y= "No. of samples") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))

plot1

#breast vs mixed vs formula
marchres_feedec <- table(marchres_wholemd$feedtype, marchres_wholemd$eczema)
chisq.test(marchres_feedec)
CMHtest(marchres_feedec)

# breast vs non-breast
marchres_bfec <- table(marchres_wholemd$breastfeed_type, marchres_wholemd$eczema)
chisq.test(marchres_bfec)

## RESONANCE
#breast vs mixed vs formula
res_feedec <- table(res_wholemd$feedtype, res_wholemd$eczema)
chisq.test(res_feedec)
CMHtest(res_feedec)

# breast vs non-breast
res_bfec <- table(res_wholemd$breastfeed_type, res_wholemd$eczema)
chisq.test(res_bfec)

## MARCH
#feedtype
march_feedec <- table(march_wholemd$feedtype, march_wholemd$eczema)
chisq.test(march_feedec)
CMHtest(march_feedec)

# breastfed vs non breastfed
march_bfec <- table(march_wholemd$breastfeed_type, march_wholemd$eczema)
chisq.test(march_bfec)

####Relative abundance of HMO metabolizers in the two cohorts ####
res_hmomet <- res_relab %>% 
  select(-NCBI_tax_id) %>%
  pivot_longer(cols = -1, names_to = "uid", values_to = "relab") %>%
  separate(uid, into = c("uid", NA), sep = "_") %>%
  pivot_wider(names_from = "clade_name", values_from = "relab") %>%
  select(uid, contains(c("s__Bifidobacterium_bifidum", "s__Bifidobacterium_breve", 
                    "s__Bifidobacterium_longum", "s__Bifidobacterium_catenulatum",
                    "s__Bifidobacterium_animalis", "s__Bacteroides_fragilis", 
                    "s__Lactobacilus_casei", "s__Lactobacilus_acidophilus", 
                    "s__Escherichia_coli", "s__Klebsiella_pneumoniae", 
                    "s__Faecalibacterium_prausnitzii", "s__Ruminococcus_gnavus", 
                    "s__Akkermansia_muciniphila"))) %>%
  pivot_longer(-1, names_to = "species", values_to = "relab") %>%
  separate(species, into = c(NA, "species"), sep = "s__") %>%
  mutate(species = str_replace(species, "_", " ")) %>%
  left_join(res_md, .)

march_hmomet <- march_relab %>% 
  select(-NCBI_tax_id) %>%
  pivot_longer(cols = -1, names_to = "uid", values_to = "relab") %>%
  separate(uid, into = c("uid", NA), sep = "_") %>%
  pivot_wider(names_from = "clade_name", values_from = "relab") %>%
  select(uid, contains(c("s__Bifidobacterium_bifidum", "s__Bifidobacterium_breve", 
                         "s__Bifidobacterium_longum", "s__Bifidobacterium_catenulatum",
                         "s__Bifidobacterium_animalis", "s__Bacteroides_fragilis", 
                         "s__Lactobacilus_casei", "s__Lactobacilus_acidophilus", 
                         "s__Escherichia_coli", "s__Klebsiella_pneumoniae", 
                         "s__Faecalibacterium_prausnitzii", "s__Ruminococcus_gnavus", 
                         "s__Akkermansia_muciniphila"))) %>%
  pivot_longer(-1, names_to = "species", values_to = "relab") %>%
  separate(species, into = c(NA, "species"), sep = "s__") %>%
  mutate(species = str_replace(species, "_", " ")) %>%
  left_join(march_md, .) %>% rename(childAgeMonths = child_age_month, childGender = SEX)

marchres_hmomet <- full_join(march_hmomet, res_hmomet)
## the common columns between the two dataframes found during the join
# Joining with `by = join_by(feedtype, eczema, feed_exposure, childAgeMonths, childGender, uid, filename,
                           # S_well, cohort, species, relab)`


plot2 <- marchres_hmomet %>%
  filter(relab>0, !is.na(eczema)) %>%
  ggplot(aes(x = species, y=relab, fill = eczema)) +
  geom_boxplot() +
  geom_point(alpha = 0.5, size = 0.5, 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             dodge.width = 0.6, seed = 511)) +
  labs(x=NULL, y = "Relative Abundance") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, face = "italic"),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 15))



plot2
#### Bray Curtis Distance and PERMANOVA ####

## RESONANCE
# have samples as rows and taxa as column
res_trans <- res_ab %>% select(-NCBI_tax_id) %>% pivot_longer(cols = -1) %>%
  pivot_wider(names_from = "clade_name", values_from = "value") %>%
  select(name, contains("s__")) %>%
  mutate(seqname = str_replace_all(name, "_profile", "")) %>% select(-name) %>%
  select(seqname, everything())

# Calculate Bray-Curtis
bc_input_res <- res_trans %>%
  column_to_rownames("seqname")

# Bray-Curtis Distance calculation
set.seed(1994) # for reproducibility 
res_bc <- vegdist(bc_input_res, method = "bray")

## PERMANOVA
##distance matrix for adonis input
set.seed(1994)
adonis_input_res <- as.matrix(res_bc)
adonis_input_res <- as_tibble(adonis_input_res, rownames = "seqname") %>%
  left_join(., res_md, by = "seqname")

dist_res <- adonis_input_res %>%
  select(all_of(.[["seqname"]])) %>%
  as.dist()

perm_res1 <- adonis2(dist_res~eczema, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm_res2 <- adonis2(dist_res~child_age_month, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm_res3 <- adonis2(dist_res~feedtype, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm_res4 <- adonis2(dist_res~feedtype_breastfed, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm_res5 <- adonis2(dist_res~feedtype_mixed, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)
perm_res6 <- adonis2(dist_res~feedtype_formulafed, data = adonis_input_res, 
                     permutations = 9999, na.action = na.omit)

## MARCH
march_trans <- march_ab %>% select(-NCBI_tax_id) %>% pivot_longer(cols = -1) %>%
  pivot_wider(names_from = "clade_name", values_from = "value") %>% 
  select(name, contains("s__")) %>%
  mutate(seqname = str_replace_all(name, "_profile", "")) %>% select(-name) %>%
  select(seqname, everything())

# Calculate Bray-Curtis
bc_input_march <- march_trans %>%
  column_to_rownames("seqname")

# Bray-Curtis Distance calculation
set.seed(1994) # for reproducibility 
march_bc <- vegdist(bc_input_march, method = "bray")

## PERMANOVA
##distance matrix for adonis input
set.seed(1994)
adonis_input_march <- as.matrix(march_bc)
adonis_input_march <- as_tibble(adonis_input_march, rownames = "seqname") %>%
  left_join(., march_md, by = "seqname")

dist_march <- adonis_input_march %>%
  select(all_of(.[["seqname"]])) %>%
  as.dist()

perm_march1 <- adonis2(dist_march~eczema, data = adonis_input_march, 
                       permutations = 9999, na.action = na.omit)
perm_march2 <- adonis2(dist_march~child_age_month, data = adonis_input_march, 
                       permutations = 9999, na.action = na.omit)
perm_march3 <- adonis2(dist_march~feedtype, data = adonis_input_march, 
                     permutations = 9999, na.action = na.omit)
perm_march4 <- adonis2(dist_march~feedtype_breastfed, data = adonis_input_march, 
                     permutations = 9999, na.action = na.omit)
perm_march5 <- adonis2(dist_march~feedtype_mixed, data = adonis_input_march, 
                     permutations = 9999, na.action = na.omit)
perm_march6 <- adonis2(dist_march~feedtype_formulafed, data = adonis_input_march, 
                     permutations = 9999, na.action = na.omit)

### Log the R^2 and p values in an excel file to make the data frame
cohort <- c("MARCH", "MARCH", "MARCH", "MARCH", "MARCH", "RESONANCE", "RESONANCE", "RESONANCE", "RESONANCE", "RESONANCE")
variable <- c("Age", "Eczema", "Breast Fed", "Formula Fed", "Mixed Feed", "Age", "Eczema", "Breast Fed", "Formula Fed", "Mixed Feed")
R_sq <- c(0.00564, 0.00789, 0.03659, 0.00749, 0.02424, 0.03673, 0.11384, 0.02138, 0.0114, 0.01851)
pvalue <- c(0.1969, 0.046, 1e-04, 0.0592, 1e-04, 1e-04, 0.0264, 0.0042, 0.2354, 0.0161)
permanova_df <- tibble(cohort, variable, R_sq, pvalue)

permanova_barplot <- permanova_df %>%
  ggplot(aes(x = variable, y = R_sq, fill = cohort)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = "R^2") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))

permanova_barplot
#### Alpha Diversity ####
## RESONANCE
res_ad <- res_ab %>% select(-NCBI_tax_id) %>%
  pivot_longer(cols = -1) %>% 
  pivot_wider(names_from = "clade_name", values_from = "value") %>%
  select(name, contains("s__")) %>%
  rename(seqname = name)  %>%
  mutate(seqname = str_replace_all(seqname, "_profile", "")) %>%
  pivot_longer(cols = -1) %>%
  group_by(seqname) %>%
  summarize(shannon = diversity(value,index = "shannon"), 
            simpson = diversity(value, index = "simpson"),
            inverse_simpson = diversity(value, index = "invsimpson"))

res_admd <- left_join(res_md, res_ad, by = "seqname") %>%
  mutate(feedtype = factor(feedtype, 
                           levels = c("BreastFed", "Mixed", "FormulaFed")),
         eczema = factor(eczema, levels = c("Yes", "No")))

resad_plot <- res_admd %>% drop_na(eczema) %>% 
  ggplot(aes(x = eczema, y = shannon, fill = eczema)) + 
  geom_boxplot(width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
  scale_fill_manual(values=c("#e00201","#8faadc")) +
  xlab(label = "") +
  theme(legend.position = "none") +
  theme_classic()

res_plot2 <- march_admd %>% drop_na(feedtype) %>% drop_na(eczema) %>%
  ggplot(aes(x = feedtype, y = shannon, fill = eczema)) + 
  geom_boxplot() +
  geom_point(alpha = 0.5, size = 0.5, 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             dodge.width = 0.6, seed = 511)) +
  scale_fill_manual(values=c("#e00201","#8faadc"))+
  xlab(label = "") +
  theme_classic()
marchad_plot2
two.way.march <- aov(shannon ~ feedtype + eczema, data = march_admd)
summary(two.way.march)

## MARCH
march_ad <- march_ab %>% select(-NCBI_tax_id) %>%
  pivot_longer(cols = -1) %>% 
  pivot_wider(names_from = "clade_name", values_from = "value") %>%
  select(name, contains("s__")) %>%
  rename(seqname = name) %>%
  mutate(seqname = str_replace_all(seqname, "_profile", "")) %>%
  pivot_longer(cols = -1) %>%
  group_by(seqname) %>%
  summarize(shannon = diversity(value,index = "shannon"), 
            simpson = diversity(value, index = "simpson"),
            inverse_simpson = diversity(value, index = "invsimpson"))

march_admd <- left_join(march_md, march_ad, by = "seqname") %>%
  mutate(feedtype = factor(feedtype, 
                           levels = c("BreastFed", "Mixed", "FormulaFed")),
         eczema = factor(eczema, levels = c("Yes", "No")))

marchad_plot <- march_admd %>% drop_na(eczema) %>% 
  ggplot(aes(x = eczema, y = shannon, fill = eczema)) + 
  geom_boxplot(width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
  scale_fill_manual(values=c("#e00201","#8faadc")) +
  xlab(label = "") +
  theme(legend.position = "none") +
  theme_classic()

marchad_plot2 <- march_admd %>% drop_na(feedtype) %>% drop_na(eczema) %>%
  ggplot(aes(x = feedtype, y = shannon, fill = eczema)) + 
  geom_boxplot() +
  geom_point(alpha = 0.5, size = 0.5, 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             dodge.width = 0.6, seed = 511)) +
  scale_fill_manual(values=c("#e00201","#8faadc"))+
  xlab(label = "") +
  theme_classic()
marchad_plot2
two.way.march <- aov(shannon ~ feedtype + eczema, data = march_admd)
summary(two.way.march)

## Breast feeding alpha diversity by cohort

res_admd <- res_admd %>% mutate(cohort = "RESONANCE")
march_admd <- march_admd %>% mutate(cohort = "MARCH")

marchres_admd <- full_join(march_admd, res_admd)

marchresad_plot2 <- marchres_admd %>% drop_na(feedtype) %>%
  ggplot(aes(x = cohort, y = shannon, fill = feedtype)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#e00201","#4472c4","#92d051")) +
  geom_point(alpha = 0.5, size = 0.5, 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             dodge.width = 0.6, seed = 511)) +
  xlab(label = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18))

marchresad_plot2
two.way <- aov(shannon ~ feedtype + cohort, data = marchres_admd)
summary(two.way)

#### Make a Biplot ####

## table with relab counts per OTU
# In Pat's CC203 relab_otu = shared_tbl

shared_tbl <- rel_ab  %>%
  filter(grepl("s__", clade_name)) %>% pivot_longer(cols = -1) %>% 
  rename(samples = name) %>%
  pivot_wider(names_from = "clade_name", values_from = "value")

## Need a distance matrix for vegdist input
# In Pat's CC203 otu_matrix = shared_ 
shared_df <- shared_tbl %>%
  column_to_rownames("samples")

## Calculate the Bray-Curtis distance
set.seed(10)
mice_dist <- vegdist(shared_df, method = "bray")

## Need to calculate nmds ordination
set.seed(10)
nmds <- metaMDS(mice_dist)
nmds_positions <- scores(nmds) %>%
  as_tibble(rownames = "samples")

## generate subsample shared file
subsample_shared_tbl <- shared_tbl %>% pivot_longer(-samples)

#correlate OTUs with x and y axis positions
nmds_positions
subsample_shared_tbl

nmds_shared <- inner_join(subsample_shared_tbl, nmds_positions)

cor_x <- nmds_shared %>%
  nest(data = -name) %>%
  mutate(cor_x = 
           map(data, ~cor.test(.x$value, .x$NMDS1, method = "spearman", 
                               exact = FALSE) %>% tidy())) %>% 
  unnest(cor_x) %>% select(name, estimate, p.value)

cor_y <- nmds_shared %>%
  nest(data = -name) %>%
  mutate(cor_y = 
           map(data, ~cor.test(.x$value, .x$NMDS2, method = "spearman", 
                               exact = FALSE) %>% tidy())) %>% 
  unnest(cor_y) %>% select(name, estimate, p.value)

correlations <- inner_join(cor_x, cor_y, by = "name")
high_corr <- correlations %>%
  filter(p.value.x < 0.05 | p.value.y <0.05) %>%
  filter(abs(estimate.x > 0.3) | abs(estimate.y > 0.3))

# plot segments from (0,0) to (Rx, Ry)

high_corr %>%
  ggplot (aes(x = 0, xend = estimate.x, y = 0, yend = estimate.y)) +
  geom_segment()

vector <- c("R. gnavus", "Flavonifactor plautii", "C. innocuum", "B. uniformis", "H. hatheway")

high_corr$name <- vector


ggplot() +
  geom_point(data = nmds_positions, aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = high_corr, aes(x = 0, xend = estimate.x, y = 0, yend = estimate.y), inherit.aes = FALSE ) +
  #geom_text(data = high_corr, aes(x = estimate.x, y = estimate.y), label = high_corr$name) +
  geom_text_repel(data=high_corr, aes(x=estimate.x, y=estimate.y, label=name), min.segment.length = 0.01, 
                  xlim=c(-1.2, 0.6), inherit.aes=FALSE)


# plot text at (Rx, Ry)

nmds_positions %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point() +
  geom_segment(data = high_corr, aes(x = 0, xend = estimate.x, y = 0, yend = estimate.y), inherit.aes = FALSE ) +
  geom_text(data = high_corr, aes(x = estimate.x, y = estimate.y), inherit.aes = FALSE, label = name)
