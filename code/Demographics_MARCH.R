# Demographic Data Scavenging for the MARCH dataset
## contains the calculation of breastfeeding and eczema + other parameters
## at the whole community and sequencing levels

## look at files in the working directory
list.files()
## clear the R environment
rm(list = ls())

library(tidyverse)
library(readxl)

#### Import the raw data/ manipulated input data ####

# Call in all the raw files

# this table has all the 3mo samples 
# for joins: 'mother_id' like 'P80004' ; 'SAMPLEID' in the format '80004M1'
metadata_3mo <- read.csv("./raw_data_march/Sarah_R15_three_month_data06162022.csv")

# for joins:`Study ID` in the format '25FE8663'
sif_data <-  read.csv("./raw_data_march/MARCH_3mos_SIF_2023MAR12.csv")
## Format the SIF table for ease of downstream use
colnames(sif_data) <- sif_data[2,] # Assign the second row as the colnames of the dataset
sif_data <- sif_data[-(1:3),] # Delete all the irrelevant rows
# rename the columns that is needed
sif_data <- sif_data[, -c(28,29,31,32)] %>% #remove the duplicate columns interfering with rename. not relevant to data of interest
  rename(specimen_ID = `Study ID`, bb_weight = `What does baby weigh now?`,
                                bb_length = `How long is baby now?`, bb_med_now = 
                                  `Is baby taking medicine(s) now?`,
                                bb_med_ever = 
                                  `If the baby is not taking medicine now, has baby taken any medicine(S) since birth?`,
                                bb_abx_birth = `Has baby had antibiotics since birth?`,
                                feed_exposure = `During the past week, my baby ate:`,
                                diet_cowmilk = `Select anything that baby ate in the past 24 hours: - Selected Choice - Cow's milk (not infant formula)`,
                                diet_other_milk = `Select anything that baby ate in the past 24 hours: - Selected Choice - Other milk: soy, rice, goat, etc. (not infant formula)`,
                                diet_nonmilk_dairy = `Select anything that baby ate in the past 24 hours: - Selected Choice - Other dairy foods: yogurt, cheese, ice cream, pudding, etc.`,
                                diet_soy_food = `Select anything that baby ate in the past 24 hours: - Selected Choice - Other soy foods: tofu, frozen soy desserts, etc.`,
                                diet_cereal = `Select anything that baby ate in the past 24 hours: - Selected Choice - Baby cereal`,
                                diet_grains = `Select anything that baby ate in the past 24 hours: - Selected Choice - Other grains: cheerios, other breakfast cereals, teething biscuits, crackers, breads, pasta, rice, etc.`,
                                diet_fruits = `Select anything that baby ate in the past 24 hours: - Selected Choice - Carrots, sweet potatoes, mangos, apricots, bell peppers`,
                                diet_veg = `Select anything that baby ate in the past 24 hours: - Selected Choice - Spinach, kale, Swiss chard, romaine lettuce`,
                                diet_other_fruits = `Select anything that baby ate in the past 24 hours: - Selected Choice - Other fruit or fruit purees`,
                                diet_other_veg = `Select anything that baby ate in the past 24 hours: - Selected Choice - Other vegetables or vegetable purees`,
                                diet_fries =`Select anything that baby ate in the past 24 hours: - Selected Choice - French fries`,
                                diet_meat =`Select anything that baby ate in the past 24 hours: - Selected Choice - Meat: chicken, beef, ham, combination dinners`,
                                diet_fish = `Select anything that baby ate in the past 24 hours: - Selected Choice - Fish or shellfish`,
                                diet_nut_products = `Select anything that baby ate in the past 24 hours: - Selected Choice - Peanut butter, other peanut foods (Bamba), other nuts`,
                                diet_liver = `Select anything that baby ate in the past 24 hours: - Selected Choice - Liver or other organ meats`,
                                diet_eggs = `Select anything that baby ate in the past 24 hours: - Selected Choice - Eggs`,
                                diet_legumes = `Select anything that baby ate in the past 24 hours: - Selected Choice - Beans, lentils, peas`,
                                diet_sweets = `Select anything that baby ate in the past 24 hours: - Selected Choice - Sweet foods: candy, chocolate, cookies, cakes, etc.`,
                                baby_prebiotic = `Select anything that baby ate in the past 24 hours: - Selected Choice - Prebiotic supplement (Gos, Fos, inulin, beta-glucan, etc.)`,
                                baby_probiotic = `Select anything that baby ate in the past 24 hours: - Selected Choice - Probiotic supplement, kefir, kimchi`,
                                baby_multivitamin = `Select anything that baby ate in the past 24 hours: - Selected Choice - Multi-vitamin supplement`)

# for joins: 'SAMPLEID' in the format 'P85821'; 'SPECIMEN_ID' in the format of '8600'
prenatal_data <- read.csv("./raw_data_march/Sarah_R15_prenatal_data06162022.csv")

# for joins: 'SAMPLEID' in the format 'P89094'; 
mod_data <- read.csv("./raw_data_march/MARCH_MOD_data_2023March10_SSC.csv")

# for joins: 'specimen_ID' in the format '25FE7014'; 'child_id' in the format '81283M1'
sample_collection_data <- read.csv("./raw_data_march/sample_collection_date_IMR3M.csv")

# for joins: 'SAMPLEID' in the format 'P80437'
birth_cert_data <- read.csv("./raw_data_march/Sarah_R15_birth_certificate_data06162022.csv")

# This table merged all outputs from metaphlan using a metaphlan subscript
# merge_metaphlan_tables.py *_profile.tsv > MARCH_merged_ab_table.tsv
# for joins: the column names correspond to the sequenced sample names 
# that correspond to the file_index in the crosswalk file
rel_ab <- read_tsv("./raw_data_march/march_mergedab.tsv", skip = 1)

## get rid of the columns with the samples sequenced in duplicate or unknown
## full reasoning outlined for the selection/deletion outlined 
# in CommunityAnalysis_MARCH.R subheading 'Resolve Tengfei's QC issue'

rel_ab <- rel_ab %>% select(-c("144_3mos_S48","166_3mos_S70",
                               "193-3mos_S1","71-3mos_S81" ))

# this table only contains the data for the sequenced samples
# for joins: 'specimen_ID' in the format '25FE8580'; 'mother_id' in the format of 'P80004'
# 'SAMPLEID' in the format of '80004M1'
# file index contains the name of the whole read file which can be modified for join with the merged abundance table
crosswalk_file <- read.csv("./raw_data_march/IMR_WGS_QC_TMA.csv")

# for joins: 'specimen_ID' in the format '25FE8580'; child_id in the format
# contains infant ages for sequenced samples
infant_ages <- read.csv("./raw_data_march/sample_collection_date_IMR3M_TMA04212023.csv")

#### Calculate the Breastfeeding and Eczema ####

# Breastfeeding and Eczema for the whole community 
## BreastFeeding Conditionals 
A <- metadata_3mo$BRST_MILK == 1 | metadata_3mo$FORMULA != 1
## Formula Conditionals
B <- metadata_3mo$FORMULA == 1 | metadata_3mo$BRST_MILK != 1

metadata_bf_ec <- metadata_3mo %>%
  mutate(feedtype = case_when(A == TRUE & (B != TRUE | is.na(B)) ~ "BreastFed",
                              (A != TRUE | is.na(A)) & B == TRUE ~ "FormulaFed",
                              A == TRUE & B == TRUE ~ "Mixed"),
         eczema = case_when(ECZEMA_CHILD == 1 ~ "Yes",
                            ECZEMA_CHILD == 5 ~ "No"))

## Eczema and breastfeeding for the sequenced community
seq_samplenames <- colnames(rel_ab)
seq_samplenames <-  str_remove(seq_samplenames, "_profile")
#get rid of the column name 'clade_name'
seq_samplenames <- seq_samplenames[!seq_samplenames == 'clade_name']
seq_samplenames <- seq_samplenames[!seq_samplenames == 'NCBI_tax_id']

metadata_seq <- crosswalk_file %>%
  mutate(seqname = str_remove(file_index, "_L001_R1_001.fastq.gz")) %>%
  filter(seqname %in% seq_samplenames) %>%
  inner_join(., metadata_bf_ec, by = c("SAMPLEID", "mother_id")) %>%
  select(SAMPLEID, seqname, specimen_ID, mother_id, feedtype, eczema)

# BreastFeeding Exposure
# Can make whole community counts directly from sif_data$feed_exposure
# Sequenced
metadata_seq <- sif_data %>% select(specimen_ID, feed_exposure) %>% 
  left_join(metadata_seq,., by = "specimen_ID") %>% 
  distinct(SAMPLEID, .keep_all=TRUE) 

# Solid Food Exposure
## Whole Community - sif
bb_diet <- sif_data %>% select(specimen_ID, starts_with("diet_")) %>% 
  select(-c("diet_cowmilk", "diet_other_milk")) %>%
  unite(col = 'diet_solid', c("diet_nonmilk_dairy", "diet_soy_food", 
                              "diet_cereal","diet_grains","diet_fruits", 
                              "diet_veg","diet_other_fruits","diet_other_veg",
                              "diet_fries","diet_meat","diet_fish",
                              "diet_nut_products","diet_liver","diet_eggs",
                              "diet_legumes","diet_sweets"), sep = '_') %>%
  mutate(solids_started = if_else(diet_solid == "_______________", "No", "Yes"))

## Sequenced
metadata_seq <- bb_diet %>% select(specimen_ID, solids_started) %>%
  left_join(metadata_seq, ., by = "specimen_ID") %>%
  distinct(SAMPLEID, .keep_all=TRUE) 


#### General Metadata Calculations - Sex, Age and Race ####

#Age
##Whole Community Data Not Available
##Sequenced
# Original Plan was to calculate age from birth date to sample collection date
# but birth date cannot be shared with me
# The age at time of sampling was calculated by Tengfei [MARCH personnel] and sent to me
metadata_seq <- infant_ages %>% rename(SAMPLEID = child_id) %>% 
  mutate(child_age_month = round(child_age_day/30.417, digit = 2)) %>%
  left_join(metadata_seq, ., by = c("SAMPLEID", "specimen_ID")) %>%
  select(-collected_date2)

# Sex
## Whole Community
sex <- birth_cert_data %>% 
  rename(mother_id = SAMPLEID) %>%
  select(SEX, mother_id) %>%
  left_join(metadata_3mo, ., by = "mother_id") %>%
  mutate(SEX = case_when(SEX == 1 ~ "Male",
                         SEX == 2 ~ "Female")) %>% 
  distinct(SAMPLEID, .keep_all=TRUE) 
# duplicate rows are being caused by the same mother_id having multiple births
#checked using unique() that all SAMPLEIDs are unique

## Sequenced
metadata_seq <- sex %>%
  select(SAMPLEID, mother_id, SEX) %>%
  left_join(metadata_seq, ., by = c("SAMPLEID", "mother_id"))


# RACE
## Whole Community

# Find out the mixed race infants
race <- metadata_3mo %>%
  select("SAMPLEID", starts_with("BABY_RACE")) %>%
  mutate(RACE_NUMBER = rowSums(.[2:9])) %>% # Add the binary coded racial columns
  mutate(MIXED_RACE = ifelse(RACE_NUMBER > 1, "Yes","No"))  # RACE_NUMBER >1 indicates more than one race column was checked 

# Find out the race of the infants after mixed race is removed
non_mixedrace <- race %>%
  select(-c("BABY_RACE_DK", "BABY_RACE_RF", "BABY_RACE_SPEC", "RACE_NUMBER")) %>%
  pivot_longer(-c("SAMPLEID", "MIXED_RACE")) %>% 
  separate(name, into = c(NA, "race"), sep = "RACE_") %>%
  filter (value == 1) %>% # only keep the rows where the race rows have been selected
  filter(MIXED_RACE == "No") # but get rid of the mixed race infants before counts

## Sequenced
# mixed race info
metadata_seq <- race %>% select(SAMPLEID, MIXED_RACE) %>%
  left_join(metadata_seq, ., by = "SAMPLEID")

# race information
metadata_seq <- non_mixedrace %>% select(SAMPLEID, race) %>%
  left_join(metadata_seq, ., by = "SAMPLEID")


#### Mother/Pregnancy Metadata Calculations ####

# Abx and steroids during pregnancy
## Whole Community
preg_meds <- prenatal_data %>%
  rename(mother_id = SAMPLEID) %>%
  mutate(abx_preg = case_when(PREGANTIB == 1 ~ "Yes",
                                   PREGANTIB == 5 ~ "No",
                                   PREGANTIB == 98 ~ NA_character_)) %>%
  mutate(steroid_preg = case_when(STEROID_PREG == 1 ~ "Yes",
                                  STEROID_PREG == 5 ~ "No")) %>%
  left_join(metadata_3mo, ., by = "mother_id") %>%
  select(SAMPLEID, mother_id, abx_preg, steroid_preg)
## Sequenced Community
metadata_seq <- preg_meds %>%
  left_join(metadata_seq, ., by = c("SAMPLEID","mother_id"))

# Abx during labor
## Whole Community
ld_abx <- birth_cert_data %>% 
  rename(mother_id = SAMPLEID) %>%
  select(mother_id, LD_ANTIBIOTICS) %>%
  mutate(Abx_Labor = case_when(LD_ANTIBIOTICS == 1 ~ "Yes",
                               LD_ANTIBIOTICS == 2 ~ "No")) %>%
  left_join(metadata_3mo, ., by = "mother_id") %>% 
  distinct(SAMPLEID, .keep_all=TRUE) #checked using unique() that all SAMPLEIDs are unique
## Sequenced
metadata_seq <- ld_abx %>%
  select(SAMPLEID, mother_id, Abx_Labor) %>%
  left_join(metadata_seq, ., by = c("SAMPLEID", "mother_id"))

# Delivery Mode
## Whole Community
delivery <- mod_data %>% rename(mother_id = SAMPLEID) %>% 
  mutate(DELIVERY = case_when(MD_FINAL_ROUTE == 1 | MD_FINAL_ROUTE == 2 
                              | MD_FINAL_ROUTE == 3 ~ "Vaginal",
                              MD_FINAL_ROUTE == 4 ~ "Cesarean"),
         DELIVERY_EXT = case_when(MD_FINAL_ROUTE == 1 ~ "Vaginal-spontaneous",
                                  MD_FINAL_ROUTE == 2 ~ "Vaginal-forceps",
                                  MD_FINAL_ROUTE == 3 ~ "Vaginal-vacuum",
                                  MD_FINAL_ROUTE == 4 ~ "Cesarean")) %>%
  left_join(metadata_3mo, ., by = "mother_id") %>%
  select(SAMPLEID, mother_id, DELIVERY, DELIVERY_EXT)
## Seqeunced
metadata_seq <- delivery %>% left_join(metadata_seq, ., by = c("SAMPLEID", "mother_id"))

# Maternal Education
## Whole Community
maternal_educ <- prenatal_data %>% select(SAMPLEID,EDUC_LVL) %>% 
  rename(mother_id = SAMPLEID) %>%
  mutate(maternal_education = case_when(EDUC_LVL == 2 ~ '=<8th Grade',
                                        EDUC_LVL == 3 ~ 'high school, no diploma',
                                        EDUC_LVL == 4 ~ 'high school diploda/GED',
                                        EDUC_LVL == 5 ~ 'college cerdit, no degree',
                                        EDUC_LVL == 6 ~ 'trade/technical/vocational training',
                                        EDUC_LVL == 7 ~ 'associate degree',
                                        EDUC_LVL == 8 ~ 'bachelors degree',
                                        EDUC_LVL == 9 ~ 'masters degree',
                                        EDUC_LVL == 10 ~ 'doctorate/professional degree')) %>%
  left_join(metadata_3mo, ., by = "mother_id")

## Sequenced 
metadata_seq <- maternal_educ %>% select(SAMPLEID, mother_id, EDUC_LVL, maternal_education) %>%
  left_join(metadata_seq, ., by = c("SAMPLEID", "mother_id"))

#### Infant Metadata Calculations ####
# Baby Length
## Whole Communtiy
# transform the data to numeric
bb_length_md <- sif_data %>% select(specimen_ID, bb_length) %>%
  mutate(bb_length = gsub('1/2', '.5', bb_length)) %>%
  mutate(bb_length = gsub('1/4', '.25', bb_length)) %>%
  mutate(bb_length = gsub('[^[:digit:].]','', bb_length))
# checked the row for manually for unusual values and deleting them
# picked out rows 48, 92, 310, 428
# row numbers are -3 for some reason [b/c sif_data has top 3 rows removed during import]
bb_length_md <- bb_length_md[-c(45,89,307,425),]
bb_length_md$bb_length <- as.numeric(bb_length_md$bb_length)
## Sequenced
test <- bb_length_md %>% left_join(metadata_seq, . , by = "specimen_ID") %>%
  distinct(SAMPLEID, .keep_all=TRUE) # duplicate specimen_IDs for some

# Baby Weight
## Whole Community
# transform the data to numeric
bb_weight_md <- sif_data %>% select(specimen_ID, bb_weight) %>%
  mutate(bb_weight = gsub('lb|lb.|lbs.|pounds|POUNDS|LB|LBS|LBS.|lbs,', 'lbs', bb_weight)) %>%
  mutate(bb_weight = gsub('ozs|ozs.|oz.|OZ|ounces|OUNCES', 'oz', bb_weight)) %>%
  separate(bb_weight, into = c("weight_lbs", "weight_oz"), sep = "lbs" ) %>%
  separate(weight_oz, into = c("weight_oz", NA), sep = "oz") %>%
  mutate(weight_lbs = gsub('[^[:digit:].]','', weight_lbs),
         weight_oz = gsub('[^[:digit:].]','', weight_oz))
bb_weight_md$weight_oz <- as.numeric(bb_weight_md$weight_oz)
bb_weight_md$weight_lbs <- as.numeric(bb_weight_md$weight_lbs)
# Checked manually for unusual numbers same as for length
# only 92 identified. similarly to length need to delete row that is 92-3 = 89
bb_weight_md <- bb_weight_md[-89,]
# total weight calculations
bb_weight_md <- bb_weight_md %>% mutate(total_weight_lbs = weight_lbs + weight_oz/16) %>%
  mutate(total_weight_grams = total_weight_lbs*453.592)
## Sequenced
metadata_seq <- bb_weight_md %>% select(specimen_ID, total_weight_lbs, total_weight_grams) %>%
  left_join(metadata_seq, ., by = "specimen_ID") %>%
  distinct(SAMPLEID, .keep_all=TRUE) # duplicate specimen_IDs for some

# Baby Meds - Any Meds Now, Antibiotics Prebiotic, Probiotic, Multivitamin
## Whole Community -Sif
bb_meds <- sif_data %>% 
  select(specimen_ID, baby_prebiotic, baby_probiotic, baby_multivitamin, 
         bb_med_ever, bb_abx_birth) %>%
  mutate(baby_prebiotic = case_when(baby_prebiotic != "" ~ "Yes"),
         baby_probiotic = case_when(baby_probiotic != "" ~ "Yes"),
         baby_multivitamin = case_when(baby_multivitamin != "" ~ "Yes"))
## Sequenced
metadata_seq <- bb_meds %>% left_join(metadata_seq, ., by = "specimen_ID") %>%
  distinct(SAMPLEID, .keep_all=TRUE)


#### Social Metadata Calculations ####
# Childcare
## Whole Community
# private care conditional
private = metadata_3mo$CHILDCARE_PLACE_BABYHOME + metadata_3mo$CHILDCARE_PLACE_OTHPRIVHOMENOCH
# social care conditional
social = metadata_3mo$CHILDCARE_PLACE_DAYCARE + metadata_3mo$CHILDCARE_PLACE_BABYHOMEOTHCHILD +
  metadata_3mo$CHILDCARE_PLACE_OTHPRIVHOMEOTHCH

childcare <- metadata_3mo %>%
  select("SAMPLEID", "mother_id", starts_with("CHILDCARE_PLACE_")) %>%
  mutate(childcare = case_when(private > 0 & social == 0 ~ "private",
                               private == 0 & social > 0 ~ "social",
                               private > 0 & social > 0 ~ "mixed"))
## Sequenced
metadata_seq <- childcare %>% select(SAMPLEID, mother_id, childcare) %>%
  left_join(metadata_seq, ., by = c("SAMPLEID", "mother_id"))

# Pets
## Whole Community
pets <- metadata_3mo %>% select("SAMPLEID", starts_with("PETS"), -PETS_SPEC) %>%
  mutate(pet_owned = case_when(PETS == 1 ~ "Yes", PETS == 5 ~ "No"),
         pet_cat = case_when(PETS_TYPE_CATS == 1 ~ "Yes"),
         pet_dog = case_when(PETS_TYPE_DOGS == 1 ~ "Yes"),
         pet_other = case_when(PETS_TYPE_OTHER == 1 ~ "Yes"))

metadata_seq <- pets %>% select(SAMPLEID, pet_owned, pet_cat, pet_dog, pet_other) %>%
  left_join(metadata_seq, ., by = "SAMPLEID")

# Sibling
## Calculated from the question to mother: how many previous births are still living
sibling <- birth_cert_data %>% select(SAMPLEID, NOWLIVING) %>% 
  rename(mother_id = SAMPLEID) %>%
  mutate(sibling = ifelse(NOWLIVING > 0, "Yes", "No")) %>%
  left_join(metadata_3mo, ., by = "mother_id") %>%
  distinct(SAMPLEID, .keep_all=TRUE) # mothers have more than one child
## Sequenced
metadata_seq <- sibling %>% select(SAMPLEID, mother_id, sibling) %>% 
  left_join(metadata_seq, ., by = c("SAMPLEID", "mother_id")) 

#### General Metadata Calculations - Miscellaneous ####

# Calculate the age of the baby when milk [direct breast or expressed] was started and stopped
metadata_bf <- metadata_3mo %>% 
  mutate(BRST_BBSTART_AGE_MONTHS = case_when(BRST_BBSTART_AGE_UNIT == 1 ~ round(BRST_BBSTART_AGE/30.417, digit = 2),
                                             BRST_BBSTART_AGE_UNIT == 2 ~ round(BRST_BBSTART_AGE/4.345, digit =2),
                                             BRST_BBSTART_AGE_UNIT == 3 ~ BRST_BBSTART_AGE), 
         BRST_BBSTOP_AGE_MONTHS = case_when(BRST_BBSTOP_AGE_UNIT == 1 ~ round(BRST_BBSTOP_AGE/30.417, digit = 2),
                                            BRST_BBSTOP_AGE_UNIT == 2 ~ round(BRST_BBSTOP_AGE/4.345, digit = 2),
                                            BRST_BBSTOP_AGE_UNIT == 3 ~ BRST_BBSTOP_AGE),
         EXPMILK_BBSTART_AGE_MONTHS = case_when(EXPMILK_BBSTART_AGE_UNIT == 1 ~ round(EXPMILK_BBSTART_AGE/30.417, digit = 2),
                                                EXPMILK_BBSTART_AGE_UNIT == 2 ~ round(EXPMILK_BBSTART_AGE/4.345, digit = 2),
                                                EXPMILK_BBSTART_AGE_UNIT == 3 ~ EXPMILK_BBSTART_AGE),
         EXPMILK_BBSTOP_AGE_MONTHS = case_when(EXPMILK_BBSTOP_AGE_UNIT == 1~ round(EXPMILK_BBSTOP_AGE/30.417, digit = 2),
                                               EXPMILK_BBSTOP_AGE_UNIT == 2 ~ round(EXPMILK_BBSTOP_AGE/4.345, digit = 2),
                                               EXPMILK_BBSTOP_AGE_UNIT == 3 ~ EXPMILK_BBSTOP_AGE))


#### Counts ####
## % calculated in excel
## counts are calculated separately here so I can keep track of the tables 
## from which they were counted

# Sex
## Whole Community
sex %>% count(SEX)
## Sequenced
metadata_seq %>% count(SEX)

# Age
# Age in Days
metadata_seq %>% filter(!is.na(child_age_day)) %>% pull(child_age_day) %>% mean()
# Age in Months
metadata_seq %>% filter(!is.na(child_age_day)) %>% pull(child_age_month) %>% mean()
# 0-3 Months : No kids this age
metadata_seq %>% filter(child_age_month < 3) %>% pull(child_age_month) 
# 3-6 Months (n = 219)
metadata_seq %>% filter(child_age_month > 3 & child_age_month < 6) %>% 
  pull(child_age_month) %>% length()
# 6-9 Months (n=11)
# one kid is 12.03 months old but just included them with the 6-12 month olds
metadata_seq %>% filter(child_age_month > 6) %>% pull(child_age_month) %>% length()

# Race
## Whole Community
### Calculate those people who are mixed race and mixed 
race %>% count(MIXED_RACE) # Yes/No Mixed Race
non_mixedrace %>% count(race) # Count the race where mixed race has been already excluded
## Sequenced
metadata_seq %>% count(MIXED_RACE)
metadata_seq %>% filter(MIXED_RACE == "No") %>% #include only non-mixed race kids
  count(race)

# Eczema
## Whole Community
metadata_bf_ec %>% count(eczema)
## Sequenced
metadata_seq %>% count(eczema)

# BreastFeeding
## Whole Community
metadata_bf_ec %>% count(feedtype)
## Sequenced
metadata_seq %>% count(feedtype)

# Feed Exposure
## Whole Community (directly from sif_data)
sif_data %>% select(specimen_ID, feed_exposure) %>% count(feed_exposure)
## Sequenced
metadata_seq %>% count(feed_exposure)

#Baby Ate Solids in the last 24 hours
## Whole Community
bb_diet %>% count(solids_started)
## Seqeunced
metadata_seq %>% count(solids_started)

# Antibiotics During Pregnancy
## Whole Community
preg_meds %>% count(abx_preg)
## Seqeunced
metadata_seq %>% count(abx_preg)

# Steroid During Pregnancy
## Whole Community
preg_meds %>% count(steroid_preg)
## Seqeunced
metadata_seq %>% count(steroid_preg)

# Abx During Labor
## Whole Community
ld_abx %>% count(Abx_Labor)
## Sequenced
metadata_seq %>% count(Abx_Labor)

# Delivery Mode
## Whole Community
delivery %>% count(DELIVERY)
delivery %>% count(DELIVERY_EXT)
## Sequenced
metadata_seq %>% count(DELIVERY)
metadata_seq %>% count(DELIVERY_EXT)

# Mean baby length
## Whole Community
vec_length <- bb_length_md %>% filter(!is.na(bb_length)) %>% pull(bb_length) 
length(vec_length)
mean(vec_length)
## Sequenced
vec_length_seq <- metadata_seq %>% filter(!is.na(bb_length)) %>% pull(bb_length) 
length(vec_length_seq)
mean(vec_length_seq)

# Mean baby weight
## Whole Community
vec_weight_grams <- bb_weight_md %>% filter(!is.na(total_weight_grams)) %>% pull(total_weight_grams) 
length(vec_weight_grams)
mean(vec_weight_grams) # in grams
vec_weight_lbs <- bb_weight_md %>% filter(!is.na(total_weight_lbs)) %>% pull(total_weight_lbs) 
length(vec_weight_lbs)
mean(vec_weight_lbs) # in pounds
## Sequenced
vec_weight_grams_seq <- metadata_seq %>% filter(!is.na(total_weight_grams)) %>% 
  pull(total_weight_grams) 
length(vec_weight_grams_seq)
mean(vec_weight_grams_seq) # in grams
vec_weight_lbs_seq <- metadata_seq %>% filter(!is.na(total_weight_lbs)) %>% pull(total_weight_lbs) 
length(vec_weight_lbs_seq)
mean(vec_weight_lbs_seq) # in pounds

# Baby Meds
## Whole Community
bb_meds %>% count(bb_med_ever)
bb_meds %>% count(bb_abx_birth)
bb_meds %>% count(baby_prebiotic)
bb_meds %>% count(baby_probiotic)
bb_meds %>% count(baby_multivitamin)
##Sequenced
metadata_seq %>% count(bb_med_ever)
metadata_seq %>% count(bb_abx_birth)
metadata_seq %>% count(baby_prebiotic)
metadata_seq %>% count(baby_probiotic)
metadata_seq %>% count(baby_multivitamin)

# Childcare
## Whole Community
childcare %>% count(childcare)
## Seqeunced
metadata_seq %>% count(childcare)

# Pets
## Whole Community
pets %>% count(pet_owned)
pets %>% count(pet_cat)
pets %>% count(pet_dog)
pets %>% count(pet_other)
##Sequenced
metadata_seq %>% count(pet_owned)
metadata_seq %>% count(pet_cat)
metadata_seq %>% count(pet_dog)
metadata_seq %>% count(pet_other)

# Siblings
## Whole Community
sibling %>% count(sibling)
## Seqeunced
metadata_seq %>% count(sibling)

#Maternal Education
## Whole Community
maternal_educ %>% count(maternal_education)
## Sequenced
metadata_seq %>% count(maternal_education)

#### Master Metadata File/ Other types of metadata file ####
## make one metadata file with relevant columns for ease of publishing raw data
# metadata_seq is the master metadata table
# it has been made with multiple joins as each metadata was collected
# march had too many raw data sheets and ID types to make one master sheet first

write_csv(metadata_seq, file = "./processed/march_seq_metadata_allvar.csv")
