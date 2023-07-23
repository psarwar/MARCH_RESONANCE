# Make one joint MARCH and RESONANCE metadata file
# This will correspond to all the demographics included in the full demographics supplement
# occasionally might need to clear environment for sanity
rm(list =ls())

# load required packages
library(tidyverse)

## Call in the two metadata files 
## generated using the Demographic_RES.R and Demographic_MARCH.R scricpts

resmd <- read_csv("./metadata_exports/seqdem_resmd.csv") %>% 
  rename(child_age_month = childAgeMonths, 
         sex = childGender, total_weight_lbs = childWeight, 
         bb_length = childHeight, eczema_timepoint = eczema, 
         maternal_education = mateduc_grp, race = simple_race,
         solids_started = feed_solid, pet_owned = pet_own) %>%
  rename(eczema = eczema_studyID) %>%
  mutate(childcare = case_when(daycare == "Daycare" ~ "social",
                               daycare == "Mixed" ~ "mixed",
                               daycare == "Private" ~ "private"),
         pet_dog = case_when(pet_type == "Dog" ~ "Yes"),
         pet_cat = case_when(pet_type == "Cat" ~ "Yes"),
         pet_other = case_when(pet_type == "Other" ~ "Yes")) %>%
    select(-c("childWei_kg", "eczema_timepoint", "feed_present", "feed_past", 
              "feed_past_studyid", "infant_race", "childHei_cm", 
              "childcare_private", "childcare_daycare", "daycare", "pet_type")) %>%
  mutate(cohort = "RES")

marchmd <- read_csv("./metadata_exports/seqdem_marchmd.csv") %>%
  rename(sample = SAMPLEID, sex = SEX, abx_labor = Abx_Labor, 
         delivery = DELIVERY, household_children = sibling) %>%
  mutate(race = case_when(MIXED_RACE == "Yes" ~ "Mixed",
                          race == "WHITE" ~ "White",
                          race == "BLACK" ~ "Black",
                          race == "OTHER" ~ "Other",
                          race == "ASIAN" ~ "Asian")) %>%
  mutate(cohort = "MARCH") %>% 
  select(-c("DELIVERY_EXT", "EDUC_LVL", "child_age_day", "MIXED_RACE", 
            "total_weight_grams", "baby_probiotic", "baby_prebiotic", 
            "baby_multivitamin", "mother_id"))

## Call in the different ID/Crosswalk files
newIDs <- read_csv("./raw_data_mixed/old_new_seqID_mgx.csv") %>%
  mutate(old_filename = paste(biospecimen,S_well, sep = "_")) %>%
  rename(new_filename = filename) %>%
  select(uid, biospecimen, S_well, new_filename, old_filename)

all_mgx <- read_csv("./raw_data_res/Prioty_allmgx.csv") %>%
  select(sample, sid_old)

#Final table joins
marchres_md <- full_join(marchmd, resmd)%>%
  left_join(. , all_mgx) %>%
  separate(seqname, into = c("biospecimen", NA), sep = "_S") %>%
  mutate(biospecimen = case_when( !is.na(biospecimen) ~ biospecimen,
                                  .default = sample),
         maternal_education = case_match(maternal_education, 
                                         "college cerdit, no degree" ~ "college, no degree",
                                         "masters degree" ~ "masters",
                                         "high school diploda/GED" ~ "high school/GED",
                                         "associate degree" ~ "associates",
                                         "high school, no diploma" ~ "high school, no degree",
                                         "bachelors degree" ~ "bachelors",
                                         "doctorate/professional degree" ~ "professional/doctorate",
                                         "masters degree" ~ "masters",
                                         "trade/technical/vocational training" ~ "partial college/special training",
                                         .default = maternal_education)) %>%
  select(-c("specimen_ID", "mother_id")) %>%
  left_join(. , newIDs) %>%
  distinct(sample, .keep_all = TRUE)

write.csv(marchres_md, file = "./metadata_exports/marchres_combinedmd_ids.csv")
