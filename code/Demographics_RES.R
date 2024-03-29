# Demographic Data Scavenging for the RESONANCE dataset
## contains the calculation of breastfeeding and eczema + other parameters
## at the whole community and sequencing levels cross-sectionally for infants < 1yr

# occasionally might need to clear environment for sanity
rm(list =ls())

# load required packages
library(tidyverse)

# RESONANCE Metadata is very long
# In addition to only calling in subsets of Jen's full export
# Need to increasr the size of the connection buffer
# connection size found by trial error just increase number until read.csv worked

Sys.setenv(VROOM_CONNECTION_SIZE=2000073)

#### Major Raw/Slightly Modified Data Inputs for Downstream Analysis ####
# Call in the full RESONANCE metadata Jan 2023 export from Jen
# change the long column names for easier downstream data wrangling
# these are the airway columns that relate to eczema
# calculate child age in months from the column age of months and days on the assessment
metadata <- read.csv("./raw_data_res/Prioty_Data_013023.csv") %>%
  mutate(studyID = str_pad(.$studyID, 4, pad = "0")) %>% #easier format of studyIDs for joins
  mutate(childAgeMonths = assessmentAgeMonths + round(assessmentAgeDays/30.417, digits = 0)) %>% # calculate child age in months
  select_all(~gsub("RedCap_Ess_CPH_Air_Inf..|Redcap_Ess_CPH_Air_EC..|Redcap_Ess_CPH_Air_MC..|Redcap_Ess_CPH_Air_Adol..|BasicFamilyAndChild..", "", .))

## Find the cross-sectional sequenced community infants < 1 year
# This table connects the sampleIDs tot he studyID and timepoint
sequenced <- read.csv("./raw_data_res/Prioty_allmgx.csv") %>%
  rename(studyID = subject) %>% #standardize the name of studyID according to Jen Export
  mutate(studyID = str_pad(.$studyID, 4, pad = "0")) %>% # make the study ID a 4 character string
  filter(childAgeMonths < 13) # only infant samples

# Call in filenames of all the RESONANCE metaphlan profile outputs
# This table was generated by saving all the RES metaphlan output filenames using ls > filenames.txt
filenames <- read_csv("./raw_data_res/filenames.csv", col_names = FALSE) %>%
  separate(X1, into = c("sample", NA), sep = "_" , remove = FALSE) 

# vector of all sequenced RESONANCE samples
sample_list <- filenames$sample

# remove the longitudinal samples
# this table with 151 samples are the infant cross-sectional samples that have been sequenced
seq_cross <- sequenced %>% 
  filter(sample %in% sample_list) %>% # select only the files for which we have metaphlan profiles
  filter(Fecal_EtOH == "F") %>% # don't include ethanol samples
  group_by(studyID) %>% arrange(timepoint, .by_group = TRUE) %>% # early to latest timepoint
  distinct(studyID, .keep_all = TRUE) %>% ungroup() # only keeps the first entry ie earliest timepoint

#cross-sectional sequenced metadata for <1 year infants
seq_metadata <- seq_cross %>% select(sample, studyID, timepoint) %>%
  left_join(., metadata, by = c("studyID","timepoint"))

## Need a cross-sectional version of the metadata

# leave out later time points per studyID to calculate 
# the problem with this table is the first timepoint collected =/= first timepoint sequenced
# therefore this is not the accurate version of the whole community
metadata_1yr <- metadata %>%
  filter(childAgeMonths < 13) %>% 
  group_by(studyID) %>% arrange(timepoint, .by_group = TRUE) %>%
  distinct(studyID, .keep_all = TRUE) %>% ungroup() # get's rid of all longitudinal data 

# Define a cross-sectional whole community that includes only one timepoint per studyID
# But that timepoint should be the timepoint sequenced
# for unsequenced studyIDs it should be the first timepoint collected

# vector of sequenced cross-sectional studyIDs
seq_studyID <- seq_cross$studyID
# filter out the studyIDs already represented by the sequenced community table
# join the table with the seqeunced metadata
comm_metadata <- metadata_1yr %>% filter(!studyID %in% seq_studyID) %>%
  full_join(. , seq_metadata)

#### Eczema incidence by study ID at any timepoint ####
## Since we are making use of longitudinal data we need to start from the whole metadata sheet ##

# Get out the eczema related columns
# Related Columns determined by reading the intake forms
# [Old Form] Diagnosis list is a checklist of diagnosis for many diseases including eczema 
# [New Forms] Geared at different ages but asks the same question
# 7 or 8: Has the child ever had eczema (also called atopic dermatitis)?
# 7a or 8a: Has a doctor or healthcare provider ever diagnosed the child with eczema or atopic dermatitis?

# Whole community
eczema <- metadata %>% 
  select(studyID, timepoint, ECHOTPCoded, assessmentAgeMonths, 
         Participants..Child_Diagnosis_List, air_inf_d7a, air_inf_d7, 
         air_ec_e8, air_ec_e8a, air_mc_e8, air_mc_e8a, air_adol_e8,
         air_adol_e8a) %>% # this step just limits the columns to make double checking easier; doesn't serve a data transformation purpose
  mutate(diagnosis = ifelse(Participants..Child_Diagnosis_List == "", NA, Participants..Child_Diagnosis_List)) %>%
  mutate(eczema = case_when(air_inf_d7 == 1 | air_inf_d7a == 1 | 
                              air_ec_e8 == 1 | air_ec_e8a == 1 | 
                              air_mc_e8 == 1 | air_mc_e8a == 1 | 
                              air_adol_e8 == 1 | air_adol_e8a == 1 ~ "Yes",
                            air_inf_d7 == 2 | air_inf_d7a == 2 | 
                              air_ec_e8 == 2 | air_ec_e8a == 2 | 
                              air_mc_e8 == 2 | air_mc_e8a == 2 | 
                              air_adol_e8 == 2 | air_adol_e8a == 2 ~ "No",
                            str_detect(diagnosis, "(?i)Eczema") ~ "Yes", #(?i) ignores case sensitivity
                            str_detect(diagnosis, "(?i)Eczema", negate = TRUE) ~ "No")) %>%
  group_by(studyID) %>% arrange(timepoint, .by_group = TRUE) %>%
  mutate(eczema_studyID = case_when(any(eczema == "Yes") ~ "Yes", #yes need to come first to override conflicting No
                                    any(eczema == "No") ~ "No")) %>% ungroup() #if even one sample has eczema they all have eczema
## eczema count
comm_metadata <- eczema %>% select(studyID, timepoint, eczema, eczema_studyID) %>%
  left_join(comm_metadata, ., by = c("studyID", "timepoint")) 
comm_metadata %>% count(eczema_studyID)

# Sequenced community
seq_metadata <- eczema %>% select(studyID, timepoint, eczema, eczema_studyID) %>%
  left_join(seq_metadata, ., by = c("studyID", "timepoint"))
## eczema count  
seq_metadata %>% count(eczema_studyID)

#### Breastfeeding Data ####
## Since we are making use of longitudinal data we need to start from the whole metadata sheet ##

#Make the columns easier to wrangle and double check steps
feeding <- metadata %>% select(studyID, timepoint, ECHOTPCoded, assessmentAgeMonths,
                               starts_with("Redcap_Ess_CHB_IFP.."), starts_with("Redcap_Ess_CHB_CFH.."), 
                               starts_with("BreastfeedingDone.."), starts_with("BreastfeedingStill..")) %>% #limit the table to make wrangling and double-checking easier
  select_if(~!(all(is.na(.)))) %>% #select the columns that are not all entirely NA
  select(-contains(c("participantid", "cohortid", "Discrepancy", "respondent",
                     "formdt", "pin", "protocolid","sequencenum", "siteid", "visitname",
                     "Redcap_Ess_CHB_IFP..StudyID", "Redcap_Ess_CHB_IFP..Timepoint",
                     "Redcap_Ess_CHB_CFH..StudyID", "Redcap_Ess_CHB_CFH..Timepoint"))) %>% # get rid of all columns that are duplicates and preventing getting rid of form name in next step
  select_all(~gsub("Redcap_Ess_CHB_IFP..|Redcap_Ess_CHB_CFH..|BreastfeedingDone..|BreastfeedingStill..", "", .))

## Identify the important columns ##
# description of the column titles for reference:
## From the Form Infant Feeding Practices
#   ifp_b01 = ever breastfed or fed pumped milk
#   ifp_b02 = how long after birth was baby put on breastmilk
#   ifp_b03 = did mother fortify breastmilk with HMF supplement
#   ifp_b04 = did mother completely stop breastfeeding and pumping milk
#   ifp_b04a = did the mother breastfeed for as long as she wanted to
#   ifp_b05 = does the child feed from both breasts while breastfeeding
#   ifp_b08_1 = last 7 days fed breast or pumped milk
#   
#   ifp_c01 = has the child ever been formulafed
#   ifp_c02 = has the child been forluma fed in the last 7 days
#   
#   ifp_d02 = In the past 7 days, has the child been fed anything other than breastmilk/formula
#   ifp_d03 = In the past 7 days, has the child been fed baby cereal
#   ifp_d04 = has the child ever been fed solid food
#   ifp_d05 = In the past 7 days, any food in addition to breastmilk/formula/cereal
#   In the past 7 days, has the child eaten [x]
#   ifp_d06h1 = cereals/starches | ifp_d06i1 = fruit| ifp_d06j1 = vegetables |
#     ifp_d06k1 = meat/chicken/combo | ifp_d06l1 = fish/shellfish | 
#     ifp_d06m1 = eggs | ifp_d06n1 = french fries | ifp_d06o1 = nuts/peanut butter | 
#     ifp_d06p1 = sweets like cake/candy
# 
# ##From the Form Complementary Feeding History [retroactively filled for infants 1st year of life]
#   cfh_01 = was the child ever breastfed [straight/pumped]
#   cfh_04 = has the biological mother completely stopped breastfeeding/pumping
#   cfh_05 = was the child ever fed someone elses breastmilk
#   cfh_06 = was the child ever fed formula
#   cfh_08 = was the child fed food/ beverage other than breastmilk/formula before 6 months
#   when was child first fed [x] (never: 12 months)
#   cfh_12a = breast milk 
#   cfh_12b = pumped breast milk 
#   cfh_12c = formula
#   cfh_12h = other dairy[yoghurt/pudding etc]
#   cfh_12i = other soy food [tofu etc]
#   cfh_12l = baby cereal
#   cfh_12m = other start/carbs
#   cfh_12n = fruit
#   cfh_12o = vegetables
#   cfh_12i = meat/chicken
#   cfh_12i = fish/shellfish
#   cfh_12i = eggs
#   cfh_12i = french fries
#   cfh_12i = nuts/peanut butter
#   cfh_12i = sweets [cake etc]

## Set up Conditionals for Feeding

#number vector representing metadata options
number_list <- c(2,3,4,5,6)

## Breastfeeding Conditionals
A_Present <- feeding$ifp_b01 == 1 | !is.na(feeding$ifp_b02) | feeding$ifp_b03 == 1 | 
  feeding$ifp_b04 == 1 | feeding$ifp_b04a == 1 | feeding$ifp_b05 == 1 |feeding$ifp_b08_1 == 1 |
  feeding$ifp_c01 == 2 | feeding$breastFedPercent == 100.00

A_Past <- feeding$cfh_01 == 1 | feeding$cfh_05 == 1 | feeding$cfh_06 == 2 | 
  feeding$cfh_12a %in% number_list | feeding$cfh_12b %in% number_list | 
  feeding$exclusivelyNursed == "Yes" | feeding$exclusivelyBottlefedBreastmilk == "Yes"


## Formula Feeding Conditionals
B_Present <- feeding$ifp_b01 == 2 | feeding$ifp_c01 == 1 | feeding$ifp_c02 == 1 | 
  feeding$breastFedPercent == 0.00
  
B_Past <- feeding$cfh_01 == 2 | feeding$cfh_06 == 1 | feeding$cfh_12c %in% number_list | 
  feeding$exclusiveFormulaFed == "Yes" | feeding$mostOftenFormulaBrand != ""


## Solids Conditionals
S_Present_Yes <- feeding$ifp_d05 == 1 | feeding$ifp_d06h1 == 1 | feeding$ifp_d06i1 == 1 | 
  feeding$ifp_d06j1 == 1 | feeding$ifp_d06k1 == 1 | feeding$ifp_d06l1 == 1 | 
  feeding$ifp_d06m1 == 1 | feeding$ifp_d06n1 == 1 | feeding$ifp_d06o1 == 1 | 
  feeding$ifp_d06p1 == 1 | feeding$ifp_d02 == 1 | feeding$ifp_d03 == 1 | feeding$ifp_d04 == 1 | 
  feeding$solidFood =="Yes"

S_Present_No <- feeding$ifp_d05 == 2 | feeding$ifp_d06h1 == 2 | feeding$ifp_d06i1 == 2 | 
  feeding$ifp_d06j1 == 2 | feeding$ifp_d06k1 == 2 | feeding$ifp_d06l1 == 2 | 
  feeding$ifp_d06m1 == 2 | feeding$ifp_d06n1 == 2 | feeding$ifp_d06o1 == 2 | 
  feeding$ifp_d06p1 == 2 | feeding$ifp_d02 == 2 | feeding$ifp_d03 == 2 | feeding$ifp_d04 == 2 | 
  feeding$solidFood =="No"
  
S_Past <- feeding$cfh_12h %in% number_list | feeding$cfh_12i %in% number_list | 
  feeding$cfh_12l %in% number_list | feeding$cfh_12m %in% number_list | 
  feeding$cfh_12n %in% number_list | feeding$cfh_12o %in% number_list |
  feeding$cfh_12p %in% number_list | feeding$cfh_12q %in% number_list | 
  feeding$cfh_12r %in% number_list | feeding$cfh_12s %in% number_list | 
  feeding$cfh_12t %in% number_list | feeding$cfh_12u %in% number_list 

## feed_present: the breastfeeding data from current info forms like infant feeding practices and breastfeeding Still
## feed_past: the breastfeeding data from older retractively filled forms like complementary feeding history and Breastfeeding Done
## feed_past_studyid: Feeding information from the retroactively filled form applied to all timepoints in the studyID
## feedtype: breast feeding information for a timepoint of a studyID compiling current information and retroactive information
## feed_solid: solid food info

feeding <- feeding %>% 
  mutate(feed_present = case_when(A_Present == TRUE & B_Present == TRUE ~ "Mixed",
                                  A_Present == TRUE ~ "BreastFed", B_Present == TRUE ~ "FormulaFed"),
         feed_past = case_when(A_Past == TRUE & B_Past == TRUE ~ "Mixed",
                               A_Past == TRUE ~ "BreastFed", B_Past == TRUE ~ "FormulaFed")) %>%
  group_by(studyID) %>% arrange(timepoint, .by_group = TRUE) %>%
  mutate(feed_past_studyid = case_when(any(feed_past == "Mixed") ~ "Mixed", #yes need to come first to override conflicting No
                                    any(feed_past == "FormulaFed") ~ "FormulaFed",
                                    any(feed_past == "BreastFed") ~ "BreastFed")) %>% 
  ungroup() %>%
  mutate(feedtype = case_when(feed_present == "BreastFed" ~ "BreastFed",
                              feed_present == "FormulaFed" ~ "FormulaFed",
                              feed_present == "Mixed" ~ "Mixed",
                              feed_past == "BreastFed" ~ "BreastFed",
                              feed_past == "FormulaFed" ~ "FormulaFed",
                              feed_past == "Mixed" ~ "Mixed",
                              feed_past_studyid == "BreastFed" ~ "BreastFed",
                              feed_past_studyid == "FormulaFed" ~ "FormulaFed",
                              feed_past_studyid == "Mixed" ~ "Mixed"),
         feed_solid = case_when(S_Present_Yes == TRUE | S_Past == TRUE ~ "Yes",
                                S_Present_No == TRUE ~ "No")) ## DO NOT add S_Past == FALSE this is a MCQ based question so false =/= no
                              

## feeding counts
## Whole Community
comm_metadata <- feeding %>% select(studyID, timepoint, feed_present, feed_past, 
                                    feed_past_studyid, feedtype, feed_solid) %>%
  left_join(comm_metadata, .)

comm_metadata %>% count(feedtype)
comm_metadata %>% count(feed_solid)

# Sequenced community
seq_metadata <- feeding %>% select(studyID, timepoint, feed_present, feed_past, 
                                   feed_past_studyid, feedtype, feed_solid) %>%
  left_join(seq_metadata, ., by = c("studyID", "timepoint"))
seq_metadata %>% count(feedtype)
seq_metadata %>% count(feed_solid)

## Breastfeeding exposure
# Whole community
comm_metadata <- comm_metadata %>% rename(bf_percent = BreastfeedingStill..breastFedPercent) %>%
  mutate(feed_exposure = case_when(bf_percent == 100 ~ "100% breastfeeding",
                                   bf_percent >= 80 ~ ">=80%", bf_percent > 50 ~ "50-80%",
                                   bf_percent == 50 ~ "50%", bf_percent > 20 ~ "20-50%",
                                   bf_percent > 0 ~ "<=20%", bf_percent == 0 ~ "100% formmulafeeding"))
comm_metadata %>% count(feed_exposure)
# Sequenced community
seq_metadata <- seq_metadata %>% rename(bf_percent = BreastfeedingStill..breastFedPercent) %>%
  mutate(feed_exposure = case_when(bf_percent == 100 ~ "100% breastfeeding",
                                   bf_percent >= 80 ~ ">=80%", bf_percent > 50 ~ "50-80%",
                                   bf_percent == 50 ~ "50%", bf_percent > 20 ~ "20-50%",
                                   bf_percent > 0 ~ "<=20%", bf_percent == 0 ~ "100% formmulafeeding"))
seq_metadata %>% count(feed_exposure)

#### Other Infant Metric ####
## Age
# Whole Community
mean(comm_metadata$childAgeMonths, na.rm = TRUE) # # mean infant age in months
mean(comm_metadata$childAgeMonths, na.rm = TRUE)*30.417 # mean infant age in days
sum(!is.na(comm_metadata$childAgeMonths)) # the number of samples for which we have ages 
comm_metadata %>% filter(childAgeMonths <= 3) %>% nrow() # no. of sample <= 3mo
comm_metadata %>% filter(childAgeMonths > 3 & childAgeMonths <= 6 ) %>% nrow() # no. of samples b/w 3 and 6 mo
comm_metadata %>% filter(childAgeMonths > 6) %>% nrow() # no. of samples b/w 6 and < 13 mo

# Sequenced
# Days
mean(seq_metadata$childAgeMonths, na.rm = TRUE) # mean infant age in months
mean(seq_metadata$childAgeMonths, na.rm = TRUE)*30.417 # mean infant age in days
sum(!is.na(seq_metadata$childAgeMonths)) # the number of samples for which we have ages 
seq_metadata %>% filter(childAgeMonths <= 3) %>% nrow() # no. of sample <= 3mo
seq_metadata %>% filter(childAgeMonths > 3 & childAgeMonths <= 6 ) %>% nrow() # no. of samples b/w 3 and 6 mo
seq_metadata %>% filter(childAgeMonths > 6) %>% nrow() # no. of samples b/w 6 and < 13 mo

## Sex
## Whole Community
comm_metadata %>% count(childGender)
## Sequenced
seq_metadata %>% count(childGender)

## Race 
## Need to simplify the various racial categories observedf by running the following
#print(metadata %>% count(Participants..Merge_Dem_Child_Race))
## Whole community
comm_metadata <- comm_metadata %>% 
  rename(infant_race = Participants..Merge_Dem_Child_Race) %>% 
  mutate(simple_race = case_when(str_detect(infant_race, "\n") | infant_race == "Mixed Race" ~ "Mixed",
                                 infant_race == "Asian "| infant_race == "Other Asian" ~ "Asian", # the space after Asian is not a typo
                                 infant_race == "White" ~ "White",
                                 infant_race == "Black or African American" ~ "Black",
                                 infant_race == "American Indian or Alaska Native" ~ "American Indian",
                                 infant_race == "Some other race"| infant_race == "Unknown" ~ "Other"))
comm_metadata %>% count(simple_race)    
## Sequenced
seq_metadata <- seq_metadata %>% 
  rename(infant_race = Participants..Merge_Dem_Child_Race) %>% 
  mutate(simple_race = case_when(str_detect(infant_race, "\n") | infant_race == "Mixed Race" ~ "Mixed",
                                 infant_race == "Asian "| infant_race == "Other Asian" ~ "Asian",
                                 infant_race == "White" ~ "White",
                                 infant_race == "Black or African American" ~ "Black",
                                 infant_race == "American Indian or Alaska Native" ~ "American Indian",
                                 infant_race == "Some other race"| infant_race == "Unknown" ~ "Other"))
seq_metadata %>% count(simple_race) 

## Height
# Whole Community
mean(comm_metadata$childHei_cm, na.rm = TRUE)
mean(comm_metadata$childHeight, na.rm = TRUE) # in inches
sum(!is.na(comm_metadata$childHei_cm)) # the number of samples for which we have height
sum(!is.na(comm_metadata$childHeight))
#Sequenced
mean(seq_metadata$childHei_cm, na.rm = TRUE)
mean(seq_metadata$childHeight, na.rm = TRUE) # in inches
sum(!is.na(seq_metadata$childHei_cm)) # the number of seq samples for which we have height
sum(!is.na(seq_metadata$childHeight))

## Weight
# Whole Community
mean(comm_metadata$childWei_kg, na.rm = TRUE)*100 # in grams
mean(comm_metadata$childWeight, na.rm = TRUE) # in pounds
sum(!is.na(comm_metadata$childWei_kg)) # the number of samples for which we have weight
sum(!is.na(comm_metadata$childWeight))
#Sequenced
mean(seq_metadata$childWei_kg, na.rm = TRUE)*100 # in grams
mean(seq_metadata$childWeight, na.rm = TRUE) # in inches
sum(!is.na(seq_metadata$childWei_kg)) # the number of seq samples for which we have weight
sum(!is.na(seq_metadata$childWeight))


#### Maternal/Pregnancy/Delivery Metrics ####

# Delivery..GBS = Did the mother receive an antibiotic during labor for GBS?
# Delivery..birthType = Was the delivery vaginal or cesarean?
# Delivery..birthAssistance = Was the delivery assisted by a vacuum and/or forceps?
# From the current PMCI forms
# Redcap_Ess_Prg_PMCI..pmci_c5 = Did you have an infection in your uterus/womb during labor for which they gave you antibiotics?
# Redcap_Ess_Prg_PMCI..pmci_c6___1 = Vaginal Birth
# Redcap_Ess_Prg_PMCI..pmci_c6___2 = C-section
# Redcap_Ess_Prg_PMCI..pmci_c6a: Was the birth assisted? [1=Forceps; 2=vacuum; 3=both]
# Redcap_Ess_Prg_PMCI..pmci_g06a : Did you receive IV antibiotics during labor?

## Antibiotics During Labor
# Whole Community
comm_metadata <- comm_metadata %>% 
  mutate(abx_labor = case_when(Delivery..GBS == "Yes" | Redcap_Ess_Prg_PMCI..pmci_c5 == 1 | Redcap_Ess_Prg_PMCI..pmci_g06a == 1   ~ "Yes",
                               Delivery..GBS == "No" | Redcap_Ess_Prg_PMCI..pmci_c5 == 2 | Redcap_Ess_Prg_PMCI..pmci_g06a == 2 ~ "No"))
comm_metadata %>% count(abx_labor)
# Sequenced
seq_metadata <- seq_metadata %>% 
  mutate(abx_labor = case_when(Delivery..GBS == "Yes" | Redcap_Ess_Prg_PMCI..pmci_c5 == 1 | Redcap_Ess_Prg_PMCI..pmci_g06a == 1   ~ "Yes",
                               Delivery..GBS == "No" | Redcap_Ess_Prg_PMCI..pmci_c5 == 2 | Redcap_Ess_Prg_PMCI..pmci_g06a == 2 ~ "No"))
seq_metadata %>% count(abx_labor)

## Delivery
# Whole Community
comm_metadata <- comm_metadata %>% 
  mutate(delivery = case_when(Delivery..birthType == "Vaginal" | Redcap_Ess_Prg_PMCI..pmci_c6___1 == 1 ~ "Vaginal",
                              Delivery..birthType == "Cesarean" | Redcap_Ess_Prg_PMCI..pmci_c6___2 == 1 ~ "Cesarean"))
comm_metadata %>% count(delivery)
# Sequenced
seq_metadata <- seq_metadata %>% 
  mutate(delivery = case_when(Delivery..birthType == "Vaginal" | Redcap_Ess_Prg_PMCI..pmci_c6___1 == 1 ~ "Vaginal",
                              Delivery..birthType == "Cesarean" | Redcap_Ess_Prg_PMCI..pmci_c6___2 == 1 ~ "Cesarean"))
seq_metadata %>% count(delivery)

## Assistance Type During Delivery
### This has been abandoned due to poor metadata availability
# Whole Community
test <- comm_metadata %>% 
  mutate(delivery_ass = case_when(Delivery..birthAssistance == "Forceps"| Redcap_Ess_Prg_PMCI..pmci_c6a == 1 ~ "Forceps",
                                  Delivery..birthAssistance == "Vacuum"| Redcap_Ess_Prg_PMCI..pmci_c6a == 2 ~ "Vacuum",
                                  Delivery..birthAssistance == "Forceps\nVacuum"| Redcap_Ess_Prg_PMCI..pmci_c6a == 3 ~ "Both"))
table(test$delivery_ass)

## Maternal Education
## Whole Community
comm_metadata <- comm_metadata %>% 
  rename(mateduc = Participants..Merge_Dem_Mom_Education) %>%
  mutate(mateduc_grp = 
           case_when(
             mateduc == "Professional or Doctorate Degree (PhD, EdD, MD, JD)" ~ "professional/doctorate",
             mateduc == "8th grade or less" ~ "8th grade or less",
             mateduc == "Associate's degree (AA, AS)" ~ "associates",
             mateduc == "Bachelor's degree (BA, BS)" ~ "bachelors",
             mateduc == "Master's degree (MA, MS, MEd, MSW, MBA, MPH)" ~ "masters",
             str_detect(mateduc, "GED or equivalent|High School Graduate|High school degree") ~ "high school/GED",
             mateduc == "Partial College or Specialized Training" ~ "partial college/special training",
             mateduc == "Some college, no degree" ~ "college, no degree",
             str_detect(mateduc, "Partial High School|Some high school, no degree") ~ "high school, no degree",
             mateduc == "College Graduate" ~ "bachelors"))
comm_metadata %>% count(mateduc_grp)

## Sequenced Community
seq_metadata <- seq_metadata %>% 
  rename(mateduc = Participants..Merge_Dem_Mom_Education) %>%
  mutate(mateduc_grp = 
           case_when(
             mateduc == "Professional or Doctorate Degree (PhD, EdD, MD, JD)" ~ "professional/doctorate",
             mateduc == "8th grade or less" ~ "8th grade or less",
             mateduc == "Associate's degree (AA, AS)" ~ "associates",
             mateduc == "Bachelor's degree (BA, BS)" ~ "bachelors",
             mateduc == "Master's degree (MA, MS, MEd, MSW, MBA, MPH)" ~ "masters",
             str_detect(mateduc, "GED or equivalent|High School Graduate|High school degree") ~ "high school/GED",
             mateduc == "Partial College or Specialized Training" ~ "partial college/special training",
             mateduc == "Some college, no degree" ~ "college, no degree",
             str_detect(mateduc, "Partial High School|Some high school, no degree") ~ "high school, no degree",
             mateduc == "College Graduate" ~ "bachelors"))
seq_metadata %>% count(mateduc_grp)
#### Social Metadata ####
## Childcare
# no relevant columns in the ChildcareStatic form
# relevant columns from Childcare Dynamic form
#ChildcareDynamic..familyMemberTakeCare # non-parental family member and nanny
#ChildcareDynamic..stayAtHomeParent # is there a stay at home parent in the household
#EarlyCare Education form logs hours of the time of childcare used
# Q1 and 2 ask about daycare and Q6-11 ask about private care by sibling/nanny etc
## Whole community
comm_metadata <- comm_metadata %>%
  mutate(childcare_daycare = case_when(EarlyCareEducation..Q1 > 0 | EarlyCareEducation..Q2 > 0 ~ "Yes"),
         childcare_private = case_when(EarlyCareEducation..Q6 > 0 | EarlyCareEducation..Q7 > 0 |
           EarlyCareEducation..Q8 > 0 | EarlyCareEducation..Q9 > 0 | 
           EarlyCareEducation..Q10 > 0 | EarlyCareEducation..Q11 > 0 ~ "Yes")) %>%
  mutate(daycare = case_when(ChildcareDynamic..daycareOutsideHome == "Yes" & 
                               ChildcareDynamic..familyMemberTakeCare == "Yes" ~ "Mixed",
                             childcare_daycare == "Yes" & childcare_private == "Yes" ~ "Mixed",
                             ChildcareDynamic..daycareOutsideHome == "Yes" | 
                               childcare_daycare == "Yes" ~ "Daycare",
                             ChildcareDynamic..familyMemberTakeCare == "Yes" |
                               childcare_private == "Yes" ~ "Private"))
comm_metadata %>% count(daycare)

## Sequenced community
seq_metadata <- seq_metadata %>%
  mutate(childcare_daycare = case_when(EarlyCareEducation..Q1 > 0 | EarlyCareEducation..Q2 > 0 ~ "Yes"),
         childcare_private = case_when(EarlyCareEducation..Q6 > 0 | EarlyCareEducation..Q7 > 0 |
                                         EarlyCareEducation..Q8 > 0 | EarlyCareEducation..Q9 > 0 | 
                                         EarlyCareEducation..Q10 > 0 | EarlyCareEducation..Q11 > 0 ~ "Yes")) %>%
  mutate(daycare = case_when(ChildcareDynamic..daycareOutsideHome == "Yes" & 
                               ChildcareDynamic..familyMemberTakeCare == "Yes" ~ "Mixed",
                             childcare_daycare == "Yes" & childcare_private == "Yes" ~ "Mixed",
                             ChildcareDynamic..daycareOutsideHome == "Yes" | 
                               childcare_daycare == "Yes" ~ "Daycare",
                             ChildcareDynamic..familyMemberTakeCare == "Yes" |
                               childcare_private == "Yes" ~ "Private"))
seq_metadata %>% count(daycare)

## Siblings
# the number of children living in the household
FamilyInformation..numberChildren > 1 # the number of children including the child of interest
HouseholdComposition..NumberChildren > 0 # the number of additional children with child of interest
Redcap_Ess_Dem_HHC_C..hhc_c_02a1 #the number of additional children currently living with the child
metadata %>% count(Redcap_Ess_Dem_HHC_C..hhc_c_02a1)
## Whole Community
comm_metadata <- comm_metadata %>%
  mutate(household_children = case_when(Redcap_Ess_Dem_HHC_C..hhc_c_02a1 == 0 |
                                          FamilyInformation..numberChildren <= 1 |
                                          HouseholdComposition..NumberChildren == 0 ~ "No",
                                        FamilyInformation..numberChildren > 1 |
                                          HouseholdComposition..NumberChildren > 0 |
                                          Redcap_Ess_Dem_HHC_C..hhc_c_02a1 > 0 ~ "Yes"))
comm_metadata %>% count(household_children)  
## Seq Community
seq_metadata <- seq_metadata %>%
  mutate(household_children = case_when(Redcap_Ess_Dem_HHC_C..hhc_c_02a1 == 0 |
                                          FamilyInformation..numberChildren <= 1 |
                                          HouseholdComposition..NumberChildren == 0 ~ "No",
                                        FamilyInformation..numberChildren > 1 |
                                          HouseholdComposition..NumberChildren > 0 |
                                          Redcap_Ess_Dem_HHC_C..hhc_c_02a1 > 0 ~ "Yes"))
seq_metadata %>% count(household_children)

## Pets
## Whole Community
comm_metadata <- comm_metadata %>% 
  mutate(pet_type = case_when(Pets..petType == "Cat" ~ "Cat",
                                Pets..petType == "Dog" ~ "Dog",
                             str_detect(Pets..petType, 
                                        "Axolotl|Bearded Dragon|Beta Fish|Bird|Bunny|Chinchilla|Fish|Reptile") ~ "Other")) %>%
  mutate(pet_own = case_when(str_detect(pet_type, "Cat|Dog|Other") ~ "Yes"))

comm_metadata %>% count(pet_own)
comm_metadata %>% count(pet_type)

## Seq Community
seq_metadata <- seq_metadata %>% 
  mutate(pet_type = case_when(Pets..petType == "Cat" ~ "Cat",
                              Pets..petType == "Dog" ~ "Dog",
                              str_detect(Pets..petType, 
                                         "Axolotl|Bearded Dragon|Beta Fish|Bird|Bunny|Chinchilla|Fish|Reptile") ~ "Other")) %>%
  mutate(pet_own = case_when(str_detect(pet_type, "Cat|Dog|Other") ~ "Yes"))

seq_metadata %>% count(pet_own)
seq_metadata %>% count(pet_type)

#### Export all the relevant metadata ####
# metadata corresponding to the demographics table for the whole community <1yr
comm_mdexport <- comm_metadata %>% 
  select(sample, studyID, timepoint, childAgeMonths, childGender, infant_race, 
         simple_race, childHeight, childHei_cm, childWeight, childWei_kg,
         eczema, eczema_studyID, feed_present, feed_past, feed_past_studyid, 
         feedtype, feed_solid, feed_exposure, delivery, abx_labor, mateduc_grp,
         childcare_daycare, childcare_private, daycare, household_children, 
         pet_own, pet_type)
write_csv(comm_mdexport, "./metadata_exports/commdem_resmd.csv")
# metadata corresponding to the demographics table for the sequenced community <1yr
seq_mdexport <- seq_metadata %>% 
  select(sample, studyID, timepoint, childAgeMonths, childGender, infant_race, 
         simple_race, childHeight, childHei_cm, childWeight, childWei_kg,
         eczema, eczema_studyID, feed_present, feed_past, feed_past_studyid, 
         feedtype, feed_solid, feed_exposure, delivery, abx_labor, mateduc_grp,
         childcare_daycare, childcare_private, daycare, household_children, 
         pet_own, pet_type)
write_csv(seq_mdexport, "./metadata_exports/seqdem_resmd.csv")

#### Checking that I have the SAMPLEIDs from all MGX runs ####
# The following confirmed I have all the samples we have sequenced as of April 2023
# it calls in its own raw data instead of building off previous wrangling sections

# Call in the 1 year RESONANCE metadata 
## Data is from the January 2023 export from Jen
## Large dataframe - to import into R I subsetted the whole dataset 
## for only infants less than 1 yo in excel first

RES_1yr <- read.csv("./raw_data_res/Prioty_data_013023_12months.csv")

## Call in Kevin's File that has limited amount of information on all our samples
# Most important contains the age in months
all_samples <- read.csv("./raw_data_res/Samples-Everything.csv")

## Call in all the mgx samples from airtable has limited amount of info
# direct csv download of the airtable environment
# Most important are the SampleIDs and age in months
old_mgx <- read.csv("./raw_data_res/Prioty_allmgx.csv")

old_mgx_1year <- old_mgx %>%
  filter(childAgeMonths < 13)
old_mgx_1year <- old_mgx_1year$sample
  ## Call in all the mgx samples from the mgx run that came in march 2023
  # list is from the samples found in IMR-SampleSubmissionSheet_v18(MCCANN_SEP06).xlsx
  # copied and pasted from clipboard using datapasta
  # the list is also saved as a separate text file: Seq_MCCANN_Sep06
  
new_mgx <- c("FG02743", "FG02709", "FG02711", "FG02713", "FG02715", "FG02717", 
              "FG02719", "FG02721", "FG02723", "FG02725", "FG02727", "FG02729", 
              "FG02731", "FG02733", "FG02735", "FG02737", "FG02739", "FG02745", 
              "FG02747", "FG02749", "FG02751", "FG02753", "FG02755", "FG02757", 
              "FG02759", "FG02761", "FG02763", "FG02765", "FG02767", "FG02769", 
              "FG02771", "FG02773", "FG02775", "FG02777", "FG02779", "FG02781", 
              "FG02783", "FG02785", "FG02787", "FG02789", "FG02791", "FG02793", 
              "FG02795", "FG02797", "FG02799", "FG02801", "FG02804", "FG02808", 
              "FG02811", "FG02809", "FG02813", "FG02815", "FG02817", "FG02819", 
              "FG02821", "FG02823", "FG02825", "FG02827", "FG02829", "FG02831", 
              "FG02833", "FG02835", "FG02837", "FG02839", "FG02907", "FG02909", 
              "FG02914", "FG02916", "FG02919", "FG02921", "FG02923", "FG02925", 
              "FG02875", "FG02877", "FG02879", "FG02881", "FG02883", "FG02885", 
              "FG02888", "FG02890", "FG02858", "FG02860", "FG02862", "FG02864", 
              "FG02866", "FG02868", "FG02870", "FG02872", "FG02841", "FG02843", 
              "FG02845", "FG02847", "FG02849", "FG02852", "FG02854", "FG02856", 
              "FG02892", "FG02894", "FG02896", "FG02898", "FG02901", "FG02903", 
              "FG02905", "FG02927")

## Find all the mgx samples

all_mgx <- all_samples %>%
  filter(sample %in% new_mgx | sample %in% old_mgx_1year) %>%
  filter(childAgeMonths < 13) 

# there was no need to do it this way - kevin's airtable mgx file 
# has all the new sequences we received on march 2023 - it is batch 19
# can just filter for childAgeMonths < 13 in the 'old_mgx' file instead

all_mgx <- old_mgx %>%
  filter(childAgeMonths < 13) %>%
  left_join(., all_samples, by = "sample")
