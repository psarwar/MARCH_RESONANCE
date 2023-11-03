##look at files inthe working directory
list.files()
# occasionally might need to clear environment for sanity
rm(list =ls())
#installation with tutorial
remotes::install_github("biobakery/anpan", build_vignettes = TRUE)
vignette("anpan_tutorial", package = 'anpan')
library(anpan)
library(dplyr)
library(readr)
library(stringr)
browseVignettes("anpan")
# to run anpan you need to download the cmdstanr package
# this package actually assists in downloading cmdstan but we need to do it in a separate step
# steps from anpan tutorial
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
# after running the above got the following output: 
# *Finished installing CmdStan to /Users/psarwar/.cmdstan/cmdstan-2.32.1
# CmdStan path set to: /Users/psarwar/.cmdstan/cmdstan-2.32.1

## Join the humann pathabundance outputs using the humann_join_tables command from humann3.7
# humann_join_tables 
# --input /Users/psarwar/VKC_Lab/CodingRelated/MARCH_RESONANCE/humann_pathabundance 
# --output humann_pathabundance.tsv --file_name pathabundance

####anpan path abundance####
## make the naming convention so that the column names match a sample_id column in the metadata file
marchres_pathab <- read_tsv("./raw_data_mixed/humann_pathabundance.tsv") %>%
  rename_with(~str_remove_all(., '_Abundance')) %>% rename(pathway = `# Pathway`)

## Check and harmonize the filenames for humann pathab outputs
meta_file <- read.csv("./metadata_exports/marchres_combinedmd_ids.csv") %>% 
  rename(sample_id = old_filename) %>%
  mutate(eczema = case_when(eczema == "Yes" ~ TRUE,
                            eczema == "No" ~ FALSE)) %>%
  mutate(feedtype = factor(feedtype, levels = c("FormulaFed","Mixed","BreastFed"))) %>%
  mutate(sid_old = str_replace_all(sid_old, '_','-')) %>%
  mutate(old_seqname_res = paste0(sid_old,"_",S_well)) 

humann_samplenames <- colnames(marchres_pathab)

## make the naming convention so that the column names match a sample_id column in the metadata file
sid_old_res <- meta_file %>%
  filter(old_seqname_res %in% humann_samplenames)
col.from <- sid_old_res$old_seqname_res
col.to <- sid_old_res$sample_id
marchres_pathab <- marchres_pathab %>% rename_at(vars(col.from), ~col.to)

## Make the input files for all the HMO metabolizing bacteria

## Bifidobacterium longum
blongum <- marchres_pathab %>%
  filter(grepl('s__Bifidobacterium_longum', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(blongum, "./HMOstrains_pathab/blongum_pathab.tsv")

## Bifidobacterium breve
bbreve <- marchres_pathab %>%
  filter(grepl('s__Bifidobacterium_breve', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(bbreve, "./HMOstrains_pathab/bbreve_pathab.tsv")

## Bifidobacterium bifidum
bbifidum <- marchres_pathab %>%
  filter(grepl('s__Bifidobacterium_bifidum', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(bbifidum, "./HMOstrains_pathab/bbifidum_pathab.tsv")

## Bifidobacterium catenulatum
bcatenulatum <- marchres_pathab %>%
  filter(grepl('s__Bifidobacterium_catenulatum', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(bcatenulatum, "./HMOstrains_pathab/bcatenulatum_pathab.tsv")

## Bifidobacterium animalis
banimalis <- marchres_pathab %>%
  filter(grepl('s__Bifidobacterium_animalis', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(banimalis, "./HMOstrains_pathab/banimalis_pathab.tsv")

## Bacteroides fragilis
bfragilis <- marchres_pathab %>%
  filter(grepl('s__Bacteroides_fragilis', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(bfragilis, "./HMOstrains_pathab/bfragilis_pathab.tsv")

## Lactobacillus casei
lcasei<- marchres_pathab %>%
  filter(grepl('s__Lactobacillus_casei', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(lcasei, "./HMOstrains_pathab/lcasei_pathab.tsv")

## Escherichia coli
ecoli<- marchres_pathab %>%
  filter(grepl('s__Escherichia_coli', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(ecoli, "./HMOstrains_pathab/ecoli_pathab.tsv")

## Staphylococcus aureus
saureus<- marchres_pathab %>%
  filter(grepl('s__Staphylococcus_aureus', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(saureus, "./HMOstrains_pathab/saureus_pathab.tsv")

## Klebsiella pneumoniae
kpneumoniae <- marchres_pathab %>%
  filter(grepl('s__Klebsiella_pneumoniae', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(kpneumoniae, "./HMOstrains_pathab/kpneumoniae_pathab.tsv")
## Faecalibacterium prausnitzii
fprausnitzii <- marchres_pathab %>%
  filter(grepl('s__Faecalibacterium_prausnitzii', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(fprausnitzii, "./HMOstrains_pathab/fprausnitzii_pathab.tsv")
## Ruminococcus gnavus
rgnavus <- marchres_pathab %>%
  filter(grepl('s__Ruminococcus_gnavus', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(rgnavus, "./HMOstrains_pathab/rgnavus_pathab.tsv")
## Akkermansia muciniphila
amuciniphila <- marchres_pathab %>%
  filter(grepl('s__Akkermansia_muciniphila', pathway)) %>%
  select(where(~ any(. != 0)))
write_tsv(amuciniphila, "./HMOstrains_pathab/amuciniphila_pathab.tsv")

hmo_anpan <- anpan_batch(bug_dir = "./HMOstrains_pathab",
                         meta_file = meta_file,
                         out_dir = "./anpan_humann",
                         filtering_method = 'kmeans',
                         model_type = 'fastglm',
                         covariates = c('child_age_month', 'feedtype'),
                         outcome = 'eczema',
                         plot_ext = 'pdf',
                         save_filter_stats = TRUE,
                         omit_na = TRUE)


hmo_anpan_nofeeding <- anpan_batch(bug_dir = "./HMOstrains_pathab",
                                                meta_file = meta_file,
                                                out_dir = "./anpan_humann",
                                                filtering_method = 'kmeans',
                                                model_type = 'fastglm',
                                                covariates = 'child_age_month',
                                                outcome = 'eczema',
                                                plot_ext = 'pdf',
                                                save_filter_stats = TRUE,
                                                omit_na = TRUE)
meta_file_bfmod <- meta_file %>%
  mutate(breastfed = case_when(feedtype == 'BreastFed' ~ TRUE,
                               feedtype == 'Mixed' | feedtype == 'FormulaFed' ~ FALSE))

hmo_anpan_breastfed <- anpan_batch(bug_dir = "./HMOstrains_pathab",
                                   meta_file = meta_file_bfmod,
                                   out_dir = "./anpan_humann",
                                   filtering_method = 'kmeans',
                                   model_type = 'fastglm',
                                   covariates = 'child_age_month',
                                   outcome = 'breastfed',
                                   plot_ext = 'pdf',
                                   save_filter_stats = TRUE,
                                   omit_na = TRUE)
### anpan test run
bbifidum_anpan <- anpan(bug_file = "./bbifidumanpantestfile.tsv",
                        meta_file = marchresmd,
                        out_dir = "./anpan_humann",
                        filtering_method = 'kmeans',
                        model_type = 'fastglm',
                        covariates = c('child_age_month', 'feedtype'),
                        outcome = 'eczema',
                        plot_ext = 'pdf',
                        save_filter_stats = TRUE,
                        omit_na = TRUE)

test.p <- plot_results(res            = bbifidum_anpan,
             model_input    = "./anpan_humann/filter_stats/filtered_bbifidumanpantestfile.tsv.gz",
             covariates     = NULL, 
             outcome        = "eczema",
             bug_name       = "Bifidobacterium_bifidum",
             cluster        = "both",
             show_trees     = TRUE,
             n_top          = 20,
             q_threshold    = NULL,
             beta_threshold = NULL)  

## anpan tutorial
meta_path = system.file("extdata", "fake_metadata.tsv", 
                        package = "anpan", mustWork = TRUE)

bug_path  = system.file("extdata", "g__Madeuppy.s__Madeuppy_Fakerii.genefamilies.tsv.gz", 
                        package = "anpan", mustWork = TRUE)

####MARCH anpan pglmm ####
#Get in the march metadatafile
meta_file <- read.csv("./processed/march_seq_metadata_allvar.csv") %>% 
  rename(sample_id = seqname) %>%
  mutate(eczema = case_when(eczema == "Yes" ~ TRUE,
                            eczema == "No" ~ FALSE)) %>% #use true/false not 1/0 binary [this confuses anpan with a binary as a continuous scale]
  mutate(feedtype = as.factor(feedtype))

#Set paths to the different treefiles
tree_file = "./treefiles/march_binfantis/RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre"
tree_file_2 = "./treefiles/march_blongum/RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre"
tree_file_3 = "./treefiles/march_bfraglis/RAxML_bestTree.s__Bacteroides_fragilis.StrainPhlAn3.tre"

# Run anpan
# the node labels of the tree has to match a sample_id column on the metadata file

##b infantis
result <- anpan_pglmm(meta_file,
            tree_file,
            out_dir = "./treefiles/anpan_march_binfantis_poster",
            outcome = "eczema",
            family = "binomial",
            covariates = "feedtype",
            omit_na = TRUE, reg_noise = FALSE, show_plot_tree = TRUE,
            run_diagnostics = TRUE)

p = plot_tree_with_post(tree_file,
                        meta_file,
                        fit        = result$pglmm_fit,
                        covariates = "feedtype",
                        outcome    = "eczema",
                        color_bars = TRUE,
                        labels = result$model_input$sample_id)
p
result$model_input
help('pareto-k-diagnostic')
?plot_tree_with_post

## b longum
result_2 <- anpan_pglmm(meta_file,
                      tree_file_2,
                      out_dir = "./treefiles/anpan_march_blongum_poster",
                      outcome = "eczema",
                      family = "binomial",
                      covariates = "feedtype",
                      omit_na = TRUE, reg_noise = FALSE, show_plot_tree = TRUE,
                      run_diagnostics = TRUE)

p_2 = plot_tree_with_post(tree_file_2,
                        meta_file,
                        fit = result_2$pglmm_fit,
                        covariates = "feedtype",
                        outcome    = "eczema",
                        color_bars = TRUE,
                        ladderize = TRUE,
                        labels = result_2$model_input$sample_id)
p_2
## b fragilis

result_3 <- anpan_pglmm(meta_file,
                        tree_file_3,
                        out_dir = "./treefiles/anpan_march_bfragilis_poster",
                        outcome = "eczema",
                        family = "binomial",
                        covariates = "feedtype",
                        omit_na = TRUE, reg_noise = FALSE, show_plot_tree = TRUE,
                        run_diagnostics = TRUE)

p_3 = plot_tree_with_post(tree_file_3,
                          meta_file,
                          fit = result_3$pglmm_fit,
                          covariates = "feedtype",
                          outcome    = "eczema",
                          color_bars = TRUE,
                          ladderize = TRUE,
                          labels = result_2$model_input$sample_id)
thing <- result_3$pglmm_fit
fileprint <- thing$summary()
write.csv(fileprint, file = "file.csv")
##Check Kevin's metadata file
#rsync -avP ada:/home/psarwar/Desktop/strainphlan/manuscript_meta.csv ./

####anpan pglmm on Deniz updated trees####
meta_file <- read.csv("./metadata_exports/marchres_combinedmd_ids.csv") %>% 
  rename(sample_id = new_filename) %>%
  mutate(eczema = case_when(eczema == "Yes" ~ TRUE,
                            eczema == "No" ~ FALSE)) %>%
  mutate(feedtype = factor(feedtype, levels = c("FormulaFed","Mixed","BreastFed")))

## fragilis
tree_file_fragilis = "./Deniz_Trees/fragilis_output/RAxML_bestTree.s__Bacteroides_fragilis.StrainPhlAn3.tre"
result_fragilis <- anpan_pglmm(meta_file,
                                 tree_file_fragilis,
                                 out_dir = "./Deniz_Trees",
                                 outcome = "eczema",
                                 family = "binomial",
                                 covariates = "feedtype",
                                 omit_na = TRUE, reg_noise = FALSE, show_plot_tree = TRUE,
                                 run_diagnostics = TRUE)

## all_infantis
tree_file_74infantis = "./Deniz_Trees/infantis_output/RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre"
result_74infantis <- anpan_pglmm(meta_file,
                            tree_file_74infantis,
                            out_dir = "./Deniz_Trees",
                            outcome = "eczema",
                            family = "binomial",
                            covariates = "feedtype",
                            omit_na = TRUE, reg_noise = FALSE, show_plot_tree = TRUE,
                            run_diagnostics = TRUE)

## 55 HMO
tree_file_55hmo = "./Deniz_Trees/55_hmo/RAxML_bestTree.s__Bifidobacterium_infantis.StrainPhlAn3.tre"

result_55hmo <- anpan_pglmm(meta_file,
                      tree_file_55hmo,
                      out_dir = "./Deniz_Trees",
                      outcome = "eczema",
                      family = "binomial",
                      covariates = "feedtype",
                      omit_na = TRUE, reg_noise = FALSE, show_plot_tree = TRUE,
                      run_diagnostics = TRUE)

id_list <- meta_file$sample_id
write.csv(id_list, file = "marchres_newseqnames.csv")

####anpan batch with kos reclassification####
marchres_kos <- read_tsv("./raw_data_mixed/humann_kos_joined.tsv") %>%
  rename_with(~str_remove_all(., '_Abundance-RPKs')) %>% rename(gene_family = `# Gene Family`)

## Check and harmonize the filenames for humann kos outputs
meta_file <- read.csv("./metadata_exports/marchres_combinedmd_ids.csv") %>% 
  rename(sample_id = old_filename) %>%
  mutate(eczema = case_when(eczema == "Yes" ~ TRUE,
                            eczema == "No" ~ FALSE)) %>%
  mutate(feedtype = factor(feedtype, levels = c("FormulaFed","Mixed","BreastFed"))) %>%
  mutate(sid_old = str_replace_all(sid_old, '_','-')) %>%
  mutate(old_seqname_res = paste0(sid_old,"_",S_well)) 

humann_samplenames <- colnames(marchres_kos)

## make the naming convention so that the column names match a sample_id column in the metadata file
sid_old_res <- meta_file %>%
  filter(old_seqname_res %in% humann_samplenames)
col.from <- sid_old_res$old_seqname_res
col.to <- sid_old_res$sample_id
marchres_kos <- marchres_kos %>% rename_at(vars(col.from), ~col.to)

## Make the input files for all the HMO metabolizing bacteria
## Bifidobacterium longum
blongum <- marchres_kos %>%
  filter(grepl('s__Bifidobacterium_longum', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(blongum, "./HMOstrains_kos/blongum_kos.tsv")

## Bifidobacterium breve
bbreve <- marchres_kos %>%
  filter(grepl('s__Bifidobacterium_breve', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(bbreve, "./HMOstrains_kos/bbreve_kos.tsv")

## Bifidobacterium bifidum
bbifidum <- marchres_kos %>%
  filter(grepl('s__Bifidobacterium_bifidum', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(bbifidum, "./HMOstrains_kos/bbifidum_kos.tsv")

## Bifidobacterium catenulatum
bcatenulatum <- marchres_kos %>%
  filter(grepl('s__Bifidobacterium_catenulatum', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(bcatenulatum, "./HMOstrains_kos/bcatenulatum_kos.tsv")

## Bifidobacterium animalis
banimalis <- marchres_kos %>%
  filter(grepl('s__Bifidobacterium_animalis', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(banimalis, "./HMOstrains_kos/banimalis_kos.tsv")

## Bacteroides fragilis
bfragilis <- marchres_kos %>%
  filter(grepl('s__Bacteroides_fragilis', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(bfragilis, "./HMOstrains_kos/bfragilis_kos.tsv")

## Lactobacillus casei
lcasei<- marchres_kos %>%
  filter(grepl('s__Lactobacillus_casei', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(lcasei, "./HMOstrains_kos/lcasei_kos.tsv")

## Escherichia coli
ecoli<- marchres_kos %>%
  filter(grepl('s__Escherichia_coli', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(ecoli, "./HMOstrains_kos/ecoli_kos.tsv")

## Staphylococcus aureus
saureus<- marchres_kos %>%
  filter(grepl('s__Staphylococcus_aureus', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(saureus, "./HMOstrains_kos/saureus_kos.tsv")

## Klebsiella pneumoniae
kpneumoniae <- marchres_kos %>%
  filter(grepl('s__Klebsiella_pneumoniae', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(kpneumoniae, "./HMOstrains_kos/kpneumoniae_kos.tsv")
## Faecalibacterium prausnitzii
fprausnitzii <- marchres_kos %>%
  filter(grepl('s__Faecalibacterium_prausnitzii', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(fprausnitzii, "./HMOstrains_kos/fprausnitzii_kos.tsv")
## Ruminococcus gnavus
rgnavus <- marchres_kos %>%
  filter(grepl('s__Ruminococcus_gnavus', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(rgnavus, "./HMOstrains_kos/rgnavus_kos.tsv")
## Akkermansia muciniphila
amuciniphila <- marchres_kos %>%
  filter(grepl('s__Akkermansia_muciniphila', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(amuciniphila, "./HMOstrains_kos/amuciniphila_kos.tsv")

hmo_anpan <- anpan_batch(bug_dir = "./HMOstrains_kos",
                         meta_file = meta_file,
                         out_dir = "./anpan_humann_kos",
                         filtering_method = 'kmeans',
                         model_type = 'fastglm',
                         covariates = c('child_age_month', 'feedtype'),
                         outcome = 'eczema',
                         plot_ext = 'pdf',
                         save_filter_stats = TRUE,
                         omit_na = TRUE)


hmo_anpan_nofeeding <- anpan_batch(bug_dir = "./HMOstrains_kos",
                                   meta_file = meta_file,
                                   out_dir = "./anpan_humann_kos",
                                   filtering_method = 'kmeans',
                                   model_type = 'fastglm',
                                   covariates = 'child_age_month',
                                   outcome = 'eczema',
                                   plot_ext = 'pdf',
                                   save_filter_stats = TRUE,
                                   omit_na = TRUE)


meta_file_bfmod <- meta_file %>%
  mutate(breastfed = case_when(feedtype == 'BreastFed' ~ TRUE,
                               feedtype == 'Mixed' | feedtype == 'FormulaFed' ~ FALSE))

hmo_anpan_breastfed <- anpan_batch(bug_dir = "./HMOstrains_kos",
                                   meta_file = meta_file_bfmod,
                                   out_dir = "./anpan_humann_kos",
                                   filtering_method = 'kmeans',
                                   model_type = 'fastglm',
                                   covariates = 'child_age_month',
                                   outcome = 'breastfed',
                                   plot_ext = 'pdf',
                                   save_filter_stats = TRUE,
                                   omit_na = TRUE)


####anpan batch with pfam reclassification####
marchres_pfams <- read_tsv("./raw_data_mixed/humann_pfams_joined.tsv") %>%
  rename_with(~str_remove_all(., '_Abundance-RPKs')) %>% rename(gene_family = `# Gene Family`)

## Check and harmonize the filenames for humann kos outputs
meta_file <- read.csv("./metadata_exports/marchres_combinedmd_ids.csv") %>% 
  rename(sample_id = old_filename) %>%
  mutate(eczema = case_when(eczema == "Yes" ~ TRUE,
                            eczema == "No" ~ FALSE)) %>%
  mutate(feedtype = factor(feedtype, levels = c("FormulaFed","Mixed","BreastFed"))) %>%
  mutate(sid_old = str_replace_all(sid_old, '_','-')) %>%
  mutate(old_seqname_res = paste0(sid_old,"_",S_well)) 

humann_samplenames <- colnames(marchres_pfams)

## make the naming convention so that the column names match a sample_id column in the metadata file
sid_old_res <- meta_file %>%
  filter(old_seqname_res %in% humann_samplenames)
col.from <- sid_old_res$old_seqname_res
col.to <- sid_old_res$sample_id
marchres_pfams <- marchres_pfams %>% rename_at(vars(col.from), ~col.to)

## Make the input files for all the HMO metabolizing bacteria
## Bifidobacterium longum
blongum <- marchres_pfams %>%
  filter(grepl('s__Bifidobacterium_longum', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(blongum, "./HMOstrains_pfams/blongum_pfams.tsv")

## Bifidobacterium breve
bbreve <- marchres_pfams %>%
  filter(grepl('s__Bifidobacterium_breve', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(bbreve, "./HMOstrains_pfams/bbreve_pfams.tsv")

## Bifidobacterium bifidum
bbifidum <- marchres_pfams %>%
  filter(grepl('s__Bifidobacterium_bifidum', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(bbifidum, "./HMOstrains_pfams/bbifidum_pfams.tsv")

## Bifidobacterium catenulatum
bcatenulatum <- marchres_pfams %>%
  filter(grepl('s__Bifidobacterium_catenulatum', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(bcatenulatum, "./HMOstrains_pfams/bcatenulatum_pfams.tsv")

## Bifidobacterium animalis
banimalis <- marchres_pfams %>%
  filter(grepl('s__Bifidobacterium_animalis', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(banimalis, "./HMOstrains_pfams/banimalis_pfams.tsv")

## Bacteroides fragilis
bfragilis <- marchres_pfams %>%
  filter(grepl('s__Bacteroides_fragilis', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(bfragilis, "./HMOstrains_pfams/bfragilis_pfams.tsv")

## Lactobacillus casei
lcasei<- marchres_pfams %>%
  filter(grepl('s__Lactobacillus_casei', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(lcasei, "./HMOstrains_pfams/lcasei_pfams.tsv")

## Escherichia coli
ecoli<- marchres_pfams %>%
  filter(grepl('s__Escherichia_coli', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(ecoli, "./HMOstrains_pfams/ecoli_pfams.tsv")

## Staphylococcus aureus
saureus<- marchres_pfams %>%
  filter(grepl('s__Staphylococcus_aureus', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(saureus, "./HMOstrains_pfams/saureus_pfams.tsv")

## Klebsiella pneumoniae
kpneumoniae <- marchres_pfams %>%
  filter(grepl('s__Klebsiella_pneumoniae', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(kpneumoniae, "./HMOstrains_pfams/kpneumoniae_pfams.tsv")
## Faecalibacterium prausnitzii
fprausnitzii <- marchres_pfams %>%
  filter(grepl('s__Faecalibacterium_prausnitzii', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(fprausnitzii, "./HMOstrains_pfams/fprausnitzii_pfams.tsv")
## Ruminococcus gnavus
rgnavus <- marchres_pfams %>%
  filter(grepl('s__Ruminococcus_gnavus', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(rgnavus, "./HMOstrains_pfams/rgnavus_pfams.tsv")
## Akkermansia muciniphila
amuciniphila <- marchres_pfams %>%
  filter(grepl('s__Akkermansia_muciniphila', gene_family)) %>%
  select(where(~ any(. != 0)))
write_tsv(amuciniphila, "./HMOstrains_pfams/amuciniphila_pfams.tsv")

hmo_anpan <- anpan_batch(bug_dir = "./HMOstrains_pfams",
                         meta_file = meta_file,
                         out_dir = "./anpan_humann_pfams",
                         filtering_method = 'kmeans',
                         model_type = 'fastglm',
                         covariates = c('child_age_month', 'feedtype'),
                         outcome = 'eczema',
                         plot_ext = 'pdf',
                         save_filter_stats = TRUE,
                         omit_na = TRUE)


hmo_anpan_nofeeding <- anpan_batch(bug_dir = "./HMOstrains_pfams",
                                   meta_file = meta_file,
                                   out_dir = "./anpan_humann_pfams",
                                   filtering_method = 'kmeans',
                                   model_type = 'fastglm',
                                   covariates = 'child_age_month',
                                   outcome = 'eczema',
                                   plot_ext = 'pdf',
                                   save_filter_stats = TRUE,
                                   omit_na = TRUE)


meta_file_bfmod <- meta_file %>%
  mutate(breastfed = case_when(feedtype == 'BreastFed' ~ TRUE,
                               feedtype == 'Mixed' | feedtype == 'FormulaFed' ~ FALSE))

hmo_anpan_breastfed <- anpan_batch(bug_dir = "./HMOstrains_pfams",
                                   meta_file = meta_file_bfmod,
                                   out_dir = "./anpan_humann_pfams",
                                   filtering_method = 'kmeans',
                                   model_type = 'fastglm',
                                   covariates = 'child_age_month',
                                   outcome = 'breastfed',
                                   plot_ext = 'pdf',
                                   save_filter_stats = TRUE,
                                   omit_na = TRUE)

