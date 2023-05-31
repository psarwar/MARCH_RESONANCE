##look at files inthe working directory
list.files()
# occasionally might need to clear environment for sanity
rm(list =ls())
#installation with tutorial
remotes::install_github("biobakery/anpan", build_vignettes = TRUE)
vignette("anpan_tutorial", package = 'anpan')
library(anpan)
library(dplyr)

# to run anpan you need to download the cmdstanr package
# this package actually assists in downloading cmdstan but we need to do it in a separate step
# steps from anpan tutorial
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
# after running the above got the following output: 
# *Finished installing CmdStan to /Users/psarwar/.cmdstan/cmdstan-2.32.1
# CmdStan path set to: /Users/psarwar/.cmdstan/cmdstan-2.32.1

#### MARCH strainphlan ####
#Get in the march metadatafile
meta_file <- read.csv("./processed/march_seq_metadata_allvar.csv") %>% 
  rename(sample_id = seqname) %>%
  mutate(eczema = case_when(eczema == "Yes" ~ TRUE,
                            eczema == "No" ~ FALSE)) %>% #use true/falso not 1/0 binary [this confuses anpan with a binary as a continuous scale]
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
                          ladderize = TRUE
                          labels = result_2$model_input$sample_id)
thing <- result_3$pglmm_fit
fileprint <- thing$summary()
write.csv(fileprint, file = "file.csv")
##Check Kevin's metadata file
rsync -avP ada:/home/psarwar/Desktop/strainphlan/manuscript_meta.csv ./



bugs = [

"Bifidobacterium_longum",
"Bifidobacterium_breve",
"Bifidobacterium_bifidum",
"Bifidobacterium_catenulatum",
"Bifidobacterium_animalis",
"Bacteroides_fragilis",
"Lactobacillus_casei",
"Lactobacillus_acidophilus",
"Escherichia_coli",
"Staphylococcus_aureus",
"Klebsiella_pneumoiae",
"Faecalibacterium_prausnitzii",
"Ruminococcus_gnavus",
"Akkermansia muciniphila"
