## to learn how to run the singularity exec for metaphlan for ada
cat /usr/local/bin/metaphlan

# the following is the output
cat /usr/local/bin/metaphlan                            │
#!/bin/sh                                                                                                                                            │
singularity exec --bind $PWD:$PWD --bind /babbage:/babbage /babbage/containers/metaphlan.sif metaphlan

## B longum longum
# to run raxml on singularity exec do
singularity exec --bind $PWD:$PWD --bind /babbage:/babbage /babbage/containers/metaphlan.sif /opt/conda/bin/raxmlHPC-PTHREADS-SSE3 -p 1989 -m GTRCAT -T 8 -t /home/psarwar/Desktop/strainphlan/raxml/march/RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre -b 12345 -N 100 -s /home/psarwar/Desktop/strainphlan/raxml/march/s__Bifidobacterium_longum.StrainPhlAn3_concatenated.aln -n bootstrap_tree
# map the bootstrap numbers onto the ML tree with singularity exec
singularity exec --bind $PWD:$PWD --bind /babbage:/babbage /babbage/containers/metaphlan.sif /opt/conda/bin/raxmlHPC-PTHREADS-SSE3 -p 1989 -m GTRCAT -f b -T 8 -t /home/psarwar/Desktop/strainphlan/raxml/march/RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre -z /home/psarwar/Desktop/strainphlan/raxml/march/RAxML_bootstrap.bootstrap_tree -n bootmapped_tree
## B longum infantis
# 100 Bootstrap
singularity exec --bind $PWD:$PWD --bind /babbage:/babbage /babbage/containers/metaphlan.sif /opt/conda/bin/raxmlHPC-PTHREADS-SSE3 -p 1989 -m GTRCAT -T 8 -t ./RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre -b 12345 -N 100 -s ./s__Bifidobacterium_longum.StrainPhlAn3_concatenated.aln -n bootstrap_tree
# map the bootstrap numbers onto the ML tree
singularity exec --bind $PWD:$PWD --bind /babbage:/babbage /babbage/containers/metaphlan.sif /opt/conda/bin/raxmlHPC-PTHREADS-SSE3 -p 1989 -m GTRCAT -f b -T 8 -t ./RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre -z ./RAxML_bootstrap.bootstrap_tree -n bootmapped_tree

#Bacteroides fragilis
singularity exec --bind $PWD:$PWD --bind /babbage:/babbage /babbage/containers/metaphlan.sif /opt/conda/bin/raxmlHPC-PTHREADS-SSE3 -p 1989 -m GTRCAT -T 8 -t ./RAxML_bestTree.s__Bacteroides_fragilis.StrainPhlAn3.tre -b 12345 -N 100 -s ./s__Bacteroides_fragilis.StrainPhlAn3_concatenated.aln -n bootstrap_tree
# map the bootstrap numbers onto the ML tree
singularity exec --bind $PWD:$PWD --bind /babbage:/babbage /babbage/containers/metaphlan.sif /opt/conda/bin/raxmlHPC-PTHREADS-SSE3 -p 1989 -m GTRCAT -f b -T 8 -t ./RAxML_bestTree.s__Bacteroides_fragilis.StrainPhlAn3.tre -z ./RAxML_bootstrap.bootstrap_tree -n bootmapped_tree

# raxml run from strainphlan
/opt/conda/bin/raxmlHPC-PTHREADS-SSE3 -p 1989 -m GTRCAT -T 8 -w /home/deniz/Repos/DenizStrains/MSA_longumlongum -s MSA_longumlongum/./s__Bifidobacterium_longum.StrainPhlAn3_concatenated.aln -n s__Bifidobacterium_longum.StrainPhlAn3.tre
##bootstraping command
/opt/conda/bin/raxmlHPC-PTHREADS-SSE3 -p 1989 -m GTRCAT -T 8 -t /home/psarwar/Desktop/strainphlan/raxml/march/RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre -b 12345 -N 100 -s /home/psarwar/Desktop/strainphlan/raxml/march/s__Bifidobacterium_longum.StrainPhlAn3_concatenated.aln -n bootstrap_tree
## for drawing bipartisions on the ML tree with the bootstrapping numbers use
raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15.
-p is the random number seed and can be replaced with -t if i give a starter/seed tree
-T is threads to use
-w unknown but could be the folder
-s is the alignment
-n is just the output names

for boostrapping need to add random seed number -b
also bootstrap replicates via -#
 
# various file transfers for bootstrapping the raxml tree
rsync -avP ./RAxML_bestTree.RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre ada:/home/psarwar/Desktop/strainphlan/raxml/march
rsync -avP ada:/home/psarwar/Desktop/strainphlan/raxml/march/march_bfraglis ./
# Transfer resonance sequenced files

# Call in filenames of all metaphlan outputs
filenames <- read_csv("./raw_data_res/filenames.csv", col_names = FALSE) %>%
  separate(X1, into = c("sample", NA), sep = "_" , remove = FALSE) 
## metaphlan profiles
filenames <- filenames %>% mutate(samfiles = gsub("_profile.tsv", ".sam.bz2", X1))
sequenced_filenames <-sequenced_pruned %>% select(sample, studyID, timepoint) %>%
  left_join(., filenames, by = "sample") 
write_csv(sequenced_filenames, file = "RES_seq_filenames.txt")

#get the metaphlan profiles to your local computer
## example: rsync -avP --dry-run ada: /lovelace/sequencing/processed/mgx/metaphlan/FG00012_S51_profile.tsv current_directory/
rsync -avP --files-from=./filetransfer.txt ada:/lovelace/sequencing/processed/mgx/metaphlan/ ./res_metaphlan_outputs_sequenced_filtered

##get samfiles to ada folder
##upload filetransfer file to ada
rsync -avP ./filetransfer.txt ada:/home/psarwar/Desktop
rsync -avP --dry-run --files-from=./filetransfer.txt /lovelace/sequencing/processed/mgx/metaphlan/ ./strainphlan/res_samfiles

##all files to transfer
all <- sequenced_filenames$samfiles
#list of files that synced during dry-run
files_available <- c("FG00048_S85.sam.bz2", "FG00058_S28.sam.bz2", 
                     "FG00065_S29.sam.bz2", "FG00077_S55.sam.bz2", "FG00079_S79.sam.bz2", 
                     "FG00084_S9.sam.bz2", "FG00085_S33.sam.bz2", "FG00103_S73.sam.bz2", 
                     "FG00147_S26.sam.bz2", "FG00155_S27.sam.bz2", "FG00167_S16.sam.bz2", 
                     "FG00170_S52.sam.bz2", "FG00186_S10.sam.bz2", "FG00211_S35.sam.bz2", 
                     "FG00228_S52.sam.bz2", "FG00297_S60.sam.bz2", "FG00306_S61.sam.bz2", 
                     "FG00315_S10.sam.bz2", "FG00329_S59.sam.bz2", "FG00331_S21.sam.bz2", 
                     "FG00349_S44.sam.bz2", "FG00374_S58.sam.bz2", "FG00375_S69.sam.bz2", 
                     "FG00389_S79.sam.bz2", "FG00395_S89.sam.bz2", "FG00397_S82.sam.bz2", 
                     "FG00398_S95.sam.bz2", "FG00403_S96.sam.bz2", "FG00413_S8.sam.bz2", 
                     "FG00416_S11.sam.bz2", "FG00423_S20.sam.bz2", "FG00434_S28.sam.bz2", 
                     "FG00436_S30.sam.bz2", "FG00445_S41.sam.bz2", "FG00450_S39.sam.bz2", 
                     "FG00453_S48.sam.bz2", "FG00457_S52.sam.bz2", "FG00502_S90.sam.bz2", 
                     "FG00508_S2.sam.bz2", "FG00528_S23.sam.bz2", "FG00533_S21.sam.bz2", 
                     "FG00534_S34.sam.bz2", "FG00538_S35.sam.bz2", "FG00557_S13.sam.bz2", 
                     "FG00561_S9.sam.bz2", "FG00565_S16.sam.bz2", "FG00578_S37.sam.bz2", 
                     "FG00583_S29.sam.bz2", "FG00590_S49.sam.bz2", "FG00599_S35.sam.bz2", 
                     "FG00606_S66.sam.bz2", "FG00622_S54.sam.bz2", "FG00648_S3.sam.bz2", 
                     "FG00660_S14.sam.bz2", "FG00665_S21.sam.bz2", "FG00666_S24.sam.bz2", 
                     "FG00670_S56.sam.bz2", "FG00671_S9.sam.bz2", "FG00673_S28.sam.bz2", 
                     "FG00686_S11.sam.bz2", "FG00706_S12.sam.bz2", "FG00714_S57.sam.bz2", 
                     "FG00719_S62.sam.bz2", "FG00721_S64.sam.bz2", "FG00723_S16.sam.bz2", "FG00738_S72.sam.bz2", "FG00748_S82.sam.bz2", "FG00750_S80.sam.bz2", "FG00765_S95.sam.bz2", "FG00767_S25.sam.bz2", "FG00771_S32.sam.bz2", "FG00775_S22.sam.bz2", "FG00786_S67.sam.bz2", "FG00794_S77.sam.bz2", "FG00795_S84.sam.bz2", "FG00797_S87.sam.bz2", "FG00811_S49.sam.bz2", "FG00839_S81.sam.bz2", "FG00843_S6.sam.bz2", "FG00857_S19.sam.bz2", "FG00864_S28.sam.bz2", "FG00878_S42.sam.bz2", "FG00885_S45.sam.bz2", "FG00896_S56.sam.bz2", "FG00900_S60.sam.bz2", "FG00901_S55.sam.bz2", "FG00905_S22.sam.bz2", "FG00953_S72.sam.bz2", "FG00992_S62.sam.bz2", "FG02272_S23.sam.bz2", "FG02278_S17.sam.bz2", "FG02318_S17.sam.bz2", "FG02342_S58.sam.bz2", "FG02390_S2.sam.bz2", "FG02401_S42.sam.bz2", "FG02406_S88.sam.bz2", "FG02458_S11.sam.bz2", "FG02520_S36.sam.bz2", "FG02555_S61.sam.bz2", "FG02608_S89.sam.bz2", "FG02612_S4.sam.bz2", "FG02622_S34.sam.bz2", "FG02626_S33.sam.bz2", "FG02634_S17.sam.bz2", "FG02644_S6.sam.bz2", "FG02650_S41.sam.bz2", "FG02689_S18.sam.bz2")

files_need <- setdiff(all,files_available)
write_lines(files_need, file = "res_missing_samfiles.txt")

##rsync for the humann files
#directory on ada
/lovelace/sequencing/processed/mgx/humann/main

rsync -avP --files-from=./new_filenames_march_res.txt ada:/lovelace/sequencing/processed/mgx/metaphlan ./
./merge_metaphlan_tablesv2.py *_profile.tsv > marchres_merged_relab.tsv

## get the pathabundance files
rsync -avP --files-from=./pathabundance_files.txt ada:/lovelace/sequencing/processed/mgx/humann/main ./
./merge_metaphlan_tablesv2.py *_profile.tsv > marchres_merged_relab.tsv
rsync -avP --files-from=./pathcoverage.txt ada:/lovelace/sequencing/processed/mgx/humann/main ./

##Join the humann output tables
humann_join_tables --input /Users/psarwar/VKC_Lab/CodingRelated/MARCH_RESONANCE/humann_pathabundance --output humann_genefamilies.tsv --file_name genefamilies_relab

humann_join_tables --input /Users/psarwar/VKC_Lab/CodingRelated/MARCH_RESONANCE/humann_pathabundance --output humann_pathabundance.tsv --file_name pathabundance

humann_join_tables --input /Users/psarwar/VKC_Lab/CodingRelated/MARCH_RESONANCE/humann_kos --output humann_kos_joined.tsv --file_name kos
humann_join_tables --input /Users/psarwar/VKC_Lab/CodingRelated/MARCH_RESONANCE/humann_pfams --output humann_pfams_joined.tsv --file_name pfams
humann_join_tables --input /Users/psarwar/VKC_Lab/CodingRelated/MARCH_RESONANCE/humann_ecs --output humann_ecs_joined.tsv --file_name ecs


cp -R

###download all of Deniz new trees
rsync -avP ada:/home/deniz/Repos/DenizStrains/55hmo_output/RAxML_bestTree.s__Bifidobacterium_infantis.StrainPhlAn3.tre ./

rsync -avP ada:/home/deniz/Repos/DenizStrains/longum_output_allstrains/RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre ./

rsync -avP ada:/home/deniz/Repos/DenizStrains/real_output/RAxML_bestTree.s__Bacteroides_fragilis.StrainPhlAn3.tre ./

rsync -avP ada:/home/deniz/Repos/DenizStrains/real_output/RAxML_bestTree.s__Bifidobacterium_longum.StrainPhlAn3.tre ./

##download regrouped samples
rsync -avP --files-from=./marchres_filenames_pfam.txt ada:/lovelace/sequencing/processed/mgx/humann/regroup ./
rsync -avP --files-from=./marchres_filenames_ecs.txt ada:/lovelace/sequencing/processed/mgx/humann/regroup ./
rsync -avP --files-from=./marchres_filenames_kos.txt ada:/lovelace/sequencing/processed/mgx/humann/regroup ./

tar czvf march.tar.gz ./MARCH_RESONANCE
