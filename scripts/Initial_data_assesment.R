########load packages
library(BSgenome)
library(BSgenome.Ecoli.NCBI.11strains)
library(MutationalPatterns)
library(Biostrings)
library(ggplot2)
library(reshape2)


##########################
#We use here all mutational profiles from Deepa, LTEE, and various papers from Foster lab (Lee, 2012; Foster, 2015 and Niccum, 2018)
###load vcfs
setwd("../vcf/")
sample_names <- c("mutH","mutL","mutS","mutY",
                  "Ara-1","Ara-2","Ara-3","Ara-4",
                  "Ara+3","Ara+6",
                  "dnaQ_N",
                  "mutMmutY_F","mutT_F","mutY_F","uvrA_F",
                  "mutL_L",
                  "polB_F","umuDC_dinB_F",
                  "nfi_F","xthA_nfo_F","nth_nei_F")
vcf_files <- c("Deepa/mutH.vcf","Deepa/mutL.vcf",
               "Deepa/mutS.vcf","Deepa/mutY.vcf",
               "LTEE/Ara-1.vcf","LTEE/Ara-2.vcf","LTEE/Ara-3.vcf",
               "LTEE/Ara-4.vcf","LTEE/Ara+3.vcf","LTEE/Ara+6.vcf",
               "Niccum_2018/dnaQ_PFM163_combined_filt.vcf",
               "Foster2015/Foster2015_mutMmutY.vcf", 
               "Foster2015/Foster2015_mutT.vcf",
               "Foster2015/Foster2015_mutY.vcf",
               "Foster2015/Foster2015_uvrA.vcf",
               "Lee2012/Lee2012_mutL.vcf",
               "Foster2015/Foster15_polB.vcf","Foster2015/Foster15_umuDC_dinB.vcf",
               "Foster2015/Foster15_nfi.vcf","Foster2015/Foster15_xthA_nfo.vcf","Foster2015/Foster15_nth_nei.vcf")
###genome
ref_genome<-"BSgenome.Ecoli.NCBI.11strains"

######work with data
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)
summary(vcfs)

##########################
#Plotting
#set path
setwd("../Figures/")

###Phase one: one-nucleotide contexts
#simple initial plot for MA lines
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
MA_spectrum <- plot_spectrum(type_occurrences,by=rownames(type_occurrences),legend = TRUE)

ggsave("Fig1_All_MA_spectrum.svg", plot= MA_spectrum, 
       width = 25, height = 19, dpi=350, units = "cm")

##Phase two: now let's switch for 3nucleotides contexts
####Build mut matrix
mut_mat <- mut_matrix(vcfs, ref_genome)

######################################
###Finally plotting

#Let's plot cosine similarities for All sample
cos_sim_all_samples<-cos_sim_matrix(mut_mat,mut_mat)

library(pheatmap)
Fig_CosineSimilarityAllMA<-pheatmap(cos_sim_all_samples,
         border_color = NA,
         display_numbers =T,
         number_format = "%.1f")
#Figure with similarities
ggsave("Fig2_CosineSimilarity_All_MA_samples.svg", plot= Fig_CosineSimilarityAllMA, 
       width = 20, height = 20, dpi=350, units = "cm")

############
####Now let's save all the mutational profiles for the groups: MMR, mutT,mutY, others
hclust=cluster_signatures(mut_mat, method = "complete")
hclust[["order"]]
###Order data frame, so it will be easier to get groups
mut_mat_reorder <-as.data.frame(mut_mat)[,hclust[["order"]]]
#Check groups after this procedure, the order can change

#mutT
part1<-as.matrix(mut_mat_reorder[,c(1:3)])
plot_mutT<-plot_96_profile(part1)
ggsave("FigS1a_MutT_mutator_samples_profile.svg", plot= plot_mutT,
       width = 18, height = 12, dpi=350, units = "cm")
#mutY
part2<-as.matrix(mut_mat_reorder[,c(4:6)])
plot_mutY<-plot_96_profile(part2)
ggsave("FigS1b_MutY_mutator_samples_profile.svg", plot= plot_mutY,
       width = 18, height = 12, dpi=350, units = "cm")
#MMR
part3<-as.matrix(mut_mat_reorder[,c(7:12)])
plot_MMR<-plot_96_profile(part3)
plot_MMR
ggsave("FigS1c_pol_BER_NER_mutator_profile.svg", plot= plot_MMR,
       width = 18, height = 23, dpi=350, units = "cm")
#Else
part4<-as.matrix(mut_mat_reorder[,c(13:21)])
plot_else<-plot_96_profile(part4)
plot_else
ggsave("FigS1d_dnaQ_MMR_mutator_samples_profile.svg", plot= plot_else,
       width = 18, height = 30, dpi=350, units = "cm")

###
#Plot example for presentation
example<-as.matrix(mut_mat_reorder[,c(3,4,7,10,13,17)])
plot_example<-plot_96_profile(example)
ggsave("FigS1e_example_mutator_samples_profile.svg", plot= plot_example,
       width = 18, height = 17, dpi=350, units = "cm")
##################################

