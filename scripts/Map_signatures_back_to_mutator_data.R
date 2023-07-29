########load genome
library(BSgenome)
library(BSgenome.Ecoli.NCBI.11strains)
library(MutationalPatterns)
library(Biostrings)
library(ggplot2)
library(reshape2)

#We load signatures
setwd("../signatures/")
#Read signatures from csv files
sig6_NMF<-read.csv("Ecoli_signatures6_NMF_no_norm.csv",
                  sep =" ")
sig6_data<-read.csv("Ecoli_signatures6_from_data.csv",
                    sep =" ")

##########
#Let's read the vcf data again
#First we download all vcfs and build the mutational profiles
###load vcfs
setwd("../vcf")
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
#######
######work with data
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)
summary(vcfs)
####Build mut matrix
mut_mat <- mut_matrix(vcfs, ref_genome)
#############################################

setwd("../Figures/")
#############################################
#Make a function that will do everything
sig_contributions<-function(df,name,...){
  sig_NMF<-as.matrix(df)
  nmf_fit_res<-fit_to_signatures(mut_mat,sig_NMF)
  select <- which(rowSums(nmf_fit_res$contribution) > 0.001)
  NMFContributionPlot<-plot_contribution(nmf_fit_res$contribution[select,], sig_NMF[,select],
                                         coord_flip = FALSE, mode ="relative") +
    theme_classic()+
    theme(axis.text.x = element_text(angle =90),
          legend.title = element_blank()) +
    scale_fill_manual(values = c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"))
  NMFContributionPlot
  ggsave(paste("Fig",name,"a_contributions_allMA_NMF_",name,"signatures.svg",sep=""), 
         plot= NMFContributionPlot,
         width = 19, height = 19, dpi=350, units = "cm")
  ##
  #Let's calculate how well reconstructed profile fit to original one
  cos_sim_ori_rec_NMF <- cos_sim_matrix(mut_mat, nmf_fit_res$reconstructed)
  # extract cosine similarities per sample between original and reconstructed
  cos_sim_ori_rec_NMF <- as.data.frame(diag(cos_sim_ori_rec_NMF))
  # Adjust data frame for plotting with gpplot
  colnames(cos_sim_ori_rec_NMF) = "cos_sim"
  cos_sim_ori_rec_NMF$sample = row.names(cos_sim_ori_rec_NMF)
  
  # Make barplot
  Explained_variance_allMA_NMFSigns<-ggplot(cos_sim_ori_rec_NMF, aes(y=cos_sim, x=sample)) +
    geom_bar(stat="identity", fill = "skyblue4") +
    coord_cartesian(ylim=c(0.4, 1)) +
    ylab("Cosine similarity\n original VS reconstructed") +
    xlab("") +
    theme_classic() +
    theme(panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          axis.text.x = element_text(angle=90, size = 12)) +
    # Add cut.off line
    geom_hline(aes(yintercept=.95))
  Explained_variance_allMA_NMFSigns
  
  ggsave(paste("Fig",name,"_b_explained_variance_allMA_",name,"signatures.svg",sep=""), 
         plot= Explained_variance_allMA_NMFSigns,
         width = 19, height = 19, dpi=350, units = "cm")
}
###
#Run
sig_contributions(sig6_NMF,"6_NMF")
sig_contributions(sig6_data,"6_data")


########
signatures_data<-as.matrix(sig6_data)

#New in this version as in Alexandrov, 2020
#Let's do it more strictly
strict_refit <- fit_to_signatures_strict(mut_mat, signatures_data, max_delta = 0.004)
fig_list <- strict_refit$sim_decay_fig
fig_list[[1]]
fit_res_strict <- strict_refit$fit_res
data_strict<-plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + scale_fill_manual(values = c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"))
data_strict
ggsave("Fig7_signatures_data_mutator_contribution_strict.svg",plot=data_strict,
       width = 19, height = 19, dpi=350, units = "cm")

#New
###Now let's do it with bootstraps
contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
                                               signatures_data,
                                               n_boots = 100,
                                               method = "strict"
)

data_boots<-plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "dotplot")
ggsave("Fig8_signatures_data_mutator_contribution_boots.svg",plot=data_boots,
       path="Figures",width = 19, height = 19, dpi=350, units = "cm")


####Check signatures orthogonality
cos_sim_sigNMF<-cos_sim_matrix(signatures_data,signatures_data)

library(pheatmap)
Fig_CosineSimilaritySIgnsNMF<-pheatmap(cos_sim_sigNMF,
                                       border_color = NA,
                                       display_numbers =T,
                                       number_format = "%.1f",
                                       fontsize = 12)
ggsave("Fig9_signatures6_from_data_cosine_similarity.svg", plot= Fig_CosineSimilaritySIgnsNMF, 
       path="Figures",
       width = 19, height = 19, dpi=350, units = "cm")