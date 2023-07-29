########load genome
library(BSgenome)
library(BSgenome.Ecoli.NCBI.11strains)
library(MutationalPatterns)
library(Biostrings)
library(ggplot2)
library(reshape2)

##########################
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
####Normalization matrix for mutators and non-mutators
genome<-getBSgenome(ref_genome)
######work with data
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)
summary(vcfs)
####Build mut matrix
mut_mat <- mut_matrix(vcfs, ref_genome)
########################
#Do the same for non-mutators
nm_sample_names<-c("non-mutator","wt_F","ED1a","IAI1")
nm_vcf_files<-c("LTEE/nonmutator_allSNPs.vcf",
                "Foster2015/Foster2015_nonmut.vcf",
                "Foster2015/Foster2015_nonmut_ED1a.vcf",
                "Foster2015/Foster2015_nonmut_IAI1.vcf")

nm_vcfs <- read_vcfs_as_granges(nm_vcf_files, nm_sample_names, genome = ref_genome)

nm_mut_mat <- mut_matrix(nm_vcfs, ref_genome)
#intergenic
##non mutators intergenic
nmi_sample_names<-c("nm_i","wt_F_i", "ED1a_i","IAI1_i")
nmi_vcf_files<-c("LTEE/nonmutator_intergenic.vcf",
                 "Foster2015/Foster2015_nonmut_intergenic.vcf",
                 "Foster2015/Foster2015_nonmut_ED1a_intergenic.vcf",
                 "Foster2015/Foster2015_nonmut_IAI1_intergenic.vcf")
nmi_vcfs <- read_vcfs_as_granges(nmi_vcf_files, nmi_sample_names, genome = ref_genome)
nmi_mut_mat <- mut_matrix(nmi_vcfs, ref_genome)
nmi<-rowSums(nmi_mut_mat)

nm_mut_mat_b<-cbind(nm_mut_mat, nmi)

##Merge two matrices
mut_mat_all<-cbind(mut_mat,nm_mut_mat_b)

setwd("../Figures")
mut_mat_s <- mut_mat_all + 0.0001
###############################
##Let's try to extract signatures de novo with NMF
#This part is commented because it is extremely slow and have to be done once
library("NMF")
estimate <- nmf(mut_mat_s, rank = 2:8, method = "brunet",
                nrun = 200, .opt = "v-p")
##plot estimations
NMFestimateplot<-plot(estimate)
NMFestimateplot
ggsave("Fig_NMF_estimate_plot_no_norm.svg", NMFestimateplot,
       width = 25, height = 25, dpi=350, units = "cm")

#We use the initial suggestion to use the number after which the cophenetic starts to decrease
#Which means that the best selection is 5 signatures
#I try 4 and 6 just to see what's going on there
##Extract signatures
nmf_res4 <- extract_signatures(mut_mat_s, rank = 4)
nmf_res5 <- extract_signatures(mut_mat_s, rank = 5)
nmf_res6 <- extract_signatures(mut_mat_s, rank = 6)

##plot 4 and 6 signatures
p4<-plot_96_profile(nmf_res4$signatures, condensed = TRUE)
p4
p5<-plot_96_profile(nmf_res5$signatures, condensed = TRUE)
p5
p6<-plot_96_profile(nmf_res6$signatures, condensed = TRUE)
p6
#Save
ggsave("Fig_NMF_signatures4.svg", p4,
       width = 25, height = 23, dpi=350, units = "cm")
ggsave("Fig_NMF_signatures5.svg", p5,
       width = 25, height = 23, dpi=350, units = "cm")

#Assign names to 6 signatures
colnames(nmf_res6$signatures) <- c("Signature A", "Signature B",
                                  "Signature C", "Signature D","Signature E", "Signature F")
rownames(nmf_res6$contribution) <- c("Signature A", "Signature B",
                                    "Signature C", "Signature D","Signature E", "Signature F")

#Plot  6 signatures
NMF_signatures6_plot<-plot_96_profile(nmf_res6$signatures, condensed = TRUE)
##Save
ggsave("Fig_NMF_signatures6.svg", NMF_signatures6_plot,
       width = 25, height = 23, dpi=350, units = "cm")


#Save signatures to file
# write.table(nmf_res4$signatures,
#             file="Ecoli_signatures4_NMF_no_norm.csv",
#             col.names = TRUE)
# write.table(nmf_res5$signatures,
#             file="Ecoli_signatures5_NMF_no_norm.csv",
#             col.names = TRUE)
write.table(nmf_res6$signatures,
            file="../signatures/Ecoli_signatures6_NMF_no_norm.csv",
            col.names = TRUE)
#

####################
#New in this version
#Let's use Bayesian NMF approach to signature extraction
library(ccfindR)
sc <- scNMFSet(count = mut_mat_s)
#set.seed(4)
estimate_bayes <- vb_factorize(sc, ranks = 1:10, nrun = 200, 
                               progress.bar = FALSE, verbose = 0)
plot(estimate_bayes)
ggsave("Fig_Bayes_NMF_estimate_plot.svg", plot(estimate_bayes),
       path="Figures_no_normalization",
       width = 25, height = 25, dpi=350, units = "cm")
#
# nmf_res_bayes_4 <- extract_signatures(mut_mat_s, rank = 4,nmf_type = "variational_bayes")
# nmf_res_bayes_5 <- extract_signatures(mut_mat_s, rank = 5,nmf_type = "variational_bayes")
nmf_res_bayes_6 <- extract_signatures(mut_mat_s, rank = 6,nmf_type = "variational_bayes")

# p4b<-plot_96_profile(nmf_res_bayes_4$signatures, condensed = TRUE)
# p4b
# p5b<-plot_96_profile(nmf_res_bayes_5$signatures, condensed = TRUE)
# p5b
p6b<-plot_96_profile(nmf_res_bayes_6$signatures, condensed = TRUE)
p6b

write.table(nmf_res_bayes_6$signatures,
            file="../signatures/Ecoli_signatures6_Bayes_NMF_no_norm.csv",
            col.names = TRUE)

####
#Let's draw beatiful ones
require(ggplot2)
require(reshape2)


forplotting<-as.data.frame(nmf_res_bayes_6$signatures)
forplotting$contexts<-rownames(forplotting)
forplotting$contexts <- factor(forplotting$contexts,levels = forplotting$contexts)

forplotting_melted<-melt(forplotting)
forplotting_melted$subst<-substr(forplotting_melted$contexts,3,5)

p<-ggplot(data=forplotting_melted, aes(x=contexts, y=value, fill = subst)) +
  geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  theme(axis.text.x = element_text(size=6, angle = 90, hjust = 1)) +
  facet_wrap( ~ variable, nrow = 6)
p
ggsave("Fig_BayesNMF_signatures6_beatiful.svg", p,
       width = 25, height = 28, dpi=350, units = "cm")

#Same for just NMF
# sig6_df<-read.csv("Ecoli_signatures6_NMF_no_norm_20201208.csv",
#                   sep =" ")

forplotting<-sig6_df
forplotting$contexts<-rownames(forplotting)
forplotting$contexts <- factor(forplotting$contexts,levels = forplotting$contexts)

forplotting_melted<-melt(forplotting)
forplotting_melted$subst<-substr(forplotting_melted$contexts,3,5)

p<-ggplot(data=forplotting_melted, aes(x=contexts, y=value, fill = subst)) +
  geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  theme(axis.text.x = element_text(size=6, angle = 90, hjust = 1)) +
  facet_wrap( ~ variable, nrow = 6)
p
ggsave("Fig_NMF_signatures6_beatiful.svg", p,
       width = 25, height = 28, dpi=350, units = "cm")



