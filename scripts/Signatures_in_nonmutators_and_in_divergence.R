########load genome
library(BSgenome)
library(BSgenome.Ecoli.NCBI.11strains)
library(MutationalPatterns)
library(Biostrings)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggpubr)

#load signatures
setwd("../signatures/")
#Read signatures from csv files
sig6_NMF<-read.csv("Ecoli_signatures6_NMF_no_norm.csv",
                  sep =" ")
sig6_data<-read.csv("Ecoli_signatures6_from_data.csv",
                    sep =" ")


ref_genome<-"BSgenome.Ecoli.NCBI.11strains"

##########
#Let's read the nonmutator data
###load vcfs
setwd("../vcf/")
nm_sample_names<-c("non-mutator","wt_F","ED1a","IAI1")#, "nm_i","wt_F_i", "ED1a_i","IAI1_i")
nm_vcf_files<-c("LTEE/nonmutator_allSNPs.vcf",
                "Foster2015/Foster2015_nonmut.vcf",
                "Foster2015/Foster2015_nonmut_ED1a.vcf",
                "Foster2015/Foster2015_nonmut_IAI1.vcf")
#
nm_vcfs <- read_vcfs_as_granges(nm_vcf_files, nm_sample_names, genome = ref_genome, predefined_dbs_mbs =T)

##non mutators intergenic
nmi_sample_names<-c("nm_i","wt_F_i", "ED1a_i","IAI1_i")
nmi_vcf_files<-c("LTEE/nonmutator_intergenic.vcf",
                 "Foster2015/Foster2015_nonmut_intergenic.vcf",
                 "Foster2015/Foster2015_nonmut_ED1a_intergenic.vcf",
                 "Foster2015/Foster2015_nonmut_IAI1_intergenic.vcf")
nmi_vcfs <- read_vcfs_as_granges(nmi_vcf_files, nmi_sample_names, genome = ref_genome, predefined_dbs_mbs =T)


#############################################

setwd("../Figures/")
###First
type_occurrences_nm <- mut_type_occurrences(nm_vcfs, ref_genome)
NM_spectrum_one <- plot_spectrum(type_occurrences_nm,legend = TRUE)

ggsave("FigS2_nonmutators_one_spectrum.svg", plot= NM_spectrum_one, 
       width = 19, height = 19, dpi=350, units = "cm")
ggsave("FigS2_nonmutators_one_spectrum.png", plot= NM_spectrum_one, 
       width = 19, height = 19, dpi=700, units = "cm")

##Only nmi
type_occurrences_nmi <- mut_type_occurrences(nmi_vcfs, ref_genome)
NMI_spectrum_int_one <- plot_spectrum(type_occurrences_nmi,legend = TRUE)

NMI_spectrum_int_sep <- plot_spectrum(type_occurrences_nmi,by=rownames(type_occurrences_nmi),
                                      legend = TRUE, condensed = T)

ggsave("FigS3_nonmutators_intergenic_one_spectrum.svg", plot= NMI_spectrum_int_one, 
       width = 19, height = 19, dpi=350, units = "cm")
ggsave("FigS3_nonmutators_intergenic_individual.svg", plot= NMI_spectrum_int_sep, 
       width = 19, height = 10, dpi=700, units = "cm")


###First all and then intergenic

nm_mut_mat <- mut_matrix(nm_vcfs, ref_genome)
nm<-rowSums(nm_mut_mat)
forplotting_nm<-as.data.frame(as.matrix(nm))
forplotting_nm$contexts<-rownames(forplotting_nm)
forplotting_nm$contexts <- factor(forplotting_nm$contexts,levels = forplotting_nm$contexts)

forplotting_nm_melted<-melt(forplotting_nm)
forplotting_nm_melted$subst<-substr(forplotting_nm_melted$contexts,3,5)

p_nm<-ggplot(data=forplotting_nm_melted, aes(x=contexts, y=value, fill = subst)) +
  geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c("#2ebaed", "#000000", "#de1c14", "#d4d2d2", "#adcc54", "#f0d0ce")) +
  theme(axis.text.x = element_text(size=6, angle = 90, hjust = 1)) +
  ylab("")+
  facet_wrap( ~ variable, nrow = 6)
p_nm
#
nmi_mut_mat <- mut_matrix(nmi_vcfs, ref_genome)
nmi<-rowSums(nmi_mut_mat)
forplotting_nmi<-as.data.frame(as.matrix(nmi))
forplotting_nmi$contexts<-rownames(forplotting_nmi)
forplotting_nmi$contexts <- factor(forplotting_nmi$contexts,levels = forplotting_nmi$contexts)

forplotting_nmi_melted<-melt(forplotting_nmi)
forplotting_nmi_melted$subst<-substr(forplotting_nmi_melted$contexts,3,5)

p_nmi<-ggplot(data=forplotting_nmi_melted, aes(x=contexts, y=value, fill = subst)) +
  geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c("#2ebaed", "#000000", "#de1c14", "#d4d2d2", "#adcc54", "#f0d0ce")) +
  theme(axis.text.x = element_text(size=6, angle = 90, hjust = 1)) +
  ylab("")+
  facet_wrap( ~ variable, nrow = 6)
p_nmi

#
NM_NMI_single<-ggarrange(NM_spectrum_one,NMI_spectrum_int_one,
                         ncol =2,
                         labels = c("a","b"),
                         align = 'hv',
                         common.legend = T,
                         legend = "right")
NM_NMI_single
#
NM_NMI_96_plot<-ggarrange(p_nm, p_nmi,
          ncol = 1,
          legend = "none",
          labels = c("c","d"),
          align ='hv')
NM_NMI_96_plot

Full_nm_nmi<-ggarrange(NM_NMI_single,
                       NM_NMI_96_plot,
                       ncol =1)
Full_nm_nmi

ggsave("Suppl_figure_nonmutators_spectra_all.svg", plot= Full_nm_nmi, 
       width = 18, height = 20, dpi=700, units = "cm")
######################################
###Now let's do divergence
###Read indata
#all intergenic
setwd("../divergence_data/")
B1_IA_all<-read.csv("B1_IA_all_phy.paths.matrix", header = TRUE, sep = "\t")
B2_IA_all<-read.csv("B2_IA_all_phy.paths.matrix", header = TRUE, sep = "\t")
E_IA_all<-read.csv("E_IA_all_phy.paths.matrix", header = TRUE, sep = "\t")
A_IA_all<-read.csv("A_IA_all_phy.paths.matrix", header = TRUE, sep = "\t")
IA_all_combined_B1_B2<-merge(B1_IA_all, B2_IA_all, by ="X")
IA_all_combined_B1_B2_E<-merge(IA_all_combined_B1_B2, E_IA_all, by ="X")
IA_all_combined<-merge(IA_all_combined_B1_B2_E, A_IA_all, by ="X")
#convergent intergenic regions
B1_IA_c<-read.csv("B1_IA_convergent.paths.matrix", header = TRUE, sep = "\t")
B2_IA_c<-read.csv("B2_IA_convergent.paths.matrix", header = TRUE, sep = "\t")
E_IA_c<-read.csv("E_IA_convergent.paths.matrix", header = TRUE, sep = "\t")
A_IA_c<-read.csv("A_IA_convergent.paths.matrix", header = TRUE, sep = "\t")
IA_c_combined_B1_B2<-merge(B1_IA_c, B2_IA_c, by ="X")
IA_c_combined_B1_B2_E<-merge(IA_c_combined_B1_B2, E_IA_c, by ="X")
IA_c_combined<-merge(IA_c_combined_B1_B2_E, A_IA_c, by ="X")

##
IA_combined<-merge(IA_all_combined,IA_c_combined)
##select columns that contain all baseml reconstructions

#this is main df with raw data
IA_combined_all_baseml<- IA_combined %>% select(contains(".all"))
rownames(IA_combined_all_baseml) <- IA_combined$X

#Save to file
#write.csv(IA_combined_all_baseml, file = "IA_combined_all_baseml.csv",quote = F)

#draw profile
#from absolute to relative
IA_combined_all_norm = apply(IA_combined_all_baseml,2,function(x){x/sum(x)})
#
forplotting<-as.data.frame(IA_combined_all_norm)
forplotting$subst<-substr(rownames(forplotting),3,5)
forplotting$contexts <-rownames(forplotting)
forplotting_o<-forplotting[order(forplotting$subst,
                                 forplotting$contexts),order(names(forplotting))]
forplotting_o$contexts <- factor(forplotting_o$contexts,levels = forplotting_o$contexts)
forplotting_melted<-melt(forplotting_o)
# forplotting_melted_o <- forplotting_melted[order(forplotting_melted$subst,
#                                                  forplotting_melted$contexts),]

p<-ggplot(data=forplotting_melted, aes(x=contexts, y=value, fill = subst)) +
  geom_bar(stat="identity") + theme_bw() + 
  scale_fill_manual(values=c("#2ebaed", "#000000", "#de1c14", "#d4d2d2", "#adcc54", "#f0d0ce")) +
  theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1)) +
  facet_wrap( ~ variable, nrow = 4)
p

################################
####Do 6 for divergence
DoSixProfiles<-IA_combined_all_baseml
DoSixProfiles$Change<-substr(rownames(DoSixProfiles),3,5)
PreforplottingSix<-DoSixProfiles[,c(1,9)] %>% group_by(Change) %>% 
  summarise(B1=sum(B1_IA_all_phy.paths.all))
PreforplottingSix$perc<-PreforplottingSix$B1/sum(PreforplottingSix$B1)
forplottingSix_m<-melt(PreforplottingSix[,c(1,3)])

pS<-ggplot(data=forplottingSix_m, aes(x=Change, y=value, fill = Change)) +
  geom_bar(stat="identity") + theme_bw() + 
  scale_fill_manual(values=c("#2ebaed", "#000000", "#de1c14", "#d4d2d2", "#adcc54", "#f0d0ce")) +
  theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1))
pS

#################################
setwd("../Figures/")
ggsave("Suppl_figure_5_All_mutational_profiles_divergence.svg", 
       plot= p,
       width = 25, height = 16, dpi=350, units = "cm")

####Draw only for one case
forplotting_all_intergenic<-forplotting_o[,c(1,7,10)]
colnames(forplotting_all_intergenic)<-c("B1","contexts","subst")
forplotting_all_intergenic_melted<-melt(forplotting_all_intergenic)

p_B1_intergenic<-ggplot(data=forplotting_all_intergenic_melted, aes(x=contexts, y=value, fill = subst)) +
  geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c("#2ebaed", "#000000", "#de1c14", "#d4d2d2", "#adcc54", "#f0d0ce")) +
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1)) +
  facet_wrap( ~ variable, nrow = 4)
p_B1_intergenic
##############################################

#Now it's time for signatures contributions
#Let's combine non-mutators and divergence
#sort IA_combined
IA_combined_all_baseml_sorted<-IA_combined_all_baseml[order(match(rownames(IA_combined_all_baseml), row.names(nmi_mut_mat))), , 
                                                      drop = FALSE]
NM_div_mat<-cbind(nmi,IA_combined_all_baseml_sorted)
############################


#########
#More strict parameters

sigs<-as.matrix(sig6_data)

strict_refit <- fit_to_signatures_strict(NM_div_mat, sigs, max_delta = 0.02)
fig_list <- strict_refit$sim_decay_fig
fig_list[[1]]
fit_res_strict <- strict_refit$fit_res
data_strict<-plot_contribution(fit_res_strict$contribution,sigs,
                               coord_flip = FALSE,
                               mode = "relative"
) + scale_fill_manual(values = c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac")) +
  theme_classic()+
  theme(axis.text.x = element_text(angle =90),
        legend.title = element_blank())
data_strict

####Figure only all intergenic and all intergenic non-mutator
#Explained variance
NM_all_intergenic <- NM_div_mat[,c(1:5)]

cos_sim_ori_rec <- cos_sim_matrix(NM_all_intergenic, fit_res_strict$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)

Explained_variance_divergence_str<-ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
  geom_bar(stat="identity", aes(fill = sample)) +
  geom_text(aes(label=formatC(cos_sim, digits = 2, format = "f")),
                vjust=0, 
            size =4)+
  coord_cartesian(ylim=c(0.4, 1)) +
  ylab("Cosine similarity\n original VS reconstructed") +
  xlab("") +
  scale_fill_manual(values = c("#ffcc00","#00ff00","#008000",
                               "#0000ff","#dd55ff"))+
  theme_classic() +
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.text.y= element_text(size = 12),
        axis.text.x = element_text(size=12, angle =90),
        legend.position = "none") +
  scale_x_discrete(labels=c("A","B1","B2","E","Non-mutator"))+
  # Add cut.off line
  geom_hline(aes(yintercept=.85), color = "#004529",
             size =1.2)
Explained_variance_divergence_str

#Draw one case with residuals
#nmi
NonMut_orvsrec<-plot_compare_profiles(NM_all_intergenic[, 1],
                      fit_res_strict$reconstructed[, 1],
                      profile_names = c("Nonmutator\noriginal", "Nonmutator\nreconstructed"),
                      condensed = TRUE
) +scale_y_continuous(limits = c()) +
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=7))
NonMut_orvsrec
#B1
B1Mut_orvsrec<-plot_compare_profiles(NM_all_intergenic[, 2],
                      fit_res_strict$reconstructed[, 2],
                      profile_names = c("B1 original", "B1 reconstructed"),
                      condensed = TRUE,

) +scale_y_continuous(limits = c()) +
  theme(axis.text.y = element_text(size=12),
         axis.text.x = element_text(size=7))
B1Mut_orvsrec

###Plot top simple
forplotting_top<-as.data.frame(NM_all_intergenic[,c(1,2)])
names(forplotting_top)<-c("Non-mutator","B1")
forplotting_top$subst<-substr(rownames(forplotting_top),3,5)
forplotting_top$contexts <-rownames(forplotting_top)
forplotting_top_o<-forplotting_top[order(forplotting_top$subst,
                                 forplotting_top$contexts),order(names(forplotting_top))]
forplotting_top_o$contexts <- factor(forplotting_top_o$contexts,levels = forplotting_top_o$contexts)
forplotting_melted<-melt(forplotting_top_o)

Topplot<-ggplot(data=forplotting_melted, aes(x=contexts, y=value, fill = subst)) +
  geom_bar(stat="identity") + theme_bw() + 
  scale_fill_manual(values=c("#2ebaed", "#000000", "#de1c14", "#d4d2d2", "#adcc54", "#f0d0ce")) +
  theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1)) +
  facet_wrap( ~ variable, nrow = 2,scales = "free_y")+
  ylab("# of mutations")
Topplot


# #######
#Boots
contri_boots <- fit_to_signatures_bootstrapped(NM_all_intergenic,
                                               sigs,
                                               n_boots = 100,
                                               max_delta = 0.02,
                                               method = "strict"
)
# 
# 
# data_boots<-plot_bootstrapped_contribution(contri_boots, 
#                                            mode = "relative", 
#                                            plot_type = "dotplot") +
#   scale_y_discrete(labels=c("A","E","B2","B1","Nonmutator"))+
#   theme(#panel.grid.minor.y=element_blank(),
#         #panel.grid.minor.x=element_blank(),
#         axis.text = element_text(size = 12))
# data_boots
##
#Plot in more human friendly way
contri_boots_df<-as.data.frame(contri_boots)
contri_boots_df$name_full<-row.names(contri_boots_df)

spl <- strsplit(as.character(contri_boots_df$name_full), "_")
contri_boots_df$name<-sapply(lapply(spl, head, 1), paste, collapse="_")
contri_boots_df$boot<-sapply(lapply(spl, tail, 1), paste, collapse="_")

contri_boots_df_perc<-as.data.frame(t(apply(t(as.matrix(contri_boots[,c(1:5)])), 2, function(i) i/sum(i))))
contri_boots_df_perc$boot<-contri_boots_df$boot
contri_boots_df_perc$name<-contri_boots_df$name
contri_boots_melted<-melt(contri_boots_df_perc)

######################
#Calculate numbers for paper
#evol
contri_boots_df_perc_evol<-subset(contri_boots_df_perc, 
                                  contri_boots_df_perc$name != "nmi")
mean(contri_boots_df_perc_evol$dnaQ)
sd(contri_boots_df_perc_evol$dnaQ)

contri_boots_df_perc_lab<-subset(contri_boots_df_perc, 
                                  contri_boots_df_perc$name == "nmi")
mean(contri_boots_df_perc_lab$dnaQ)
sd(contri_boots_df_perc_lab$dnaQ)
######################

Signatures_contribution_evolution<-ggplot(data = contri_boots_melted,
       aes(x = variable,y=value, fill = name))+
  geom_point(pch = 21, position = position_jitterdodge(), size =.5,
             alpha = .5,
             aes(color =name))+
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(-0.02,0.7), breaks=seq(0.00, 1.00, by =0.05),
                     name = "Relative contbution")+
  xlab("")+
  labs(fill="", color = "")+
  scale_fill_manual(values = c("#ffcc00","#00ff00","#008000",
                                   "#0000ff","#dd55ff"))+
  scale_colour_manual(values = c("#ffcc00","#00ff00","#008000",
                                   "#0000ff","#dd55ff"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12),
        legend.position = "none")
Signatures_contribution_evolution

# ggsave("Fig3_sign_boxplot.svg", 
#        plot= Signatures_contribution_evolution,
#        path="Figures",
#        width = 26, height = 14, dpi=350, units = "cm")

#Combine
Bottom<-ggarrange(Explained_variance_divergence_str,
                  Signatures_contribution_evolution,
          ncol = 2,
          labels = c("b","c"),
          widths = c(0.4,1))
Bottom
##TopPlot simple

Top<-ggarrange(Topplot,
               labels=c("a"))
  

SupplFig7<-ggarrange(NonMut_orvsrec,
               B1Mut_orvsrec,
               labels = c("a", "b"),
               ncol =1,
               align = 'hv')


Full_fig3<-ggarrange(Top,Bottom,
                     ncol =1,
                     heights = c(1,1),
                     align = 'hv')
Full_fig3

#Save
ggsave("SupplFig7_Real_vs_recontructed_nmi_B1_20221211.svg", 
       plot= SupplFig7,
       width = 28, height = 35, dpi=300, units = "cm")
ggsave("SupplFig7_Real_vs_recontructed_nmi_B1_20221211.png", 
       plot= SupplFig7,
       width = 28, height = 35, dpi=300, units = "cm")

ggsave("Fig3_NM_Divergence_signatures_contributions_20221211.svg", 
              plot= Full_fig3,
              width = 28, height = 25, dpi=300, units = "cm")

###Supplementary Figure 6
NM_all_con <- NM_div_mat[,c(1,6:9)]

cos_sim_con <- cos_sim_matrix(NM_all_con, fit_res_strict$reconstructed[,c(1,6:9)])
cos_sim_con_df <- as.data.frame(diag(cos_sim_con))
colnames(cos_sim_con_df) = "cos_sim"
cos_sim_con_df$sample = row.names(cos_sim_con)

Explained_variance_divergence_str_con<-ggplot(cos_sim_con_df, aes(y=cos_sim, x=sample)) +
  geom_bar(stat="identity", aes(fill = sample)) +
  geom_text(aes(label=formatC(cos_sim, digits = 2, format = "f")),
            vjust=0, 
            size =4)+
  coord_cartesian(ylim=c(0.4, 1)) +
  ylab("Cosine similarity\n original VS reconstructed") +
  xlab("") +
  scale_fill_manual(values = c("#ffcc00","#00ff00","#008000",
                               "#0000ff","#dd55ff"))+
  theme_classic() +
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.text= element_text(size = 12),
        legend.position = "none") +
  scale_x_discrete(labels=c("A","B1","B2","E","Nonmutator"))+
  # Add cut.off line
  geom_hline(aes(yintercept=.85), color = "#004529",
             size =1.2)
Explained_variance_divergence_str_con
#######
#Boots
contri_boots_con <- fit_to_signatures_bootstrapped(NM_all_con,
                                               sigs,
                                               n_boots = 100,
                                               max_delta = 0.02,
                                               method = "strict"
)


data_boots_con<-plot_bootstrapped_contribution(contri_boots_con, 
                                           mode = "relative", 
                                           plot_type = "dotplot") +
  scale_y_discrete(labels=c("A","E","B2","B1","Nonmutator"))+
  theme(#panel.grid.minor.y=element_blank(),
    #panel.grid.minor.x=element_blank(),
    axis.text = element_text(size = 12))
data_boots_con
#Combine
SupplementaryFig6<-ggarrange(Explained_variance_divergence_str_con,
                  data_boots_con,
                  ncol = 2,
                  labels = c("a","b"),
                  widths = c(0.6,1))
#save
ggsave("Suppl_Fig7_NM_div_explained_variance_convergent_sign6_data_strict.svg", 
       plot= SupplementaryFig6,
       width = 28, height = 12, dpi=350, units = "cm")
#
ggsave("Suppl_Fig7_NM_div_explained_variance_convergent_sign6_data_strict.png", 
       plot= SupplementaryFig6,
       width = 28, height = 12, dpi=350, units = "cm")




# ######################################
#Exclude parallelisms
IA_combined_all_baseml_no_par<- IA_combined %>% select(contains(".no_parallelism"))
rownames(IA_combined_all_baseml_no_par) <- IA_combined$X

#draw profile
#from absolute to relative
IA_combined_all_norm_no_par = apply(IA_combined_all_baseml_no_par,2,function(x){x/sum(x)})
#
forplotting<-as.data.frame(IA_combined_all_norm_no_par)
forplotting$subst<-substr(rownames(forplotting),3,5)
forplotting$contexts <-rownames(forplotting)
forplotting_o<-forplotting[order(forplotting$subst,
                                 forplotting$contexts),order(names(forplotting))]
forplotting_o$contexts <- factor(forplotting_o$contexts,levels = forplotting_o$contexts)
forplotting_melted<-melt(forplotting_o)
# forplotting_melted_o <- forplotting_melted[order(forplotting_melted$subst,
#                                                  forplotting_melted$contexts),]

p_no_par<-ggplot(data=forplotting_melted, aes(x=contexts, y=value, fill = subst)) +
  geom_bar(stat="identity") + theme_bw() + 
  scale_fill_manual(values=c("#2ebaed", "#000000", "#de1c14", "#d4d2d2", "#adcc54", "#f0d0ce")) +
  theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1)) +
  facet_wrap( ~ variable, nrow = 4)
p_no_par

ggsave("Suppl_figure_5_All_mutational_profiles_divergence_no_par.svg", 
       plot= p_no_par,
       width = 25, height = 16, dpi=350, units = "cm")

###
IA_combined_all_baseml_no_par_sorted<-IA_combined_all_baseml_no_par[order(match(rownames(IA_combined_all_baseml_no_par), 
                                                                                row.names(nmi_mut_mat))), , 
                                                      drop = FALSE]
NM_div_mat_no_par<-cbind(nmi,IA_combined_all_baseml_no_par_sorted)
#########
#More strict parameters
strict_refit_no_par <- fit_to_signatures_strict(NM_div_mat_no_par, sigs, max_delta = 0.02)
fig_list_no_par <- strict_refit_no_par$sim_decay_fig
fig_list_no_par[[1]]
fit_res_strict_no_par <- strict_refit_no_par$fit_res

# data_strict_no_par<-plot_contribution(fit_res_strict_no_par$contribution,sig_NMF,
#                                coord_flip = FALSE,
#                                mode = "relative"
# ) + scale_fill_manual(values = c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac")) +
#   theme_classic()+
#   theme(axis.text.x = element_text(angle =90),
#         legend.title = element_blank())
# data_strict_no_par

#Explained variance
NM_all_intergenic_no_par <- NM_div_mat_no_par[,c(1:5)]

cos_sim_ori_rec <- cos_sim_matrix(NM_all_intergenic_no_par, 
                                  fit_res_strict_no_par$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)

Explained_variance_divergence_str_no_par<-ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
  geom_bar(stat="identity", aes(fill = sample)) +
  geom_text(aes(label=formatC(cos_sim, digits = 2, format = "f")),
            vjust=0, 
            size =4)+
  coord_cartesian(ylim=c(0.4, 1)) +
  ylab("Cosine similarity\n original VS reconstructed") +
  xlab("") +
  scale_fill_manual(values = c("#ffcc00","#00ff00","#008000",
                               "#0000ff","#dd55ff"))+
  theme_classic() +
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.text.y= element_text(size = 12),
        axis.text.x = element_text(size=12, angle =90),
        legend.position = "none") +
  scale_x_discrete(labels=c("A","B1","B2","E","Non-mutator"))+
  # Add cut.off line
  geom_hline(aes(yintercept=.85), color = "#004529",
             size =1.2)
Explained_variance_divergence_str_no_par
###
#Plot in more human friendly way

contri_boots_no_par <- fit_to_signatures_bootstrapped(NM_all_intergenic_no_par,
                                               sigs,
                                               n_boots = 100,
                                               max_delta = 0.02,
                                               method = "strict"
)

############Complicated plot
contri_boots_df<-as.data.frame(contri_boots_no_par)
contri_boots_df$name_full<-row.names(contri_boots_df)

spl <- strsplit(as.character(contri_boots_df$name_full), "_")
contri_boots_df$name<-sapply(lapply(spl, head, 1), paste, collapse="_")
contri_boots_df$boot<-sapply(lapply(spl, tail, 1), paste, collapse="_")

contri_boots_df_perc<-as.data.frame(t(apply(t(as.matrix(contri_boots[,c(1:5)])), 2, function(i) i/sum(i))))
contri_boots_df_perc$boot<-contri_boots_df$boot
contri_boots_df_perc$name<-contri_boots_df$name
contri_boots_melted<-melt(contri_boots_df_perc)
# 
Signatures_contribution_evolution_no_par<-ggplot(data = contri_boots_melted,
                                          aes(x = variable,y=value, fill = name))+
  geom_point(pch = 21, position = position_jitterdodge(), size =.5,
             alpha = .5,
             aes(color =name))+
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(-0.02,0.7), breaks=seq(0.00, 1.00, by =0.05),
                     name = "Relative contbution")+
  xlab("")+
  labs(fill="", color = "")+
  scale_fill_manual(values = c("#ffcc00","#00ff00","#008000",
                               "#0000ff","#dd55ff"))+
  scale_colour_manual(values = c("#ffcc00","#00ff00","#008000",
                                 "#0000ff","#dd55ff"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12),
        legend.position = "none")
Signatures_contribution_evolution_no_par

Signatures_contribution_evolution_no_par_simp<-plot_bootstrapped_contribution(contri_boots_no_par, 
                                               mode = "relative", 
                                               plot_type = "dotplot") +
  scale_y_discrete(labels=c("A","E","B2","B1","Nonmutator"))+
  theme(#panel.grid.minor.y=element_blank(),
    #panel.grid.minor.x=element_blank(),
    axis.text = element_text(size = 12))
Signatures_contribution_evolution_no_par_simp

####Draw supplementary figure 8
SupplFig9<-ggarrange(Explained_variance_divergence_str_no_par,
          Signatures_contribution_evolution_no_par,
          widths = c(0.5,1),
          labels = c("a","b"),
          align = 'hv')

ggsave("Suppl_figure_9_Mutation_contributions_no_par.svg", 
       plot= SupplFig9,
       width = 25, height = 14, dpi=350, units = "cm")

ggsave("Suppl_figure_9_Mutation_contributions_no_par.png", 
       plot= SupplFig9,
       width = 25, height = 14, dpi=350, units = "cm")
# ######################################