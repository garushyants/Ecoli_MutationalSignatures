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
#######
####Normalization matrix for mutators and non-mutators
genome<-getBSgenome(ref_genome)
#MG1655
NC_000913_trinucfreqs_all<-as.matrix(trinucleotideFrequency(genome$NC_000913)+
                                       trinucleotideFrequency(reverseComplement(genome$NC_000913)))
#MG1655 Old for Foster data
NC_000913_F_trinucfreqs_all<-as.matrix(trinucleotideFrequency(genome$NC_000913_F)+
                                         trinucleotideFrequency(reverseComplement(genome$NC_000913_F)))
#REL606
NC_012967_trinucfreqs_all<-as.matrix(trinucleotideFrequency(genome$NC_012967)+
                                       trinucleotideFrequency(reverseComplement(genome$NC_012967)))
#ED1a
NC_011741_trinucfreqs_all<-as.matrix(trinucleotideFrequency(genome$NC_011741)+
                                       trinucleotideFrequency(reverseComplement(genome$NC_011741)))
#IAI1
NC_011745_trinucfreqs_all<-as.matrix(trinucleotideFrequency(genome$NC_011745)+
                                       trinucleotideFrequency(reverseComplement(genome$NC_011745)))
include_list<-c(rep(c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT",
                      "TCA","TCC","TCG","TCT"),3),
                rep(c("ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT",
                      "TTA","TTC","TTG","TTT"),3))
control_freqs_matrix<-matrix(c(rep(NC_000913_trinucfreqs_all[include_list,],4),
                               rep(NC_012967_trinucfreqs_all[include_list,],6),
                               rep(NC_000913_F_trinucfreqs_all[include_list,],11)),
                             96,21)
rownames(control_freqs_matrix)<-include_list

######work with data
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)
summary(vcfs)
####Build mut matrix
mut_mat <- mut_matrix(vcfs, ref_genome)
mut_mat_norm<-mut_mat / control_freqs_matrix #normalization by contexts
###Normalize
mut_mat_contr = apply(mut_mat_norm,2,function(x){x/sum(x)})

MMR<-rowSums(mut_mat_contr[,c("Ara-2","Ara-4","Ara+3","mutL",
                                  "mutH","mutS","Ara-3","mutL_L")])/8
mutT<-rowSums(mut_mat_contr[,c("Ara-1","Ara+6","mutT_F")])/3
mutY<-rowSums(mut_mat_contr[,c("mutY","mutMmutY_F","mutY_F")])/3
dnaQ<-mut_mat_contr[,"dnaQ_N"]
nth_nei<-mut_mat_contr[,"nth_nei_F"]
polB<-rowSums(mut_mat_contr[,c("polB_F","umuDC_dinB_F")])/2

mutational_signatures_no_norm<-cbind(MMR,mutT,mutY, dnaQ,nth_nei,polB)
mutational_signatures = apply(mutational_signatures_no_norm,2,function(x){x/sum(x)})

setwd("../signatures/")
##Write signatures to the file
write.table(mutational_signatures,
            file="Ecoli_signatures6_from_data.csv",
            col.names = TRUE)

##############Plot these
#Draw beatiful figure with 6 signatures from data
setwd("../Figures/")
forplotting<-as.data.frame(mutational_signatures)
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
ggsave("Fig_signatures6_from_data_beatiful.svg", p,
       width = 25, height = 28, dpi=350, units = "cm")