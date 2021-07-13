library("stringr")
library("dplyr")
setwd("~/Documents/Research/cosmic 38/sql")

# load the dataframe
comsic_mut <- read.csv("combine_n_cosmic_filtered.csv")
comsic_mut<- comsic_mut[,-1]

comsic_mut<- comsic_mut[comsic_mut$Mutation.Description != "Unknown", ]

#Choose only canonical genes (to avoid duplicated mutations)
is_split<- function(x){ 
  len_x<- length(unlist(strsplit(x, "_")))
  if (len_x == 1) {return(TRUE)}
  else {{return(FALSE)}}
}
comsic_mut<- comsic_mut[,-19]
canonical<- sapply(comsic_mut$Gene_name, function(x) is_split(x))
comsic_mut_f<- comsic_mut[canonical, ]

#remove Y chromosome ----
comsic_mut_f<- comsic_mut_f[comsic_mut_f$chrom != "Y",]
#filter out sex specific tissues
comsic_mut_f<- comsic_mut_f[!comsic_mut_f$Primary.Site %in% c("testis", "prostate", "placenta", "NS", "ovary", "breast", "cervix",
                                                              "endometrium", "genital_tract", "penis"),]


# Organize the dataset ---- 
#mut<- sapply(comsic_mut_f$HGVSG, function(x) unlist(strsplit(as.character(x), ">"))[2])
#head(mut)
#comsic_mut_f<- cbind.data.frame(comsic_mut_f, mut)
#     
ref<- sapply(comsic_mut_f$HGVSG, function(x) unlist(strsplit(as.character(x), ">"))[1])
ref<- sapply(ref, function(x) str_sub(x,-1,-1))
head(ref)
comsic_mut_f<- cbind.data.frame(comsic_mut_f, ref)


pos<- sapply(comsic_mut_f$'Mutation.genome.position', function(x) unlist(strsplit(x, "-"))[2])
comsic_mut_f<- cbind.data.frame(comsic_mut_f, pos)
table(comsic_mut_f$Mutation.Description)

#run the dndscv analysis on males and females speratly
library(remotes)
remotes::install_github("im3sanger/dndscv")

library("dndscv")


df_male<- comsic_mut_f[comsic_mut_f$gender == "m", c("sample_id", "chrom","pos", "ref", "mut")]
df_female<- comsic_mut_f[comsic_mut_f$gender == "f", c("sample_id", "chrom","pos", "ref", "mut")]

prepare_for_dndnscv<- function(df){
  colnames(df)<- c("sampleID", "chr","pos", "ref", "mut")
  df$chr <- sub("^", "chr", df$chr )
  return(df)
}


df_male<- prepare_for_dndnscv(df_male)
df_female<- prepare_for_dndnscv(df_female)


dndsout_female = dndscv(df_female,refdb="RefCDS_human_GRCh38.p12.rda",cv=NULL,
                        max_muts_per_gene_per_sample = 1000,max_coding_muts_per_sample = 40e3,
                        outmats=T)
print(dndsout_female$sel_cv) 

sel_cv_female  = dndsout_female$sel_cv
setwd("/home/ethel/Documents/Research/cosmic 38/final_code")
write.csv(sel_cv_female, "sel_cv_female_new.csv")

dndsout_male = dndscv(df_male,refdb="RefCDS_human_GRCh38.p12.rda",cv=NULL,
                      max_muts_per_gene_per_sample = 1000,max_coding_muts_per_sample = 40e3,
                      outmats=T)
print(dndsout_male$sel_cv) 

sel_cv_male  = dndsout_male$sel_cv
setwd("/home/ethel/Documents/Research/cosmic 38/final_code")
write.csv(sel_cv_male, "sel_cv_male.csv")

#read the data
setwd("/home/ethel/Documents/Research/cosmic 38/final_code")
sel_cv_male<- read.csv("sel_cv_male.csv")
sel_cv_female<- read.csv("sel_cv_female.csv")

#match the rows of both tables
table(sel_cv_female$gene_name == sel_cv_male$gene_name)
m<- match(sel_cv_female$gene_name, sel_cv_male$gene_name)

sel_cv_male<- sel_cv_male[m,]

# plot the results
plot((sel_cv_male$wmis_cv)~(sel_cv_female$wmis_cv))
text(sel_cv_male$wmis_cv~sel_cv_female$wmis_cv, labels=sel_cv_female$gene_name, cex=0.6, font=1, pos = 4)

plot(log(sel_cv_female$wmis_cv), log(sel_cv_male$wmis_cv))

plot((sel_cv_female$wnon_cv), (sel_cv_male$wnon_cv))
plot(log(sel_cv_female$wnon_cv), log(sel_cv_male$wnon_cv))

t.test((sel_cv_female$wmis_cv), (sel_cv_male$wmis_cv))
t.test((sel_cv_female$wnon_cv), (sel_cv_male$wnon_cv))

cor.test((sel_cv_female$wmis_cv), (sel_cv_male$wmis_cv))
cor.test((sel_cv_female$wnon_cv), (sel_cv_male$wnon_cv))

effect<- (sel_cv_female$wmis_cv)- (sel_cv_male$wmis_cv)

# calculation of the confidance interval of the dN/dS ratio
setwd("/home/ethel/Documents/Research/cosmic 38/final_code")

ci_female<-geneci(dndsout_female)
write.csv(ci_female, "ci_female.csv")

ci_male<-geneci(dndsout_male)
write.csv(ci_male, "ci_male.csv")

# interaction test to determine genes with different dN/dS ratio

sd_calc<- function(low, high) {
  df<- cbind(high, low)
  sd <- apply(df, 1, function(x) (x[1]-x[2])/3.92)
  return(sd)
}

sd_male_mis<- sd_calc(ci_male$mis_low,ci_male$mis_high)
sd_female_mis<- sd_calc(ci_female$mis_low,ci_female$mis_high) 

d_mis<- ci_male$mis_mle - ci_female$mis_mle

d_mis_sd<- sqrt((sd_male_mis)^2 + (sd_female_mis)^2)

z_mis <- d_mis/ d_mis_sd

p_mis = 2*pnorm(-abs(z_mis))

plot((d_mis), -log10(p_mis), cex = 0.5)
##nonsense
sd_male_non<- sd_calc(ci_male$tru_low,ci_male$tru_high)
sd_female_non<- sd_calc(ci_female$tru_low,ci_female$tru_high) 

d_non<- ci_male$tru_mle - ci_female$tru_mle

d_non_sd<- sqrt((sd_male_non)^2 + (sd_female_non)^2)

z_non <- d_non/ d_non_sd

p_non = 2*pnorm(-abs(z_non))

interaction_p<- cbind.data.frame(ci_male$gene,d_mis, d_mis_sd, p_mis,d_non, d_non_sd, p_non )
colnames(interaction_p)[c(1,3,6)]<- c("gene", "mis_sd", "non_sd")

#choose only genes that passed filtration
fisher_df<- read.csv("fisher test results.csv")
m<- interaction_p$gene %in% fisher_df$gene
interaction_p<- interaction_p[m,]


fdr_mis<- p.adjust(interaction_p$p_mis, method = "fdr")
interaction_p<- cbind(interaction_p,fdr_mis )

fdr_non<- p.adjust(interaction_p$p_non, method = "fdr")
interaction_p<- cbind(interaction_p,fdr_non )


write.csv(interaction_p, "interaction_p_val.csv")

plot(log10(ci_female$mis_mle[m]/ci_male$mis_mle[m]), -log10(interaction_p$p_mis), cex = 0.3, main = "intercation missense", 
     xlab= "log(dnds_female/dnds_male)", ylab = "-log10(interaction p_value)")
text(log10(ci_female$mis_mle[m]/ci_male$mis_mle[m]), -log10(interaction_p$p_mis), ci_female$gene[m], cex=0.66)

plot(log(ci_female$tru_mle[m]/ci_male$tru_mle[m]), -log(interaction_p$p_non), cex = 0.3, main = "intercation nonsense", 
     xlab= "log(dnds_female/dnds_male)", ylab = "log(interaction p_value")
text(log(ci_female$tru_mle[m]/ci_male$tru_mle[m]), -log(interaction_p$p_non), ci_female$gene[m], cex=0.66)

