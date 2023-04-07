### Code to calculate the ICC value of different diversity index and species abundance in each cohort

### Load packages
library(vegan)
library(phyloseq)
library(ape)
library(bestNormalize)
library(rptR)

### Set wroking directory
setwd("~/Microbiome_response/ICC_Calculation/")

### Import the species abundance data and phylogenetic information of each cohort
species_tree <- read.tree("species.nwk")
species_tax <- tax_table(as.matrix(read.table("tax_ranks.csv",sep= ",",header = TRUE,quote = "",row.names = 1)))

H1_species <- read.csv("H1_species.csv", row.names = 1)
H2_species <- read.csv("H2_species.csv", row.names = 1)
A1_species <- read.csv("A1_species.csv", row.names = 1)
A2_species <- read.csv("A2_species.csv", row.names = 1)
A3_species <- read.csv("A3_species.csv", row.names = 1)
A4_species <- read.csv("A4_species.csv", row.names = 1)
I1_species <- read.csv("I1_species.csv", row.names = 1)
I2_species <- read.csv("I2_species.csv", row.names = 1)
I3_species <- read.csv("I3_species.csv", row.names = 1)
I4_species <- read.csv("I4_species.csv", row.names = 1)
I5_species <- read.csv("I5_species.csv", row.names = 1)

### Diversity calculation
diversity_calculation <- function(input_1, input_2){
  input_1<-input_1[,colSums(input_1 != 0)>= nrow(input_1)*0.1]
  
  ### Alpha diversity calculating
  alpha_diversity<-as.data.frame(diversity(input_1))
  
  colnames(alpha_diversity)<-"species_shannon"
  alpha_diversity$species_simpson<-diversity(input_1,index = "simpson")
  alpha_diversity$observed_OTU_species<-apply(input_1, 1, function(c)sum(c!=0))
  
  ### Beta_diveristy calculating
  bray <- vegdist(input_1, method="bray")
  
  ### PCoA analysis
  pcoaresult<-pcoa(bray)
  vector<-pcoaresult$vectors
  vector<-as.data.frame(vector)
  
  alpha_diversity$PCoA1_species <- vector$Axis.1
  alpha_diversity$PCoA2_species <- vector$Axis.2
  alpha_diversity$PCoA3_species <- vector$Axis.3
  alpha_diversity$PCoA4_species <- vector$Axis.4
  alpha_diversity$PCoA5_species <- vector$Axis.5
  
  OTU = otu_table(t(input_1), taxa_are_rows = TRUE)
  physeq2 = phyloseq(OTU,species_tax,species_tree)
  
  unweighted_unifrac<-UniFrac(physeq2,weighted=FALSE)
  pcoaresult<-pcoa(unweighted_unifrac)
  vector<-pcoaresult$vectors
  vector<-as.data.frame(vector)
  
  alpha_diversity$Unweighted_Unifrac_PCoA1_species<-vector$Axis.1
  alpha_diversity$Unweighted_Unifrac_PCoA2_species<-vector$Axis.2
  alpha_diversity$Unweighted_Unifrac_PCoA3_species<-vector$Axis.3
  alpha_diversity$Unweighted_Unifrac_PCoA4_species<-vector$Axis.4
  alpha_diversity$Unweighted_Unifrac_PCoA5_species<-vector$Axis.5
  
  weighted_unifrac<-UniFrac(physeq2,weighted=TRUE)
  pcoaresult<-pcoa(weighted_unifrac)
  vector<-pcoaresult$vectors
  vector<-as.data.frame(vector)
  
  alpha_diversity$weighted_Unifrac_PCoA1_species<-vector$Axis.1
  alpha_diversity$weighted_Unifrac_PCoA2_species<-vector$Axis.2
  alpha_diversity$weighted_Unifrac_PCoA3_species<-vector$Axis.3
  alpha_diversity$weighted_Unifrac_PCoA4_species<-vector$Axis.4
  alpha_diversity$weighted_Unifrac_PCoA5_species<-vector$Axis.5
  
  return(alpha_diversity)
}

H1_diversity <- diversity_calculation(H1_species, "species")
H2_diversity <- diversity_calculation(H2_species, "species")
A1_diversity <- diversity_calculation(A1_species, "species")
A2_diversity <- diversity_calculation(A2_species, "species")
A3_diversity <- diversity_calculation(A3_species, "species")
A4_diversity <- diversity_calculation(A4_species, "species")
I1_diversity <- diversity_calculation(I1_species, "species")
I2_diversity <- diversity_calculation(I2_species, "species")
I3_diversity <- diversity_calculation(I3_species, "species")
I4_diversity <- diversity_calculation(I4_species, "species")
I5_diversity <- diversity_calculation(I5_species, "species")

### Data transformation function
Diversity_Transform <- function(p) { 
  a<-bestNormalize(p, loo = TRUE)
  a<-a$x.t}
Taxon_Transform <- function(p) { asin(sqrt(0.01*p)) }

### ICC calculation
ICC_Calculation <- function(input_1, input_2){
  study_name <- deparse(substitute(input_1))
  if (input_2 == "diversity") {
    input_1_transformed<-as.data.frame(lapply(input_1, Diversity_Transform))
  } else  {
    input_1 <- input_1[, colSums(input_1 != 0) >= nrow(input_1)*.1]
    input_1_transformed<-as.data.frame(lapply(input_1, Taxon_Transform))
  }
  row.names(input_1_transformed)<-row.names(input_1)
  if (study_name == "H1_diversity" | study_name == "H1_species"){
    study_subject <- sapply(strsplit(row.names(input_1),split = "_"), "[[", 2)
    study_subject_2 <- sapply(strsplit(row.names(input_1),split = "_"), "[[", 3)
    study_subject_2 <- ifelse(study_subject_2=="SF05"|study_subject_2=="SF06","_1","_2")
    study_subject <- paste(study_subject,study_subject_2,sep = "")
  } else (
    study_subject <- sapply(strsplit(row.names(input_1),split = "_"), "[[", 2)
  )
  
  input_1_transformed$Subject <- study_subject
  formula_list <- lapply(colnames(input_1_transformed[1:(length(input_1_transformed)-1)]),function(x) paste(x, " ~ (1|Subject)", sep = ""))
  icc_value <- lapply(formula_list, function(x) rptGaussian(as.formula(x), grname = "Subject", data = input_1_transformed, nboot = 10, npermut = 10)$R)
  #icc_value <- lapply(formula_list, function(x) rptGaussian(as.formula(x), grname = "Subject", data = input_1_transformed, nboot = 10, npermut = 10)$se)
  icc_value <-as.data.frame( unlist(icc_value))
  row.names(icc_value) <- colnames(input_1_transformed[1:(length(input_1_transformed)-1)])
  colnames(icc_value) <- c(study_name)
  return(icc_value)
}

H1_Diversity_ICC <- ICC_Calculation(H1_diversity, "diversity")
H1_Species_ICC <- ICC_Calculation(H1_species, "species")
H2_Diversity_ICC <- ICC_Calculation(H2_diversity, "diversity")
H2_Species_ICC <- ICC_Calculation(H2_species, "species")
A1_Diversity_ICC <- ICC_Calculation(A1_diversity, "diversity")
A1_Species_ICC <- ICC_Calculation(A1_species, "species")
A2_Diversity_ICC <- ICC_Calculation(A2_diversity, "diversity")
A2_Species_ICC <- ICC_Calculation(A2_species, "species")
A3_Diversity_ICC <- ICC_Calculation(A3_diversity, "diversity")
A3_Species_ICC <- ICC_Calculation(A3_species, "species")
A4_Diversity_ICC <- ICC_Calculation(A4_diversity, "diversity")
A4_Species_ICC <- ICC_Calculation(A4_species, "species")
I1_Diversity_ICC <- ICC_Calculation(I1_diversity, "diversity")
I1_Species_ICC <- ICC_Calculation(I1_species, "species")
I2_Diversity_ICC <- ICC_Calculation(I2_diversity, "diversity")
I2_Species_ICC <- ICC_Calculation(I2_species, "species")
I3_Diversity_ICC <- ICC_Calculation(I3_diversity, "diversity")
I3_Species_ICC <- ICC_Calculation(I3_species, "species")
I4_Diversity_ICC <- ICC_Calculation(I4_diversity, "diversity")
I4_Species_ICC <- ICC_Calculation(I4_species, "species")
I5_Diversity_ICC <- ICC_Calculation(I5_diversity, "diversity")
I5_Species_ICC <- ICC_Calculation(I5_species, "species")

### Combine and store output from each cohort 
merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by, all = T)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}

diversity_ICC_all <- merge.all(H1_Diversity_ICC, H2_Diversity_ICC, A1_Diversity_ICC, A2_Diversity_ICC, A3_Diversity_ICC, A4_Diversity_ICC,
                               I1_Diversity_ICC, I2_Diversity_ICC, I3_Diversity_ICC, I4_Diversity_ICC, I5_Diversity_ICC)
write.csv(diversity_ICC_all,"species_diversity_ICC.csv")

species_ICC_all <- merge.all(H1_Species_ICC, H2_Species_ICC, A1_Species_ICC, A2_Species_ICC, A3_Species_ICC, A4_Species_ICC,
                             I1_Species_ICC, I2_Species_ICC, I3_Species_ICC, I4_Species_ICC, I5_Species_ICC)
write.csv(species_ICC_all,"species_ICC.csv")

