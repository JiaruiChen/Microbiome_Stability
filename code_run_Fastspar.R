### Code to obtain the SparCC network using Fastspar as implemented in https://rdrr.io/github/danlooo/coabundance/

### Load packages
library(dplyr)
library(tidyr)

### set the working directory
setwd("~/Microbiome_response/network_analysis/")

### Setting Environments
# import the needed functions from https://rdrr.io/github/danlooo/coabundance/src/R/correlate.R
# download the code in a .R file (here as 'correlate.R') and source it as shown
source("/Microbiome_response/network_analysis/correlate.R")

### Create a conda environment and call it fastspar, you can do it very easily as shown here https://anaconda.org/bioconda/fastspar

### Then run the command below to tell the system to use this environment
Sys.setenv(PATH = paste0("/home/user/anaconda3/envs/fastspar/bin", ":", Sys.getenv("PATH")))
Sys.setenv(LD_LIBRARY_PATH = paste0("/home/user/anaconda3/pkgs/mkl-2019.1-144/lib/", ":", Sys.getenv("LD_LIBRARY_PATH")))

### Load the following function to format the data
sparcc_format_species <- function(data) {
  
  data = data[ ,colSums(data!= 0) >= nrow(data)*0.2]
  data = round(data*10^6)
  
  return(data)
}


### Generate the networks 
# 'data' is a dataframe with rows = subjects and columns = species 

# Responders network
data <- read.csv('data_species_R.csv', row.names = 1)
FastSpar_spe_R <- 
  data %>% 
  sparcc_format_species(.) %>% 
  as.matrix(.) %>% 
  correlate_fastspar(., iterations = 20, exclude_iterations = 10, bootstraps = 100, threads = 1) %>%
  filter(p.value < 0.05) %>% 
  select(from,to,estimate) %>%
  rename(FP_Cor = "estimate") %>% 
  filter(abs(FP_Cor) > 0)
write.csv(FastSpar_spe_R, 'FastSpar_spe_R.csv')

# Non responders network
data <- read.csv('data_species_NR.csv', row.names = 1)
FastSpar_spe_NR <- 
  data %>% 
  sparcc_format_species(.) %>% 
  as.matrix(.) %>% 
  correlate_fastspar(., iterations = 20, exclude_iterations = 10, bootstraps = 100, threads = 1) %>%
  filter(p.value < 0.05) %>% 
  select(from,to,estimate) %>%
  rename(FP_Cor = "estimate") %>% 
  filter(abs(FP_Cor) > 0)
write.csv(FastSpar_spe_NR, 'FastSpar_spe_NR.csv')

# The .csv files can be imported to Cytoscape to analyze the networks with networkAnalyzer (https://apps.cytoscape.org/apps/networkanalyzer)

