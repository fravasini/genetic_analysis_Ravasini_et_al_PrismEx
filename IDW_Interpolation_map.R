library(admixtools)
library(tidyverse)
library(ape)
library(purrr)
library(pheatmap)
library(viridis)
library(plotly)

# load f2 matrix
f2_mean <- read.table("f2_matrix.txt")

# remove outgroup
mat2 <- f2_mean[rownames(f2_mean) != 'Mbuti',
                colnames(f2_mean) != 'Mbuti']

