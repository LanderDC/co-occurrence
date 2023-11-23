#!/usr/bin/env Rscript
if( ! require( "argparse", character.only = TRUE, quietly = T) ){
  #  If package was not able to be loaded then re-install
  install.packages( "argparse" , dependencies = TRUE, repo="http://cran.rstudio.com/" )
  #  Load package after installing
  library( "argparse" , character.only = TRUE )
} else {
  library( "argparse" , character.only = TRUE )
}

#' Specify CLI options
#option_list = list(
#  make_option(c("-i", "--input"), type="character", default=NULL, 
#              help="Abundance table.", metavar="file"),
#  make_option(c("-p", "--prevalence"), type="double", default=0.1, 
#              help="Minimum percentage of samples that the contig should be present in to be considered for correlation analysis. [default = %default]", 
#              metavar="float"),
#  make_option(c("-c", "--correlation"), type="double", default=0.3, 
#              help="Minimum correlation to keep pairs in the pairwise correlation dataframe. [default = %default]", 
#              metavar="float"),
#  make_option(c("-t", "--threads"), type="integer", default=4, 
#              help="Number of threads to use. [default = %default]", 
#              metavar="integer")
#  #make_option(c("-o", "--out"), type="character", default="rotapies.pdf", 
#  #            help="output file name [default= %default]", metavar="output")
#)

#opt_parser = OptionParser(option_list=option_list)
#opt = parse_args(opt_parser)

#if (opt$threads <= 0) {
#  stop("You can not supply 0 or negative threads.")
#}


parser <- ArgumentParser()

# Add arguments
parser$add_argument("-i", "--input", type="character", default=NULL,
                    help="Abundance table.", metavar="file", required=T)
parser$add_argument("-p", "--prevalence", type="double", default=0.1,
                    help="Minimum percentage of samples that the contig should be present in to be considered for correlation analysis. [default = %(default)s]",
                    metavar="float")
parser$add_argument("-c", "--correlation", type="double", default=0.3,
                    help="Minimum correlation to keep pairs in the pairwise correlation dataframe. [default = %(default)s]",
                    metavar="float")
parser$add_argument("-t", "--threads", type="integer", default=4,
                    help="Number of threads to use. [default = %(default)s]",
                    metavar="integer")

# Uncomment the following lines if you want to include an output option
# parser$add_argument("-o", "--out", type="character", default="rotapies.pdf",
#                     help="output file name [default= %(default)s]", metavar="output")

# Parse the command-line arguments
opt <- parser$parse_args()


if (opt$threads <= 0) {
  stop("You can not supply 0 or negative threads.")
}

#' Check if required packages are installed and install/load them
#' function that evaluates installed packages
packages <- c("readr", "dplyr", "purrr", "foreach", "doSNOW")

for( pkg in packages ){
  #  require returns TRUE invisibly if it was able to load package
  if( ! require( pkg , character.only = TRUE, quietly = T) ){
    #  If package was not able to be loaded then re-install
    install.packages( pkg , dependencies = TRUE, repo="http://cran.rstudio.com/" )
    #  Load package after installing
    library( pkg , character.only = TRUE )
  } else {
    library( pkg , character.only = TRUE )
  }
}

message("Read in abundance table.")
OTU <- read.table(opt$input, header=TRUE, sep="\t", dec=".")

df <- OTU %>%
  select(-contains("NC")) %>% 
  mutate(sample_count = rowSums(. != 0)-1,
         proportion_samples = sample_count / (ncol(.) - 1))  # Subtract 2 for the ID and sample count columns

# Define the threshold
threshold <- opt$prevalence

# Filter rows where the proportion of 0s is less than or equal to the threshold
filtered_df <- df %>%
  filter(proportion_samples >= threshold) %>% 
  select(-sample_count, -proportion_samples)

integer_columns <- sapply(filtered_df, is.integer)

prevalence_df <- filtered_df %>%
  mutate_if(integer_columns, ~ifelse(. > 0, 1, .))

n <- length(prevalence_df$Contig)

message(paste("Make empty matrix of", n, "dimensions (contig prevalence in samples =", threshold*100, "%)."))
Cor <-  data.frame(matrix(NA, ncol = n+1, nrow = n))

colnames(Cor) <- c("Contig", as.vector(prevalence_df[,1]))
Cor$Contig <- prevalence_df[,1]

message("Start correlation analysis, this might take some time...")

## Set up a parallel backend with a specified number of cores
num_cores <- opt$threads

cl <- makeCluster(num_cores)
registerDoSNOW(cl)

# Define a function for calculating correlations
correlation_function <- function(i, j, prevalence_df) {
  cor(as.numeric(prevalence_df[i, -1]), as.numeric(prevalence_df[j, -1]))
}

# Add progressbar
iterations <- length(Cor$Contig)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# Parallelize the nested loop using foreach
results <- foreach(i = 1:length(Cor$Contig), .combine = rbind, 
.packages = 'foreach', .inorder=T, .options.snow=opts) %dopar% {
  foreach(j = 1:length(Cor$Contig), .combine = c) %do% {
    if (i >= j) {
      correlation_function(i, j, prevalence_df)
    } else {
      NA
    }
  }
}

close(pb)

# Stop the parallel backend
stopCluster(cl)

# Update Cor with the correlation values
Cor[, 2:(length(prevalence_df[,1]) + 1)] <- results

# Original from Lore
#for (i in 1:length(Cor$Contig)){
#  for (j in 1:length(Cor$Contig)){
#    if (i < j){ 
#      next
#    }  
#    Cor[j,i+1] <- cor(as.numeric(prevalence_df[i,-1]), as.numeric(prevalence_df[j,-1]))
#  }
#}

# Not parallelized
# Reworked
#for (i in 1:length(Cor$Contig)){
#  if (i %% 500 == 0){
#   message(paste("Starting calculation", i)) 
#  }
#  for (j in 1:length(Cor$Contig)){
#    if (i <= j){ 
#      Cor[j,i+1] <- cor(as.numeric(prevalence_df[i,-1]), as.numeric(prevalence_df[j,-1]))
#    }  else {
#      NA
#    }
#  }
#}

message("Write correlation matrix.")
write_tsv(Cor, "./correlation_matrix.tsv", col_names = T)

rownames(Cor) <- Cor$Contig
Cor <- Cor[-1]

threshold <- opt$correlation
correlation_matrix <- as.matrix(Cor)

# Find the row and column indices of entries with correlation above the threshold
message("Create pairwise dataframe.")
indices <- which(correlation_matrix >= threshold, arr.ind = TRUE)

related_contigs_df <- map_df(1:nrow(indices), function(i) {
  row_idx <- indices[i, 1]
  col_idx <- indices[i, 2]
  
  data.frame(
    Contig1 = rownames(correlation_matrix)[row_idx],
    Contig2 = colnames(correlation_matrix)[col_idx],
    Correlation = correlation_matrix[row_idx, col_idx]
  )
})

related_contigs_df_filtered <- related_contigs_df %>% 
  filter(Contig1!=Contig2) %>% 
  arrange("Contig2")

message("Write pairwise dataframe.")
write_tsv(related_contigs_df_filtered, "./related_contigs.tsv", col_names = T)

message("Finished.")
