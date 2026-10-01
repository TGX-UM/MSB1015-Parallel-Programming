###########################################
# -----_- Install required packages -----_#
###########################################
install.packages("data.table")
install.packages("microbenchmark")
install.packages("doParallel")
install.packages("httr")
install.packages("jsonlite")
install.packages("randomForest")
install.packages("caret")
install.packages("BiocManager")
BiocManager::install("EBImage")

# if you got an error when installing EBImage and you are using Linux OS
# then run the following command to install the required packages
# then try installing EBImage again

# apt update && apt-get install -y libfftw3-dev libtiff-dev

######################################################
# Check the number of logical cores (threads) you have
######################################################

library(doParallel)
library(microbenchmark)

detectCores()
detectCores(logical = FALSE)


###########################################
# ------------ Side note ----------------_#
###########################################

# In order to measure the time needed to execute code, we will use the "microbenchmark" package
# which tests the code execution time through multiple runs and reports the average time.
# Documentation: https://cran.r-project.org/web/packages/microbenchmark/microbenchmark.pdf

# However, if you want to have a similar functionality in your R scripts aside from "microbenchmark", you can use the following code snippet.

time_runs <- c()

n_cores <- c(1,2,3,4,5,6,7,8)

set.seed(123)

for(n in n_cores){
  
  print(paste0("Starting run ",n,": ",0," seconds"))
  
  start_time <- Sys.time()
  
  ### Start of code
  
  
  ### End of code
  
  end_time <- Sys.time()
  
  time_taken <- round(end_time - start_time, digits = 2)
  
  print(paste0("Time taken for run ",n,": ",time_taken," seconds"))
  
  time_runs <- c(time_runs, time_taken)
}

###############################################################
# Exercise 1: Read files in parallel, using the "camog" package
###############################################################

# Download the augmented gene expression file from the following link, unzip the file and place the resulting CSV inside the data/ folder
# https://filesender.surf.nl/?s=download&token=f33f51d8-9037-4dd7-beb5-e10d8de445f3 (until 2026-10-22)
# 
# This is a made-up large file simulating gene expression data for 3000 sample and ~43000 human genes
# 
# The file size is ~2.4GB and the large size is considered on purpose to allow measuring reasonable difference in performance when reading it in a parallel way

library(doParallel)
library(data.table)

setwd('/home/rstudio/workspace/parallel-practical')


# The following 4 code snippets read the file "augmented_gene_expression.csv" using 1,2,4,8 threads respectively and print the execution time underneath each cell.
# Run those cells and compare the obtained times.
# 
# What are your observations?


setDTthreads(1)

tm <- microbenchmark(fread("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)

setDTthreads(2)

tm <- microbenchmark(fread("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)


setDTthreads(4)

tm <- microbenchmark(fread("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)


setDTthreads(8)

tm <- microbenchmark(fread("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)

# Lets check now the run time when using the built-in read.csv function to read the file

tm <- microbenchmark(read.csv("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)

# The "data.table" package is an open source project available on GitHub https://github.com/Rdatatable/data.table
# 
# Check the source code file reponsible for reading the csv file in parallel and explain what type of parallel programming is used in this package, which programming language and what library is used.
#   
# Hint: start with this file https://github.com/Rdatatable/data.table/blob/master/src/openmp-utils.c

###############################################################################################
# Exercise 2: Populate a dataframe with API calls (enrich UniProt IDs with protein information)
###############################################################################################

library(doParallel)
library(httr)
library(jsonlite)

# The code cell bellow define a function that will be applied on each dataframe chunk processed in parallel. 
# It takes a dataframe as input (in the parallel way, the dataframe is split into multiple chunks, each to be handled by a different thread) 
# and returns the same dataframe after filling the empty column.

# Explain the function code below using comments showing the role of each piece of code

process_df_row <- function(df){
  
  for(i in 1:nrow(df)) {
    row <- df[i,]
    
    uniprot_id <- row[["uniprot_id"]]
    
    uniprot_url <- paste0("https://rest.uniprot.org/uniprotkb/",uniprot_id,".json")
    
    response <- GET(url = uniprot_url)
    
    response_text <- content(response, "text", encoding = "UTF-8")
    # Parsing data in JSON
    response_json <- fromJSON(response_text)
    
    protein_name <- response_json$proteinDescription$recommendedName$fullName$value
    
    if(is.null(protein_name)){
      protein_name <- ""
    }
    
    protein_length <- response_json$sequence$length
    
    if(is.null(protein_length)){
      protein_length <- ""
    }

    protein_organism <- response_json$organism$scientificName
    
    if(is.null(protein_organism)){
      protein_organism <- ""
    } 
    
    protein_sequence <- response_json$sequence$value
    
    if(is.null(protein_sequence)){
      protein_sequence <- ""
    } 
    
    df[i, 'protein_name'] <- protein_name
    df[i, 'protein_length'] <- protein_length
    df[i, 'protein_organism'] <- protein_organism
    df[i, 'protein_sequence'] <- protein_sequence
    
  }
  
  return(df)
}

# The following code line reads a dataframe of 100 rows and five columns. The first column contains UniProt IDs for 100 proteins and the remaing columns are empty.
# The purpose of this exercise is to apply parallel computing on a dataframe and make external API calls to retrieve information about the protiens in the first column and fill the rest of empty column with relevant information about those proteins.

uniprot_df <- read.csv('data/uniprot_ids_df.csv')

# this line is meant to show how the output look like and not meant to assess performance
uniprot_df_enriched <- process_df_row(uniprot_df)

# The following two statements compare the sequential and parallel processing of the dataframe.
# Run the code and compare the results.
# 
# You can try the parallel part with different number of cores and see how does that effect the execution time

microbenchmark(process_df_row(uniprot_df), times=10L)

microbenchmark( "exp" = {
  
    NUM_OF_CORES <- 4
    
    cl <- makeCluster(NUM_OF_CORES)
    registerDoParallel(cl)
    
    uniprot_df <- read.csv('data/uniprot_ids_df.csv')
    
    df_splits <- split(uniprot_df, (as.numeric(rownames(uniprot_df))-1) %/% 25)
    
    uniprot_df_enriched <- foreach (i = 1:length(df_splits), .combine = rbind, .packages = c("httr", "jsonlite")) %dopar% {
    
        tmp_df = df_splits[[i]]
        
        tmp_df = process_df_row(tmp_df)
    }
    
    stopCluster(cl)
  
  }, times = 10L
)

#########################################################################
# Exercise 3: Apply edge detection to blood cell image files in parallel
#########################################################################

# In this exercise, we will work with applying an image processing function on a large number of blood cell images (simulating what you would do in a similar research project) and we will compare this process with and without parallel computing.
# The dataset was originally obtained from Kaggle (https://www.kaggle.com/datasets/paultimothymooney/blood-cells/). However, the images were copied and multiplied a couple of times to increase their number in order the observe reasonable differences between sequential and parallel approaches.
# Therefore, download the image dataset from the following URL: https://drive.google.com/file/d/1a5EPJPSrrpaKTu6tvIY37sdtqoPdTMNd/view?usp=sharing
# Unzip the folder into the data/ folder and make sure that you have the images directly under data/original_images/

library(doParallel)
library(EBImage)

# Get a list of all image file names in the specified path
list_of_files <- list.files(path="data/original_images", pattern=".jpg", all.files=FALSE, full.names=FALSE)

length(list_of_files)

process_image <- function(image_file){
  
  original_folder <- "data/original_images"
  processed_folder <- "data/processed_images"
  
  img <- readImage(paste0(original_folder,'/',image_file))
  
  #display(img)
  
  gray_img = channel(img, "gray")
  thresh_img = gray_img < 0.7
  binary_img = bwlabel(thresh_img)
  
  #display(binary_img)

  writeImage(binary_img, paste0(processed_folder,'/',image_file))
}

microbenchmark( "exp" = {
  
    for(image_file in list_of_files){
      process_image(image_file)
    }
    
  }, times = 3L
)


microbenchmark( "exp" = {
  
    NUM_OF_CORES <- 4
    
    cl <- makeCluster(NUM_OF_CORES)
    registerDoParallel(cl)
    
    foreach (i = 1:length(list_of_files), .combine = rbind, .packages = c("EBImage")) %dopar% {
      
      image_file = list_of_files[i]
      
      process_image(image_file)
    }
    
    stopCluster(cl)
    
  }, times = 3L
)

# What is the difference in using pool.map() in this excercise compared to excercise 2?

#######################################################################
# Exercise 4: machine learning model training using parallel computing
#######################################################################

library(doParallel)
library(caret)

# make up a regression dataset sample to use it for machine learning model training
dataset <- matrix(rnorm(13000),nrow=500)

time_runs <- c()

n_cores <- c(2,3,4,5,6,7,8)

set.seed(825)

for(n in n_cores){

  print(paste0("Starting run ",n,": ",0," seconds"))
  
  start_time <- Sys.time()
  
  cl <- makeCluster(n)
  registerDoParallel(cl)
  
  model <- train(V26 ~ ., data = as.data.frame(dataset), method = "rf", verbose = FALSE)  
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  
  time_taken <- round(end_time - start_time, digits = 2)
  
  print(paste0("Time taken for run ",n,": ",time_taken," seconds"))
  
  time_runs <- c(time_runs, time_taken)
}

plot(n_cores, time_runs)


# To learn more about the caret package and parallel processing in caret, check this link: https://topepo.github.io/caret/parallel-processing.html

#################################################################################
## Exercise 5: Solve the following problem using a parallel programming approach
#################################################################################

# 1. Load the CSV file in the "data" folder, named "pdb_ids.csv" using read.csv function
# 2. Create a function to process each row of the dataframe. The function should take one argument of type DataFrame and return a dataframe object.
# The function should iterate through the dataframe rows and perform the following steps:
#   * Get the PDB ID from the relevant column and make an HTTP call to download the protein image
#       from PDB (use the following URL template: https://cdn.rcsb.org/images/structures/{}_assembly-1.jpeg
#       In this template the {} is where you put the PDB ID (in lower case), e.g. 1crn
#   * Save the content of the response (binary content) to an image file stored in the folder "data/pdb_images" named with the PDB id and the extension "jpeg"
#   * read the image file from the folder using OpenCV and extract the size of the image (i.e. width and height)
#   * store the width and the height of the image in the relevant columns in the dataframe
# 3. Save the dataframe to a file
# 4. Use the timing template provided at the beginning of this notebook to time your code and test it using 2, 4 and 8 cores


