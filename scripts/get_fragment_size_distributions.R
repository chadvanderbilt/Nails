#!/opt/R/4.1.3/bin/Rscript
#Description. This script pulls data related to sequencing of nails and allows for the comparison of the the two method
#for extracting DNA from nails.  The data is subsequently deidentified and shared as "fragment_quartiles.csv" in data directory. 

## Load required packages
# Data manipulation
library(dplyr)            

# Excel file reading
library(readxl)     




nails <- read_excel("all_cv3_nail_variants_13Dec2021.xlsx") %>%
select(Nail_SampleID, Run, Mnumber_Nail) %>%
distinct(Nail_SampleID, .keep_all = TRUE) %>% mutate(mean_coverage = NA)

nails_frag <- nails %>% select(-mean_coverage) %>%
mutate(frag_mean = NA, frag_median = NA, frag_1Q = NA, frag_3q = NA, frag_sd = NA, frag_peak = NA)

# Define the base directory path
base_dir <- "/dmp/dms/qc-data/IMPACT"



i <- 5
# Iterate through each value of "variable"
for(i in 1:length(nails_frag$Nail_SampleID)) {
  # Create the file pattern with the variable as a wildcard
  variable <- nails_frag$Nail_SampleID[i]
  Run <- nails$Run[i]
  file_pattern_frag <- paste0("ls ",base_dir, "/*/", Run, "/DEFAULT/", paste0(Run, "_ALL_insertsizemetrics.txt"))
  file_pattern <- paste0("ls ",base_dir, "/*/", Run, "/DEFAULT/", paste0(variable, "_bc*.target.covg"))
  print(file_pattern)
  tryCatch({
    matching_files <- system(file_pattern, intern = TRUE)
    matching_files_frag <- system(file_pattern_frag, intern = TRUE)
    barcode <- sub(".*_(bc\\d+)_.*", "\\1", matching_files)
    print(matching_files)
    print(barcode)
    print(matching_files_frag)
    print(i)

    if (length(matching_files) > 0) {
      coverage_file <- read.delim(matching_files[1], header = TRUE, sep = "\t")
      frag_file <- read.delim(matching_files_frag[1], header = TRUE, sep = "\t")
      print(head(frag_file))
      #make new dataframe with columns "frag_file" and barcode 
      new_frag <- frag_file[,c(1, which(colnames(frag_file) == barcode))]
      cumulative <- cumsum(new_frag[,2])
      x=as.numeric(new_frag[1:500,1])
      y=as.numeric(new_frag[1:500,2])
      #mean_value
      nails_frag[i,4] <- weighted_mean <- sum(x * y) / sum(y)
      #median_value 
      nails_frag[i,5]<- as.numeric(new_frag[,1][min(which(cumulative >= 0.5))])
      #first_quartile 
      nails_frag[i,6]<- as.numeric(new_frag[,1][min(which(cumulative >= 0.25))])
      #third_quartile 
      nails_frag[i,7]<- as.numeric(new_frag[,1][min(which(cumulative >= 0.75))])
      #sd_value
      nails_frag[i,8]<- sqrt(sum(y * (x - mean_value)^2))
      #peak 
      nails_frag[i,9]<- (new_frag %>% filter(insert_size=="Peak"))[1,2]
      # nails_frag[i, 4] <- mean(coverage_file$mean_coverage)
    } else {
      cat("No matching files found for variable:", variable, "\n")
    }
  }, error = function(e) {
    cat("Error occurred while executing system command:", conditionMessage(e), "\n")
  })
 
}

nails_frag %>% write.table("frag_nails.csv", sep = ",", row.names = F)




for(i in 1:length(nails$Nail_SampleID)) {
  # Create the file pattern with the variable as a wildcard
  variable <- nails$Nail_SampleID[i]
  Run <- nails$Run[i]
  file_pattern <- paste0("ls ",base_dir, "/*/", Run, "/DEFAULT/", paste0(variable, "_bc*.target.covg"))
  print(file_pattern)
  tryCatch({
    matching_files <- system(file_pattern, intern = TRUE)
    print(matching_files)

    if (length(matching_files) > 0) {
      coverage_file <- read.delim(matching_files[1], header = TRUE, sep = "\t")
      nails[i, 4] <- mean(coverage_file$mean_coverage)
    } else {
      cat("No matching files found for variable:", variable, "\n")
    }
  }, error = function(e) {
    cat("Error occurred while executing system command:", conditionMessage(e), "\n")
  })
 
}

nails %>% write.table("coverage_nails.csv", sep = ",", row.names = F)
