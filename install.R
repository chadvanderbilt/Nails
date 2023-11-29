#!/opt/R/4.1.3/bin/Rscript


#This repository was built in below environment
#               _                           
# platform       x86_64-pc-linux-gnu         
# arch           x86_64                      
# os             linux-gnu                   
# system         x86_64, linux-gnu           
# status                                     
# major          4                           
# minor          1.3                         
# year           2022                        
# month          03                          
# day            10                          
# svn rev        81868                       
# language       R                           
# version.string R version 4.1.3 (2022-03-10)
# nickname       One Push-Up    





# Below is a list of required packages:
required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "varhandle",
  "ComplexHeatmap",
  "circlize",
  "ggforce",
  "readxl",
  "ggpubr",
  "beyonce",
  "gmodels",
  "tibble",
  "forcats",
  "grid",
  "gtable",
  "gridExtra"
)

# Install packages if not already installed
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# Install each package
lapply(required_packages, install_if_missing)
