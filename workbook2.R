library("MetaboAnalystR")
library("KEGGREST")
library("tidyverse")
library("httr")
library("RColorBrewer")
library("XML")
library("methods")
library("rvest")
library("data.table")


source("/Users/don/Documents/pathway-mapping/maca-to-ipath/myfunc.R")

#  ==================== INPUT FIELDS ==================== 
setwd("/Users/don/Documents/pathway-mapping/maca-to-ipath/")

fn_in <- "sample-input.csv"
kegg_species_id <- "dme"

tbl0 <- read_csv(fn_in)
