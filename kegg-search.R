library("MetaboAnalystR")
library("KEGGREST")
#library("org.Hs.eg.db")
#library("limma")
library("tidyverse")
library("httr")
library("xml2")
library("RColorBrewer")
source("/Users/don/Documents/bin/sprida-src.R")

setwd("/Users/don/Documents/pathway-mapping/maca-to-ipath/data")

fn_in <- ("drossie_minus_pbqc.csv")
kegg_species_id <- "dme"



# Get KEGG name-to-map.id tibble
kl.dt <- as.data.table(keggList("pathway", kegg_species_id), keep.rownames = T)
kl.dt[, 2] <- as.vector(unlist(kl.dt[ , lapply(V2, function(x) strsplit(x, " - ")[[1]][1])]))
kl.dt[, 1] <- as.vector(unlist(kl.dt[ , lapply(V1, function(x) gsub("path:", "", x))]))
colnames(kl.dt) <- c("pw_id", "pw_name")

# ========== ma.ca PW analysis - compound list only # ========== 
#mSet<-InitDataObjects("conc", "pathora", FALSE)
#cmpd.vec<-c("Acetoacetic acid","Beta-Alanine","Creatine","Dimethylglycine","Fumaric acid","Glycine","Homocysteine","L-Cysteine","L-Isolucine","L-Phenylalanine","L-Serine","L-Threonine","L-Tyrosine","L-Valine","Phenylpyruvic acid","Propionic acid","Pyruvic acid","Sarcosine")
#cmpd.vec <- scan("/Users/don/OneDrive\ -\ The\ University\ of\ Melbourne/kenji/data/cpd_list.txt", character(), quote = "")
#mSet<-Setup.MapData(mSet, cmpd.vec);
#mSet<-CrossReferencing(mSet, "name");
#mSet<-CreateMappingResultTable(mSet)
#mSet<-SetKEGG.PathLib(mSet, "hsa")
#mSet<-SetMetabolomeFilter(mSet, F);
#mSet<-CalculateOraScore(mSet, "rbc", "hyperg")


# # ==========  ma.ca PW analysis - from run summary with 2 groups# ========== 
mSet<-InitDataObjects("conc", "pathqea", FALSE)
mSet<-Read.TextData(mSet, fn_in, "rowu", "disc");
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-SetKEGG.PathLib(mSet, kegg_species_id)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateQeaScore(mSet, "rbc", "gt")


# filter out FDR < 0.01
my_fdr <- 0.01
tbl.out <- as_tibble(mSet$analSet$qea.mat, rownames="map")
tbl.out2 <- tbl.out %>% filter(FDR < my_fdr)
#all.pws <- names(mSet$analSet$ora.hits)

# Get the significant pathways
# But not all pw names can be found in KEGG
available.pws.in.kegg <- c()
for (pw in significant.pws) {
  d_t <- kl.dt %>% filter(pw_name==pw)
  if (nrow(d_t)==0) {
    print(paste0(pw, " not found in KEGG!"))
  } else {
    available.pws.in.kegg <- c(available.pws.in.kegg, pw)
  }
}

significant.pws <- as.vector(unlist(tbl.out2[, "map"]))
significant.pw.ids <- c()
significant.metabs <- c()
for (pw in significant.pws) {
  pw_id <- kl.dt[pw_name==pw, pw_id]
  significant.pw.ids <- c(significant.pw.ids, pw_id)
  x <- as.vector(unlist(mSet$analSet$qea.hits[pw_id]))
  significant.metabs <- c(significant.metabs, x)
}
# drop duplicates
significant.metabs <- unique(significant.metabs)

print(paste0(length(significant.pws), " significant PWs found at FDR < 0.01"))
print(paste0(length(significant.metabs), " significant unique metabolites found"))

selection.str <- ""
for (x in significant.metabs) {
  print(paste0(x, " W20 #ff0000"))
  selection.str <- paste0(selection.str, paste0(x, " W20 #ff0000"), "\n")
}

#palette <- brewer.pal(length(significant.pws), name="Set3")
url <- "https://pathways.embl.de/mapping.cgi"
# First POST call to iPath where `whole_modules`=0
call1 <- list(selection = selection.str, 
              export_type="svg", 
              default_opacity="0.7",
              default_width="2",
              default_radius="5",
              whole_modules="0")
r <- POST(url, body = call1, encode = "form", verbose())
print(http_status(r))
xml_raw1 <- content(r, "text")
write(xml_raw1, "ipath_nodes_only.svg")


# POST request to ipath where `whole_modules`=1
body <- list(selection = selection.str, 
             export_type="svg", 
             default_opacity="0.7",
             default_width="2",
             default_radius="5",
             whole_modules="1")
r <- POST(url, body = body, encode = "form", verbose())
print(http_status(r))
write(content(r, "text"), "ipath_whole_module.svg")







# Test GET call to KEGG
url <- "http://rest.kegg.jp/get/hsa01100/kgml"
r <- GET(url, verbose())

kxml <- content(r)

write(content(r, "text"), "hsa01100.xml")








page<-read_xml('<data>
    <result name="r">
               <item>
               <str name="id">123</str>
               <str name="xxx">aaa</str>
               </item>
               <item>
               <str name="id">456</str>
               <str name="xxx">aaa</str>
               </item>
               </result>
               </data>')


nodes<-xml_find_all(page, ".//str")
#filter out the nodes where the attribute name=id
nodes<-nodes[xml_attr(nodes, "name")=="id"]
#get values (as character strings)
xml_text(nodes)


