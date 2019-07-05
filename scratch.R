library("httr")
library("XML")
library("methods")
library("rvest")
#library("dplyr")
library("data.table")



# Get (cx, cy) coords of significant ellipses, i.e. where fill=<red>
x0 <- read_html("ipath_nodes_only.svg")
ellipse_ls <- xml_find_all(x0, '//ellipse')
x0.attr.dt <- as.data.table(t(as.data.table(xml_attrs(ellipse_ls))))
colnames(x0.attr.dt) <- names(xml_attrs(ellipse_ls)[[1]])
x0.attr.dt <- x0.attr.dt[fill=="#ff0000"]
coords_ls <- c()
for (i in 1:nrow(x0.attr.dt)) {
  coords_ls <- rbind(coords_ls, c(x0.attr.dt[i, cx], x0.attr.dt[i, cy]))
}

# Modify whole_module.svg
x1 <- scan("ipath_whole_module.svg", what="", sep="\n")

for (j in 1:length(x1)) {
  # Change colour of significant ellipses to blue
  if (grepl("ellipse", x1[j])) {
    attrs_ls <- xml_attrs(read_xml(x1[j]))

    for (i in 1:nrow(coords_ls)) {
      if (attrs_ls["cx"] == coords_ls[i,1] & attrs_ls["cy"] == coords_ls[i,2]) {
        x1[j] <- gsub("#ff0000", "#00A1E1", x1[j])
      }
    }
  }
}

# Set highlighted path width to 7, and red ellipse radius to 5
highlight_path_width <- 7
module_ellipse_radius <- 5
for (j in 1:length(x1)) {
  if (grepl("path", x1[j]) & grepl("stroke-width='20'", x1[j])) {
    x1[j] <- gsub("stroke-width='20'", paste0("stroke-width='", highlight_path_width, "'"), x1[j])
  }
  if (grepl("ellipse", x1[j]) & grepl("fill='#ff0000'", x1[j])) {
    x1[j] <- gsub("rx='20'", paste0("rx='", module_ellipse_radius, "'"), x1[j])
    x1[j] <- gsub("ry='20'", paste0("ry='", module_ellipse_radius, "'"), x1[j])
  }
}

# Write out 
write(x1, file = "test2.html")
