library("MetaboAnalystR")

call_maca_pw_analysis <- function(fn_auc_csv, kegg_species_id) {
	"Calls the pathway enrichment analysis module from MetaboAnalystR. 

	PARAMS
	------
	fn_auc_csv: filename of input run summary table as a csv file, with AUCs as values. 
	rownames are the sample names, column names are the metabolite names. Column 1 are the 
	experimental groupings. Because of the way this module works, only 2 groups are supported.


	RETURNS
	-------
	tbl.out: output tibble of pathway analysis enrichment. columns:
	map (compound common name), total cmpd, Hits, raw p (raw p value), -log p, 
	Holm adjust(ed p value), FDR, Impact. 
	"

	mSet<-InitDataObjects("conc", "pathqea", FALSE)
	mSet<-Read.TextData(mSet, fn_auc_csv, "rowu", "disc");
	mSet<-CrossReferencing(mSet, "name");
	mSet<-CreateMappingResultTable(mSet)
	mSet<-SanityCheckData(mSet)
	mSet<-ReplaceMin(mSet);
	mSet<-PreparePrenormData(mSet)
	mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
	#mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
	#mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
	mSet<-SetKEGG.PathLib(mSet, kegg_species_id)
	mSet<-SetMetabolomeFilter(mSet, F);
	mSet<-CalculateQeaScore(mSet, "rbc", "gt")

	tbl.out <- as_tibble(mSet$analSet$qea.mat, rownames="map")
	pw.dict <- mSet$analSet$qea.hits

	return(list(tbl.out, pw.dict))
}


get_significant_ellipse_coords <- function(url, selection.str) {
	"Makes a POST req to iPath with `selection.str` (which should comprise all
	metabolite hits), and `whole_modules` set to 0. Scrapes the resulting output
	to retrieve the `cx` and `cy` attributes of each `ellipse` html object. 
	This function got modularized for containment. 

	PARAMS
	------
	url: string; API url.
	selection.str: vector of characters, comprising the KEGG Ids of metabs of interest

	RETURNS
	-------
	coords_ls: mtrix of coordinates, shape (n_significant_ellipses, 2)
	"
	call1 <- list(selection = selection.str, 
              export_type="svg", 
              default_opacity="0.7",
              default_width="1",
              default_radius="5",
              whole_modules="0")
	r <- POST(url, body = call1, encode = "form")
	print(http_status(r)$message)
	# Get (cx, cy) coords of significant ellipses, i.e. where fill=<red>
	ellipse_ls <- xml_find_all(read_html(content(r, "text")), '//ellipse')
	x0.attr.dt <- as.data.table(t(as.data.table(xml_attrs(ellipse_ls))))
	colnames(x0.attr.dt) <- names(xml_attrs(ellipse_ls)[[1]])
	x0.attr.dt <- x0.attr.dt[fill=="#ff0000"]
	coords_ls <- c()
	for (i in 1:nrow(x0.attr.dt)) {
	  coords_ls <- rbind(coords_ls, c(x0.attr.dt[i, cx], x0.attr.dt[i, cy]))
	}
	return(coords_ls)
}