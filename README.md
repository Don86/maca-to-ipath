# MA.ca to iPath

Call the pathway (enrichment) analysis module in `MetaboAnalystR`, then POST to iPath for network viz.

Main workbook: `workbook.R`. `myfunc.R` is a library of functions which are sourced within `workbook.R`. There are POST requests involved, so an active internet connection is required. 

### Libraries

```
library("MetaboAnalystR")
library("KEGGREST")
library("tidyverse")
library("httr")
library("RColorBrewer")
library("XML")
library("methods")
library("rvest")
library("data.table")
```