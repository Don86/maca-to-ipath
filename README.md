# A Collection of Metabolomic Network Analysis Scripts

Hopefully, all this will one day coalesce in some organized fashion into a more robust piece of software, like a Shiny app. 

* `workbook.R` calls the pathway (enrichment) analysis module in `MetaboAnalystR`, then makes POST requests to iPath to call network `svg` files. `myfunc.R` is a library of functions which are sourced within `workbook.R`. The methodology behind pathway enrichment is the `globaltest` algorithm, which determines whether some pre-specified *group of metabolites* is related to some outcome. It does not matter to the test whether the group of metabolites consists of increased or decreased metabolite intensities (or both), all that matters is the aggregate behaviour of the group. 

* `workbook2.R` (under construction) comprises the traditional method of pathway analysis, i.e.:
1. For each metabolite, conduct a t-test between pairs of experimental groups to determine statistical significance of change. Each metabolite thus has a p-value associated with it. Adjust this p-value to account for multiple corrections. 
2. For each metabolite, compute all pairwise fold changes between the means of each experimental group. 
3. From (1), filter out non-significant metabolite changes. 
4. From (3), threshold to determine which fold changes can be considered increased/decrease/no change. 
5. Plot the results of (4) on iPath `svg` files. 

### Notes
* iPath doesn't explicitly allow species selection (whereas `kegg` does). I'm using only `map011000`, but I'm not sure if this global metabolomic map differs between species (I'm pretty sure it does). 

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