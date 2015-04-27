# PLAYRDesign

## Install required R packages

You need to install the devtools package, available from CRAN, and a number of packages from Bioconductor. The rest of the dependencies for PLAYRDesign will be automatically installed

#### Devtools

Open an R session, type the following command and select a CRAN mirror when prompted.

`install.packages("devtools")`

#### FlowCore

Open an R session and type the following commands

```
source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "AnnotationFuncs", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", 
"Biostrings", "GenomicFeatures", "GenomicRanges", "IRanges")

```

## Install PLAYRDesign

Once you have succesfully completed the steps above, you have to create a Github token by following [these instructions.](https://help.github.com/articles/creating-an-access-token-for-command-line-use/) (This won't be necessary anymore when the repository goes public).
Copy the token, start an R session and type the following commands, substituing your Github token

```
library(devtools)
install_github("nolanlab/PLAYRDesign", auth_token = "YOUR TOKEN HERE")
```

This will install the PLAYRDesign R package together with all the required dependencies. If evertyhing was successful you should be able to start PLAYRDesign by typing the following commands

```
library(PLAYRDesign)
PLAYRDesign.run()
```
to stop PLAYRDesign simply hit the "ESC" key in your R session.

