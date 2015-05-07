# PLAYRDesign

## Install required R packages

You need to install the devtools package, available from CRAN, and a number of packages from Bioconductor. The rest of the dependencies for PLAYRDesign will be automatically installed

#### Devtools

Open an R session, type the following command and select a CRAN mirror when prompted.

`install.packages("devtools")`

#### Bioconductor packages

Open an R session and type the following commands

```
source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "AnnotationFuncs", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", 
"Biostrings", "GenomicFeatures", "GenomicRanges", "IRanges", "org.Hs.eg.db"))
```

## Primer3

Download the latest version of Primer3 from [here](http://primer3.sourceforge.net/releases.php) and follow the [installation instructions](http://primer3.sourceforge.net/primer3_manual.htm). At the end of the installation the **primer3_core** executable must be in your PATH. There are several ways to accomplish this. On OSX or Linux the easiest way is probably to copy or link the executable in the /usr/local/bin directory. Please note that on OSX GUI applications do not necessarily have the same PATH as the shell, so using other install locations might be cumbersome if you are running the default R GUI. Whichever approach you choose, at the end of the installation you need to be able to run the following command from within R without errors

```
system("primer3_core")
```

## Install BLAST+

The installation of BLAST+ could be the subject of an entire book. Only minimal instructions are given here for the purpose of setting up a barebones BLAST+ environment that will interact with PLAYRDesign. First download BLAST+ [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). You will then need to specify the location where you want your BLAST+ sequence database to be stored. The setup procedure differs according to the platform you are using, please refer to the BLAST+ [manual](http://www.ncbi.nlm.nih.gov/books/NBK1762/) for details. PLAYRDesign makes use of two sequence databases, which need to have these **exact** names:






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

