# PLAYRDesign

## Installation

### Install required R packages

You need to install the devtools package, available from CRAN, and a number of packages from Bioconductor. The rest of the dependencies for PLAYRDesign will be automatically installed

### Devtools

Open an R session, type the following command and select a CRAN mirror when prompted.

`install.packages("devtools")`

### Bioconductor packages

Open an R session and type the following commands

```
source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "AnnotationFuncs", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", 
"Biostrings", "GenomicFeatures", "GenomicRanges", "IRanges", "org.Hs.eg.db"))
```
### Install PLAYRDesign

Start an R session and type the following commands

```
library(devtools)
install_github("nolanlab/PLAYRDesign")
```

This will install the PLAYRDesign R package together with all the required dependencies. 

However before you can use PLAYRDesign you will need to install two additional programs: **Primer3** and **BLAST+**



### Installing Primer3

Download the latest version of Primer3 from [here](http://primer3.sourceforge.net/releases.php) and follow the [installation instructions](http://primer3.sourceforge.net/primer3_manual.htm). At the end of the installation the **primer3_core** executable must be in your PATH. There are several ways to accomplish this. On OSX or Linux the easiest way is probably to copy or link the executable in the /usr/local/bin directory. Please note that on OSX GUI applications do not necessarily have the same PATH as the shell, so using other install locations might be cumbersome if you are running the default R GUI. Whichever approach you choose, at the end of the installation you need to be able to run the following command from within R without errors

```
system("primer3_core")
```

Primer3 also requires a directory that contains thermodynamic parameters for primer design. When you download the primer3 package these files are located in the *primer3_config* subfolder. This directory has to reside in one of the following locations:

- the */opt* folder (OSX and Linux)
- the working directory from which the **primer3_core** executable is run (OSX, Linux, Windows)

If you choose the first option, simply copy the entire *primer3_config* folder in */opt*. If you go with the second option, the *primer3_config* directory needs to be located in your PLAYRDesign **working directory** (see Below *Usage*).

### Installing BLAST+

The installation of BLAST+ could be the subject of an entire book. Only minimal instructions are given here for the purpose of setting up a barebones BLAST+ environment that will interact with PLAYRDesign. First download BLAST+ [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and make sure that the *blastn* and *makeblastdb* executables are in your PATH (see the above considerations for primer3 regarding the best way to do this). You will then need to specify the location where you want your BLAST+ sequence database to be stored. The setup procedure differs according to the platform you are using, please refer to the BLAST+ [manual](http://www.ncbi.nlm.nih.gov/books/NBK1762/) for details. PLAYRDesign makes use of two sequence databases, which need to have these **exact** names:

```
repbase.fa: this is a database of repetitive sequences
rna_human_high_qual.fa: this is a database of reference human RNA sequences
```
The repetitive sequences can be downloaded from [Repbase](http://www.girinst.org/repbase/). Download the *humrep.ref* and *simple.ref* in FASTA format and concatenate them to generate the *repbase.fa* file.

The Human RefSeq RNA sequences can be downloaded by visiting the NCBI ftp server (ftp://ftp.ncbi.nlm.nih.gov/), navigating to the **refseq -> H_sapiens -> H_sapiens -> RNA** folder and selecting the *rna.fa.gz* file. We reccomend filtering this file to only contain *NR* and *NM* records. You can do so by using the *filter_refseq_file* function included in the PLAYRDesign R package (unpack the *rna.fa.gz* file first).

```
library(PLAYRDesign)
PLAYRDesign.filter_refseq_file("PUT PATH TO INPUT FILE HERE", "PATH TO OUTPUT FILE HERE")
```

At this point you should have two sequence files in FASTA format, called *repbase.fa* and *rna_human_high_qual.fa*. Move both files to the BLAST+ database directory and convert them to BLAST+ databases by typing these commands (the *makeblastdb* program must be in your PATH, see above). Refer to the BLAST+ [manual](http://www.ncbi.nlm.nih.gov/books/NBK279688/) for additional details on how to use *makeblastdb*.

```
makeblastdb -in repbase.fa -parse_seqids -dbtype nucl
makeblastdb -in rna_human_high_qual.fa -parse_seqids -dbtype nucl

```

If everything was successful you should be able to type the following commands in an R session

```
system("blastn -db repbase.fa")
system("blastn -db rna_human_high_qual.fa")
```

and get the following error message, which indicates that R can find your *blastn* executable and your database file, even though not input sequence is provided

```
Command line argument error: Query is Empty!
```


If evertyhing was successful you should be able to start PLAYRDesign by typing the following commands

```
library(PLAYRDesign)
PLAYRDesign.run()
```
to stop PLAYRDesign simply hit the "ESC" key in your R session.

### Regenerating EST and exon information

PLAYRDesign depends on two pieces of data to determine the exon structure of a gene and its overlap with ESTs. The data files are provided in the **inst/** subdirectory of the package and are named **spliced_est_hg19.RData** and **UCSC_Refseq_transcripts.sqlite**. When you install the R package these files are automatically copied in the installation directory. This section explains how to regenerate these data files from scratch, in case you wanted to upgrade EST or exon definitions. After you have generated new versions of these files you will have to move them in their installation directory. To determine the location of the installation directory on your machine, type the following command in an R session

```
system.file("spliced_est_hg19.RData", package = "PLAYRDesign")
```

#### Regenerating EST information

Access the UCSC Table Browser [here](https://genome.ucsc.edu/cgi-bin/hgTables) and select the **intronEst** table from the **Spliced ESTs** track in the **mRNA and EST** group. To shorten the download and processing times, after you press the *get output* button, you can select only the following fields, which are used by PLAYRDesign.

```
strand, tName, tStart, tEnd, blockSizes, tStarts
```
Save the file on your computer and use the following function from the PLAYRDesign R package to convert the txt file into the format that will be used by PLAYRDesign

```
PLAYRDesign.convert_est_to_RData("PUT THE PATH TO THE INPUT FILE HERE")
```

The command will take a while to run and will  generate a file called **spliced_est_hg19.RData** which you can then move to the PLAYRDesign installation folder (see above).

#### Regenerating Exon information

Type the following commands in R (in older version of Bioconductor the **makeTxDbFromUCSC** was called **makeTranscriptDbFromUCSC**, but it works the same)

```
library(GenomicFeatures)
txdb <-  makeTxDbFromUCSC(genome = "hg19", tablename = "refGene")
saveDb(txdb, "UCSC_Refseq_transcripts.sqlite")
```

This will create a file named **UCSC_Refseq_transcripts.sqlite** in your working directory, which you can then move in the PLAYRDesign installation directory.

## Usage

### Starting the analysis

First download the sequence of the transcript for which you want to design probes in FASTA format and save it in a plain text file with a .fasta extension. Because of the various databases that PLAYRDesign uses during the analysis, at present the software only works with human sequences. Support for additional species may be added in the future.

PLAYRDesign parses the FASTA descritpion line to extract the accession number of the transcript. The accession is used to:
- eliminate BLAST matches to the same transcript
- eliminate BLAST matches to different isoforms of the same transcript
- retrieve the exon structure of the gene

The FASTA line has to use the standard NCBI format which looks similar to this

```
>gi|61676094|ref|NM_006137.6|
```

in practice the best option is to download RefSeq transcripts from the NCBI [nucleotide](http://www.ncbi.nlm.nih.gov/nuccore/) database, preferably choosing *NR* and *NM* records. We recommend choosing the longest isoform of the transcript because the software will show which exons can undergo alternative splicing. 

When you start the PLAYRDesign software you will be prompted to select a file: you can choose *any* .fasta file that is located in the directory which contains your transcript sequences. This directory will be your **working directory**.

Your R window will then show the message 

```
"Loading EST data..."
```

this will take a couple of minutes. When the process is completed the PLAYRDesign controls will appear in your browser window. In the GUI use the "Select input file" dropdown to select the fasta file you want to design probes for. The boxes with the numeric values are for setting parameteres for the primer3 software. The defaults are the values used in the paper.
Once you are ready hit the "Start analysis" button. Several messages should appear in your R window as the software is running. Once the analysis is completed a number of plots will appear in the browser window

### Selecting probes

The candidate probes are displayed as red rectangles at the bottom of the interface. Each pair is identified by a unique number on the rectangle. If you click on a probe both oligos in a pair will be selected. Selected oligos appear in the "Select oligos" box, and can be removed from there if desired. It is also possible to generate a probe pair by selecting individual oligos from two different primer3 pairs. To do so ALT+Click on the first and then ALT+Click on the second (to clear the working selection ALT+Click on any blank region of the plot). When you combine oligos from two different pairs, the "Select oligos" box will display their unique ids (which are different from the ids you see in the plot, and can be visualized by hovering over an oligo) separated by a **"_"** character.

Once you have selected the probe pairs use the "Select PLAYR system" dropdown to select an insert system and enter an id for the first oligo. Our standard is for the 5' oligo (on the transcript) of a pair to be the first one and to have an odd number. Hit the "Write oligos" button and a text file with the **.playrdesign_out.txt** extension will appear in your working directory (the same directory where your fasta files are located). The directory will also contain two **.blast_out.txt** files which contain the BLAST+ result and which are your free to delete after the analysis is completed.






