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

Download the latest version of Primer3 from [here](http://primer3.sourceforge.net/releases.php), follow the [installation instructions](http://primer3.sourceforge.net/primer3_manual.htm) and copy the **primer3_core** executable in a location of your choice (you will need to specify the location later, see *Configuring PLAYRDesign*).

Primer3 also requires a directory that contains thermodynamic parameters for primer design. When you download the primer3 package these files are located in the *primer3_config* subfolder. Copy the entire folder to a location of your choice.

### Installing BLAST+

The installation of BLAST+ could be the subject of an entire book. Only minimal instructions are given here for the purpose of setting up a barebones BLAST+ environment that will interact with PLAYRDesign. First download BLAST+ [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and copy the *blastn* executable to a location of your choice. You will also have to select a directory where you want your BLAST+ sequence database to be stored (see below, *Configuring PLAYRDesign*). 

PLAYRDesign makes use of two sequence databases, one contains repetitive sequences for the organism of interest, the other is a database of all the transcripts in the organism of interest. These database are used to avoid as much as possible the selection of probes that either match repetive sequences, or transcripts for different genes. (All FASTA sequence database files have to be saved with the *.fa* extension).

Instructions are given here for designing probes for human transcripts. The repetitive sequences can be downloaded from [Repbase](http://www.girinst.org/repbase/). Download the *humrep.ref* and *simple.ref* in FASTA format and concatenate them to generate the *repbase.fa* file. The Human RefSeq RNA sequences can be downloaded by visiting the NCBI ftp server (ftp://ftp.ncbi.nlm.nih.gov/), navigating to the **refseq -> H_sapiens -> H_sapiens -> RNA** folder and selecting the *rna.fa.gz* file. We reccomend filtering this file to only contain *NR* and *NM* records. You can do so by using the *filter_refseq_file* function included in the PLAYRDesign R package (unpack the *rna.fa.gz* file first). For the purpose of this example use *rna_human_high_qual.fa* as the output file name.

```
library(PLAYRDesign)
PLAYRDesign.filter_refseq_file("PUT PATH TO INPUT FILE HERE", "PATH TO OUTPUT FILE HERE")
```

Move both files to the BLAST+ database directory and convert them to BLAST+ databases by typing these commands (the *makeblastdb* program must be in your PATH, or you have to call it by specifying the full path to the executable). Refer to the BLAST+ [manual](http://www.ncbi.nlm.nih.gov/books/NBK279688/) for additional details on how to use *makeblastdb*.

```
makeblastdb -in repbase.fa -parse_seqids -dbtype nucl
makeblastdb -in rna_human_high_qual.fa -parse_seqids -dbtype nucl

```


### Generating EST and exon information (Optional)

The following section is optional. If you don't want to take EST information into account skip this.

PLAYRDesign depends on two pieces of data to determine the exon structure of a gene and its overlap with ESTs. This section explains how to generate this data, see the *Configuring PLAYRDesign* section for details regarding where these files should be stored.

#### Generating EST information

Access the UCSC Table Browser [here](https://genome.ucsc.edu/cgi-bin/hgTables) and select the **intronEst** table from the **Spliced ESTs** track in the **mRNA and EST** group. To shorten the download and processing times, after you press the *get output* button, you can select only the following fields, which are used by PLAYRDesign.

```
strand, tName, tStart, tEnd, blockSizes, tStarts
```
Save the file on your computer and use the following function from the PLAYRDesign R package to convert the txt file into the format that will be used by PLAYRDesign. The output file must have the extension *.RData* (The command will probably take a while to run).

```
PLAYRDesign.convert_est_to_RData("PUT THE PATH TO THE INPUT FILE HERE", "PUT THE PATH TO THE INPUT FILE HERE")
```

#### Generating Exon information

Type the following commands in R (in older version of Bioconductor the **makeTxDbFromUCSC** was called **makeTranscriptDbFromUCSC**, but it works the same). The output file must have the extension *.sqlite*. If you want to design probes for a different organism, change the *genome* parameter (refer to the GenomicFeatures documentation for details).

```
library(GenomicFeatures)
txdb <-  makeTxDbFromUCSC(genome = "hg19", tablename = "refGene")
saveDb(txdb, "PUT THE PATH TO THE OUTPUT FILE HERE")
```

## Usage

If evertyhing was successful you should be able to start PLAYRDesign by typing the following commands

```
library(PLAYRDesign)
PLAYRDesign.run()
```

When you start PLAYRDesign you will be prompted to select a file: select *any* file that is located in the directory you want to be used as the PLAYRDesign **working directory**. (Ideally we would have you select the directory itself instead of a file, but R does not allow to do this in a platform-independent way). Note that the content of the directory is only read at startup so all the necessary files need to be present when you start the GUI, in order for them to be accessible.To stop PLAYRDesign simply hit the "ESC" key in your R session. 

### Configuring PLAYRDesign

PLAYRDesign needs to know the location of external programs and data files to run. These locations are specified in a file that must be named *playrdesign_opt.txt" and must be located in your PLAYRDesign working directory. An example of the format of the file is given below, substitute the relevant paths for your specific installation.

```
BLASTN_EXEC=/usr/bin/blastn           (The full path to the blastn executable)
BLASTN_DB=/opt/BLAST/                 (The directory containing your BLAST database files)
PRIMER3_EXEC=/usr/bin/primer3_core    (The full path to the primer3_core executable)
PRIMER3_CONFIG=/opt/primer3_config/   (The primer3_config directory that is found in the primer3 distribution, see above)
PLAYRDESIGN_DATA=/opt/PLAYRDesign_data   (The directory containing the EST and exon data, if you are using them, see above)
```

### Starting the analysis

First download the sequence of the transcript for which you want to design probes in FASTA format and save it in a plain text file with a .fasta extension.

PLAYRDesign parses the FASTA descritpion line to extract the accession number of the transcript. The accession is used to:
- eliminate BLAST matches to the same transcript
- eliminate BLAST matches to different isoforms of the same transcript
- retrieve the exon structure of the gene

The FASTA line has to use the standard NCBI format which looks similar to this

```
>gi|61676094|ref|NM_006137.6|
```

in practice the best option is to download RefSeq transcripts from the NCBI [nucleotide](http://www.ncbi.nlm.nih.gov/nuccore/) database, preferably choosing *NR* and *NM* records. We recommend choosing the longest isoform of the transcript because the software will show which exons can undergo alternative splicing. 

In the GUI first select the different EST, exon, repetitive sequence and transcript databases and then use the "Select input file" dropdown to select the fasta file you want to design probes for. The boxes with the numeric values are for setting parameteres for the primer3 software. The defaults are the values used in the paper. Once you are ready hit the "Start analysis" button. Several messages should appear in your R window as the software is running. Once the analysis is completed a number of plots will appear in the browser window

### Selecting probes

The candidate probes are displayed as red rectangles at the bottom of the interface. Each pair is identified by a unique number on the rectangle. If you click on a probe both oligos in a pair will be selected. Selected oligos appear in the "Select oligos" box, and can be removed from there if desired. It is also possible to generate a probe pair by selecting individual oligos from two different primer3 pairs. To do so ALT+Click on the first and then ALT+Click on the second (to clear the working selection ALT+Click on any blank region of the plot). When you combine oligos from two different pairs, the "Select oligos" box will display their unique ids (which are different from the ids you see in the plot, and can be visualized by hovering over an oligo) separated by a **"_"** character.

Once you have selected the probe pairs use the "Select PLAYR system" dropdown to select an insert system and enter an id for the first oligo. Our standard is for the 5' oligo (on the transcript) of a pair to be the first one and to have an odd number. Hit the "Write oligos" button and a text file with the **.playrdesign_out.txt** extension will appear in your working directory (the same directory where your fasta files are located). The directory will also contain BLAST and Primer3 output files which are your free to delete after the analysis is completed.
