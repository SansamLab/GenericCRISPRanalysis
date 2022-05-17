Vignette3: Run MAGeCK MLE
================
Kevin Boyd
2022-05-17

## ![Alt text](vignette3Data/Magic.jpg)

[MAGeCK Website](https://sourceforge.net/p/mageck/wiki/Home/)

## Step 1: Check Count Tables are Complete with All Samples

### Counts Table Example

Counts tables should have the header “sgRNA”, “Gene”, and a column for
each sample. Make sure that you have included your control in this
counts table (possibly Day0 or a pDNA sample). Also sometime we get
“extra” reads a few weeks after the sequencing run. This also shows how
to combine those additional reads into one counts table.

[CountsTable](vignette3Data/All_Day9_CDC45_Counts_pDNA.txt)

``` r
# read in data
# counts_1 <- read.table(file="~/Dropbox (OMRF)/Github/CRISPRScreenAnalysis/Vignettes/vignette3Data/All_Day9_CDC45_Counts_pDNA.txt", header = T)
counts_2 <- read.table(file="~/Dropbox (OMRF)/Github/CRISPRScreenAnalysis/Vignettes/vignette3Data/sample1.count.txt", header = T)

# order the guides in column 2 in the same order as column 1
# counts_2 <- counts_2[match(counts_1$sgRNA,counts_2$sgRNA),]

# check to see if columns are equal
# all.equal(counts_1$sgRNA,counts_2$sgRNA)

# remove first 2 columns from both dataframes and add the similar values together
# new.counts <- counts_1[c(8:11)] + counts_2[c(4:7)]

# combine counts tables together. Use when not adding additional counts.
# new.counts <- cbind(counts_1,counts_2[,-c(1,2)])
# head(new.counts)

# add back the first two columns (with sgRNA and gene name). Use when adding additional counts
# new.counts <- cbind(counts_1[c(1:7)],new.counts)

# make a text file for the combined counts tables
# write.table(new.counts, file = "AllCounts_Day9_pDNA_AdditionalReads.txt", quote=F,row.names = F,sep="\t")
```

## Step 2: Gather Count Tables and Design Matrix

### Design Matrix Example

The Design Matrix should have a header row with “sample”, “baseline”,
and +/- of what you are interested in measuring “Auxin” and a replicate
identifier “May5”. The rows of the design matrix are named the same as
the columns of the Count Table (ie: the rows are named for each sample)

[Design matrix](vignette3Data/DesignMatrixAll.txt)

| Samples                | baseline | DMSO | Gilteritinib | Midosaturin |
|------------------------|----------|------|--------------|-------------|
| Day0                   | 1        | 0    | 0            | 0           |
| MOLM-13.DMSO_1         | 1        | 1    | 0            | 0           |
| MOLM-13.DMSO_2         | 1        | 1    | 0            | 0           |
| MOLM-13.DMSO_3         | 1        | 1    | 0            | 0           |
| MOLM-13.Midostaurin_3  | 1        | 0    | 0            | 1           |
| MOLM-13.Midostaurin_4  | 1        | 0    | 0            | 1           |
| MOLM-13.Midostaurin_5  | 1        | 0    | 0            | 1           |
| MOLM-13.Gilteritinib_1 | 1        | 0    | 1            | 0           |
| MOLM-13.Gilteritinib_2 | 1        | 0    | 1            | 0           |
| MOLM-13.Gilteritinib_3 | 1        | 0    | 1            | 0           |

## Step 3: Write and execute entire sbatch –wrap command.

Be sure to load the slurm module.

### Entire sbatch –wrap command (all steps through MAGeCK mle):

Run Mageck MLE using “DesignMatrixPaired.txt on new.counts.txt”

``` bash
sbatch --output slurm-%x.%A.%a.log \
--mail-type END,FAIL,ARRAY_TASKS \
--mem 32G \
--cpus-per-task 4 \
-p serial \
--wrap="\
mageck mle \
-k sample1.count.txt \
-d MolmDesignMatrix.txt \
-n NameOfExperiment"
```

#### For the “-k” (counts_table) option:

``` bash
-k sample1.count.txt
```

#### For the “-d” (design_matrix_file) option:

``` bash
-d MolmDesignMatrix.txt
```

#### For the “-n” (naming_file) option:

``` bash
-n NameOfExperiment
```

### Inputs:

1.  [CRISPRScreenAnalysisWrapperPairedEnd_ver01.sh](../Scripts/CRISPRScreenAnalysisWrapperPairedEnd_ver01.sh)  
2.  [DesignMatrixAll.txt](vignette3Data/DesignMatrixAll.txt)

### Outputs:

1.  MAGeCK MLE  

-   [2022March23_All_CDC45_Day9.gene_summary.txt](vignette3Data/sample1.gene_summary.txt)
-   [2022March23_All_CDC45_Day9.sgrna_summary.txt](vignette3Data/sample1.sgrna_summary.txt)
