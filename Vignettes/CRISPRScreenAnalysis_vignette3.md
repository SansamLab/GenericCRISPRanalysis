Vignette3: Run MAGeCK MLE
================
Kevin Boyd
2022-04-26

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
counts_1 <- read.table(file="~/Dropbox (OMRF)/Github/CRISPRScreenAnalysis/Vignettes/vignette3Data/All_Day9_CDC45_Counts_pDNA.txt", header = T)
counts_2 <- read.table(file="~/Dropbox (OMRF)/Github/CRISPRScreenAnalysis/Vignettes/vignette3Data/sample1.count.txt", header = T)

# order the guides in column 2 in the same order as column 1
counts_2 <- counts_2[match(counts_1$sgRNA,counts_2$sgRNA),]

# check to see if columns are equal
all.equal(counts_1$sgRNA,counts_2$sgRNA)
#> [1] TRUE

# remove first 2 columns from both dataframes and add the similar values together
# new.counts <- counts_1[c(8:11)] + counts_2[c(4:7)]

# combine counts tables together. Use when not adding additional counts.
new.counts <- cbind(counts_1,counts_2[,-c(1,2)])
head(new.counts)
#>                                    sgRNA    Gene Sansam_Brunello_lentiCRISPRv2
#> 354   CABLES2_52045_TGTGTGACCAGGGCGTTGGT CABLES2                          1630
#> 39152  OR52E6_71598_ATGATACGATGTCCACAGAA  OR52E6                          1453
#> 33795   ZC3H8_53739_AAGGATCAAATGCTTTGCTG   ZC3H8                          1441
#> 48468    CLDN4_3791_CAAGTGTACCAACTGCCTGG   CLDN4                          1767
#> 75177    PPAN_44242_TGGCTCAGTGAGCAAGACGG    PPAN                           421
#> 9663  ZDHHC12_54363_AGGGATGGCTGGAGGAACCA ZDHHC12                          1524
#>       CDC45_Day9_May5 CDC45_Day9_PlusAux_May5 CDC45_Day9_May19
#> 354               884                     574             1110
#> 39152             564                     301              551
#> 33795             528                     310              484
#> 48468             739                     535              928
#> 75177              51                      16               56
#> 9663             1139                     542              977
#>       CDC45_Day9_PlusAux_May19 CDC45_Day9_PlusAux_Rep1 CDC45_Day9_PlusAux_Rep2
#> 354                       1861                     755                     810
#> 39152                      892                     328                     581
#> 33795                     2058                     328                     503
#> 48468                     1558                     216                     814
#> 75177                      209                      39                      88
#> 9663                      1504                     707                    1286
#>       CDC45_Day9_Rep1 CDC45_Day9_Rep2 CDC45_Day12_Jan21_Rep1
#> 354               997             609                     86
#> 39152             484             234                     32
#> 33795             419             265                     31
#> 48468             729             545                     60
#> 75177              57              15                      2
#> 9663              706             535                     63
#>       CDC45_Day21_.Aux_Jan21_Rep1 CDC45_Day12_Jan21_Rep2
#> 354                           136                     54
#> 39152                          57                     17
#> 33795                          24                     18
#> 48468                           1                     53
#> 75177                           0                      0
#> 9663                          324                     36
#>       CDC45_Day21_.Aux_Jan21_Rep2
#> 354                            63
#> 39152                          17
#> 33795                          32
#> 48468                          80
#> 75177                           1
#> 9663                          183

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

| Samples                  | baseline | NoAuxin | Auxin | May5 | May19 | Nov23\_Rep1 | Nov23\_Rep2 |
|--------------------------|----------|---------|-------|------|-------|-------------|-------------|
| lentiCRISPRv2\_pDNA      | 1        | 1       | 0     | 0    | 0     | 0           | 0           |
| Day9 May5                | 1        | 1       | 0     | 1    | 0     | 0           | 0           |
| Day9 PlusAux May5        | 1        | 0       | 1     | 1    | 0     | 0           | 0           |
| Day9 May19               | 1        | 1       | 0     | 0    | 1     | 0           | 0           |
| Day9 PlusAux May19       | 1        | 0       | 1     | 0    | 1     | 0           | 0           |
| Day9 Nov23\_Rep1         | 1        | 1       | 0     | 0    | 0     | 1           | 0           |
| Day9 PlusAux Nov23\_Rep1 | 1        | 0       | 1     | 0    | 0     | 1           | 0           |
| Day9 Nov23\_Rep1         | 1        | 1       | 0     | 0    | 0     | 0           | 1           |
| Day9 PlusAux Nov23\_Rep2 | 1        | 0       | 1     | 0    | 0     | 0           | 1           |

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
-k AllCounts_Day9_pDNA_AdditionalReads.txt \
-d DesignMatrixAll.txt \
-n 2022March23_All_CDC45_Day9"
```

#### For the “-k” (counts\_table) option:

``` bash
-k AllCounts_Day9_pDNA_AdditionalReads.txt
```

#### For the “-d” (design\_matrix\_file) option:

``` bash
-d DesignMatrixAll.txt
```

#### For the “-n” (naming\_file) option:

``` bash
-n 2022March23_All_CDC45_Day9
```

### Inputs:

1.  [CRISPRScreenAnalysisWrapperPairedEnd\_ver01.sh](../Scripts/CRISPRScreenAnalysisWrapperPairedEnd_ver01.sh)  
2.  [DesignMatrixAll.txt](vignette3Data/DesignMatrixAll.txt)

### Outputs:

1.  MAGeCK MLE  

-   [2022March23\_All\_CDC45\_Day9.gene\_summary.txt](vignette3Data/2022March23_All_CDC45_Day9.gene_summary.txt)
-   [2022March23\_All\_CDC45\_Day9.sgrna\_summary.txt](vignette3Data/2022March23_All_CDC45_Day9.sgrna_summary.txt)
-   [2022March23\_All\_CDC45\_Day9.log](vignette3Data/2022March23_All_CDC45_Day9.log)
