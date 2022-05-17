Vignette1: Download fastq files from SRA and run MAGeCK count and mle
================
Chris Sansam
2020-06-16

## Step 1: Gather SRA IDs, Sample labels, and fastq filenames

### Sample labels and fastq filenames

The sample labels should be specified by you. The SRA IDs can be found
at the sequence read archive. They begin with SRR… The fastq filenames
for SRA files will be the SRA IDs with
.fastq.

| **sample\_labels**      | **fastq\_filenames**                                  |
| ----------------------- | ----------------------------------------------------- |
| Day0                    | SRR10852434.fastq,SRR10852437.fastq,SRR10852440.fastq |
| MOLM-13.DMSO\_1         | SRR10852435.fastq                                     |
| MOLM-13.DMSO\_2         | SRR10852438.fastq                                     |
| MOLM-13.DMSO\_3         | SRR10852441.fastq                                     |
| MOLM-13.Midostaurin\_3  | SRR10852443.fastq                                     |
| MOLM-13.Midostaurin\_4  | SRR10852447.fastq                                     |
| MOLM-13.Midostaurin\_5  | SRR10852450.fastq                                     |
| MOLM-13.Gilteritinib\_1 | SRR10852436.fastq                                     |
| MOLM-13.Gilteritinib\_2 | SRR10852439.fastq                                     |
| MOLM-13.Gilteritinib\_3 | SRR10852442.fastq                                     |

#### For the “-s” (SRA IDs) option:

``` bash
-s 'SRR10852434 SRR10852437 SRR10852440 SRR10852435 SRR10852438 SRR10852441 \
SRR10852443 SRR10852447 SRR10852450 SRR10852436 SRR10852439 SRR10852442'
```

#### For the “-b” (sample labels) option:

``` bash
-b 'Day0,MOLM-13.DMSO_1,MOLM-13.DMSO_2,MOLM-13.DMSO_3,MOLM-13.Midostaurin_3,\
MOLM-13.Midostaurin_4,MOLM-13.Midostaurin_5,MOLM-13.Gilteritinib_1,\
MOLM-13.Gilteritinib_2,MOLM-13.Gilteritinib_3'
```

#### For the “-q” (fastq filenames) option:

``` bash
-q 'SRR10852434.fastq,SRR10852437.fastq,SRR10852440.fastq SRR10852435.fastq \
SRR10852438.fastq SRR10852441.fastq SRR10852443.fastq SRR10852447.fastq \
SRR10852450.fastq SRR10852436.fastq SRR10852439.fastq SRR10852442.fastq'
```

## Step 2: Put the sgRNA library files in the working folder.

#### For the “-l” (sgRNA\_library\_file) option:

[broadgpp-brunello-library-corrected.txt](vignette1Data/broadgpp-brunello-library-corrected.txt)

``` bash
-l broadgpp-brunello-library-corrected.txt
```

#### For the “-n” (negative\_control\_sgRNA\_file) option:

[brunello\_nonTargeting.txt](vignette1Data/brunello_nonTargeting.txt)

``` bash
-n brunello_nonTargeting.txt
```

## Step 3: Make the design matrix and place in working folder.

#### Design matrix

| Samples                 | baseline | DMSO | Gilteritinib | Midostaurin |
| ----------------------- | -------- | ---- | ------------ | ----------- |
| Day0                    | 1        | 0    | 0            | 0           |
| MOLM-13.DMSO\_1         | 1        | 1    | 0            | 0           |
| MOLM-13.DMSO\_2         | 1        | 1    | 0            | 0           |
| MOLM-13.DMSO\_3         | 1        | 1    | 0            | 0           |
| MOLM-13.Midostaurin\_3  | 1        | 0    | 0            | 1           |
| MOLM-13.Midostaurin\_4  | 1        | 0    | 0            | 1           |
| MOLM-13.Midostaurin\_5  | 1        | 0    | 0            | 1           |
| MOLM-13.Gilteritinib\_1 | 1        | 0    | 1            | 0           |
| MOLM-13.Gilteritinib\_2 | 1        | 0    | 1            | 0           |
| MOLM-13.Gilteritinib\_3 | 1        | 0    | 1            | 0           |

#### For the “-d” (design\_matrix\_file) option:

[MolmDesignMatrix.txt](data/MolmDesignMatrix.txt)

``` bash
-d MolmDesignMatrix.txt
```

## Step 4: Write and execute entire sbatch –wrap command.

Be sure to load the slurm module.

### Entire sbatch –wrap command (all steps through MAGeCK mle):

``` bash
sbatch \
--output slurm-%x.%A.%a.log \
--mail-type END,FAIL,ARRAY_TASKS \
--mem 96 \
--cpus-per-task 12 \
-p serial \
--wrap=\
"/Volumes/Sansam/hpc-nobackup/scripts/CRISPRScreenAnalysisWrapper_ver01.sh \
-s 'SRR10852434 SRR10852437 SRR10852440 SRR10852435 SRR10852438 SRR10852441 \
SRR10852443 SRR10852447 SRR10852450 SRR10852436 SRR10852439 SRR10852442' \
-b 'Day0,MOLM-13.DMSO_1,MOLM-13.DMSO_2,MOLM-13.DMSO_3,MOLM-13.Midostaurin_3,\
MOLM-13.Midostaurin_4,MOLM-13.Midostaurin_5,MOLM-13.Gilteritinib_1,\
MOLM-13.Gilteritinib_2,MOLM-13.Gilteritinib_3' \
-q 'SRR10852434.fastq,SRR10852437.fastq,SRR10852440.fastq SRR10852435.fastq \
SRR10852438.fastq SRR10852441.fastq SRR10852443.fastq SRR10852447.fastq \
SRR10852450.fastq SRR10852436.fastq SRR10852439.fastq SRR10852442.fastq' \
-l broadgpp-brunello-library-corrected.txt \
-d MolmDesignMatrix.txt\
-n brunello_nonTargeting.txt"
```

### Inputs:

1.  [CRISPRScreenAnalysisWrapper\_ver01.sh](../CRISPRScreenAnalysisWrapper_ver01.sh)  
2.  [broadgpp-brunello-library-corrected.txt](vignette1Data/broadgpp-brunello-library-corrected.txt)  
3.  [brunello\_nonTargeting.txt](vignette1Data/brunello_nonTargeting.txt)  
4.  [MolmDesignMatrix.txt](vignette1Data/MolmDesignMatrix.txt)

### Outputs:

1.  MAGeCK
    Count  

<!-- end list -->

  - [sample1\_countsummary.R](vignette1Data/sample1_countsummary.R)  
  - [sample1\_countsummary.Rnw](vignette1Data/sample1_countsummary.Rnw)  
  - [sample1.count\_normalized.txt](vignette1Data/sample1.count_normalized.txt)  
  - [sample1.count.txt](vignette1Data/sample1.count.txt)  
  - [sample1.countsummary.txt](vignette1Data/sample1.countsummary.txt)  

<!-- end list -->

2.  MAGeCK
    mle  

<!-- end list -->

  - [sample1.gene\_summary.txt](vignette1Data/sample1.gene_summary.txt)  
  - [sample1.sgrna\_summary.txt](vignette1Data/sample1.sgrna_summary.txt)  

<!-- end list -->

3.  [multiQC
    Report](vignette1Data/fastqc_results/multiqc_results/multiqc_report.html)
