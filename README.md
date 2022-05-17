CRISPR Screen Analysis Usage and Results
================
Kevin Boyd & Chris Sansam
2022-05-17

![Alt text](Images/crispr_diagram.jpg)

------------------------------------------------------------------------

# General Flow

![Alt text](Images/flowChart.svg)

## Example input files (see vignette 1):

1.  [sgRNA library file](Data/broadgpp-brunello-library-corrected.txt)  
2.  [Negative control sgRNA file](Data/brunello_nonTargeting.txt)  
3.  [Design matrix](vignettes/MolmDesignMatrix.txt)

## Vignettes

-   [Vignette 1: Download fastq files from SRA and run MAGeCK Counts and
    MLE](Vignettes/CRISPRScreenAnalysis_vignette1.md)
-   [Vignette 2: Get fastq files from object store and SRA and make
    Counts table](Vignettes/CRISPRScreenAnalysis_vignette2.md)
-   [Vignette 3: Run MAGeCK MLE with Counts and Design
    Matrix](Vignettes/CRISPRScreenAnalysis_vignette3.md)
-   [Vignette 4: Plot Beta Scores from Gene
    Summary](Vignettes/CRISPRScreenAnalysis_vignette4.md)
-   [Vignette 5: Individual Guide
    Assay](Vignettes/CRISPRScreenAnalysis_vignette5.md)

------------------------------------------------------------------------

# Pipeline Bash Script

The script CRISPRScreenAnalysisWrapper_verXX.sh is to run the MAGeCK
pipeline on specified fastq files. If specified, it will also download
fastq files from SRA or the object store to a specified directory. The
usage is describe below.

### Usage:

      Usage: ./CRISPRScreenAnalysisWrapper_verXX.sh 
                [ -o Object_store_paths ]
                [ -s SRA_IDs ]
                [ -b sample_labels ]
                [ -q fastq_filenames ] 
                [ -f fastq_directory ] 
                [ -l sgRNA_library_file ]
                [ -i sample_index_file ]
                [ -d design_matrix_file ]
                [ -n negative_control_sgRNA_file ]

      Note:  
            When multiple arguments are provided for a single option, they must
            all be enclosed in a set of single quotes. For example:
            -s 'SRR4880321 SRR4880322'

[CRISPRScreenAnalysisWrapper_ver01.sh](Scripts/CRISPRScreenAnalysisWrapperPairedEnd_ver01.sh)

#### Arguments

| **option flag** | **variable name**           | **description**                                                                                                                                                                                                                                                                                                                    |
|-----------------|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -o              | Object_store_paths          | Full path of any fastq.gz files that must be transferred from the object store                                                                                                                                                                                                                                                     |
| -s              | SRA_IDs                     | SRR IDs for any fastq files to transfer from SRA using fastq-dump                                                                                                                                                                                                                                                                  |
| -b              | sample_labels               | Labels of samples. These will be the column headers for the counts table. Sample labels must be separated by commas, and they must be in the same order as the fastq_filenames. A single sample label should be given for technical replicates (ie fastq files separated by commas in “-q fastq_filenames”).                       |
| -q              | fastq_filenames             | Filenames of fastq files. These should be listed in the same oder as the sample_labels. Fastq filenames must be separated by spaces unless they are technical replicates. Technical replicates, which will be analyzed as one sample, should be separated by commas. A single sample label will be given for technical replicates. |
| -f              | fastq_directory             | Optional directory where fastq files will be stored. This option should not be used unless the user wants to store the files locally. The default and preferred is to store the fastq files in scratch space.                                                                                                                      |
| -l              | sgRNA_library_file          | Name of the required sgRNA library file. See the MAGeCK documentation for details. This file must be stored in the working directory.                                                                                                                                                                                              |
| -d              | design_matrix_file          | Name of the required design matrix file. See the MAGeCK documentation for details. This file must be stored in the working directory.                                                                                                                                                                                              |
| -n              | negative_control_sgRNA_file | Name of the required negative control file. See the MAGeCK documentation for details. This file must be stored in the working directory.                                                                                                                                                                                           |

# Possible Outcomes from Hits

## Supressors

![Alt text](Images/Supressors.jpg)

## Enhancers

![Alt text](Images/Enhancers.jpg)

# Follow up Hits

## Individual Guide Assays

Measure the IC50 of cells expressing one guide and see if the IC50
shifts. Also we can measure how much the cells grow when in different
concentrations of Auxin. This will allow to to determine if hits are
“real” or not.
