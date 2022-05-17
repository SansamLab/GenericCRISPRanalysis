Vignette 2: Get fastq files from object store and SRA and make count
table
================
Chris Sansam
2020-06-16

## Step 1: Gather SRA IDs, Object\_store\_paths, Sample labels, and fastq filenames

| **sample label**                 | **fastq filename**                     | **Object Store Path or SRA ID**                                                                         |
| -------------------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Sansam\_Brunello\_lentiguide     | brunello\_TNscreen\_S32\_R1\_001.fastq | o3://LDAP\_o3-sansamc/2020/2020Mar19\_TDN\_CRISPRLibraryPreps/brunello\_TNscreen\_S32\_R1\_001.fastq.gz |
| SRR8983579\_Brunello\_lentiguide | SRR8983579.fastq                       | SRR8983579                                                                                              |
| SRR6737232\_Brunello\_lentiguide | SRR6737232.fastq                       | SRR6737232                                                                                              |
| SRR8297998\_Brunello\_lentiguide | SRR8297998.fastq                       | SRR8297998                                                                                              |

#### For the “-o” (Object\_store\_paths) option:

``` bash
-o 'o3://LDAP_o3-sansamc/2020/2020Mar19_TDN_CRISPRLibraryPreps/brunello_TNscreen_S32_R1_001.fastq.gz'
```

#### For the “-s” (SRA IDs) option:

``` bash
-s 'SRR8983579 SRR6737232 SRR8297998'
```

#### For the “-b” (sample labels) option:

``` bash
-b 'Sansam_Brunello_lentiguide,SRR8983579_Brunello_lentiguide,\
SRR6737232_Brunello_lentiguide,SRR8297998_Brunello_lentiguide'
```

#### For the “-q” (fastq filenames) option:

``` bash
-q 'brunello_TNscreen_S32_R1_001.fastq SRR8983579.fastq \
SRR6737232.fastq SRR8297998.fastq'
```

## Step 2: Put the sgRNA library files in the working folder.

#### For the “-l” (sgRNA\_library\_file) option:

[broadgpp-brunello-library-corrected.txt](vignette1Data/broadgpp-brunello-library-corrected.txt)

``` bash
-l broadgpp-brunello-library-corrected.txt
```

## Step 3: Write and execute entire sbatch –wrap command.

### Sbatch -wrap command

``` bash
sbatch \
--output slurm-%x.%A.%a.log \
--mail-type END,FAIL,ARRAY_TASKS \
--mem 96G \
--cpus-per-task 12 \
-p serial \
--wrap=\
"/Volumes/Sansam/hpc-nobackup/scripts/CRISPRScreenAnalysisWrapper_ver01.sh \
-o 'o3://LDAP_o3-sansamc/2020/2020Mar19_TDN_CRISPRLibraryPreps/brunello_TNscreen_S32_R1_001.fastq.gz' \
-s 'SRR8983579 SRR6737232 SRR8297998' \
-b 'Sansam_Brunello_lentiguide,SRR8983579_Brunello_lentiguide,\
SRR6737232_Brunello_lentiguide,SRR8297998_Brunello_lentiguide' \
-q 'brunello_TNscreen_S32_R1_001.fastq SRR8983579.fastq \
SRR6737232.fastq SRR8297998.fastq' \
-l broadgpp-brunello-library-corrected.txt"
```

### Output

[output files folder](vignette2Data/vignette2_output/)
