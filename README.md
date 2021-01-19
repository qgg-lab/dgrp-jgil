# DGRP-JGIL

Pipeline to call variants in DGRP using BWA-GATK-JGIL (Eric Stone, 2012)

## Prerequisites

The pipeline requires a set of softwares (e.g. BWA, GATK, R, etc.) and the genome assembly etc. They are all specified in [resource/resource.env](resource/resource.env). This file is sourced in most scripts that require the softwares and files.

## Main steps

### 1. Prepare FastQ files

There must be one directory for each sample, within each sample, there can ben multiple sets of fastq files. For example, the directory structure should look like:

```bash
[huangw53@dev-intel18 379]$ ls -rlt
total 2007297
-r--r----- 1 huangw53 qgg 435080631 Jan  1 00:00 pe.USI-EAS034_1_PE_FC30C50AAXX.4_1.fastq.gz
-r--r----- 1 huangw53 qgg 390866335 Jan  1 00:00 se.USI-EAS034_1_PE_FC30MHJAAXX.4_1.fastq.gz
-r--r----- 1 huangw53 qgg 434065487 Jan  1 00:00 pe.USI-EAS034_1_PE_FC30C50AAXX.4_2.fastq.gz
-r--r----- 1 huangw53 qgg 795328105 Jan  1 00:00 se.USI-EAS376_2_PE1_FC30HHMAAXX.7_1.fastq.gz
-rw-r----- 1 huangw53 qgg       191 Jan  1 00:00 file.info
[huangw53@dev-intel18 379]$ cat file.info 
se.USI-EAS034_1_PE_FC30MHJAAXX.4 USI-EAS034_1_PE_FC30MHJAAXX 4
se.USI-EAS376_2_PE1_FC30HHMAAXX.7 USI-EAS376_2_PE1_FC30HHMAAXX 7
pe.USI-EAS034_1_PE_FC30C50AAXX.4 USI-EAS034_1_PE_FC30C50AAXX 4
```

There is a file `file.info` that stores the fasqt information. This file has three columns, first is an identifier for the file, second is the flowcell identifier (unique), third is the lane number within the flowcell. This file is sufficient to determine that there are four `fastq` files associated with this sample. Two single end files and one paired end samples. This file is used to loop through the `fastq` files in the pipeline.

### 2. Run `fq2gvcf.bash`

This is the script to convert `fastq` to `gvcf`. It is a wrapper around a bunch of `sbatch` scripts that are in the directory [sbatch/](sbatch/). A basic command is like this:

```bash
# 
bash /mnt/research/qgg/resource/dgrp/fq2gvcf.bash --resource /mnt/research/qgg/resource/dgrp --dir /mnt/gs18/scratch/users/tansuxu/dgrp --tmp /mnt/gs18/scratch/users/tansuxu/tmp --out /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs --sample 379 > /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/379.out 2>&1 &
```

You can either loop through all samples, which will queue everything and is perfectly fine. Alternatively, this `bash` script is designed to be used in combination with the `GNU parallel` that can specify the maximum number of jobs being run. Refer to the `bash` script [bash/main.bash](bash/main.bash).

### 3. Run `combineGVCFs.sbatch` and `genotypeGVCFs.sbatch`

These two `sbatch` scripts take the GVCFs