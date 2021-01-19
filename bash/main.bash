# ==============
# = main steps =
# ==============

# 1. convert DGRP freeze 2 bam files to fastq's
# retain old quality score
# ============================================================

# get list of bam files and their sample names to process
for file in `find /mnt/research/qgg/dgrp/bwa/ -type f -name "*.bam"`
do
  sample=`echo $file | sed 's/\/mnt\/research\/qgg\/dgrp\/bwa\///' | sed 's/\.slx.*//'`
  echo $file $sample
done > bam.file.list

# make a new directory to store logs
mkdir log

# use parallel to submit jobs, this make sure that there are 
# at most -j 16 jobs running simultaneously
# adjust this number accordingly
# it's not good to use too many because it slows down I/O
# I usually use 16
cat bam.file.list | tr ' ' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 2 -j 16 "bash /mnt/research/qgg/resource/dgrp/bam2fqwait.bash --resource /mnt/research/qgg/resource/dgrp --bam {1} --wdir /mnt/gs18/scratch/users/tansuxu/dgrp --sample {2} > log/bam2fq.{2}.log 2>&1" > log/bam2fq.parallel.12102020.log &

# 2. get gvcf from fastq files
# follow basic GATK pipeline
# for base quality recalibration
# use DGRP sites lifted from R5 genome
# this is by no means optimal but the best we could do
# ============================================================

# need to also create a temporary file location, for example
# this is where the temporary processing will go
# if the job runs successfully, the temporary files will be deleted
# otherwise, you need to clean them up periodically when no job is running
mkdir /mnt/gs18/scratch/users/tansuxu/tmp

# create a directory to store the gvcfs for example
mkdir /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs

# in the directory where you have the bam.file.list
# --dir tells the program where to look for the fastq files
# --tmp is the temporary directory created above
# --out is where to store the results, each gvcf will create an additional
# directory within --out with the name specificied by --sample
cut -d " " -f 2 bam.file.list | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 1 -j 16 "bash /mnt/research/qgg/resource/dgrp/fq2gvcf.bash --resource /mnt/research/qgg/resource/dgrp --dir /mnt/gs18/scratch/users/tansuxu/dgrp --tmp /mnt/gs18/scratch/users/tansuxu/tmp --out /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs --sample {} > /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/{}.out 2>&1" &

# 2.1 discovery that 374 has singletons that are not properly handled by bam2fq
# modify bam2fq, re-run this line
# for both bam2fq and fq2gvcf
# ============================================================

awk '$2 == 374' bam.file.list | tr ' ' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 2 -j 16 "bash /mnt/research/qgg/resource/dgrp/bam2fqwait.bash --resource /mnt/research/qgg/resource/dgrp --bam {1} --wdir /mnt/gs18/scratch/users/tansuxu/dgrp --sample {2} > log/bam2fq.{2}.log 2>&1" > log/bam2fq.parallel.12222020.log &

# after bam2fq finishes
awk '$2 == 374' bam.file.list | cut -d " " -f 2 | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 1 -j 16 "bash /mnt/research/qgg/resource/dgrp/fq2gvcf.bash --resource /mnt/research/qgg/resource/dgrp --dir /mnt/gs18/scratch/users/tansuxu/dgrp --tmp /mnt/gs18/scratch/users/tansuxu/tmp --out /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs --sample {} > /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/{}.out 2>&1" > log/fq2gvcf.parallel.12222020.log &

# 3. combine gvcfs and genotype gvcfs
# in gvcfs directory
# make a directory to store combined gvcf and vcf files
# ============================================================

mkdir -p combine/log # this will make both combine and combine/log
cd /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/

# make sure all samples have been correctly called
# grep "completed haplotype" *.out in gvcfs
# must be 205 of them
# then make a list of gvcf files
# ============================================================

awk '{print "/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/"$2"/hc.g.vcf.gz"}' /mnt/gs18/scratch/users/tansuxu/dgrp/bam.file.list > /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/dgrp.freeze2.gvcf.list

# get interval list (whole genome)
awk '{print $1"\t0\t"$2}' ~/qgg/resource/flybase/r6.36/genome/genome.fa.fai > /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/genome.bed

# run combine gvcfs
sbatch --export=env=/mnt/research/qgg/resource/dgrp/resource.env,mem=32G,tmp=/mnt/gs18/scratch/users/tansuxu/tmp,list=dgrp.freeze2.gvcf.list,prefix=all,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine,int=/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/genome.bed --output=/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/log/combineGVCFs.out --error=/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/log/combineGVCFs.err /mnt/research/qgg/resource/dgrp/sbatch/combineGVCFs.sbatch

# after making sure that the combine GVCFs works
sbatch --export=env=/mnt/research/qgg/resource/dgrp/resource.env,mem=64G,tmp=/mnt/gs18/scratch/users/tansuxu/tmp,input=/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/all.g.vcf.gz,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine,prefix=all --mem=64G --output=/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/log/genotypeGVCFs.out --error=/mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/log/genotypeGVCFs.err /mnt/research/qgg/resource/dgrp/sbatch/genotypeGVCFs.sbatch

# after genotype GVCF is done
# ============================================================

mkdir /mnt/gs18/scratch/users/tansuxu/dgrp/jgil
gunzip -c /mnt/gs18/scratch/users/tansuxu/dgrp/gvcfs/combine/all.vcf.gz | perl /mnt/research/qgg/resource/dgrp/alleleCounts.pl > /mnt/gs18/scratch/users/tansuxu/dgrp/jgil/all.counts.out 2> /mnt/gs18/scratch/users/tansuxu/dgrp/jgil/allele.counts.log &

# when this step is done
srun --mem=4G --nodes=1 --ntasks-per-node=1 --cpus-per-task=2 --time=48:00:00 /mnt/research/qgg/software/R-3.6.0/bin/Rscript /mnt/research/qgg/resource/dgrp/jgil.R /mnt/gs18/scratch/users/tansuxu/dgrp/jgil/all.counts.out /mnt/gs18/scratch/users/tansuxu/dgrp/jgil/all.jgil.vcf > /mnt/gs18/scratch/users/tansuxu/dgrp/jgil/jgil.Rout 2>&1 &
