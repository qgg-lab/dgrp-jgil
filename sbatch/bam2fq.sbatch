#!/bin/bash

# Job name:
#SBATCH --job-name=bam2fq
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=2
#
# Memory per node:
#SBATCH --mem=5G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
# ============================================================

source "$env"

# get fastq files, collated, may contain multiple read groups
# ============================================================

$SAMTOOLS collate -O $bam $dir/tmp | $SAMTOOLS fastq -O -1 $dir/pe_1.fastq.gz -2 $dir/pe_2.fastq.gz -0 $dir/se.fastq.gz -s $dir/singleton.fastq.gz -c 6 2> $dir/bam2fq.log

# consolidate se and singleton reads
# ============================================================

zcat $dir/se.fastq.gz $dir/singleton.fastq.gz | gzip > $dir/tmp.fastq.gz
mv $dir/tmp.fastq.gz $dir/se.fastq.gz
rm $dir/singleton.fastq.gz

# get run info for se reads if there are any
# ============================================================

nse=`gunzip -c $dir/se.fastq.gz | head | wc -l`
if [ $nse -eq 10 ]
then
	gunzip -c $dir/se.fastq.gz | paste - - - - | cut -f 1 | shuf | head -n 10000 | sed 's/^@//' | awk -F ":" '{print "se."$1"."$2" "$1" "$2}' | uniq | sort | uniq | awk '$3 != ""' > $dir/file.info
fi

# get run info for pe reads if there are any
# ============================================================

npe=`gunzip -c $dir/pe_1.fastq.gz | head | wc -l`
if [ $npe -eq 10 ]
then
	gunzip -c $dir/pe_1.fastq.gz | paste - - - - | cut -f 1 | shuf | head -n 10000 | sed 's/^@//' | awk -F ":" '{print "pe."$1"."$2" "$1" "$2}' | uniq | sort | uniq | awk '$3 != ""' >> $dir/file.info
fi

# line 360 is the only one with only a small number of reads not conforming to the format
# in read names machineID:laneNumber:tile...

# split fastq files
# ============================================================

while read line           
do           
	ins=$(echo $line | cut -d " " -f 2)
	lane=$(echo $line | cut -d " " -f 3)
	rg=$(echo $line | cut -d " " -f 1)

	if [[ $rg == "pe"* ]]
  then
    gunzip -c $dir/pe_1.fastq.gz | grep -A 3 ^"@$ins:$lane" --no-group-separator | gzip > $dir/"pe.$ins.$lane"_1.fastq.gz
    gunzip -c $dir/pe_2.fastq.gz | grep -A 3 ^"@$ins:$lane" --no-group-separator | gzip > $dir/"pe.$ins.$lane"_2.fastq.gz
  fi
  
  if [[ $rg == "se"* ]]
  then
    gunzip -c $dir/se.fastq.gz | grep -A 3 ^"@$ins:$lane" --no-group-separator | gzip > $dir/"se.$ins.$lane"_1.fastq.gz
  fi
done < $dir/file.info

# change fastq files to read only
# ============================================================

rm $dir/se.fastq.gz $dir/pe_1.fastq.gz $dir/pe_2.fastq.gz
chmod 440 "$dir"/*.fastq.gz

echo $(date +"%m-%d-%Y-%T"):info: done.main.process
