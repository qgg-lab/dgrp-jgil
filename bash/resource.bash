# ==============================
# = prepare resources for dgrp =
# ==============================

# 1. genome sequence has been downloaded
# in /mnt/research/qgg/resource/flybase/r6.36/genome
# ============================================================

java -jar /mnt/research/qgg/software/picard-2.23.3/picard.jar CreateSequenceDictionary R=/mnt/research/qgg/resource/flybase/r6.36/genome/genome.fa O=/mnt/research/qgg/resource/flybase/r6.36/genome/genome.dict > /mnt/research/qgg/resource/flybase/r6.36/genome/genome.dict.log 2>&1 &

# 2. index fa
# ============================================================

~/qgg/software/samtools-1.10/samtools faidx /mnt/research/qgg/resource/flybase/r6.36/genome/genome.fa > /mnt/research/qgg/resource/flybase/r6.36/genome/genome.fa.faidx.log 2>&1 &

# 3. variant database
# use the one I created for flybase, this is interim
# and only used for BQSR
# ============================================================

gunzip /mnt/research/qgg/resource/dgrp/dgrp.r6.vcf.gz
echo "##fileformat=VCFv4.2" > /mnt/research/qgg/resource/dgrp/dgrp.freeze2.r6.lift.vcf
sed -n '15p' /mnt/research/qgg/resource/dgrp/dgrp.r6.vcf | cut -f 1-8 >> /mnt/research/qgg/resource/dgrp/dgrp.freeze2.r6.lift.vcf
tail -n+16 /mnt/research/qgg/resource/dgrp/dgrp.r6.vcf | cut -f 1-8 >> /mnt/research/qgg/resource/dgrp/dgrp.freeze2.r6.lift.vcf

java -jar /mnt/research/qgg/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar IndexFeatureFile --input /mnt/research/qgg/resource/dgrp/dgrp.freeze2.r6.lift.vcf > /mnt/research/qgg/resource/dgrp/dgrp.freeze2.r6.lift.vcf.index.log 2>&1 &

# 4. index for bwa
# ============================================================

srun --mem=4G --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 --time=4:00:00 ~/qgg/software/bwa-0.7.17/bwa index -p /mnt/research/qgg/resource/flybase/r6.36/genome/bwa0717 /mnt/research/qgg/resource/flybase/r6.36/genome/genome.fa > /mnt/research/qgg/resource/flybase/r6.36/genome/bwa0717.index.log 2>&1 &
