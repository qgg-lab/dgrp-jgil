#! /bin/bash
# ============================================================
# bash to call sbatch bam2fq and work with parallel
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,w:,s:,b:" -l "resource:,wdir:,sample:,bam:" -- "$@"`
>&2 echo "command arguments given: $args"
eval set -- "$args"

# parse arguments
# ============================================================

while true;
do
  case $1 in

    -r|--resource)
      resource=$2
      shift 2;;

    -b|--bam)
      bam=$2
      shift 2;;

    -w|--wdir)
			wdir=$2
			shift 2;; 

    -s|--sample)
		  sample=$2
			shift 2;;

    --)
      shift
      break;;
      
  esac
done

# check for file availability
# ============================================================

if [[ -e $resource/resource.env && -e $resource/sbatch/bam2fq.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi

if [[ -e $bam ]]
then
  >&2 echo "converting $bam."
else
  >&2 echo "cannot find bam file."
  exit 1
fi

# run sbatch
# ============================================================

mkdir "$wdir"/"$sample"

jobID=$(sbatch --export=env=$resource/resource.env,dir=$wdir/$sample,bam=$bam,tmp=/mnt/ufs18/scratch/huangw53/tmp/$sample --output="$wdir"/"$sample"/bam2fq.out --error="$wdir"/"$sample"/bam2fq.err $resource/sbatch/bam2fq.sbatch | cut -d " " -f 4)

sleep 1m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
  sleep 5m
  jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

reRunCounter=0
until [ `wc -l "$wdir"/"$sample"/bam2fq.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$wdir"/"$sample"/bam2fq.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting bam2fq job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> $wdir/$sample/bam2fq.fail.out
	bamTime=$(($reRunCounter * 4 + 8))
	bamMem=$(($reRunCounter * 4 + 8))G
  if [[ -e "$wdir"/"$sample"/file.info ]]
  then
    rm "$wdir"/"$sample"/file.info
    rm -f "$wdir"/"$sample"/*.fastq.gz
    rm "$wdir"/tmp.*.bam
  fi
  jobID=$(sbatch --export=env=$resource/resource.env,dir=$wdir/$sample,bam=$bam,tmp=/mnt/ufs18/scratch/huangw53/tmp/$sample --time=$bamTime:00:00 --mem=$bamMem --output="$wdir"/"$sample"/bam2fq.out --error="$wdir"/"$sample"/bam2fq.err $resource/sbatch/bam2fq.sbatch | cut -d " " -f 4)  
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l "$wdir"/"$sample"/bam2fq.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$wdir"/"$sample"/bam2fq.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed bam2fq conversion."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> "$wdir"/"$sample"/bam2fq.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):info: error for bam2fq conversion."
	exit 1
fi
