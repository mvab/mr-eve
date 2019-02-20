#!/bin/bash

#SBATCH --job-name=mr
#SBATCH --array=0-1000
#SBATCH --nodes=1 --cpus-per-task=10 --time=0-12:00:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --mem=20G

echo "Running on ${HOSTNAME}"
module add languages/r/3.4.4

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}
# i=$((i + 1000))


GWASDIR='../gwas-files'
OUTPUT="$GWASDIR/$id/derived/instruments/mr.rdata"
IDLIST="resources/idlist.txt"
IDINFO="resources/idinfo.rdata"
RF="resources/rf.rdata"
THREADS="10"

cd ${HOME}/mr-eve/mr-eve
ids=($(cat $IDLIST))

echo "Total number of ids: ${#ids[@]}"
id=`echo ${ids[$i]}`
echo "Analysing $id"
if [ -z "$id" ]
then
	echo "outside range"
	exit
fi



time Rscript scripts/mr.r --idlist $IDLIST --gwasdir $GWASDIR --id $id --rf $RF --what eve --out $OUTPUT --threads $THREADS --idinfo $IDINFO
