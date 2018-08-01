#!/bin/bash 
#SBATCH -p batch 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --time=05:00:00 
#SBATCH --mem=2GB 
#SBATCH --array=0-5400
#SBATCH --err="output/Driving_force_%a.err" 
#SBATCH --output="output/Driving_force_%a.out" 
#SBATCH --job-name="Driving_force"
# NOTIFICATIONS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=phillip.j.brown@adelaide.edu.au

mkdir -p output
module load CMake
module unload OpenMPI

echo "array_job_index: $SLURM_ARRAY_TASK_ID" 
i=1
found=0
for N in $(seq 1 1 10)
do
	for MS in 0.1 0.2 0.5
	do
		for EMS in 5 10 15 20
		do
			for EES in 5 10 15 20
			do
				for YF in $(seq 1 1 20)
				do

					if [ $i = $SLURM_ARRAY_TASK_ID ]; then 
				        echo "Compression: MS" ${MS} ", EMS" ${EMS} ", EES" ${EES} ", N" ${N} ", YF" ${YF};
				        found=1 
				        break 5
				    fi
					i=$((i + 1))
				done
			done
		done
	done
done
if [ $found = 1 ]; then 
	echo "Driving force: YF" ${YF} ", MS" ${MS} ", EMS" ${EMS} ", EES" ${EES} ", N" ${N};
	echo "/home/a1738927/fastdir/chaste_build/projects/ChasteMembrane/test/TestMultipleCellDrivingForce -x .7 -y 1 -ms" ${MS} "-ems "${EMS} "-ees "${EES} "-mir 2 -eir 2 -mpr .25 -epr .75 -xf 0 -yf "${YF} "-n "${N} "-t 10 -wh 30";
	/home/a1738927/fastdir/chaste_build/projects/ChasteMembrane/test/TestMultipleCellDrivingForce -x .7 -y 1 -ms ${MS} -ems ${EMS} -ees ${EES} -mir 2 -eir 2 -mpr .25 -epr .75 -xf 0 -yf ${YF} -n ${N} -t 10 -wh 30
else 
  echo "args.csv does not have enough parameters for $SLURM_ARRAY_TASK_ID index" 
fi