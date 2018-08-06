#!/bin/bash 
#SBATCH -p batch 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --time=10:00:00 
#SBATCH --mem=2GB 
#SBATCH --array=0-10000
#SBATCH --err="output/Driving_force_%a.err" 
#SBATCH --output="output/Driving_force_%a.out" 
#SBATCH --job-name="Driving_force"
# NOTIFICATIONS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=phillip.j.brown@adelaide.edu.au

mkdir -p output
module load CMake
module unload OpenMPI
module load numpy

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
				        echo "Driving force: MS" ${MS} ", EMS" ${EMS} ", EES" ${EES} ", N" ${N} ", YF" ${YF};
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
	/home/a1738927/fastdir/chaste_build/projects/ChasteMembrane/test/TestMultipleCellDrivingForce -x .7 -y 1 -ms ${MS} -ems ${EMS} -ees ${EES} -mir 2 -eir 2 -mpr .25 -epr .75 -xf 0 -yf ${YF} -n ${N} -t 10 -wh 60
	python /home/a1738927/fastdir/Chaste/projects/ChasteMembrane/driving_force.py "n_"${N}"_MS_"${MS}"_EMS_"${EMS}"_MIR_2_MPR_0.25_EES_"${EES}"_EIR_2_EPR_0.75_YF_"${YF}
else 
  echo "args.csv does not have enough parameters for $SLURM_ARRAY_TASK_ID index" 
fi