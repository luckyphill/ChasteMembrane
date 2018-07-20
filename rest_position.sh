



# A parameter sweep for the resting position of an epithelial cell on the basement membrane
# under varying parameters

for Y in $(seq 2.51 0.01 3.01)
do
	for MS in 0.2 0.5 .1 .2 .5
	do
		for EMS in 10
		do
			for MIR in 1.0 1.25 1.5 2.0 2.5
			do
				for MPR in .05 .25 .45 .65
				do
				echo "Rest position with: Y" ${Y} ", MS" ${MS} ", EMS" ${EMS} ", MIR" ${MIR} ", MPR" ${MPR};

				/Users/phillipbrown/chaste_build/projects/ChasteMembrane/test/TestMembraneCellWall -x 0.45 -y ${Y} -ms ${MS} -ems ${EMS} -mir ${MIR} -mpr ${MPR}

				done
			done
		done
	done
done

echo "Done"