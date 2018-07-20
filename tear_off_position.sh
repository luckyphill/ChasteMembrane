



# A parameter sweep for the resting position of an epithelial cell on the basement membrane
# under varying parameters

for EMS in 10 20
	do
	for Y in $(seq 2 0.01 3.01)
		do
	
		for MIR in 2
		do
			for MPR in 1
			do
			echo "Rest position with: Y" ${Y} ", MS" ${MS} ", EMS" ${EMS} ", MIR" ${MIR} ", MPR" ${MPR};

			/Users/phillipbrown/chaste_build/projects/ChasteMembrane/test/TestMembraneCellWall -x 1 -y ${Y} -ms .5 -ems ${EMS} -mir ${MIR} -mpr ${MPR} -xf 10

			done
		done
	done
done

echo "Done"