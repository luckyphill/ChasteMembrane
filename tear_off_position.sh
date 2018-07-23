



# A parameter sweep for the resting position of an epithelial cell on the basement membrane
# under varying parameters

for EMS in 10 20
do
	for MS in .1 .2 .5
	do
		for MPR in 1 1.5 2
		do
			for MIR in 1 1.5 2 2.5
			do
				for Y in $(seq 2 .01 3.01)
				do

					echo "Rest position with: Y" ${Y} ", MS" ${MS} ", EMS" ${EMS} ", MIR" ${MIR} ", MPR" ${MPR};

					/Users/phillipbrown/chaste_build/projects/ChasteMembrane/test/TestMembraneCellWall -x .8 -y ${Y} -ms ${MS} -ems ${EMS} -mir ${MIR} -mpr ${MPR} -xf 10

				done
			done
		done
	done
done

echo "Done"