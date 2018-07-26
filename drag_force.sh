for EMS in 10 20
do
	for MS in .1 .2 .5
	do
		for MPR in 1 1.5 2
		do
			for MIR in 1 1.5 2 2.5
			do
				for YF in 1 2 3 4 5 7 10 12 15 18 20
				do

					echo "Drag force with: YF" ${YF} ", MS" ${MS} ", EMS" ${EMS} ", MIR" ${MIR} ", MPR" ${MPR};

					/Users/phillipbrown/chaste_build/projects/ChasteMembrane/test/TestMembraneCellWall -x .8 -y 1 -ms ${MS} -ems ${EMS} -mir ${MIR} -mpr ${MPR} -yf ${YF}
					python /Users/phillipbrown/Chaste/projects/ChasteMembrane/plot_drag.py "Drag_force_MS_"${MS}"_EMS_"${EMS}"_MIR_"${MIR}"_MPR_"${MPR}".txt" ${YF}
				done
			done
		done
	done
done