#!/bin/bash
#
# Script to run Wnt Wall simultations in series




CYCLE = 1
E_S = 1
EM_S = 1

for CYCLE in 1 2 5 10 20
do
	for E_S in 1 2 5 10 20
	do
		for EM_S in 1 2 5 10 20 # need to put 1 2 5 back in after three simulations
		do
		echo "Wnt Wall simulation with Cycle time " ${CYCLE} " Epithelial Stiffness " ${E_S} " and Epithleial-Membrane Stiffness " ${EM_S};
		# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
		# ">" directs std::cout to the file.
		# "2>&1" directs std::cerr to the same place.
		# "&" on the end lets the script carry on and not wait until this has finished.
		/Users/phillipbrown/chaste_build/projects/ChasteMembrane/test/TestMembraneCellCrypt -e ${E_S} -em ${EM_S} -ct ${CYCLE}
		python cell_birth.py E_${E_S}EM_${EM_S}CCT_${CYCLE}
		python cell_death.py E_${E_S}EM_${EM_S}CCT_${CYCLE}
		python cell_speed.py E_${E_S}EM_${EM_S}CCT_${CYCLE}

		done
	done
done

echo "Done"