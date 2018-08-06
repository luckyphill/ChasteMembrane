for MS in 0.1 0.2 0.5
do
	for EMS in 5 10 15 20
	do
		for EES in 5 10 15 20
		do
			python /Users/phillipbrown/Chaste/projects/ChasteMembrane/Membrane_Cell_Characterisation/make_vel_vs_num_for_plots.py "MS_"${MS}"_EMS_"${EMS}"_MIR_2_MPR_0.25_EES_"${EES}"_EIR_2_EPR_0.75";
		done
	done
done