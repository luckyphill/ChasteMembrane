import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

## Expects irst input argument to be the file name without the _YF_x.EXT at the end
## The specific force will be adde in this script.

file_part = str(sys.argv[1])
path = "/Users/phillipbrown/Chaste/projects/ChasteMembrane/Membrane_Cell_Characterisation/driving_force/" + file_part

forces = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
velocities = np.array([0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,])

for YF in forces:
	file = path + "_YF_" + str(YF) + ".vel"
	if os.path.isfile(file):
		with open(file , 'r') as posfile:
			data_list = posfile.readlines()

		if data_list:
			[dummy, speed] = data_list[0].split(" ")
			velocities[YF-1] = float(speed.rstrip("\n"))
		else:
			velocities[YF-1] = -1

		velocities[YF-1] = float(speed.rstrip("\n"))
	else:
		velocities[YF-1] = -1

plt.plot(forces[velocities > -1],velocities[velocities > -1])
plt.title(file_part)
plt.xlabel("Force")
plt.ylabel("Velocity")
axes = plt.gca()
axes.set_xlim([-0.5, 20.5])
axes.set_ylim([min(velocities) - 0.5, max(velocities) + 0.5])

fig_name ="/Users/phillipbrown/Chaste/projects/ChasteMembrane/Membrane_Cell_Characterisation/driving_force_plots/VvsF/" + file_part + ".png"
plt.savefig(fig_name)
plt.close('all')
#plt.show()