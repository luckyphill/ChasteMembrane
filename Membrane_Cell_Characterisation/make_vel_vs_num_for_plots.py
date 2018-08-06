import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

## Expects irst input argument to be the file name without the n_x part at the front or the YF part at the end
## all the values of force are included on this plot
## The specific force will be adde in this script.
file_part = str(sys.argv[1])

path = "/Users/phillipbrown/Chaste/projects/ChasteMembrane/Membrane_Cell_Characterisation/driving_force/"

forces = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
n_cells = np.array([1,2,3,4,5,6,7,8,9,10])
n = len(n_cells)
velocities = np.zeros([20,10], dtype = float)


for f in forces:
	for n in n_cells:
		file = path + "n_" + str(n) + "_" + file_part + "_YF_" + str(f) + ".vel"
		#print file
		if os.path.isfile(file):
			with open(file , 'r') as posfile:
				data_list = posfile.readlines()
			if data_list:
				[dummy, speed] = data_list[0].split(" ")
				velocities[f-1][n-1] = float(speed.rstrip("\n"))
			else:
				velocities[f-1][n-1] = -1
		else:
			# print 'failed ' + file
			velocities[f-1][n-1] = -1

for i in forces:
	y_points = velocities[i-1][velocities[i-1]> -1]
	x_points = n_cells[velocities[i-1]> -1]
	plt.plot(x_points,y_points)

plt.title(file_part)
plt.xlabel("Cells")
plt.ylabel("Velocity")
axes = plt.gca()
axes.set_xlim([0.5, 10.5])

fig_name ="/Users/phillipbrown/Chaste/projects/ChasteMembrane/Membrane_Cell_Characterisation/driving_force_plots/VvsNF/" + file_part + ".png"
plt.savefig(fig_name)
plt.close('all')
# plt.show()