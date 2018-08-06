import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

## Expects irst input argument to be the file name without the n_x part at the front
## the force is fixed for each plot
## The specific force will be adde in this script.
file_part = str(sys.argv[1])

path = "/Users/phillipbrown/Chaste/projects/ChasteMembrane/Membrane_Cell_Characterisation/driving_force/"

n_cells = np.array([1,2,3,4,5,6,7,8,9,10])
velocities = np.array([0.0,0,0,0,0,0,0,0,0,0])

for n in n_cells:
	file = path + "n_" + str(n) + "_" + file_part + ".vel"
	#print file
	if os.path.isfile(file):
		with open(file , 'r') as posfile:
			data_list = posfile.readlines()

		[dummy, speed] = data_list[0].split(" ")
		#print dummy, speed.rstrip("\n")

		velocities[n-1] = float(speed.rstrip("\n"))
	else:
		print 'failed'

end_index = 20
for i in n_cells[:-1]: # not strictly good form, but it is convenient because n_cells has the indices I need in order
	if velocities[i] == 0:
		end_index = i
		break

plt.plot(n_cells[:end_index],velocities[:end_index])
plt.title(file_part)
plt.xlabel("Cells")
plt.ylabel("Velocity")
axes = plt.gca()
axes.set_xlim([-0.5, 10.5])
axes.set_ylim([min(velocities) - 0.5, max(velocities) + 0.5])

fig_name ="/Users/phillipbrown/Chaste/projects/ChasteMembrane/Membrane_Cell_Characterisation/driving_force_plots/VvsN/" + file_part + ".png"
plt.savefig(fig_name)
plt.close('all')
#plt.show()

from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

popt, pcov = curve_fit(func, n_cells[velocities>0.01], velocities[velocities>0.01])

plt.figure()
plt.plot(n_cells, velocities, 'ko', label="Original Noised Data")
plt.plot(n_cells, func(n_cells, *popt), 'r-', label="Fitted Curve")
plt.legend()
plt.show()