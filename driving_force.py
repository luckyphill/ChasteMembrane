import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

## Expects irst input argument to be the folder location of the data
## Expects the second input argument to be the specific driving force....

folder = str(sys.argv[1])
# path = "/home/a1738927/fastdir/Chaste/testoutput/TestMultipleCellDrivingForce/" + folder + "/results_from_time_0/cell_force.dat"
path = "../../testoutput/TestMultipleCellDrivingForce/" + folder + "/results_from_time_0/cell_force.dat"

with open(path, 'r') as posfile:
	data_list = posfile.readlines()

cell_data = {}
cell_list = []
time_steps = []

# data in each row is stored as:   time | cell_id_1, x_pos, y_pos, x_force, y_force | cell_id_2 ... 
for row in data_list:
	temp_row = row.split(" | ")
	time_steps.append(float(temp_row[0].rstrip("\t")))
	for cell_info in temp_row[1:]:
		temp_cell_data = cell_info.split(",")
		current_cell = temp_cell_data[0]
		if current_cell not in cell_list:
			cell_list.append(current_cell)
			cell_data[current_cell] = {}
			cell_data[current_cell]['x_pos'] = []
			cell_data[current_cell]['y_pos'] = []
			cell_data[current_cell]['x_for'] = []
			cell_data[current_cell]['y_for'] = []
			
		cell_data[current_cell]['x_pos'].append(float(temp_cell_data[1]))
		cell_data[current_cell]['y_pos'].append(float(temp_cell_data[2]))
		cell_data[current_cell]['x_for'].append(float(temp_cell_data[3]))
		cell_data[current_cell]['y_for'].append(float(temp_cell_data[4]))
#### NOTE: ALL THIS DATA IS STORED AS STRINGS NOT INTS OR FLOATS ####

dt = time_steps[1]

slipped = False
slip_time = -1
# cell_list[0] is the driving cell
# cell_list[1] is the first cell being driven. There will always be a value in cell_list[1]
# if the driving cell moves past the first driven cell, then the simulation won't give the results we're after and should be recorded as 'slipped'
for [pos_driving, pos_driven, time] in zip(cell_data[cell_list[0]]['y_pos'], cell_data[cell_list[1]]['y_pos'], time_steps):
	if pos_driving > pos_driven:
		# the driving cell has slipped past and the simulation is done
		slipped = True
		slip_time = time
		break


if slipped:
	filename = folder + ".FAILED"
	with open(filename, 'w') as output_file:
		output_file.write("Slip occurred at " + str(slip_time) + "s")
	# code to identify the failed simulation
else:
	for cell in cell_list:
		y_positions = np.array(cell_data[cell]['y_pos'])
		#y_speed = (y_positions[1:] - y_positions[:-1])/dt
		# plt.figure()
		# plt.subplot(311)
		# plt.scatter(time_steps,y_positions)
		# plt.subplot(312)
		# plt.scatter(time_steps[1:], y_speed)
		# plt.subplot(313)
		# plt.scatter(time_steps, cell_data[cell]['y_for'])
		slope, intercept, r_value, p_value, std_err = stats.linregress(time_steps[3000:6000],y_positions[3000:6000])
		# print slope
		filename = folder + ".vel"
		with open(filename, 'a') as output_file:
			output_file.write(cell + " " + str(slope) + "\n")
		# slope, intercept, r_value, p_value, std_err = stats.linregress(time_steps[1000:3000],y_speed[1000:3000])
		# print intercept
		# slope, intercept, r_value, p_value, std_err = stats.linregress(time_steps[1000:3000],cell_data[cell]['y_for'][1000:3000])
		# print intercept
		# print "\n"
	# plt.show()

	# work out the velocity/force we're interested in
	# Write a file with a summary of the important info at the top and all the raw data thereafter

	# IMPORTANT DATA
	# Cell velocity at cruising speed given number of cells to push and force applied
	



