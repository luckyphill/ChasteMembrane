import csv
import sys

## Expects irst input argument to be the folder location of the data
## Expects the second input argument to be the specific driving force....

folder = str(sys.argv[1])
# path = "/home/a1738927/fastdir/Chaste/testoutput/TestMultipleCellDrivingForce/" + folder + "/results_from_time_0/cell_force.dat"
path = "../../testoutput/TestMultipleCellDrivingForce/" + folder + "/results_from_time_0/cell_force.dat"

with open(path, 'r') as posfile:
	data_list = posfile.readlines()

cell_position = {}
cell_force = {}
cell_list = []
time_steps = []

# data in each row is stored as:   time | cell_id_1, x_pos, y_pos, x_force, y_force | cell_id_2 ... 
for row in data_list:
	temp_row = row.split(" | ")
	time_steps.append(float(temp_row[0].rstrip("\t")))
	for cell_data in temp_row[1:]:
		temp_cell_data = cell_data.split(",")
		if temp_cell_data[0] not in cell_list:
			cell_list.append(temp_cell_data[0])
			cell_position[temp_cell_data[0]] = []
			cell_force[temp_cell_data[0]] = []
		cell_position[temp_cell_data[0]].append(temp_cell_data[1:3])
		cell_force[temp_cell_data[0]].append(temp_cell_data[3:])
#### NOTE: ALL THIS DATA IS STORED AS STRINGS NOT INTS OR FLOATS ####

slipped = True
slip_time = -1
# cell_list[0] is the driving cell
# cell_list[1] is the first cell being driven. There will always be a value in cell_list[1]
# if the driving cell moves past the first driven cell, then the simulation won't give the results we're after and should be recorded as 'slipped'
for [pos_driving, pos_driven, time] in zip(cell_position[cell_list[0]], cell_position[cell_list[1]], time_steps):
	if pos_driving[1] > pos_driven[1]:
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
	pass
	# work out the velocity/force we're interested in
	# Write a file with a summary of the important info at the top and all the raw data thereafter



