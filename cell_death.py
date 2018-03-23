import numpy as np
import matplotlib.pyplot as plt
import math
import sys

folder = str(sys.argv[1])
path = "../../testoutput/WntWallTests/" + folder + "/results_from_time_0/cell_positions.dat"


with open(path, 'r') as posfile:
	data_list = posfile.readlines()

cell_positions = {}
cell_ids = []
times = []
for t_step in data_list:
	times.append({})
	cell_separated  = t_step.split(" | ")
	time = cell_separated[0]
	for cell in cell_separated[1:]:
		details = cell.split(", ")
		cell_name = details[0]
		x_pos = float(details[1])
		y_pos = float(details[2])
		times[-1][cell_name] = [cell_name, x_pos, y_pos]

death_position =[]
for i,time_point in zip(xrange(len(times)-1),times[:-1]):
	for cell in time_point.values():
		cell_name = cell[0]
		if cell_name not in times[i+1].keys():
			death_position.append(cell[2])


bottom = min(death_position)
top = max(death_position)
height = top - bottom
n_bins = int(height) + 1
size = height/n_bins
bins = []

b_l = bottom
b_r = b_l + size

while b_l < top:
	bins.append([b_l, b_r])
	b_l = b_l + size
	b_r = b_r + size
	if b_r > top:
		b_r = top

if bins[-1][1] - bins[-1][0] < size/2:
	bins[-2][1] = top
	del bins[-1]


plt.hist(death_position)
fig_name =  "testoutput/" + folder + "_deaths.png"
plt.savefig(fig_name)
plt.close('all')

