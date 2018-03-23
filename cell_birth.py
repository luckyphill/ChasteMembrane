import numpy as np
import matplotlib.pyplot as plt
import math
import sys

folder = str(sys.argv[1])
path = "../../testoutput/WntWallTests/" + folder + "/results_from_time_0/cell_birth.dat"


with open(path, 'r') as posfile:
	data_list = posfile.readlines()

cell_positions = []
times = []
for t_step in data_list:
	cell_separated  = t_step.split(" | ")
	times.append(cell_separated[0])
	for cell in cell_separated[1:]:
		details = cell.split(", ")
		cell_positions.append([float(details[0]), float(details[1])])


y_pos = [x[1] for x in cell_positions]

bottom = min(y_pos)
top = max(y_pos)
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

plt.hist(y_pos)
fig_name =  "testoutput/" + folder + "_births.png"
plt.savefig(fig_name)
plt.close('all')
