import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import pickle
import time

plt.ion()

# this could be helpful too:
# https://stackoverflow.com/questions/11140163/plotting-a-3d-cube-a-sphere-and-a-vector-in-matplotlib

def convert_positions_to_circle(num_positions, circle_radius):
	"""
	this function converts a list of positions of length N
	to surround a circle of radius X
	"""
	theta = float(2 * np.pi) / float(num_positions)

	circle_coords = []
	for i in range(num_positions):
		angle = theta * float(i)

		x = float(circle_radius) * math.cos(angle)
		y = float(circle_radius) * math.sin(angle)

		circle_coords.append([x, y])

	return circle_coords


def calculate_arrowdelta(initial_x, initial_y, final_x, final_y, initial_z, final_z):
	"""
	Find the dx, dy for the arrow from start, final
	"""
	dx = final_x - initial_x
	dy = final_y - initial_y
	dz = final_z - initial_z

	return dx, dy, dz


def generate_color_dict(unique_symbols):
	"""
	this function returns a color map for each unique symbol via index
	"""
	color_map = {}
	"""
	for i in range(len(unique_symbols)):
		color = "C" + str(i)
		color_map[i] = color
	"""
	colors = cm.rainbow(np.linspace(0, 1, len(unique_symbols)))
	for i in range(len(unique_symbols)):
		color = colors[i]
		color_map[i] = color

	return color_map


def set_pointcircle_maxradius(cylinder_height, unique_symbols, symbol_points_sizeshrink):
	"""
	this function determines the maximum radius that a symbol point can have
	by limiting the total diameter to fit Y * the number of unqiue symbols
	"""
	cicle_diameter = float(cylinder_height) / float(len(unique_symbols) * symbol_points_sizeshrink)

	return cicle_diameter / 2


def draw_symbolpoint_circle(ax, x, y, z, radius):
	"""
	This function draws a symbol point-circle at point (x, y, z) with radius
	"""
	num_pnts = 20
	u = np.linspace(0, 2* np.pi, num_pnts)
	v = np.linspace(0, np.pi, num_pnts) 

	ax.scatter(x, y, z, radius)

	return


def generate_single_plot(symbol_position_tracking, symbol_edge_tracking, cylinder_height, \
													unique_symbols, vector_color, color_map, circle_coords, max_symbolpoint_radius, 
													p_thresh_max, p_thresh_min, pkl_plot_folder, vertical_extension, vector_width_binom, radius_mult, \
													pfam_family_id):
	"""
	Writes out a single plot, with a p-value threshold
	"""
	fig = plt.figure()
	ax = fig.add_subplot(111, projection = '3d')

	# add vectors
	num_kept, num_thrown = 0, 0
	for position_initial, symbol_index_initial, position_final, symbol_index_final, pvalue in symbol_edge_tracking:
		if pvalue <= p_thresh_max and pvalue >= p_thresh_min:
			num_kept += 1

			initial_x, initial_y = circle_coords[position_initial]
			final_x, final_y = circle_coords[position_final]

			initial_z = float(symbol_index_initial) / float(len(unique_symbols)) * float(cylinder_height)
			final_z = float(symbol_index_final) / float(len(unique_symbols)) * float(cylinder_height)

			dx, dy, dz = calculate_arrowdelta(initial_x, initial_y, final_x, final_y, initial_z, final_z)

			# set pvalue to vector weight (e.g., 1.0 - pvalue)
			vector_weight = (1.0 - pvalue) ** vector_width_binom
			
			vector_length = math.sqrt((final_x - initial_x) ** 2 + (final_y - initial_y) ** 2 + (final_z - initial_z) ** 2)
			ax.quiver(initial_x, initial_y, initial_z, dx, dy, dz, linewidths = vector_weight, color = vector_color, \
									arrow_length_ratio = 0)
									#scale_units='xy', scale = 1)

		else:
			num_thrown += 1

	print("# of vectors kept:", num_kept, " & thrown out", num_thrown)

	# Add symbols
	for position, symbol_index, fract in symbol_position_tracking:
		x, y = circle_coords[position]
		z = float(symbol_index) / float(len(unique_symbols)) * float(cylinder_height)

		radius = float(fract) * max_symbolpoint_radius
		color = color_map[symbol_index]

		ax.scatter(x, y, z, s = radius * radius_mult, c = [color], alpha = 1)

	# Draw each position's vertical, add position text
	for i in range(len(circle_coords)):
		x, y = circle_coords[i]

		initial_x, final_x = x, x
		initial_y, final_y = y, y
		initial_z, final_z = 0, vertical_extension * cylinder_height
		ax.plot(xs = [initial_x, final_x], ys = [initial_y, final_y], zs = [initial_z, final_z], color = "grey")

		position_text = str(i + 1)
		ax.text(initial_x, initial_y, final_z, position_text, color = "yellow")

	# Draw the z-axis labels (e.g., symbols)
	for i in range(len(unique_symbols)):
		x = min(circle_coords, key = lambda l: l[0])[0]
		y = min(circle_coords, key = lambda l: l[1])[1]

		z = float(i) / float(len(unique_symbols)) * float(cylinder_height)

		symbol = unique_symbols[i]
		color = color_map[i]

		ax.text(x, y, z, symbol, color = color, fontsize = 8)


	# graph parameters
	ax.set_facecolor("black")
	ax.grid(False)
	ax.axis('off')
	ax.set_title("PFAM Family ID:" + pfam_family_id + "\n" + \
									str(num_kept) + "Vectors In Correlation Plot P-Value Cutoffs:" + \
									str(p_thresh_min) + "-" + str(p_thresh_max), color = "blue")

	# dump interactive images into PKL files to re-open later
	# NOTE - this doesnt' currently work. MatPlotLib developers ackowledge and supposedly are working on a fix
	#pickle.dump(fig, open(pkl_plot_folder + "pcut=" + str(pvalue_threshold) + ".pickle", 'wb'))

	#plt.show()
	plt.draw()
	#plt.close()

	return


def generate_vector_plots(sequence_length, symbol_position_tracking, symbol_edge_tracking, circle_radius, cylinder_height, \
													unique_symbols, vector_color, symbol_points_sizeshrink, pkl_plot_folder, pvalue_thresholds, \
													vertical_extension, vector_width_binom, radius_mult, pfam_family_id):
	"""
	this function plots a series of vectors (e.g., correlation edges) with weights (p-values)
	and next plots initial symbol positions (x, y) with weights (fract)
	"""
	time_start = time.time()

	color_map = generate_color_dict(unique_symbols)
	circle_coords = convert_positions_to_circle(sequence_length, circle_radius)
	max_symbolpoint_radius = set_pointcircle_maxradius(cylinder_height, unique_symbols, symbol_points_sizeshrink)

	for p_thresh_max, p_thresh_min in pvalue_thresholds:
		print("pthresh=", p_thresh_min, p_thresh_max)

		generate_single_plot(symbol_position_tracking, symbol_edge_tracking, cylinder_height, \
													unique_symbols, vector_color, color_map, circle_coords, max_symbolpoint_radius, 
													p_thresh_max, p_thresh_min, pkl_plot_folder, vertical_extension, vector_width_binom, radius_mult, \
													pfam_family_id)

	time_end = time.time()
	print("Plotting time for", len(pvalue_thresholds), "p-value thresholds (s)", round(time_end - time_start, 2))

	print("\nEither cntrl-C in command window or close every plot to end program")
	plt.show(block = True)

	return