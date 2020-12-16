import numpy as np
"""
Contains variables, paths, etc. necessary for code execution
"""

pfam_family_id = "PF00626__23"

# pipeline choices
execute_correlations = False
generate_plots = False
execute_pdb_analysis = True

put_seqs_in_frame = False
pvalue_calc_type = "Fisher" # Binomial, Fisher


# P-Value thresholds for plotting
pvalue_thresholds = [
	[0.000001, 0.0], \
	[0.0, 0.0]
	]

# Paths
sequence_location = "seq_gather/qual_pfam/"
structure_locations = "seq_gather/pdb_files/"
pkl_plot_folder = "pkl_plot_dumps/"
vector_writeout_loc = "results/vector_results/"
vector_cor_w_dist_loc = "results/vector_cor_wdist/"

# Key Variables

circle_radius = 300.0
cylinder_height = 150.0
vertical_extension = 1.2 # factor for extending the z-axes above the cylinder

vector_color = "white"
vector_width_binom = 1 # width ^ x <-- e.g., lower x = larger vectors at p_value = 1

symbol_points_sizeshrink = 1 # >1 = smaller circles, <1 = bigger circles for symbols
radius_mult = 100 # scaling for scatter point size
