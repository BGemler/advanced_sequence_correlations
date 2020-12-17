import numpy as np
"""
Contains variables, paths, etc. necessary for code execution
"""

pfam_family_id = "PF00626__23"

# pipeline choices
execute_correlations = False
execute_multi_corr = True
generate_plots = False
execute_pdb_analysis = False

put_seqs_in_frame = False
pvalue_calc_type = "Fisher" # Binomial, Fisher
minimum_count_fract = 0.1 # note: some low counts probably aren't that meaningful
num_combinations = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] # list of num_combos

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
impact_vector_loc = "results/multi_vector_results/"
multi_vector_writeout_loc = "results/multi_vector_results/bulk_csvs/"

# Key Variables

circle_radius = 300.0
cylinder_height = 150.0
vertical_extension = 1.2 # factor for extending the z-axes above the cylinder

vector_color = "white"
vector_width_binom = 1 # width ^ x <-- e.g., lower x = larger vectors at p_value = 1

symbol_points_sizeshrink = 1 # >1 = smaller circles, <1 = bigger circles for symbols
radius_mult = 100 # scaling for scatter point size
