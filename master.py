from config import *
from load_in_seqs import load_fasta
from get_position_distributions import get_position_symbol_dists, get_position_symbols
from get_correlation_matrices import get_correlation_matrix
from plotting import generate_vector_plots
from utils import write_out_vectors
from cor_to_dists import get_atom_dists
from analyze_pvalue_dist_cor import analyze_p_dist
from analyze_pvalue_dist_symbols import analyze_dist_pvalue
import time
 
def master_pipeline(group_similar_aas):
	"""
	"""
	if execute_correlations == True:
		time_start = time.time()
		print("Load in a FASTA file of test sequences")
		sequences, sequence_length = load_fasta(pfam_family_id, sequence_location, put_seqs_in_frame, group_similar_aas)

		print("Finding Symbol Distribution at Each Position")
		unique_symbols, position_symbol_distribution = get_position_symbol_dists(sequences)

		print("Calculating Correlation Matrices")
		# This is a list of vectors - i.e., initial (x,z) and end (x,z) with a weight (p-value)
		symbol_edge_tracking = get_correlation_matrix(sequences, unique_symbols, position_symbol_distribution, pvalue_calc_type)

		print("Generating Initial Symbol Positions")
		# This is a list of symbols, also the starting point weights of the vectors - i.e., initial (x,z) and weight (fraction symbol at position)
		symbol_position_tracking = get_position_symbols(len(sequences), position_symbol_distribution, unique_symbols)

		print("Writing out vectors")
		write_out_vectors(pfam_family_id, vector_writeout_loc, symbol_edge_tracking, unique_symbols, group_similar_aas)

		time_end = time.time()
		print("Computation time for processing", len(sequences), "sequences (s)", round(time_end - time_start, 2))

	if generate_plots == True:
		print("Plotting Results")
		generate_vector_plots(sequence_length, symbol_position_tracking, symbol_edge_tracking, circle_radius, cylinder_height, \
														unique_symbols, vector_color, symbol_points_sizeshrink, pkl_plot_folder, pvalue_thresholds, \
														vertical_extension, vector_width_binom, radius_mult, \
														pfam_family_id)
	
	if execute_pdb_analysis == True:
		print("Determining atom distances from PDB files")
		#get_atom_dists(structure_locations, sequence_location, pfam_family_id, vector_cor_w_dist_loc, vector_writeout_loc, group_similar_aas)
		analyze_p_dist(pfam_family_id, vector_cor_w_dist_loc, group_similar_aas)
		#analyze_dist_pvalue(pfam_family_id, vector_cor_w_dist_loc, group_similar_aas)

	return


group_similar_aas_list = [False, True]
for group_similar_aas in group_similar_aas_list:
	master_pipeline(group_similar_aas)
