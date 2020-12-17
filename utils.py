from scipy.stats import binom, fisher_exact
import csv


def convert_seqlist_to_dists(sequences, unique_symbols):
	"""
	This function takes any list of sequences and converts them into
	a list of symbol:count for every position in the sequence
	"""
	sequence_length = len(sequences[0])
	total_num_sequences = len(sequences)

	position_symbol_distribution = []
	for position in range(sequence_length):
		symbol_dict = {}
		for seq in sequences:
			symbol = seq[position]

			if symbol != "-":
				if symbol not in symbol_dict:
					symbol_dict[symbol] = 0

				symbol_dict[symbol] += 1

		position_out_row = []
		for symbol in unique_symbols:
			if symbol not in symbol_dict:
				position_out_row.append(0)
			else:
				symbol_count = symbol_dict[symbol]
				#symbol_fract = float(symbol_cout) / float(total_num_sequences)
				position_out_row.append(symbol_count)

		position_symbol_distribution.append(position_out_row)

	return position_symbol_distribution


def calc_cdf(x, n, p):
	"""
	when x > n *p:
    use survival function (1 - cdf)
    subtract 1 from x to account for differences in definition of
    cdf and significance	
	"""
	if x > n * p:
		return binom.sf(x - 1, n, p)
	else:
		return binom.cdf(x, n, p)


def calc_fisher(table):
	"""
	scipy fisher_exact function
	assuming two-sided alternative hypothesis 
	"""
	_, pvalue = fisher_exact(table, alternative = 'greater')

	return pvalue


def generate_pvalue(num_at_proceeding_pos, count, position_symbol_distribution, position_index, position_symbol_index, \
										total_num_sequences, cont_table, pvalue_calc_type):
	"""
	This is the function where a p-value is generated
	Call scipy binomial distribution to make this estimation
	or fisher exact probability test
	"""
	# scipy binom (x, n, p)
	# x = number of observed target when downselected
	# n = number of downselected
	# p = expected probability of target (population level)

	if pvalue_calc_type == "Binomial":
		population_num_target = position_symbol_distribution[position_index][position_symbol_index]
		population_target_fract = float(population_num_target) / float(total_num_sequences)

		pvalue = calc_cdf(num_at_proceeding_pos, count, population_target_fract)

	elif pvalue_calc_type == "Fisher":
		pvalue = calc_fisher(cont_table)

	return pvalue


def write_out_vectors(pfam_family_id, vector_writeout_loc, symbol_edge_tracking, unique_symbols, group_similar_aas):
	"""
	This function writes out all of the vectors and their p-values into a CSV.
	Symbols and positions for the beginning and end of each vector are written out
	"""
	if group_similar_aas == False:
		out_loc = vector_writeout_loc + pfam_family_id + "-vector_params_aa.csv"
	elif group_similar_aas == True:
		out_loc = vector_writeout_loc + pfam_family_id + "-vector_params_groups.csv"

	with open(out_loc, "w", newline='') as f:
		out = csv.writer(f)
		out.writerow(["Starting Position on Sequence", "Starting Symbol", \
									"Ending Position on Sequence", "Ending Symbol", \
									"P-Value of Vector"])

		for position_initial, symbol_index_initial, position_final, symbol_index_final, pvalue in symbol_edge_tracking:
			initial_symbol = unique_symbols[symbol_index_initial]
			final_symbol = unique_symbols[symbol_index_final]

			out.writerow([position_initial, initial_symbol, position_final, final_symbol, pvalue])
	f.close()

	return

def write_out_multi_vectors(pfam_family_id, multi_vector_writeout_loc, multi_symbol_edge_tracking, num_combos, group_similar_aas, \
																sequences, minimum_count_fract):
	"""
	"""
	if group_similar_aas == False:
		out_loc = multi_vector_writeout_loc + pfam_family_id + "-multi_" + str(num_combos) + "-aa.csv"
	elif group_similar_aas == True:
		out_loc = multi_vector_writeout_loc + pfam_family_id + "-multi_" + str(num_combos) + "-groups.csv"

	num_seqs = len(sequences)

	with open(out_loc, "w", newline='') as f:
		out = csv.writer(f)
		out.writerow(["Start Vector Position Begin", "Start Sequence String", \
										"Number of Sequences w/ Start Sequence at Start Vector", \
										"End Vector Position Begin", "End Sequence String", \
										"P-Value of Vector"])

		for position_i, seq_por_i, count_i, position_j, seq_por_j, _, pvalue in multi_symbol_edge_tracking:
			#if count_i >= int(num_seqs * minimum_count_fract):
			if len(seq_por_i.replace(".","")) > 0.5 * len(seq_por_i) and len(seq_por_j.replace(".","")) > 0.5 * len(seq_por_j):
				out.writerow([position_i, seq_por_i, count_i, position_j, seq_por_j, pvalue])
	f.close()

	return