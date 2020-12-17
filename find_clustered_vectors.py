from get_correlation_matrices import narrow_seqs_to_position_info
from utils import calc_fisher


def narrow_seqs_to_multiposition_info(sequences, position_i, num_combos, seq_por):
	"""
	"""
	downselected_sequences = []
	for seq in sequences:
		candidate_seq_por = seq[position_i : position_i + num_combos]

		if candidate_seq_por == seq_por:
			downselected_sequences.append(seq)

	return downselected_sequences


def generate_seqstring_combos(downselected_sequences, position_i, num_combos):
	"""
	"""
	uniq_seq_combos = {}
	for seq in downselected_sequences:
		seq_por = seq[position_i:position_i + num_combos]

		if seq_por not in uniq_seq_combos:
			uniq_seq_combos[seq_por] = 0
		uniq_seq_combos[seq_por] += 1

	uniq_combo_list = []
	for seq_por in uniq_seq_combos:
		count = uniq_seq_combos[seq_por]
		uniq_combo_list.append([seq_por, count])

	return uniq_combo_list


def gen_multi_contingency_table(sequences, position_i, seq_por_i, position_j, seq_por_j, num_combos):
	"""
	"""
	num_i_not_j = 0
	num_j_not_i = 0
	num_i_j = 0
	num_not_i_not_j = 0

	for seq in sequences:
		seq_string_i = seq[position_i : position_i + num_combos]
		seq_string_j = seq[position_j : position_j + num_combos]

		if seq_string_i == seq_por_i and seq_string_j == seq_por_j:
			num_i_j += 1
		elif seq_string_i == seq_por_i and seq_string_j != seq_por_j:
			num_i_not_j += 1
		elif seq_string_i != seq_por_i and seq_string_j == seq_por_j:
			num_j_not_i += 1
		elif seq_string_i != seq_por_i and seq_string_j != seq_por_j:
			num_not_i_not_j += 1

	cont_table = [[num_not_i_not_j, num_i_not_j], [num_j_not_i, num_i_j]]

	return cont_table


def gen_multi_symbol_edge(multi_symbol_edge_tracking, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, \
														sequences, pvalue_calc_type, num_combos):
	"""
	"""

	if pvalue_calc_type == "Fisher":
		cont_table = gen_multi_contingency_table(sequences, position_i, seq_por_i, position_j, seq_por_j, num_combos)
		pvalue = calc_fisher(cont_table)
	elif pvalue_calc_type == "Binomial":
		cont_table = []
		print("Binomial for multi not done yet")
		raise Exception 

	multi_symbol_edge_tracking.append([position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, pvalue])

	return multi_symbol_edge_tracking


def multipoint_vectors(sequences, unique_symbols, position_symbol_distribution, pvalue_calc_type, num_combos):
	"""
	Calculate p-values of combinations of query symbols (# of combinations to anahyze = num_combos)
	"""
	sequence_length = len(sequences[0])
	total_num_sequences = len(sequences)
	effective_sequence_length = sequence_length - num_combos - 1

	# generate sequences indexed by combo


	multi_symbol_edge_tracking = []
	for position_i in range(effective_sequence_length):
		# At each initial position + num_combos, find symbols that occur at that position
		expressed_symbols_at_position = []

		for symbol_index_i in range(len(unique_symbols)):
			count = position_symbol_distribution[position_i][symbol_index_i]
			if count != 0:
				downselected_sequences = narrow_seqs_to_position_info(sequences, position_i, symbol_index_i, unique_symbols)
				uniq_combo_list = generate_seqstring_combos(downselected_sequences, position_i, num_combos)
				expressed_symbols_at_position.append(uniq_combo_list)


		for uniq_combo_list in expressed_symbols_at_position:
			for seq_por_i, count_i in uniq_combo_list:
				downselected_sequences = narrow_seqs_to_multiposition_info(sequences, position_i, num_combos, seq_por_i)

				# check other positions (if room)
				# preceeding
				if position_i >= num_combos:
					for position_j in range(position_i - num_combos):
						# counts of seq portions at position j from sequences downselected to position_i's seq portion
						position_j_combo_counts = generate_seqstring_combos(downselected_sequences, position_j, num_combos)

						for seq_por_j, count_j in position_j_combo_counts:
							multi_symbol_edge_tracking = gen_multi_symbol_edge(multi_symbol_edge_tracking, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, \
																																			sequences, pvalue_calc_type, num_combos)

				# proceeding
				if effective_sequence_length - position_i - num_combos > 0:
					for other_position in range(0, effective_sequence_length - position_i - num_combos):
						position_j = position_i + other_position + num_combos
						position_j_combo_counts = generate_seqstring_combos(downselected_sequences, position_j, num_combos)

						for seq_por_j, count_j in position_j_combo_counts:
							multi_symbol_edge_tracking = gen_multi_symbol_edge(multi_symbol_edge_tracking, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, \
																																			sequences, pvalue_calc_type, num_combos)

	return multi_symbol_edge_tracking	
