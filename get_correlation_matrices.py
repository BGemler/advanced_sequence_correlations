from utils import convert_seqlist_to_dists, generate_pvalue


def narrow_seqs_to_position_info(sequences, position, symbol_index, unique_symbols):
	"""
	This function finds all sequences without a list that contain a certain 
	symbol at a certain position
	"""
	downselected_sequences = []

	symbol_of_interest = unique_symbols[symbol_index]

	for seq in sequences:
		seq_symbol = seq[position]

		if seq_symbol == symbol_of_interest:
			downselected_sequences.append(seq)

	return downselected_sequences


def generate_contingency_table(sequences, position_i, symbol_index_i, position_j, symbol_index_j, unique_symbols):
	"""
	"""
	symbol_i = unique_symbols[symbol_index_i]
	symbol_j = unique_symbols[symbol_index_j]

	num_i_not_j = 0
	num_j_not_i = 0
	num_i_j = 0
	num_not_i_not_j = 0

	for seq in sequences:
		seq_symbol_i = seq[position_i]
		seq_symbol_j = seq[position_j]

		if seq_symbol_i == symbol_i and seq_symbol_j == symbol_j:
			num_i_j += 1
		elif seq_symbol_i == symbol_i and seq_symbol_j != symbol_j:
			num_i_not_j += 1
		elif seq_symbol_i != symbol_i and seq_symbol_j == symbol_j:
			num_j_not_i += 1
		elif seq_symbol_i != symbol_i and seq_symbol_j != symbol_j:
			num_not_i_not_j += 1

	cont_table = [[num_not_i_not_j, num_i_not_j], [num_j_not_i, num_i_j]]

	return cont_table



def generate_symbol_edge_entry(symbol_edge_tracking, symbol_distribution, count, position, symbol_index, position_index, total_num_sequences, \
																	position_symbol_distribution, sequences, unique_symbols, pvalue_calc_type):
	"""
	For every symbol-symbol relationship, calculate the vector's p-value and add to the 
	vector generation list (symbol_edge_tracking) to be plotted
	"""
	for position_symbol_index in range(len(symbol_distribution)):
		if pvalue_calc_type == "Fisher":
			cont_table = generate_contingency_table(sequences, position, symbol_index, position_index, position_symbol_index, unique_symbols)
		elif pvalue_calc_type == "Binomial":
			cont_table = []
		else:
			print("invalid p-value calculation method")
			raise Exception

		num_at_proceeding_pos = symbol_distribution[position_symbol_index]

		if num_at_proceeding_pos > 0.0:
			pvalue = generate_pvalue(num_at_proceeding_pos, count, position_symbol_distribution, position_index, position_symbol_index, \
																	total_num_sequences, cont_table, pvalue_calc_type)

			symbol_edge_tracking.append([position, symbol_index, position_index, position_symbol_index, pvalue])

	return symbol_edge_tracking


def get_correlation_matrix(sequences, unique_symbols, position_symbol_distribution, pvalue_calc_type):
	"""
	For each position, calculates the probability for every proceeding position and unique
	symbol at all other position of encountering another symbol.
	Note: both ends of the vector matter, e.g.,: "What is the probability of having a G at position 2 
	when I have a A at position 1" is a diff question than the converse.
	Only returns non-0 probabilities (e.g, the list is a list of vectors, with a p-value)
	"""
	sequence_length = len(sequences[0])
	total_num_sequences = len(sequences)

	symbol_edge_tracking = []
	for position in range(sequence_length):
		# At each initial position, find symbols that occur at that position
		expressed_symbols_at_position = []
		for i in range(len(unique_symbols)):
			count = position_symbol_distribution[position][i]
			if count != 0:
				expressed_symbols_at_position.append([i, count])

		# For each unique symbol at the initial position, find distribution of other symbols at ALL OTHER positions
		for symbol_index, count in expressed_symbols_at_position:
			downselected_sequences = narrow_seqs_to_position_info(sequences, position, symbol_index, unique_symbols)
			downselected_position_symbol_distribution = convert_seqlist_to_dists(downselected_sequences, unique_symbols)

			# proceding positions
			for other_position in range(1, sequence_length - position):
				position_index = other_position + position
				symbol_distribution = downselected_position_symbol_distribution[position_index]
				symbol_edge_tracking = generate_symbol_edge_entry(symbol_edge_tracking, symbol_distribution, count, position, symbol_index, position_index, \
																														total_num_sequences, position_symbol_distribution, sequences, unique_symbols, pvalue_calc_type)
				
			# preceeding positions
			for other_position in range(position):
				position_index = other_position
				symbol_distribution = downselected_position_symbol_distribution[other_position]
				symbol_edge_tracking = generate_symbol_edge_entry(symbol_edge_tracking, symbol_distribution, count, position, symbol_index, position_index, \
																														total_num_sequences, position_symbol_distribution, sequences, unique_symbols, pvalue_calc_type)

	return symbol_edge_tracking
