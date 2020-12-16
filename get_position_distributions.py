from utils import convert_seqlist_to_dists


def get_position_symbol_dists(sequences):
	"""
	This function takes in a list of equal length sequences
	and calculates the distribution of unique symbols at each position.
	
	Returns list where each line is a sequence position starting with the 
	first position in the sequnece. Each line contains the number of 
	each symbols expressed at that position (sans gaps)

	Further, returns ordered list of unique symbols in the sequences 
	which will be used to order all symbol lists. Currently, the ordered
	symbol list is generated arbitrarily (could add option for user-input)
	in the future
	"""
	unique_symbols = []
	for seq in sequences:
		for s in seq:
			if s not in unique_symbols and s is not "-":
				unique_symbols.append(s)

	position_symbol_distribution = convert_seqlist_to_dists(sequences, unique_symbols)

	return unique_symbols, position_symbol_distribution


def get_position_symbols(tot_num_seqs, position_symbol_distribution, unique_symbols):
	"""
	This function generates a list of symbol fractions at each position
	Returns a list of non-zero fractions (e.g., the size of each initial point
	of each vector)
	"""
	symbol_position_tracking = []
	for position in range(len(position_symbol_distribution)):
		for symbol_index in range(len(unique_symbols)):
			count = position_symbol_distribution[position][symbol_index]

			if count > 0:
				fract = float(count) / float(tot_num_seqs)

				# artifically make p-values
				symbol_position_tracking.append([position, symbol_index, fract])

	return symbol_position_tracking
