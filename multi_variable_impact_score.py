import csv


def eval_impact_score(count_i, num_seqs, pvalue, l, l_max, l_min, pvalue_max, pvalue_min):
	"""
	"""
	impact_score = (int(count_i) / num_seqs) * ((pvalue_max - pvalue) / (pvalue_max - pvalue_min)) * ((l - l_min) / (l_max - l_min))

	return impact_score


def calc_max_min_pvalue(num_combinations, group_similar_aas, multi_vector_writeout_loc, pfam_family_id):
	"""
	"""
	pvalue_min, pvalue_max = 0.1, 0.1

	for num_combos in num_combinations:
		if group_similar_aas == False:
			out_loc = multi_vector_writeout_loc + pfam_family_id + "-multi_" + str(num_combos) + "-aa.csv"
		elif group_similar_aas == True:
			out_loc = multi_vector_writeout_loc + pfam_family_id + "-multi_" + str(num_combos) + "-groups.csv"

		with open(out_loc, "r") as f:
			reader = csv.reader(f)
			next(reader, None)

			for row in reader:
				position_i, seq_por_i, count_i, position_j, seq_por_j, pvalue = row
				pvalue = float(pvalue)

				if pvalue > pvalue_max:
					pvalue_max = pvalue
				if pvalue < pvalue_min:
					pvalue_min = pvalue
		f.close()

	return pvalue_max, pvalue_min


def calc_impact_score(pfam_family_id, multi_vector_writeout_loc, group_similar_aas, sequences, \
												num_combinations, impact_vector_loc):
	"""
	"""
	impact_vector_all = []

	num_seqs = len(sequences)
	l_max, l_min = max(num_combinations), min(num_combinations)
	pvalue_max, pvalue_min = calc_max_min_pvalue(num_combinations, group_similar_aas, multi_vector_writeout_loc, pfam_family_id)

	for num_combos in num_combinations:
		if group_similar_aas == False:
			out_loc = multi_vector_writeout_loc + pfam_family_id + "-multi_" + str(num_combos) + "-aa.csv"
		elif group_similar_aas == True:
			out_loc = multi_vector_writeout_loc + pfam_family_id + "-multi_" + str(num_combos) + "-groups.csv"

		with open(out_loc, "r") as f:
			reader = csv.reader(f)
			next(reader, None)

			for row in reader:
				position_i, seq_por_i, count_i, position_j, seq_por_j, pvalue = row
				pvalue = float(pvalue)

				impact_score = eval_impact_score(count_i, num_seqs, pvalue, num_combos, l_max, l_min, pvalue_max, pvalue_min)
				impact_vector_all.append([impact_score, position_i, seq_por_i, count_i, position_j, seq_por_j, pvalue])
		f.close()

	impact_vector_all = sorted(impact_vector_all, key = lambda x: x[0], reverse = True)

	if group_similar_aas == False:
		out_loc = impact_vector_loc + pfam_family_id + "impact_scores-aa.tsv"
	elif group_similar_aas == True:
		out_loc = impact_vector_loc + pfam_family_id + "impact_scores-groups.tsv"

	with open(out_loc, "w") as f:
		f.write("Impact Score\tStart Vector Position Begin\tStart Sequence String\tNumber of Sequences w/ Start Sequence at Start Vector" + \
							"\tEnd Vector Position Begin\tEnd Sequence String\tP-Value of Vector\n")

		for impact_score, position_i, seq_por_i, count_i, position_j, seq_por_j, pvalue in impact_vector_all:
			f.write(str(impact_score) + "\t" + str(position_i) + "\t" + seq_por_i + "\t" + str(count_i) + "\t" + \
								str(position_j) + "\t" + seq_por_j + "\t" + str(pvalue) + "\n")
	f.close()

	return