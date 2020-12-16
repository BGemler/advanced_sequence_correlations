import csv
import numpy as np
import matplotlib.pyplot as plt


def analyze_dist_pvalue(pfam_family_id, vector_cor_w_dist_loc, group_similar_aas):
	"""
	"""
	if group_similar_aas == False:
		out_loc = vector_cor_w_dist_loc + pfam_family_id + "-vector-dist_aa.csv"
	elif group_similar_aas == True:
		out_loc = vector_cor_w_dist_loc + pfam_family_id + "-vector-dist_groups.csv"

	same_dict, diff_dict = {}, {}
	any_dict = {}
	with open(out_loc, "r") as f:
		reader = csv.reader(f)
		next(reader, None)

		for row in reader:
			posi, aai, posj, aaj, pvalue, _, dist = row

			if len(dist) > 0:
				posi, posj = int(posi), int(posj)
				position_length = abs(posj - posi)

				pvalue = float(row[4])
				dist = float(row[6])

				if aai == aaj:
					if position_length not in same_dict:
						same_dict[position_length] = []
					same_dict[position_length].append(dist)
				else:
					if position_length not in diff_dict:
						diff_dict[position_length] = []
					diff_dict[position_length].append(dist)	
				if position_length not in any_dict:
					any_dict[position_length] = []
				any_dict[position_length].append(dist)
	f.close()

	same_position_lengths, same_distances = [], []
	diff_position_lengths, diff_distances = [], []
	any_position_lengths, any_distances = [], []

	for position_length in same_dict:
		dists = same_dict[position_length]
		avg_dist = np.average(dists)

		same_position_lengths.append(position_length)
		same_distances.append(avg_dist)

	for position_length in diff_dict:
		dists = diff_dict[position_length]
		avg_dist = np.average(dists)

		diff_position_lengths.append(position_length)
		diff_distances.append(avg_dist)

	for position_length in any_dict:
		dists = any_dict[position_length]
		avg_dist = np.average(dists)

		any_position_lengths.append(position_length)
		any_distances.append(avg_dist)

	fig = plt.figure()
	ax1 = fig.add_subplot(111)

	ax1.scatter(same_position_lengths, same_distances, label = "Symbols Match")
	ax1.scatter(diff_position_lengths, diff_distances, label = "Symbols Don't Match")
	ax1.scatter(any_position_lengths, any_distances, label = "No Symbol Restrictions")
	plt.xlabel("Distance Btwn Positions (bps)")
	plt.ylabel("Distance (Angstroms)")
	plt.title(pfam_family_id + " - Distance Btwn Sequence Positions vs PDB Distance")
	plt.legend()
	plt.savefig("test.png")
	plt.close()


	return