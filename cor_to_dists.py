import csv
import math
import numpy as np
from load_in_seqs import aa_mapping


aa_triplet_tosymbol_map = {
	"GLY": "G",
	"ALA": "A",
	"LEU": "L",
	"MET": "M",
	"PHE": "F",
	"TRP": "W",
	"LYS": "K",
	"GLN": "Q",
	"GLU": "E",
	"SER": "S",
	"PRO": "P",
	"VAL": "V",
	"ILE": "I",
	"CYS": "C",
	"TYR": "Y",
	"HIS": "H",
	"ARG": "R",
	"ASN": "N",
	"ASP": "D",
	"THR": "T",
	".": "."
}


def calc_minimum_dist_twopos(position_coords_i, position_coords_j):
	"""
	"""
	minimum_dist = ""

	for _, x, y, z in position_coords_i:
		for _, xj, yj, zj in position_coords_j:
			dist = math.sqrt(((x - xj) ** 2) + \
													((y - yj) ** 2) + \
													((z - zj) ** 2))

			if minimum_dist == "":
				minimum_dist = dist

			if dist <= minimum_dist:
				minimum_dist = dist

	return minimum_dist


def get_position_dists(pdf_file_loc, position_range):
	"""
	For the positions in position_range, calculate smallest distance matrix
	"""
	# dict where key = position, value = list of [aa, x, y, z] coords of atoms at position
	position_coords = {}
	with open(pdf_file_loc, "r") as f:
		for row in f:
			aa = ""

			if row.startswith("ATOM   "):
				# Split by 5 spaces - very stupid .pdb file organization, if string >= 1000, splitting by spaces breaks row
				row = row.split("ATOM")[1].strip()
				row = row.split("     ")
				first_part, second_part = row[0], row[1]

				first_part = first_part.strip()
				first_part = ' '.join(first_part.split()) # collapse all repeating regions of spaces into 1 space
				first_part = first_part.split()

				position = int(first_part[0])
				aa_trip = first_part[2]

				if len(aa_trip) == 4:
					aa_trip = aa_trip[1:]

				if aa_trip not in aa_triplet_tosymbol_map:
					aa = ""
				else:
					aa = aa_triplet_tosymbol_map[aa_trip]

				second_part = second_part.strip()
				second_part = ' '.join(second_part.split()) # collapse all repeating regions of spaces into 1 space
				second_part = second_part.split()

				if len(second_part) == 5 and len(aa) > 0:
					x = float(second_part[0])
					y = float(second_part[1])
					z = float(second_part[2])

					# see if position is in the position range
					if position >= position_range[0] and position <= position_range[1]:
						if position not in position_coords:
							position_coords[position] = []
						position_coords[position].append([aa, x, y, z])


				"""
				row = ' '.join(row.split()) # collapse all repeating regions of spaces into 1 space
				row = row.split()

				if len(row) != 11:
					print(row)
					raise Exception
				if len(row) == 11:
					position = int(row[1])

					aa_trip = row[3]
					aa = aa_triplet_tosymbol_map[aa_trip]

					x = float(row[6])
					y = float(row[7])
					z = float(row[8])
					"""
	f.close()

	aa_seq = ""
	for position in position_coords:
		aa = position_coords[position][0][0]
		aa_seq = aa_seq + aa

	if len(aa_seq) != position_range[1] - position_range[0] + 1:
		# signifies not every required position could be read in the .PDB file (corrupted rows)
		return "", []

	# intialize position matrix - we want lowest dist of position(i) to position(i +/- n)
	position_matrix = []
	for i in range(position_range[1] - position_range[0] + 1):
		position_row = []
		for j in range(position_range[1] - position_range[0] + 1):
			position_row.append("")
		position_matrix.append(position_row)

	for i in range(position_range[1] - position_range[0] + 1):
		for j in range(position_range[1] - position_range[0] + 1):
			position_i = position_range[0] + i
			position_j = position_range[0] + j

			minimum_dist = calc_minimum_dist_twopos(position_coords[position_i], position_coords[position_j])

			position_matrix[i][j] = minimum_dist

	return aa_seq, position_matrix


def get_atom_dists(structure_locations, sequence_location, pfam_family_id, vector_cor_w_dist_loc, vector_writeout_loc, group_similar_aas):
	"""
	"""
	seqid_dict = {}
	with open(sequence_location + pfam_family_id + "-seqs.fa", "r") as f:
		subj, seq = "", ""
		for row in f:
			row = row.replace("\n","")
			if row.startswith(">"):
				if len(subj) != 0:
					seqid_dict[subj] = seq
				subj = row[1:]
				seq = ""
			else:
				seq = seq + row
		seqid_dict[subj] = seq
	f.close()

	structure_dict = []
	with open(sequence_location + pfam_family_id + "-pdb_structs.tsv", "r") as f:
		for row in f:
			row = row.replace("\n","").split("\t")
			seqid, pdb_info = row

			seq = seqid_dict[seqid]
			_, pdf_id, positions, _ = pdb_info.split(";")
			
			pdf_id = pdf_id.strip().split()[0] # strip letter code
			positions = positions.strip().split("-") 
			positions = [int(x) for x in positions] # int(pos_start), int(pos_end)

			structure_dict.append([pdf_id, positions, seq])
	f.close()

	structure_distance_dict = {}
	total_count, good_count = 0, 0
	already_done_pdbs = []
	for pdf_id, positions, norm_seq in structure_dict:
		pdb_seq, position_matrix = get_position_dists(structure_locations + pdf_id + ".pdb", positions)

		# create a dict mapping pdb_seq positions to norm_seq positions
		pdb_pos_to_norm_pos = {}
		norm_count, count = 0, 0
		for i in range(len(norm_seq)):
			residue = norm_seq[i:i+1]
			if residue != ".":
				pdb_pos_to_norm_pos[count] = norm_count
				count += 1
			norm_count += 1

		proceed = True
		if len(position_matrix) == 0:
			# indicates PDB file is incorrectly formatted
			#print(pdf_id, "'s PDB file data rows are corrupted - skipping")
			proceed = False
		elif len(norm_seq.replace(".","")) != len(pdb_seq):
			# indicates error with positional information cataloged in PFAM
			#print(pdf_id, "'s PFAM file has incorrect position or normalized seq")
			proceed = False
		if pdf_id in already_done_pdbs:
			proceed = False
		else:
			already_done_pdbs.append(pdf_id)

		total_count += 1
		if proceed == True:
			good_count += 1
			# parse over position_matrix, adding position-position & AA-AA maps to the structure_distance_dict
			for i in range(len(pdb_seq)):
				for j in range(len(pdb_seq)):

					# Distances don't care about direction - only grab 1 side of the matrix mirror
					#if i <= j:
					pdb_aa_i = pdb_seq[i:i+1]
					pdb_aa_j = pdb_seq[j:j+1]

					if group_similar_aas == True:
						pdb_aa_i = aa_mapping[pdb_aa_i]
						pdb_aa_j = aa_mapping[pdb_aa_j]

					norm_position_i = pdb_pos_to_norm_pos[i]
					norm_position_j = pdb_pos_to_norm_pos[j]

					dist = position_matrix[i][j]

					dict_key = str(norm_position_i) + "|" + str(norm_position_j) + "|" + \
											pdb_aa_i + "|" + pdb_aa_j

					if dict_key not in structure_distance_dict:
						structure_distance_dict[dict_key] = []
					structure_distance_dict[dict_key].append(dist)

	# find average distance at each dict_key
	for dict_key in structure_distance_dict:
		dists = structure_distance_dict[dict_key]
		avg_dist = np.average(dists)

		structure_distance_dict[dict_key] = [len(dists), avg_dist]

	# load in p-value matrix, write out average distance
	if group_similar_aas == False:
		out_loc = vector_cor_w_dist_loc + pfam_family_id + "-vector-dist_aa.csv"
		in_loc = vector_writeout_loc + pfam_family_id + "-vector_params_aa.csv"
	elif group_similar_aas == True:
		out_loc = vector_cor_w_dist_loc + pfam_family_id + "-vector-dist_groups.csv"
		in_loc = vector_writeout_loc + pfam_family_id + "-vector_params_groups.csv"

	with open(out_loc, "w", newline = '') as f:
		out = csv.writer(f)

		with open(in_loc, "r") as r:
			reader = csv.reader(r)

			row_count = 0
			for row in reader:
				if len(row) == 5:
					row_count += 1
					if row_count == 1:
						header_row = row + ["Number of PDB Structures with Matching AAs at Start/End Positions", \
																	"Average Minimum Distance Across All PDB Structures (Angstroms)"]
						out.writerow(header_row)
					else:
						position_i, aa_i, position_j, aa_j, _ = row

						dict_key = position_i + "|" + position_j + "|" + aa_i + "|" + aa_j
						if dict_key in structure_distance_dict:
							num_pdbs, avg_dist = structure_distance_dict[dict_key]
						else:
							num_pdbs = 0
							avg_dist = ""

						if aa_i != "." and aa_j != ".":
							out.writerow(row + [num_pdbs, avg_dist])
		r.close()
	f.close()

	print("Wrote out distances for:", pfam_family_id, " group similar aas:", group_similar_aas)
	print("Total # of PDB File References:", total_count)
	print("# of Good PDB File References:", good_count)

	return