import matplotlib.pyplot as plt
import os
from urllib.request import urlretrieve
import time


def write_out_family(qualifying_pfam_outfolder, acc, seq_out_rows, structure_out_rows):
	"""
	"""
	# write out seq overlap regions
	with open(qualifying_pfam_outfolder + acc.replace(".","__") + "-seqs.fa", "w") as f:
		for seqid, seq_por in seq_out_rows:
			f.write(">" + seqid + "\n")
			f.write(seq_por + "\n")
	f.close()

	# write out PDB structures
	with open(qualifying_pfam_outfolder + acc.replace(".","__") + "-pdb_structs.tsv", "w") as f:
		for seqid, structure_id in structure_out_rows:
			f.write(seqid + "\t" + structure_id + "\n")
	f.close()

	return


def parse_pfam_file(pfam_seq_dump_loc, pfam_seqnum_thresh, qualifying_pfam_outfolder):
	"""
	We want PFAM seq families that have # of seqs > pfam_seqnum_thresh
	and at least 1 seq has a PDB structure ID associated with it
	"""
	all_pfam_fam_lens = []
	with open(pfam_seq_dump_loc, "r", encoding='utf-8', errors='replace') as f:
		acc, num_seqs,  = "", 0
		seq_out_rows, structure_out_rows = [], []

		for row in f:
			row = row.replace("\n","")

			if row.startswith("#=GF AC"):
				acc = row.split("   ")[1]

			elif row.startswith("#=GF SQ"):
				num_seqs = int(row.split("   ")[1])
				all_pfam_fam_lens.append(num_seqs)

			elif row.startswith("#=GS") and " DR " in row:
				structure_row = row.split("#=GS ")[1]
				structure_row = ' '.join(structure_row.split())# collapse all repeating regions of spaces into 1 space
				structure_row = structure_row.split()

				seq_id = structure_row[0]
				structure_id = ""
				for i in range(len(structure_row) - 2):
					v = structure_row[i + 2]
					structure_id = structure_id + v + " "
				structure_id = structure_id[:-1]

				if "PDB" in structure_id:
					structure_out_rows.append([seq_id, structure_id])

			elif row.startswith("//"):
				if num_seqs >= pfam_seqnum_thresh and len(structure_out_rows) > 0:
					write_out_family(qualifying_pfam_outfolder, acc, seq_out_rows, structure_out_rows)
				acc, num_seqs = "", 0
				seq_out_rows, structure_out_rows = [], []

			elif row.startswith("# STOCKHOLM 1.0"):
				version = "# STOCKHOLM 1.0"

			elif row[:2] != "#=":
				# seq align regions
				seq_row = ' '.join(row.split()) # collapse all repeating regions of spaces into 1 space
				seq_id, seq_out_row = seq_row.split() 
				seq_out_rows.append([seq_id, seq_out_row])
	f.close()

	print("total # of PFAM families: ", len(all_pfam_fam_lens))
	print("total # of PFAM seqs across all families: ", sum(all_pfam_fam_lens))

	return all_pfam_fam_lens


def gen_hist(all_pfam_fam_lens, pfam_seqnum_thresh):
	"""
	"""
	list_size = len(all_pfam_fam_lens)
	narrowed_list_size = len([x for x in all_pfam_fam_lens if x >= pfam_seqnum_thresh])

	plt.hist(all_pfam_fam_lens, bins = range(min(all_pfam_fam_lens), max(all_pfam_fam_lens) + pfam_seqnum_thresh, pfam_seqnum_thresh))
	plt.xlabel("# of Family Seed Sequences Bin")
	plt.ylabel("Number of PFAM Families in Bin")
	plt.title("PFAM Family Seed Sequence Distribution" + "\n" + \
							"Fraction of Families with > " + str(pfam_seqnum_thresh) + " Sequences: " + str(round(narrowed_list_size / list_size, 2)))
	plt.savefig("results/pfam-seq-hist.png")
	plt.close()

	plt.boxplot(all_pfam_fam_lens)
	plt.ylabel("Num of PFAM Seed Sequences Per Family")
	plt.title("Box & Whisker Plot: PFAM Family Seed Sequence Distribution" + "\n" + \
							"Fraction of Families with > " + str(pfam_seqnum_thresh) + " Sequences: " + str(round(narrowed_list_size / list_size, 2)))
	plt.savefig("results/pfam-seq-bw.png")
	plt.close()

	return


def grab_pdb_files(qualifying_pfam_outfolder, pdb_file_folder):
	"""
	"""
	pdb_dl_root = "https://files.rcsb.org/download/"

	all_pdb_ids = set()
	for filename in os.listdir(qualifying_pfam_outfolder):
		if "structs.tsv" in filename:
			with open(qualifying_pfam_outfolder + filename, "r") as f:
				for row in f:
					row = row.replace("\n","").split("\t")
					_, pdb_row = row
					pdb_id = pdb_row.split()[1]
					all_pdb_ids.update(set([pdb_id]))
			f.close()

	pdb_ids_not_work = []
	for pdb_id in all_pdb_ids:
		#time.sleep(1)
		url = pdb_dl_root + pdb_id + ".pdb"
		out_loc = pdb_file_folder + pdb_id + ".pdb"
		print(url)
		try:
			urlretrieve(url, out_loc)
		except:
			pdb_ids_not_work.append(url)

	print("# of unique PDB IDs:", len(all_pdb_ids))
	print("# of PDB IDs with timeout:", len(pdb_ids_not_work))

	with open("pdb_not_work.out", "w") as f:
		for pdb_id in pdb_ids_not_work:
			f.write(pdb_id + "\n")
	f.close()

	return


def main(pfam_seq_dump_loc, pfam_seqnum_thresh, qualifying_pfam_outfolder, pdb_file_folder):
	"""
	"""
	all_pfam_fam_lens = parse_pfam_file(pfam_seq_dump_loc, pfam_seqnum_thresh, qualifying_pfam_outfolder)
	#gen_hist(all_pfam_fam_lens, pfam_seqnum_thresh)

	grab_pdb_files(qualifying_pfam_outfolder, pdb_file_folder)

	return


pfam_seq_dump_loc = "resources/Pfam-A.seed" #ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
pfam_seqnum_thresh = 10
qualifying_pfam_outfolder = "qual_pfam/"
pdb_file_folder = "pdb_files/"

main(pfam_seq_dump_loc, pfam_seqnum_thresh, qualifying_pfam_outfolder, pdb_file_folder)