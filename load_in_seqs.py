
# charged: C
# aromatic: A
# Hydrophillic: L
# Very hydrophobic: V
# Hydrophobic: H
aa_mapping = {
	"D": "C",
	"E": "C", 
	"K": "C", 
	"R": "C", 
	"W": "A",
	"Y": "A", 
	"F": "A", 
	"N": "L", 
	"Q": "L", 
	"G": "L", 
	"S": "L", 
	"V": "V", 
	"C": "V", 
	"I": "V", 
	"L": "V", 
	"M": "V", 
	"P": "V", 
	"A": "H", 
	"H": "H", 
	"T": "H", 
	".": "."
}


def load_fasta(pfam_family_id, sequence_location, put_seqs_in_frame, group_similar_aas):
	"""
	Loads in a FASTA file
	If required, puts all sequences in the same frame
	Returns list of sequences all in the same frame
	"""
	fasta_loc = sequence_location + pfam_family_id + "-seqs.fa"

	seqid_dict = {}
	with open(fasta_loc, "r") as f:
		subj, seq = "", ""
		for row in f:
			row = row.replace("\n","").replace("\r","")
			if row.startswith(">"):
				if len(subj) != 0:
					seqid_dict[subj] = seq
				subj = row[1:]
				seq = ""
			else:
				seq = seq + row.upper()
		seqid_dict[subj] = seq
	f.close()
	
	if put_seqs_in_frame == True:
		print("Sequence frame adjustment not yet available")
		raise Exception

	elif put_seqs_in_frame == False:
		sequences = []
		for seqid in seqid_dict:
			seq = seqid_dict[seqid]
			
			if group_similar_aas == False:
				sequences.append(seq)
			elif group_similar_aas == True:
				group_seq = ""
				for s in seq:
					group = aa_mapping[s]
					group_seq = group_seq + group
				sequences.append(group_seq)

	sequence_length = len(sequences[0])

	return sequences, sequence_length