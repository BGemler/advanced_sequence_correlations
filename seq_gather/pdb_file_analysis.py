import os
import csv


def convert_list_to_out_str(list_in):
	"""
	"""
	str_out = ""
	for l in list_in:
		str_out = str_out + l + ";"

	if len(str_out) > 0:
		str_out = str_out[:-1]

	return str_out


def parse_pdb_file(pdb_file):
	"""
	https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html
	"""
	pdb_id = pdb_file.split("/")[1].split(".pdb")[0]
	print(pdb_id)

	header, title, organism, taxid = "", "", "", ""
	contains_sheet, contains_helix = False, False
	sheet_ids, sheet_aa_positions = [], []
	helix_ids, helix_aa_positions = [], []

	with open(pdb_file, "r") as f:
		for row in f:
			if row.startswith("HEADER"):
				entry = row.split("HEADER")[1]
				for i in range(len(entry)):
					if entry[i:i+2] != "  ":
						header = header + entry[i:i+1]

			elif row.startswith("TITLE"):
				entry = row.split("TITLE")[1]
				for i in range(len(entry)):
					if entry[i:i+2] != "  ":
						title = title + entry[i:i+1]

			elif "ORGANISM_SCIENTIFIC" in row:
				organism = row.split("ORGANISM_SCIENTIFIC: ")[1]

			elif "ORGANISM_TAXID" in row:
				taxid = row.split("ORGANISM_TAXID: ")[1]

			elif row.startswith("SHEET"):
				contains_sheet = True

				sheet_id = row[11:11 + 3]
				init_residue = row[22 : 22 + 4].replace(" ","")
				term_residue = row[33 : 33 + 4].replace(" ","")

				sheet_ids.append(sheet_id)
				sheet_aa_positions.append(init_residue + "-" + term_residue)

			elif row.startswith("HELIX"):
				contains_helix = True

				helix_id = row[11:11 + 3]
				init_residue = row[21 : 21 + 4].replace(" ","")
				term_residue = row[33 : 33 + 4].replace(" ","")

				helix_ids.append(helix_id)
				helix_aa_positions.append(init_residue + "-" + term_residue)
	f.close()

	out_data_row = [
		pdb_id, header, title, organism, taxid, \
		contains_sheet, len(set(sheet_ids)), convert_list_to_out_str(sheet_aa_positions), \
		contains_helix, len(set(helix_ids)), convert_list_to_out_str(helix_aa_positions)
	]

	return out_data_row


def process_pdbs(pdb_file_bulk_loc, pdb_analysis_output_loc):
	"""
	"""
	with open(pdb_analysis_output_loc, "w", newline = '') as f:
		out = csv.writer(f)
		out.writerow(["PDB ID", "Header", "Title", "Organism(s)", "TaxID(s)", \
										"Contains Sheet?", "Number of Sheets", "Initial & Terminal Residues of Sheet Strands", \
										"Contains Helix?", "Number of Helices", "Initial & Terminal Residues of Helices"])

		for filename in os.listdir(pdb_file_bulk_loc):
			pdb_data_row = parse_pdb_file(pdb_file_bulk_loc + filename)
			out.writerow(pdb_data_row)
	f.close()

	return


pdb_file_bulk_loc = "pdb_files/"
pdb_analysis_output_loc = "pdb_analysis_results/pdb_summary.csv"

process_pdbs(pdb_file_bulk_loc, pdb_analysis_output_loc)