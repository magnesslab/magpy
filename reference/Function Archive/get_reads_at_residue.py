# Find the sequence of all reads covering a given residue
# Outputs a csv file with one line for each barcode in the format: cell_barcode,(residue,quality),(residue,quality)

import pysam

#KRAS-G12D @ chromosome 6, residue 145246770

def get_reads_at_residue(chromosome,residue,read_file,write_file):
	with pysam.AlignmentFile(read_file, "rb") as bamfile:
		output_dict = {}
		#Get all reads that align to the desired region
		for read in bamfile.fetch(chromosome, residue, residue+1): 
			#Get a list of all residues in the aligned read
			aligned_residues = read.get_reference_positions()

			#If our residue is present, get all the relevant information
			if residue in aligned_residues:
				cell_barcode = read.get_tag('CR')
				sequence = read.query_alignment_sequence
				qualities = read.query_alignment_qualities

				i = aligned_residues.index(residue)
				res12_seq = sequence[i]
				res12_qual = qualities[i]

				if cell_barcode in output_dict.keys():
					output_dict[cell_barcode].append((res12_seq,res12_qual))
				else: 
					output_dict[cell_barcode] = [(res12_seq,res12_qual)]

		with open(write_file,"w") as file:
			for cell_barcode in output_dict.keys():
				line = cell_barcode
				for read in output_dict[cell_barcode]:
					line += f",{read[0]}:{read[1]}"
				line += "\n"
				file.write(line)

if __name__ == "__main__":
	file_path = "/proj/smngsslb/users/jbliton/mist1kras/Mist1-Kras-scRNAseq"
	get_reads_at_residue('6', 145246770, file_path+"/1mo-Tam/11_ctl/outs/possorted_genome_bam.bam", "kras_1mo_ctl_counts.csv")
	get_reads_at_residue('6', 145246770, file_path+"/1mo-Tam/11_exp/outs/possorted_genome_bam.bam", "kras_1mo_exp_counts.csv")
	get_reads_at_residue('6', 145246770, file_path+"/4mo-Tam/10_ctl/outs/possorted_genome_bam.bam", "kras_4mo_ctl_counts.csv")
	get_reads_at_residue('6', 145246770, file_path+"/4mo-Tam/10_ctl/outs/possorted_genome_bam.bam", "kras_4mo_exp_counts.csv")
