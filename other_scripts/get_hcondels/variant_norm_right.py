import sys

arg=sys.argv

#adapted from
#https://genome.sph.umich.edu/wiki/Variant_Normalization
#and var_norm_left.y
VCFFileName=arg[1]
seqFileName=arg[2]
chrLengthFileName=arg[3]
outFileName=arg[4]

vcf_df=[]

#read in VCF file
with open(VCFFileName, "r") as VCFFile:
	for ind, l1 in enumerate(VCFFile):
		#skip header lines
		if "#"==l1[0]:
			continue
		l1=l1.split()
		chromosome=l1[0]
		ref_start_pos=int(l1[1])
		seq_id=l1[2]
		ref_allele=l1[3]
		alt_allele=l1[4]
		vcf_df.append([seq_id, chromosome, ref_start_pos, ref_allele, alt_allele])

#read in sequence table
#first column sequence id, second column sequence, third column sequence start pos, fourth column sequence end pos
with open(seqFileName, "r") as seqFile:
	for ind, l1 in enumerate(seqFile):
		l1=l1.split()
		seq_id=l1[0]
		#check if seq_id matches VCF file
		if seq_id!=vcf_df[ind][0]:
			sys.exit("discordant seq_id: "+ seq_id)
		ref_sequence=l1[1]
		ref_sequence_chr=l1[2]
		ref_sequence_start=int(l1[3])
		ref_sequence_end=int(l1[4])
		vcf_df[ind].extend([ref_sequence, ref_sequence_start, ref_sequence_end])

chr_length_map=dict()
with open(chrLengthFileName, "r") as chrLengthFile:
	for ind, l1 in enumerate(chrLengthFile):
		l1=l1.split()
		chr_name=l1[0]
		chr_length=int(l1[1])
		chr_length_map[chr_name]=chr_length

def var_norm_right(ref_allele, alt_allele, ref_sequence, ref_end_pos_on_seq, ref_sequence_end, chr_length):
	num_changes=0
	change_in_allele=True
	allele_changed="0"
	while change_in_allele:
							
		ref_allele_orig=ref_allele
		alt_allele_orig=alt_allele	
		ref_end_pos_on_seq_orig=ref_end_pos_on_seq
		
		#if ref_sequence_end is already at beginning of chr, ref_end_pos_on_seq already at the beginning, and length of the ref or alt alleles are 1 (and thus can't trim)
		if (ref_sequence_end==chr_length) and (ref_end_pos_on_seq==len(ref_sequence)-1) and (len(ref_allele)==1 or len(alt_allele)==1):
			break
		
		#if alleles start in same nucleotide, truncate leftmost nucleotide of each allele
		if ref_allele[0]==alt_allele[0]:
			ref_allele=ref_allele[1:len(ref_allele)]
			alt_allele=alt_allele[1:len(alt_allele)]
		
		#if any empty allele, extend 1 bp to the right
		if (ref_allele=="") or (alt_allele==""):
			ref_end_pos_on_seq=ref_end_pos_on_seq+1
			
			if ref_end_pos_on_seq==len(ref_sequence)-1:
				sys.exit("greater than chr length position  encountered\nref_allele: " + ref_allele + "\nalt_allele: " + alt_allele + "\nsequence: " + ref_sequence)
			
			ref_allele=ref_allele+ref_sequence[ref_end_pos_on_seq]
			alt_allele=alt_allele+ref_sequence[ref_end_pos_on_seq]
		
		#if right most nucleotide are the same, and all alleles have length 2 or more, truncate rightmost nucleotide of each allele
		while (ref_allele[len(ref_allele)-1]==alt_allele[len(alt_allele)-1]) and ((len(ref_allele)>=2) and (len(alt_allele)>=2)):
			ref_allele=ref_allele[0:len(ref_allele)-1]
			alt_allele=alt_allele[0:len(alt_allele)-1]
		
		if (ref_end_pos_on_seq_orig==ref_end_pos_on_seq) and (ref_allele_orig==ref_allele) and (alt_allele_orig==alt_allele):
			change_in_allele=False
		else:
			num_changes=num_changes+1
		
		if num_changes!=0:
			allele_changed="1"
	
	normalized_alleles=[ref_allele, alt_allele, ref_end_pos_on_seq, allele_changed]
	
	return(normalized_alleles)

outFile=open(outFileName, "w")

#iterate through all sequences and normalize
for ind in range(len(vcf_df)):
	seq_id=vcf_df[ind][0]
	chromosome=vcf_df[ind][1]
	ref_end_pos=vcf_df[ind][2]
	ref_allele=vcf_df[ind][3]
	alt_allele=vcf_df[ind][4]
	ref_sequence=vcf_df[ind][5]
	ref_sequence_start=vcf_df[ind][6]
	ref_sequence_end=vcf_df[ind][7]
	
	ref_end_pos_on_seq=ref_end_pos-ref_sequence_start
	chr_length=chr_length_map[chromosome]	 
	
	normalized_alleles=var_norm_right(ref_allele, alt_allele, ref_sequence, ref_end_pos_on_seq, ref_sequence_end, chr_length)
	
	#write out new VCF
	normalized_ref_allele=normalized_alleles[0]
	normalized_alt_allele=normalized_alleles[1]
	normalized_ref_start_pos=normalized_alleles[2]+ref_sequence_start
	allele_changed=normalized_alleles[3]
		
	outFile.write(chromosome+"\t"+str(normalized_ref_start_pos)+"\t"+seq_id+"\t"+normalized_ref_allele+"\t"+normalized_alt_allele+"\t"+allele_changed+"\n")

