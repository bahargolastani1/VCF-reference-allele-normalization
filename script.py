#!/usr/bin/env python3

# VCF reference allele normalization 
# -----------------------------------
# Takes a VCF and a reference FASTA, swaps REF/ALT alleles when the
# reference base in the VCF does not match the FASTA, and rewrites
# genotypes accordingly.

import argparse    # for clean command‑line argument parsing

parser = argparse.ArgumentParser(description='VCF reference normalisation using FASTA reference')
parser.add_argument('--vcf_file', metavar='', required=True, help='VCF input file')
parser.add_argument('--fasta_file', metavar='', required=True, help='FASTA input file')
parser.add_argument('--output_file', metavar='', required=True, help='Output VCF file')
args = parser.parse_args()

def load_fasta(fasta_file):
    # dictionary that will hold {chromosome: sequence}
    fasta = {}
    current_chr = None                    # chromosome currently being processed 
    current_seq = []                      # list to accumulate sequence lines for that chromosome
    print("Reading FASTA file...")
    # open FASTA and read line by line
    with open(fasta_file) as f:
        for line in f:                    # for every line in FASTA file
            line = line.rstrip()          # remove newline
            if line.startswith('>'):      # header line
                current_chr = line[4:]    # drop '>chr'
                current_seq = []          # reset sequence list
            else:                         # sequence line
                current_seq.append(line)  # collect the sequence chunk
                # join the sequence chunks 'so far' and store/update in the dict
                fasta[current_chr] = ''.join(current_seq)
                
    return fasta

def load_vcf(vcf_file):
    # list that will hold every parsed variant record
    vcf_data = []
    print("Reading VCF file...")
    # open the VCF and read line by line
    with open(vcf_file) as v:
        for line in v:                  # for every line in VCF file
            line = line.rstrip()        # remove trailing newline
            if line.startswith('#'):    # header/meta
                continue                # skip to next line
            # split the record by whitespace (VCF is tab‑delimited)
            fields = line.split()
            # extract the VCF columns
            chrom = fields[0]           # chromosome
            pos = int(fields[1]) - 1    # convert to 0 base position for python
            ID = fields[2]              # variant ID
            ref = fields[3]             # reference allele
            alt = fields[4]             # alternate allele
            qual = fields[5]            # quality score 
            filt = fields[6]            # PASS / filter status
            info = fields[7]            # INFO
            form = fields[8]            # format
            samples = fields[9:]        # genotype columns (one per sample)

            # Store everything in a list
            vcf_data.append([chrom, pos, ID, ref, alt, qual, filt, info, form] + samples)
    return vcf_data

def save_vcf_header(vcf_file):
    header = ''
    print("Saving header line in VCF file...")
    with open(vcf_file) as v:
        for line in v:
            # keep only the header line that starts with '#'
            if line.startswith('#') and not line.startswith('##'):
                header = '\t'.join(line.strip().split()) + '\n'    # final result is '#CHROM\tPOS\tID...' 
    return header

def extract_fasta_nts(vcf_data, fasta_data):
    # list that will contain the nucleotides from the FASTA sequence associated with each variant in the VCF file
    fasta_nt = []
    print("Extracting the corresponding nucleotide from the FASTA file for each position in the VCF file...")
    for variant in vcf_data:                      # variant = [chrom, pos, ...]
        for chromo in fasta_data:                 # scan for matching chromosomes
            if chromo == variant[0]:              # if chromosome no. in the FASTA dict matches the chromosome no. in VCF file 
                fasta_nts = fasta_data[chromo]    # then search the FASTA dict value using chromosome no. as the key, and assign to fasta_nts 
                index = variant[1]                # variant[1] is the POS in the VCF file
                # ensure the position in the VCF is in range of the FASTA seq
                if index < 0 or index >= len(fasta_nts):    # if out of range then print error
                    raise ValueError(f"Position {variant[1] + 1} is outside the length " f"of '{chromo}' (length {len(fasta_nts)})."
                                    )
                fasta_nt.append(fasta_nts[index])
    return fasta_nt

def normalise_vcf_data(vcf_data, fasta_nt):
    x = 0
    # list that will contain the normalised VCF data
    vcf_data_v2 = []
    print("Normalising the VCF file...")
    while x < len(fasta_nt):    # length will equal the no. of entries/positions in the vcf file
        for line in vcf_data:
            # convert pos back to 1-based before we output
            line[1] = line[1] + 1
            if fasta_nt[x] == line[3]:
                # FASTA agrees with REF, keep record unchanged - **named tuple or data class
                vcf_data_v2.append(line)
            elif fasta_nt[x] != line[3]:
                # FASTA disagrees, swap REF/ALT
                line[3], line[4] = line[4], line[3]
                # Adjust every genotype field
                n = 9    # define first genotype column
                while n < len(line): # while loop will ensure that each sample has genotype corrected accordingly, will terminate after the last genotype column
                    geno = line[n]
                    g, t = geno[0], geno[2]    # two alleles (diploid)
                    if '/' in geno:
                        sep = '/'    # '/' = unphased
                    else:
                        sep = '|'    # '|' = phased
                    if g == '0' and t == '1'and sep == '/':
                        swap = '0' + sep + '1'
                    elif g == '1' and t == '1' and sep == '/':
                        swap = '0' + sep + '0'
                    elif g == '0' and t =='0' and sep == '/':
                        swap = '1' + sep + '1'
                    elif g == '0' and t =='0' and sep == '/':
                        swap = '1' + sep + '1'
                    elif g == '0' and t =='1' and sep == '|':
                        swap = '1' + sep + '0'
                    elif g == '1' and t =='0' and sep == '|':
                        swap = '0' + sep + '1'
                    elif g == '0' and t =='0' and sep == '|':
                        swap = '1' + sep + '1'
                    elif g == '1' and t =='1' and sep == '|':
                        swap = '0' + sep + '0'
                    else: 
                        print('Error in formatting')
                    line[n] = swap
                    n += 1 # n+1 will shift to the next genotype column
                vcf_data_v2.append(line)
            x += 1     # keep fasta_nt in sync with vcf_data
    return vcf_data_v2

def write_vcf(vcf_v2, vcf_header, output_file):
    print("Writing output to new VCF file...")
    with open(output_file, 'w') as file:
        # writes the header line saved previously
        file.write(vcf_header)
        # loop over every normalised variant record we just produced…
        for row in vcf_v2:
            # …convert each column to a string, join them with tabs, append newline, write it out
            file.write('\t'.join(str(col) for col in row) + '\n')

# Check if this script is being run directly (not imported as a module)
if __name__ == '__main__':
    
    # Load the FASTA file into memory (contains reference sequences)
    fasta_file = load_fasta(args.fasta_file)

    # Load the VCF file into memory (contains variant information)
    vcf_file = load_vcf(args.vcf_file)

    # Extract and save the header lines from the VCF file (metadata and column names)
    vcf_header = save_vcf_header(args.vcf_file)

    # Extract the reference nucleotides from the FASTA file for each variant position
    fasta_nt = extract_fasta_nts(vcf_file, fasta_file)

    # Normalize the VCF data using the extracted nucleotides (e.g., adjust variants for consistency)
    output_data = normalise_vcf_data(vcf_file, fasta_nt)

    # Write the normalized VCF data and header to a new output file
    write_vcf(output_data, vcf_header, args.output_file)
    print("Success")

