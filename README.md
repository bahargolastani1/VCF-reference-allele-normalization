# VCF-reference-allele-normalization

Reads a VCF file (variants + sample genotypes) and a FASTA file (reference sequences).
  
Checks each variant:
- If the FASTA base matches the VCF reference allele → leave unchanged
- If it matches the alternate allele → swap REF/ALT and update genotypes accordingly
  
Ensures genotype meaning is preserved (with consistent handling of heterozygous 0/1)

Outputs a corrected, normalized VCF file
