import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import re
import random

# Initialize codon table for standard genetic code
codon_table = CodonTable.unambiguous_dna_by_name["Standard"]

# Helper function to find synonymous codons
def get_synonymous_codons(codon):
    """Return list of synonymous codons for a given codon."""
    amino_acid = str(Seq(codon).translate())
    return [c for c in codon_table.forward_table if codon_table.forward_table[c] == amino_acid and c != codon]

# Function to apply silent mutations (updated for flexibility)
def apply_silent_mutations(sequence, guide):
    """Find CRISPR sites and introduce silent mutations."""
    modified_seq = list(sequence)
    modified_count = 0
    debug_log = []
    
    # Find all occurrences of guide RNA in sequence
    for match in re.finditer(guide, sequence):
        start, end = match.start(), match.end()
        debug_log.append(f"Match found at {start}-{end}")
        
        # Adjust start to the nearest codon boundary
        codon_start = start - (start % 3) if start % 3 != 0 else start
        if codon_start < start:
            codon_start += 3  # Move to next codon if start is mid-codon
        debug_log.append(f"Adjusted codon start: {codon_start}")
        
        # Process each codon in the guide region
        for i in range(codon_start, min(end, len(sequence)), 3):
            codon = sequence[i:i+3]
            if len(codon) != 3:
                debug_log.append(f"Skipping incomplete codon at {i}")
                continue
                
            # Get synonymous codons
            syn_codons = get_synonymous_codons(codon)
            if not syn_codons:
                debug_log.append(f"No synonymous codons for {codon} at {i}")
                continue
                
            # Choose a random synonymous codon
            new_codon = random.choice(syn_codons)
            
            # Verify protein translation remains unchanged
            if str(Seq(codon).translate()) == str(Seq(new_codon).translate()):
                modified_seq[i:i+3] = list(new_codon)
                modified_count += 1
                debug_log.append(f"Mutated {codon} to {new_codon} at {i}")
            else:
                debug_log.append(f"Translation mismatch for {codon} to {new_codon} at {i}")
    
    return "".join(modified_seq), modified_count, debug_log

# Test sequences and guide RNAs (same 33 as before)
test_cases = [
    ("ATGCGTAAAGCTGGACACTATAGCATAGACCGTTAGCTAGCTAGCTAGCTAGC", "GGCCATTACTCGATCGATCG"),
    ("ATGTTGGCTATCACACGGAGTATTGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "CTGGCCATTACTCGATCGAT"),
    ("ATGCGTAAAGCTGGCCATTATTCTATAGACCGGTAGCTCGCTAGCTAGCTAGC", "TACTCGATCGATCGCTAGCT"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "ATCGATCGCTAGCTAGCTAG"),
    ("ATGTCTAGCGCTGGGCACTATAGTATAGACCGATAGCTAGCTAGCTAGCTAGC", "GGCCATTACTCGATCGATCG"),
    ("ATGCGTAAAGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "CGATCGCTAGCTAGCTAGCT"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "TTACTCGATCGATCGCTAGC"),
    ("ATGTCTAGCGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "ATTACTCGATCGATCGCTAG"),
    ("ATGCGTAAAGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "CTAGCTAGCTAGCTAGCTAG"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "GGCCATTACTCGATCGATCG"),
    ("ATGTCTAGCGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "ATCGATCGCTAGCTAGCTAG"),
    ("ATGCGTAAAGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "TTACTCGATCGATCGCTAGC"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "CGATCGCTAGCTAGCTAGCT"),
    ("ATGTCTAGCGCTGGACACTATAGCATAGACCGTTAGCTAGCTAGCTAGCTAGC", "CTAGCTAGCTAGCTAGCTAG"),
    ("ATGCGTAAAGCTGGCCATTATTCTATAGACCGGTAGCTCGCTAGCTAGCTAGC", "GGCCATTACTCGATCGATCG"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "ATTACTCGATCGATCGCTAG"),
    ("ATGTCTAGCGCTGGGCACTATAGTATAGACCGATAGCTAGCTAGCTAGCTAGC", "TACTCGATCGATCGCTAGCT"),
    ("ATGCGTAAAGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "CTGGCCATTACTCGATCGAT"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "GGCCATTACTCGATCGATCG"),
    ("ATGTCTAGCGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "CGATCGCTAGCTAGCTAGCT"),
    ("ATGCGTAAAGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "TTACTCGATCGATCGCTAGC"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "ATTACTCGATCGATCGCTAG"),
    ("ATGTCTAGCGCTGGACACTATAGCATAGACCGTTAGCTAGCTAGCTAGCTAGC", "CTAGCTAGCTAGCTAGCTAG"),
    ("ATGCGTAAAGCTGGCCATTATTCTATAGACCGGTAGCTCGCTAGCTAGCTAGC", "GGCCATTACTCGATCGATCG"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "TACTCGATCGATCGCTAGCT"),
    ("ATGTCTAGCGCTGGGCACTATAGTATAGACCGATAGCTAGCTAGCTAGCTAGC", "ATCGATCGCTAGCTAGCTAG"),
    ("ATGCGTAAAGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "CTGGCCATTACTCGATCGAT"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "GGCCATTACTCGATCGATCG"),
    ("ATGTCTAGCGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "TTACTCGATCGATCGCTAGC"),
    ("ATGCGTAAAGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGC", "ATTACTCGATCGATCGCTAG"),
    ("ATGCTGGCCATTACTCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "CGATCGCTAGCTAGCTAGCT"),
    ("ATGTCTAGCGCTGGACACTATAGCATAGACCGTTAGCTAGCTAGCTAGCTAGC", "CTAGCTAGCTAGCTAGCTAG"),
    ("ATGCGTAAAGCTGGCCATTATTCTATAGACCGGTAGCTCGCTAGCTAGCTAGC", "GGCCATTACTCGATCGATCG"),
]

# Main CLI function to test sequences and save results
def run_tests():
    """Run tests on 33 sequences and save results to result.txt."""
    with open("result.txt", "w") as f:
        for i, (sequence, guide) in enumerate(test_cases, 1):
            # Validate inputs
            if not sequence or not guide:
                result = "Invalid: Missing sequence or guide RNA"
                debug_log = ["No sequence or guide RNA provided"]
            elif len(guide) != 20:
                result = "Invalid: Guide RNA must be 20 bp"
                debug_log = ["Guide RNA length invalid"]
            elif not all(c in "ATCG" for c in sequence + guide):
                result = "Invalid: Invalid DNA sequence or guide RNA"
                debug_log = ["Invalid characters in sequence or guide RNA"]
            else:
                # Process sequence
                modified_seq, modified_count, debug_log = apply_silent_mutations(sequence, guide)
                result = f"Modified Sequence: {modified_seq}\nCRISPR Sites Evaded: {modified_count}"
            
            # Write to file
            f.write(f"Test Case {i}:\n")
            f.write(f"Original Sequence: {sequence}\n")
            f.write(f"Guide RNA: {guide}\n")
            f.write(f"{result}\n")
            f.write("Debug Log:\n" + "\n".join(debug_log) + "\n")
            f.write("-" * 50 + "\n")

# Run the CLI
if __name__ == "__main__":
    run_tests()
    print("Testing complete. Results saved to result.txt")