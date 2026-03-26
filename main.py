"""
DNA Sequence Analysis - Learning Project
A beginner-friendly toolkit for analyzing DNA sequences.
"""

from dna import DNASequence, read_fasta


def main():
    print("=== DNA Sequence Analysis ===\n")

    # ------------------------------------------------------------------
    # 1. Basic analysis on a raw sequence
    # ------------------------------------------------------------------
    raw_sequence = "ATGCGTAACGTTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGC"
    dna = DNASequence(raw_sequence)

    print(f"Sequence:  {dna.sequence}")
    print(f"Length:    {len(dna)} nucleotides\n")

    print("--- Nucleotide Counts ---")
    for base, count in dna.count_nucleotides().items():
        print(f"  {base}: {count}")

    print(f"\n--- GC Content ---")
    print(f"  GC%: {dna.gc_content():.1f}%")

    print("\n--- Complement Strands ---")
    print(f"  Complement:        {dna.complement()}")
    print(f"  Reverse complement:{dna.reverse_complement()}")

    print("\n--- Transcription (DNA → mRNA) ---")
    print(f"  mRNA: {dna.transcribe()}")

    print("\n--- Translation (mRNA → Protein) ---")
    print(f"  Protein: {dna.translate()}")

    # ------------------------------------------------------------------
    # 2. Motif search
    # ------------------------------------------------------------------
    print("\n--- Motif Search ---")
    motif = "ATCG"
    positions = dna.find_motif(motif)
    print(f"  '{motif}' found at positions: {positions}")

    # ------------------------------------------------------------------
    # 3. Open Reading Frames
    # ------------------------------------------------------------------
    print("\n--- Open Reading Frames (ORFs) ---")
    orfs = dna.find_orfs(min_length=9)
    if orfs:
        for orf in orfs:
            print(f"  Frame {orf['frame']} | pos {orf['start']}–{orf['end']} "
                  f"| {orf['length']} bp | protein: {orf['protein']}")
    else:
        print("  No ORFs found.")

    # ------------------------------------------------------------------
    # 4. Mutation detection
    # ------------------------------------------------------------------
    print("\n--- Mutation Detection ---")
    other_seq = "ATGCGTAACGTTAGCATCGATCGATCGCAGCTAGCTAGCATCGATCGATCGTAGC"
    other = DNASequence(other_seq)
    mutations = dna.mutations(other)
    similarity = dna.similarity(other)

    print(f"  Sequences differ at {len(mutations)} position(s) | Similarity: {similarity:.1f}%")
    for m in mutations:
        print(f"  Position {m['position']}: {m['original']} → {m['mutated']}")

    # ------------------------------------------------------------------
    # 5. Read from a FASTA file
    # ------------------------------------------------------------------
    print("\n--- Reading from FASTA file ---")
    sequences = read_fasta("sample.fasta")
    for seq in sequences:
        print(f"  [{seq.name}]")
        print(f"    Length: {len(seq)} bp | GC: {seq.gc_content():.1f}% | Protein: {seq.translate()}")


if __name__ == "__main__":
    main()