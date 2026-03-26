"""
dna.py — Core DNA sequence class.

Concepts covered:
  - Nucleotides: A, T, G, C
  - Complement: A↔T, G↔C
  - GC content: measure of sequence stability
  - Transcription: DNA → mRNA (T becomes U)
  - Translation: mRNA codons → amino acids (simplified)
  - Motif search: finding a sub-sequence within a sequence
  - Open Reading Frames (ORFs): regions that could encode a protein
  - Mutation detection: finding differences between two sequences
  - FASTA parsing: reading the standard bioinformatics file format
"""

# The genetic code: maps each 3-letter mRNA codon to an amino acid (1-letter code).
# '*' means stop codon — translation ends here.
CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",  # AUG = Start codon
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT_MAP = str.maketrans("ATGC", "TACG")


class DNASequence:
    """Represents a single-stranded DNA sequence and provides analysis methods."""

    def __init__(self, sequence: str, name: str = "unnamed"):
        self.name = name
        self.sequence = sequence.upper().strip()
        self._validate()

    def _validate(self):
        valid = set("ATGC")
        invalid = set(self.sequence) - valid
        if invalid:
            raise ValueError(f"Invalid nucleotides found: {invalid}. Only A, T, G, C are allowed.")

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        preview = self.sequence[:20] + "..." if len(self.sequence) > 20 else self.sequence
        return f"DNASequence('{preview}', length={len(self)})"

    # -------------------------------------------------------------------------
    # Basic stats
    # -------------------------------------------------------------------------

    def count_nucleotides(self) -> dict:
        """Count how many of each nucleotide (A, T, G, C) are in the sequence."""
        return {base: self.sequence.count(base) for base in "ATGC"}

    def gc_content(self) -> float:
        """
        Return the GC content as a percentage.
        GC content matters because G-C pairs have 3 hydrogen bonds (vs 2 for A-T),
        making high-GC sequences more thermally stable.
        """
        if len(self) == 0:
            return 0.0
        gc = self.sequence.count("G") + self.sequence.count("C")
        return (gc / len(self)) * 100

    # -------------------------------------------------------------------------
    # Strand operations
    # -------------------------------------------------------------------------

    def complement(self) -> str:
        """
        Return the complementary strand (5'→3' direction preserved).
        A pairs with T, G pairs with C.
        """
        return self.sequence.translate(COMPLEMENT_MAP)

    def reverse_complement(self) -> str:
        """
        Return the reverse complement — the actual antiparallel partner strand.
        DNA is double-stranded and antiparallel, so the partner strand runs 3'→5'
        relative to our sequence, which means we reverse it.
        """
        return self.complement()[::-1]

    # -------------------------------------------------------------------------
    # Central dogma: DNA → RNA → Protein
    # -------------------------------------------------------------------------

    def transcribe(self) -> str:
        """
        Transcribe DNA to mRNA.
        In transcription, the DNA template is read and mRNA is built:
        every T in the coding strand becomes U in mRNA.
        """
        return self.sequence.replace("T", "U")

    def translate(self) -> str:
        """
        Translate mRNA into a protein sequence (amino acid chain).
        - Reads the mRNA in 3-base windows called codons.
        - Each codon maps to one amino acid.
        - Translation starts at the first AUG (Met) codon and stops at a stop codon (*).
        - Returns the protein as a string of 1-letter amino acid codes.
        """
        mrna = self.transcribe()

        # Find the start codon
        start = mrna.find("AUG")
        if start == -1:
            return "(no start codon found)"

        protein = []
        for i in range(start, len(mrna) - 2, 3):
            codon = mrna[i:i + 3]
            if len(codon) < 3:
                break
            amino_acid = CODON_TABLE.get(codon, "?")
            if amino_acid == "*":
                break  # Stop codon reached
            protein.append(amino_acid)

        return "".join(protein) if protein else "(empty protein)"

    # -------------------------------------------------------------------------
    # Motif / pattern searching
    # -------------------------------------------------------------------------

    def find_motif(self, motif: str) -> list[int]:
        """
        Find all positions where a motif (sub-sequence) appears.
        Returns a list of 1-based positions (like biology conventions).
        """
        motif = motif.upper()
        positions = []
        start = 0
        while True:
            pos = self.sequence.find(motif, start)
            if pos == -1:
                break
            positions.append(pos + 1)  # Convert to 1-based index
            start = pos + 1
        return positions

    # -------------------------------------------------------------------------
    # Comparison
    # -------------------------------------------------------------------------

    def similarity(self, other: "DNASequence") -> float:
        """
        Compare two sequences position by position.
        Returns the percentage of positions that match.
        Uses the shorter sequence length to avoid index errors.
        """
        min_len = min(len(self), len(other))
        if min_len == 0:
            return 0.0
        matches = sum(a == b for a, b in zip(self.sequence, other.sequence))
        return (matches / min_len) * 100

    def mutations(self, other: "DNASequence") -> list[dict]:
        """
        Find point mutations (substitutions) between this sequence and another.
        A point mutation is a single nucleotide that differs between two sequences.
        Returns a list of dicts: {position, original, mutated}
        Position is 1-based (biology convention).
        """
        min_len = min(len(self), len(other))
        return [
            {"position": i + 1, "original": self.sequence[i], "mutated": other.sequence[i]}
            for i in range(min_len)
            if self.sequence[i] != other.sequence[i]
        ]

    # -------------------------------------------------------------------------
    # Open Reading Frames (ORFs)
    # -------------------------------------------------------------------------

    def find_orfs(self, min_length: int = 10) -> list[dict]:
        """
        Find all Open Reading Frames (ORFs) in the sequence.

        An ORF is a stretch of DNA that:
          1. Starts with a start codon (ATG)
          2. Ends with a stop codon (TAA, TAG, or TGA)
          3. Is read in a consistent reading frame (every 3 bases)

        ORFs are candidates for protein-coding genes. Real gene finders are
        more complex, but this gives you the basic idea.

        Returns a list of dicts with: start, end, length, protein.
        min_length filters out very short ORFs (likely not real genes).
        """
        stop_codons = {"TAA", "TAG", "TGA"}
        orfs = []

        # Check all 3 reading frames (positions 0, 1, 2)
        for frame in range(3):
            i = frame
            while i < len(self.sequence) - 2:
                codon = self.sequence[i:i + 3]
                if codon == "ATG":  # Start codon found
                    # Now scan forward for a stop codon in this frame
                    for j in range(i, len(self.sequence) - 2, 3):
                        stop = self.sequence[j:j + 3]
                        if stop in stop_codons:
                            orf_seq = self.sequence[i:j + 3]
                            if len(orf_seq) >= min_length:
                                orf = DNASequence(orf_seq)
                                orfs.append({
                                    "start": i + 1,       # 1-based
                                    "end": j + 3,          # 1-based, inclusive
                                    "length": len(orf_seq),
                                    "frame": frame + 1,
                                    "protein": orf.translate(),
                                })
                            break
                i += 3  # Stay in this reading frame

        return orfs


# -----------------------------------------------------------------------------
# FASTA file reader
# -----------------------------------------------------------------------------

def read_fasta(filepath: str) -> list["DNASequence"]:
    """
    Read a FASTA file and return a list of DNASequence objects.

    FASTA is the most common format for storing biological sequences.
    Each record looks like:
        >sequence_name optional description
        ATGCATGCATGC...
        ATGCATGC...       (sequence can span multiple lines)

    Lines starting with '>' are headers; everything else is sequence data.
    """
    sequences = []
    current_name = None
    current_bases = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save the previous record before starting a new one
                if current_name is not None:
                    seq = DNASequence("".join(current_bases))
                    seq.name = current_name
                    sequences.append(seq)
                current_name = line[1:]  # Strip the '>'
                current_bases = []
            else:
                current_bases.append(line.upper())

        # Don't forget the last record
        if current_name is not None:
            seq = DNASequence("".join(current_bases))
            seq.name = current_name
            sequences.append(seq)

    return sequences