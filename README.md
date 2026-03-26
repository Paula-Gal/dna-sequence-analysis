# DNA Sequence Analysis

A beginner-friendly Python toolkit for learning how DNA works through code.
Built to grow — start with the basics, add more as you learn.

## What it does

| Feature | Description |
|---|---|
| Nucleotide counting | How many A, T, G, C are in a sequence |
| GC content | % of G+C bases (higher = more thermally stable) |
| Complement | Flip each base to its pair (A↔T, G↔C) |
| Reverse complement | The antiparallel partner strand |
| Transcription | DNA → mRNA (T becomes U) |
| Translation | mRNA codons → protein (amino acid chain) |
| Motif search | Find all positions of a sub-sequence |
| ORF finder | Detect regions that could encode a protein |
| Mutation detection | Find point differences between two sequences |
| FASTA reader | Load sequences from the standard `.fasta` file format |

## Project structure

```
dna-sequence-analysis/
├── dna.py          # Core DNASequence class + FASTA reader
├── main.py         # Demo — runs through all features
├── sample.fasta    # Example sequences in FASTA format
└── README.md
```

## How to run

```bash
python3 main.py
```

No external libraries needed — pure Python standard library.

## Key concepts (quick reference)

**Nucleotides** — the 4 building blocks of DNA: Adenine (A), Thymine (T), Guanine (G), Cytosine (C).

**Base pairing** — A always pairs with T; G always pairs with C. This is how the double helix stays together.

**GC content** — G-C pairs have 3 hydrogen bonds vs 2 for A-T, so sequences with more G/C are harder to "melt" apart. Important for lab techniques like PCR.

**Transcription** — DNA is copied into messenger RNA (mRNA) inside the nucleus. Every T becomes U in RNA.

**Translation** — Ribosomes read mRNA in 3-base windows called *codons*. Each codon maps to one amino acid. A chain of amino acids folds into a protein.

**Open Reading Frame (ORF)** — A stretch of DNA starting at ATG (start codon) and ending at a stop codon (TAA, TAG, TGA). ORFs are candidate protein-coding genes.

**FASTA format** — The most common format for storing biological sequences:
```
>sequence_name optional description
ATGCATGCATGCATGC
ATGCATGC
```

## Ideas for what to add next

- [ ] Sliding-window GC content (see how GC% varies across a long sequence)
- [ ] Multiple sequence alignment (line up two sequences to find shared regions)
- [ ] Codon usage analysis (which codons are used most often?)
- [ ] Visualize nucleotide distribution with `matplotlib`
- [ ] Search for restriction enzyme cut sites
- [ ] Load real sequences from NCBI databases