# CRISPR Assay Phylogeny

Visualize gene detection and RPA primer/gRNA sequence variation alongside a midpoint-rooted phylogeny. Evaluates all qualifying BLAST loci per genome and displays the best assay match.

## Install

```bash
conda create -n assay-phylogeny -c conda-forge -c bioconda \
  python biopython matplotlib blast mafft
conda activate assay-phylogeny
```

## Inputs

- Genome assemblies in FASTA format (`.fa`, `.fna`, `.fasta`, optionally gzipped).
- One reference gene sequence in FASTA format.
- DNA oligo sequences named `fwdPrimer`, `revPrimer` and `gRNA`, supplied 5′→3′.
- A Newick tree with branch lengths and tip identifiers matching genome filenames.

## Run

```bash
python assay_phylogeny.py \
  --genomes genomes/ \
  --gene-reference reference.fasta \
  --assay-oligos oligos.fasta \
  --tree tree.nwk \
  --output-dir results/ \
  --threads 8
```

Defaults: ≥60% gene identity, ≥50% reference coverage, ≤6 substitutions per oligo and a 100–250 bp product. Sites containing indels or ambiguous bases are excluded.

## Outputs

PDF/PNG figure, selected gene and assay-site FASTAs, candidate-locus CSV reports and run metadata.

Use `python assay_phylogeny.py --help` for all options.
