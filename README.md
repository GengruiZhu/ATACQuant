# ATACQuant

ATACQuant is a *k*-mer-based Bayesian expectation-maximization (EM) framework designed to precisely quantify subgenome-specific chromatin accessibility, effectively resolving the issue of multi-mapping data loss in polyploid ATAC-seq analyses.

It resolves reads that multi-map across homeologous (syntenic) peak regions by combining:
- **Topological anchoring** via GraphAllele allele clusters
- **Dual-layer k-mer immune dictionaries** (transposon masking + genome uniqueness filtering)
- **N-way lineage P-value scoring** from species-specific diagnostic k-mers
- **Bayesian EM deconvolution** with physical length normalization

ATACQuant supports any number of subgenome species (N ≥ 2), making it applicable to allopolyploids of any ploidy level.

---

## Pipeline Overview

```
ATAC-seq FASTQs
      │
      ▼
[00] Bowtie2 alignment + Tn5 shift + MACS2 peak calling
      │
      ▼
[01] Topological anchoring: assign peaks to GraphAllele syntenic clusters
      │
      ├─── Syntenic (multi-subgenome) peaks ──────────────────────────────┐
      │                                                                    │
      ▼                                                                    │
[02] K-mer immune dictionary construction                                  │
      │   (hard-masked FASTA → transposon filter → Bowtie2 uniqueness)    │
      ▼                                                                    │
[02b] K-mer FASTQ scan → per-peak prior theta (opening probability)       │
      │                                                                    │
      ▼                                                                    │
[03] Per-read lineage P-value assignment                                   │
      │   (species-specific k-mer FASTA libraries, N-way)                 │
      ▼                                                                    │
[04] Bayesian EM convergence                                               │
      │   (fuses theta prior + read P-values + physical length)           │
      └──────────────┬────────────────────────────────────────────────────┘
                     │
[05] Final count matrix assembly
      │   (EM-deconvolved syntenic peaks + raw-counted lineage-specific peaks)
      ▼
[06] Local gene annotation via GraphAllele matrix
      │
      ▼
     Final annotated count matrix (DESeq2-ready)
```

---

## Input Files

| Argument | Description |
|---|---|
| `--r1` / `--r2` | Cleaned paired-end FASTQ files |
| `-f` / `--fasta` | Reference genome FASTA |
| `-m` / `--masked` | Hard-masked reference genome FASTA (RepeatMasker output) |
| `-g` / `--gff` | Genome annotation GFF3 |
| `-c` / `--clusters` | GraphAllele `my_clusters.tsv` — allele cluster assignments |
| `--matrix` | GraphAllele `PolyAlleler_Global_Matrix_Cleaned.tsv` |
| `-a` / `--ancestry` | Gene ancestry mapping TSV (columns: Gene_ID, …, …, Ancestry_Label) |
| `--kmers` | Species-specific diagnostic k-mer FASTA files, one per subgenome |
| `--labels` | Short species labels corresponding to `--kmers` (must match Ancestry_Label values) |

### Ancestry mapping file format

Tab-separated, with a header row. The **4th column** (index 3) must contain the ancestry label for each gene, and its values must **exactly match** the `--labels` you supply.

```
Gene_ID    Col2    Col3    Ancestry_Label
Gene0001   ...     ...     SP1
Gene0002   ...     ...     SP2
```

### K-mer FASTA files

Each file should contain species-specific diagnostic k-mers (e.g., produced by KmerGenoPhaser or similar). One sequence per entry; header lines are ignored. Only the first k bases of each sequence are used.

---

## Usage

```bash
cd workflow/

python ATACQuant.py \
  --r1 /path/to/clean_R1.fq \
  --r2 /path/to/clean_R2.fq \
  -f /path/to/genome.fasta \
  -m /path/to/genome.fasta.masked \
  -g /path/to/genome.gff3 \
  -c /path/to/my_clusters.tsv \
  --matrix /path/to/PolyAlleler_Global_Matrix_Cleaned.tsv \
  -a /path/to/Gene_Ancestry_Mapping.tsv \
  --kmers /path/to/SP1_k15.fasta /path/to/SP2_k15.fasta /path/to/SP3_k15.fasta \
  --labels SP1 SP2 SP3 \
  -t 48 \
  -k 15 \
  --prefix SAMPLE_NAME \
  -o /path/to/output_dir
```

See `workflow/run_ATACQuant.sh` for a ready-to-edit cluster submission template.

### Key parameters

| Flag | Default | Description |
|---|---|---|
| `-k` / `--kmer` | `15` | K-mer length for lineage P-value scoring (step 03); must match your k-mer FASTA files |
| `-t` / `--threads` | `48` | CPU threads for Bowtie2 / samtools |
| `--prefix` | `sample` | Sample name used in output filenames and count matrix column header |
| `-o` / `--outdir` | `workflow/ATACQuant_Out` | Output directory |
| `--force` | off | Re-run from scratch, ignoring checkpoints |
| `--up` | `3000` | TSS upstream window (bp) for peak anchoring |
| `--down` | `1000` | TSS downstream window (bp) for peak anchoring |

---

## Output Files

All outputs are written to `--outdir`.

| File | Description |
|---|---|
| `Syntenic_Peak_Clusters.tsv` | Peaks grouped into syntenic allele clusters |
| `EM_Target_Peaks.saf` | SAF-format coordinates for syntenic peaks |
| `Peak_Ancestry.txt` | Per-peak subgenome ancestry label |
| `Lineage_Specific_Peaks.tsv` | Peaks unique to one subgenome (not deconvolved) |
| `ATAC_Master_Immune_Dict.pkl` | Binary immune k-mer dictionary |
| `Global_Theta_Matrix.tsv` | Prior opening probability (theta) per peak |
| `Reads_P_values.tsv` | Per-read lineage likelihood scores |
| `EM_Final_Expression.tsv` | EM-deconvolved read counts + density per peak |
| `{prefix}_Final_Count_Matrix.tsv` | Merged count matrix (syntenic + lineage-specific) |
| `{prefix}_Local_Annotated_Matrix.tsv` | Final matrix annotated with local gene IDs |

---

## N-species Support

ATACQuant generalizes to any number of subgenomes (N ≥ 2). Simply provide one k-mer FASTA and one label per subgenome:

```bash
# Tetraploid (2 subgenomes)
--kmers SP1.fasta SP2.fasta --labels SP1 SP2

# Hexaploid (3 subgenomes)
--kmers SP1.fasta SP2.fasta SP3.fasta --labels SP1 SP2 SP3

# Octoploid (4 subgenomes)
--kmers A.fasta B.fasta C.fasta D.fasta --labels A B C D
```

**Important:** The label values in `--labels` must exactly match the ancestry labels (4th column, uppercase) in the gene ancestry mapping file (`-a`).

---

## Citation

Manuscript in preparation. Please check back for the citation once the paper is published.

---

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/GengruiZhu/ATACQuant/blob/main/LICENSE) file for details.

Copyright (c) 2026 Yi Chen & Gengrui Zhu
