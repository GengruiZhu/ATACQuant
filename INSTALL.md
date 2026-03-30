# Installation Guide

## Requirements

ATACQuant requires a Linux/HPC environment and the following external tools:

| Tool | Version | Purpose |
|---|---|---|
| Python | ≥ 3.8 | Core pipeline |
| Bowtie2 | ≥ 2.4 | Read alignment |
| SAMtools | ≥ 1.15 | BAM processing |
| BEDTools | ≥ 2.30 | Sequence extraction |
| Sambamba | ≥ 0.8 | PCR deduplication |
| deepTools | ≥ 3.5 | Tn5 shift correction |
| MACS2 | ≥ 2.2 | Peak calling |
| BioPython | ≥ 1.79 | FASTA parsing |
| pysam | ≥ 0.19 | BAM access |

---

## Step 1 — Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/ATACQuant.git
cd ATACQuant
```

---

## Step 2 — Create and activate the conda environment

```bash
conda env create -f environment.yml
conda activate atacquant
```

> **Note:** If `macs2` fails to install via conda on newer Python versions, install it via pip inside the environment:
> ```bash
> pip install macs2
> ```

---

## Step 3 — Make scripts executable

```bash
chmod +x bin/*.sh bin/*.py workflow/ATACQuant.py
```

---

## Step 4 — Verify installation

```bash
python workflow/ATACQuant.py --help
```

You should see the full argument list without errors.

---

## Troubleshooting

- **Bowtie2 index not found:** Make sure the index prefix you supply via `-f` has a matching `.bt2` or `.bt2l` index built in the same directory.
- **BAM not indexed:** ATACQuant expects BAM files to have a `.bai` index alongside them.
- **Module not found:** Ensure you have activated the `atacquant` conda environment before running.
