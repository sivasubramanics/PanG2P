# PanG2P

**PanG2P**: A pangenome-informed genome-to-phenome framework for trait mapping and prediction in crop plants.

PanG2P constructs a nonredundant haplotype block catalog from pangenome assemblies and enables scalable, alignment-free, reference-bias-free genotyping of pangenome-derived haplotypes across large populations.

---

## Table of Contents

- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Limitations](#limitations)
- [AI Assistance](#ai-assistance)
- [License](#license)
- [Contact](#contact)

---

## Installation

PanG2P requires [Snakemake](https://snakemake.readthedocs.io/) >=8.0

```bash
# Create and activate a dedicated environment
conda create -n pang2p -c bioconda snakemake
conda activate pang2p

# Clone the repository
git clone git@github.com:sivasubramanics/PanG2P.git
cd PanG2P
```

---

## Configuration

### Config file

The workflow is configured via `config/config.yaml`. You can edit this file directly or supply a custom config using the `--configfile` flag at runtime.

| Parameter | Default | Description |
|---|---|---|
| `genomes` | — | Path to a CSV file listing input genomes. Columns: `id`, `path` |
| `samples` | — | Path to a CSV file listing sequencing samples. Columns: `id`, `path` |
| `max_mem` | `4` | Maximum memory (GB) allocated per rule |
| `threads` | `2` | Threads allocated per rule |
| `kmer_size` | `31` | *k*-mer size for building the single-copy reference |
| `min_kmer_count` | `2` | Minimum *k*-mer occurrences required to mask a genome region |
| `min_length` | `100` | Minimum length (bp) of a single-copy region |
| `merge_gap` | `50` | Maximum gap (bp) between regions to merge into one |
| `min_seq_id` | `0.9` | Minimum sequence identity for MMseqs2 clustering |
| `min_coverage` | `0.8` | Minimum coverage for MMseqs2 clustering |
| `coverage_mode` | `1` | Coverage mode for MMseqs2 clustering |
| `pan_kmer_size` | `51` | *k*-mer size for the pan-*k*-mer reference catalogue |
| `pan_kmer_step` | `5` | *k*-mer step size for the pan-*k*-mer reference catalogue |

### Runtime environment

PanG2P supports Conda, Docker, and Singularity/Apptainer environments.

For air-gapped systems without internet access, download the pre-built Apptainer image and point to it with the `container` key in `config.yaml`:

```bash
singularity pull pangat.sif oras://ghcr.io/sivasubramanics/pangat:f58cfca10bda
```

The image is available from [GitHub Packages](https://github.com/sivasubramanics/PanG2P/pkgs/container/pangat).

---

## Usage

**Dry run** — validate the workflow configuration without executing any jobs:

```bash
snakemake --dry-run --configfile config.yaml
```

**Run with Conda:**

```bash
# Using conda
snakemake --use-conda --conda-frontend conda --configfile config.yaml --cores 8

# Using mamba (faster dependency resolution)
snakemake --use-conda --conda-frontend mamba --configfile config.yaml --cores 8
```

**Run with Apptainer/Singularity:**

```bash
snakemake --use-singularity --configfile config.yaml --cores 8
```

**Run under a SLURM:**

```bash
sbatch --time 08:00:00 -N 1 -c 42 --mem 200G -J pang2p_run -o pang2p_run.log -e pang2p_run.err snakemake --use-singularity --configfile config.yaml --cores 42
```

---

## Limitations

- **Resource allocation** — The workflow assigns most available resources to BBTools steps. Set `max_mem` and `threads` conservatively to avoid resource contention, especially on shared systems.
- **Intermediate files** — The workflow does not automatically clean up intermediate files on completion. Large runs can consume significant disk space; manual cleanup may be needed.

---

## AI Assistance

AI tools, including GitHub Copilot, were used to assist with code completion, refactoring, and documentation. All AI-generated content was reviewed by the authors.

---

## License

Licensed under the [GPL-3.0 License](LICENSE).

---

## Contact

For bugs and feature requests, please [open an issue](https://github.com/sivasubramanics/PanG2P/issues) on GitHub.

For questions or feedback, reach out to:
- [muru.jayakodi@ag.tamu.edu](mailto:muru.jayakodi@ag.tamu.edu)
- [c.s.sivasubramani@gmail.com](mailto:c.s.sivasubramani@gmail.com)
