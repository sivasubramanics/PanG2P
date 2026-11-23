# PanG2P

**PanG2P: A pangenome-informed genome-to-phenome framework for trait mapping and prediction in crop plants.**

PanG2P constructs a nonredundant haplotype block catalog from pangenome assemblies and enables scalable, alignment-free, reference-bias-free genotyping of pangenome-derived haplotypes across large populations.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input files](#input-files)
- [Usage](#usage)
- [Configuration](#configuration)
- [Output](#output)
- [Resuming and recovering runs](#resuming-and-recovering-runs)
- [Running on SLURM](#running-on-slurm)
- [Limitations](#limitations)
- [Development](#development)
- [AI Assistance](#ai-assistance)
- [License](#license)
- [Contact](#contact)

## Overview

PanG2P is a Python package that wraps a [Snakemake](https://snakemake.readthedocs.io/) workflow behind a single command, `pang2p`. The workflow has four stages:

1. **Single-copy regions:** mask repeated k-mers in each genome and extract single-copy regions.
2. **Clustering:** cluster single-copy regions across genomes with MMseqs2 and pick representatives.
3. **Pan-k-mers:** build the pan-k-mer reference catalogue from the representatives.
4. **Genotyping:** count pan-k-mers in each sample's reads (BBTools), then normalize and export genotype tables.

## Installation

PanG2P supports **Linux**. It needs Python >= 3.9, Snakemake >= 9, and one of two ways of providing the bioinformatics tools (BBTools, samtools, bedtools, MMseqs2, R and others):

- **Apptainer/Singularity** (recommended on clusters): rules run inside a pre-built container image.
- **Conda/Mamba/Micromamba**: each rule creates its own environment.

Apptainer and conda are not pip-installable, so install one of them separately (for example `module load apptainer` on a cluster, or follow the [Apptainer](https://apptainer.org/docs/admin/main/installation.html) / [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) instructions).

```bash
# Create and activate a dedicated environment
micromamba create -n pang2p -c conda-forge -c bioconda python=3.12 "snakemake>=9" pandas pyyaml apptainer
micromamba activate pang2p

# Install PanG2P
pip install git+https://github.com/sivasubramanics/PanG2P.git
# or from a local clone:
#   git clone https://github.com/sivasubramanics/PanG2P.git
#   pip install ./PanG2P

pang2p --version
```

`conda` or `mamba` can be used instead of `micromamba` in the commands above.

### macOS

Apptainer does not run natively on macOS, so use `--use-conda`. Some bioconda packages are not built for Apple Silicon (`osx-arm64`). If an environment fails to solve, try Intel builds with `CONDA_SUBDIR=osx-64 pang2p ... --use-conda`.

### Air-gapped systems

Compute nodes often have no internet access. Pull the container image once on a machine that does, then point PanG2P at the local file:

```bash
apptainer pull pangat.sif oras://ghcr.io/sivasubramanics/pangat:f58cfca10bda
pang2p ... --use-container --container /path/to/pangat.sif
```

The image is available from [GitHub Packages](https://github.com/sivasubramanics/PanG2P/pkgs/container/pangat). For conda mode on such systems, build the environments beforehand on a login node:

```bash
pang2p -g genomes.csv -s samples.csv --use-conda --snakemake-args --conda-create-envs-only
```

## Quick start

`examples/run_test.sh` downloads a small test dataset, builds the input CSVs and runs the pipeline:

```bash
./examples/run_test.sh --use-conda          # or --use-container on a cluster
THREADS=8 MEM=16g ./examples/run_test.sh    # change resources
DRY_RUN=1 ./examples/run_test.sh            # show what would run
```

The script verifies the download against the published SHA-256 checksum, and writes everything to `pang2p_test/`.

## Input files

PanG2P takes two CSV files, each with a header row.

**Genomes** (`-g`): one row per assembly.

```csv
id,path
genomeA,/data/genomes/genomeA.fa
genomeB,/data/genomes/genomeB.fa
```

**Samples** (`-s`): one row per sequenced sample. `path` is either a single `.fq.gz`/`.fastq.gz` file or a directory containing such files (searched recursively).

```csv
sample_id,path
sampleA,/data/reads/sampleA.fq.gz
sampleB,/data/reads/sampleB_dir/
```

Notes:

- Relative paths are resolved against the folder containing the CSV, not the working directory.
- IDs must be unique. Use only letters, digits, `_` and `-` in genome IDs, because they become Snakemake wildcards.
- The CSVs are validated before anything runs, and copies with absolute paths are written to `<outdir>/inputs/`.

## Usage

```bash
pang2p -g genomes.csv -s samples.csv -t 32 -m 200g -o results
```

If neither `--use-container` nor `--use-conda` is given, PanG2P uses Apptainer when it finds it, otherwise conda. Use `--skip-checks` to disable this and run rules with whatever is already on `PATH` (for example inside a Docker image that contains all tools).

Dry run to validate the workflow without executing any jobs:

```bash
pang2p -g genomes.csv -s samples.csv -n
```

Run with Apptainer or conda explicitly:

```bash
pang2p -g genomes.csv -s samples.csv -t 32 -m 200g --use-container
pang2p -g genomes.csv -s samples.csv -t 8 -m 16g --use-conda
```

### Command-line options

| Option | Description |
|---|---|
| `-g, --genomes` | Genomes CSV (`id,path`) |
| `-s, --samples` | Samples CSV (`sample_id,path`) |
| `-o, --outdir` | Output directory [`pang2p_out`] |
| `-t, --threads` | Threads per job |
| `-m, --memory` | Memory for Java tools, e.g. `200g` |
| `--use-conda` / `--use-container` | Software environment |
| `--skip-checks` | Skip dependency checks and environment auto-selection |
| `-n, --dry-run` | Show what would run |
| `--unlock` | Remove a stale lock after a crashed run, then exit |
| `--force` | Rerun everything from scratch |
| `-k, --keep-going` | Continue independent jobs after a failure |
| `--retries N` | Retry failed jobs N times |
| `--ignore-changes` | On resume, rerun only when files are missing or newer |
| `-v, --verbose` / `-q, --quiet` | Debug / warnings-only console output |
| `--snakemake-args ...` | Pass everything after it to Snakemake unchanged (must be last) |

```bash
pang2p -g genomes.csv -s samples.csv -c my_config.yaml -t 32 -m 200g
```

## Output

Everything is written below the output directory:

```
results/
├── inputs/                  # validated copies of the CSVs with absolute paths
├── pang2p_config.yaml       # effective settings of the latest run
├── logs/
│   ├── pang2p/run_<time>.log   # full log of each run (including Snakemake output)
│   └── <step>/...              # per-rule logs
├── benchmarks/              # per-rule runtime and memory
└── work/
    ├── 01.single_copy/
    ├── 02.clustering/
    ├── 03.pankmers/
    └── 04.genotyping/
        ├── final.genotype.tsv       # normalized genotype table
        ├── final.GAPIT.RData        # input for GAPIT
        ├── overall_counts.tsv.gz
        └── overall_mapped_kmers.tsv.gz
```

## Resuming and recovering runs

Re-running the same command with the same `-o` resumes the run. Finished steps are skipped and incomplete outputs are redone.

| Situation | What to do |
|---|---|
| Job failed or was interrupted | Re-run the same command |
| "Directory cannot be locked" after a crash | Run the command once with `--unlock`, then re-run it |
| One bad sample or genome | Use `-k` so independent jobs finish; add `--retries 2` for flaky jobs |
| Changed a parameter such as `--pan-kmer-size` | Affected steps rerun, and a warning lists the changed settings |
| Changed a parameter but want to keep existing results | Use `--ignore-changes` |
| Start over | Use `--force`, or delete `<outdir>/work` |

## Running on SLURM

### Single-node job (simplest)

Install PanG2P once on a login node, pre-pull the container (see [Air-gapped systems](#air-gapped-systems)), then submit a script like this:

```bash
#!/bin/bash
#SBATCH --job-name=pang2p
#SBATCH --cpus-per-task=32
#SBATCH --mem=220G
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

module load apptainer                      # site specific
eval "$(micromamba shell hook -s bash)"
micromamba activate pang2p

export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_BIND=/scratch,/project    # filesystems holding your data

pang2p -g genomes.csv -s samples.csv \
    -t "$SLURM_CPUS_PER_TASK" -m 200g -o results \
    --use-container --container /project/$USER/pang2p/pangat.sif \
    --retries 1 --keep-going
```

```bash
sbatch run_pang2p.sh
tail -f results/logs/pang2p/run_*.log
```

Request more `--mem` than `-m`: `-m` sets only the Java heap, and the JVM, samtools and R need additional memory. If the job is killed or reaches its time limit, resubmit the same script (add `--unlock` once if the directory is locked).

### One SLURM job per rule instance

Snakemake can submit each rule as its own job through the SLURM executor plugin:

```bash
pip install snakemake-executor-plugin-slurm

pang2p -g genomes.csv -s samples.csv -t 16 -m 64g -o results \
    --use-container --keep-going \
    --snakemake-args --executor slurm --jobs 20 \
        --default-resources slurm_account=<account> slurm_partition=<partition> runtime=240
```

This mode is experimental. In it, `-t` means threads per job, and memory requests come from `mem_mb`, which the workflow does not set, so provide it with `--default-resources mem_mb=...`. Run a dry run and a small test first.

## Limitations

- **Resource allocation:** the workflow assigns most available resources to BBTools steps. Set `max_mem` and `threads` conservatively to avoid contention, especially on shared systems.
- **BBTools / BBDuk JVM memory:** BBTools (including `bbduk`) uses strict JVM memory management (`-m` / `max_mem`), and the JVM itself plus non-Java tools (`samtools`, `R`, etc.) require extra memory beyond the heap. Keep `-m` (`max_mem`) within your allocated memory limits (and request higher total job memory on SLURM) so jobs are not killed with out-of-memory errors.
- **Conda environment file counts (`.snakemake/`):** using `--use-conda` creates separate Conda environments per rule and can generate many thousands of small files under `.snakemake/`. On shared filesystems or HPC clusters with strict inode or small-file quotas, prefer `--use-container` or place your work/conda cache on a filesystem with generous inode limits.
- **Per-job parallelism:** every multithreaded rule requests the full `threads` value, so jobs generally run one after another on a single node.
- **Intermediate files:** the workflow does not clean up intermediate files on completion. Large runs can use a lot of disk space, and manual cleanup may be needed.
- **Platform:** Linux only. On macOS use `--use-conda`.

## Development

```bash
git clone https://github.com/sivasubramanics/PanG2P.git
cd PanG2P
pip install -e .
```

The workflow lives in `src/pang2p/workflow/` (`Snakefile`, `rules/`, `scripts/`, `envs/`, `config.yaml`). The Apptainer recipe is in `container/singularity/`, and a GitHub Actions workflow builds it and pushes it to GitHub Packages. If you change the recipe, update the `container` tag in `config.yaml`.

## AI Assistance

AI tools, including GitHub Copilot, were used to assist with code completion, refactoring, and documentation. All AI-generated content was reviewed by the authors.

## License

Licensed under the GPL-3.0 License.

## Contact

For bugs and feature requests, please [open an issue](https://github.com/sivasubramanics/PanG2P/issues) on GitHub.

For questions or feedback, reach out to:

- [muru.jayakodi@ag.tamu.edu](mailto:muru.jayakodi@ag.tamu.edu)
- [c.s.sivasubramani@gmail.com](mailto:c.s.sivasubramani@gmail.com)
