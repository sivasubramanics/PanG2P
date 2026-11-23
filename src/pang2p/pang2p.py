#!/usr/bin/env python3
# Script        : pang2p.py
# Description   : Wrapper script to run the pang2p workflow via snakemake.
# Usage         : pang2p -g genomes.csv -s samples.csv -o outdir [options]
# Version       : 0.1.0
# Date          : 2026-11-23
# Contact       : c.s.sivasubramani@gmail.com

import argparse
import csv
import datetime as dt
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

from . import __version__

WORKFLOW_DIR = Path(__file__).parent / "workflow"
SNAKEFILE = WORKFLOW_DIR / "Snakefile"
DEFAULT_CONFIG = WORKFLOW_DIR / "config.yaml"

log = logging.getLogger("pang2p")
sm_log = logging.getLogger("pang2p.snakemake")
sm_log.propagate = False

# defaults from the config.yaml; overridden by --configfile and CLI flags
PARAMS = {
    # "kmer_size":      ("--kmer-size", int, "k-mer size for single-copy region detection"),
    # "min_kmer_count": ("--min-kmer-count", int, "minimum k-mer count for masking repeats"),
    # "min_length":     ("--min-length", int, "minimum length of a single-copy region"),
    # "merge_gap":      ("--merge-gap", int, "max gap between single-copy regions to merge"),
    # "min_seq_id":     ("--min-seq-id", float, "mmseqs minimum sequence identity"),
    # "min_coverage":   ("--min-coverage", float, "mmseqs minimum coverage"),
    # "coverage_mode":  ("--coverage-mode", int, "mmseqs coverage mode"),
    # "pan_kmer_size":  ("--pan-kmer-size", int, "pan-kmer size"),
    # "pan_kmer_step":  ("--pan-kmer-step", int, "pan-kmer step size"),
    # "container":      ("--container", str, "container image (.sif path or oras:// URI)"),
}


# Utility Function: Parse Memory
def parse_mem_gb(value) -> int:
    """'200g' / '200G' / '200' -> 200 (GB); also accepts k, m, t suffixes."""
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([kmgt]?)b?", str(value).strip().lower())
    if not m:
        raise argparse.ArgumentTypeError(f"invalid memory value: {value!r} (e.g. 200g)")
    gb = float(m.group(1)) * {"k": 1e-6, "m": 1e-3, "g": 1, "t": 1000}[m.group(2) or "g"]
    return max(1, int(gb))


# Utility Function: Load Default Config from config.yaml
def load_defaults() -> dict:
    cfg = yaml.safe_load(DEFAULT_CONFIG.read_text()) or {}
    cfg.pop("genomes", None)
    cfg.pop("samples", None)
    return cfg

# Utility Function: Check for Apptainer/Singularity
def has_apptainer() -> bool:
    return bool(shutil.which("apptainer") or shutil.which("singularity"))

# Utility Function: Check for Conda/Mamba
def has_conda() -> bool:
    return any(shutil.which(x) for x in ("mamba", "conda", "micromamba"))

# Utility Function: Argparse Parser for CLI
def build_parser(d: dict) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pang2p",
        description="Pan-genome single-copy k-mer genotyping pipeline",
        epilog="Resume: re-run the same command with the same -o; finished steps are skipped.",
    )
    p.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")

    g = p.add_argument_group("input / output")
    g.add_argument("-g", "--genomes", type=Path, help="CSV with columns: id,path")
    g.add_argument("-s", "--samples", type=Path, help="CSV with columns: sample_id,dir_path")
    g.add_argument("-o", "--outdir", type=Path, default=Path("pang2p_out"),
                   help="output directory [pang2p_out]")
    # g.add_argument("-c", "--configfile", type=Path,
    #                help="YAML overriding the packaged defaults (CLI flags still win)")

    g = p.add_argument_group("resources")
    g.add_argument("-t", "--threads", type=int, help=f"threads [{d['threads']}]")
    g.add_argument("-m", "--memory", type=parse_mem_gb,
                   help=f"memory in GB for Java tools, e.g. 200g [{d['max_mem']}g]")

    g = p.add_argument_group("workflow parameters (defaults from config.yaml)")
    for key, (flag, typ, text) in PARAMS.items():
        g.add_argument(flag, dest=key, type=typ, default=None,
                       help=f"{text} [{d.get(key)}]")

    g = p.add_argument_group("software environment")
    g.add_argument("--use-conda", action="store_true", help="use per-rule conda envs")
    g.add_argument("--use-container", action="store_true",
                   help="run rules in the Apptainer/Singularity image")
    g.add_argument("--skip-checks", action="store_true",
                   help="skip dependency checks and environment auto-selection "
                        "(use when the tools are already on PATH, e.g. inside Docker)")

    g = p.add_argument_group("run control")
    g.add_argument("-n", "--dry-run", action="store_true", help="show what would run")
    g.add_argument("--unlock", action="store_true",
                   help="remove a stale lock after a crashed run, then exit")
    g.add_argument("--force", action="store_true", help="rerun everything from scratch")
    g.add_argument("-k", "--keep-going", action="store_true",
                   help="continue independent jobs after a failure")
    g.add_argument("--retries", type=int, default=0, help="retry failed jobs N times [0]")
    g.add_argument("--ignore-changes", action="store_true",
                   help="on resume, rerun only when files are missing/newer "
                        "(ignore parameter/code/env changes)")

    g = p.add_argument_group("logging")
    g.add_argument("-v", "--verbose", action="store_true", help="debug messages")
    g.add_argument("-q", "--quiet", action="store_true",
                   help="only warnings on console (everything still goes to the log file)")

    p.add_argument("--snakemake-args", nargs=argparse.REMAINDER, default=[],
                   help="everything after this is passed to snakemake as-is")
    return p

# Logging Setup
def setup_logging(logdir: Path, verbose: bool, quiet: bool) -> Path:
    """
    Set up logging for the pang2p workflow.
    """
    logdir.mkdir(parents=True, exist_ok=True)
    path = logdir / f"run_{dt.datetime.now():%Y%m%d_%H%M%S}.log"
    log.setLevel(logging.DEBUG)
    sm_log.setLevel(logging.DEBUG)

    console = logging.StreamHandler(sys.stderr)
    console.setLevel(logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO)
    console.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S"))
    log.addHandler(console)

    fh = logging.FileHandler(path)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S"))
    log.addHandler(fh)

    raw = logging.FileHandler(path)          # snakemake lines go in verbatim
    raw.setFormatter(logging.Formatter("%(message)s"))
    sm_log.addHandler(raw)
    return path

# Utility Function: Get Snakemake Info
def snakemake_info():
    exe = shutil.which("snakemake")
    if not exe:
        return None, "", ""
    ver = subprocess.run([exe, "--version"], capture_output=True, text=True).stdout.strip()
    h = subprocess.run([exe, "--help"], capture_output=True, text=True)
    return exe, ver, h.stdout + h.stderr

# Utility Function: Preflight Checks
def preflight(use_conda: bool, use_container: bool, have_snakemake: bool):
    errors = []
    if not have_snakemake:
        errors.append("snakemake not found on PATH (pip install 'snakemake>=8').")

    if use_container and not has_apptainer():
        msg = "--use-container needs apptainer (or singularity) on PATH."
        if sys.platform == "darwin":
            msg += " Apptainer does not run natively on macOS; use --use-conda instead."
        errors.append(msg)
    if use_conda and not has_conda():
        errors.append("--use-conda needs conda, mamba or micromamba on PATH.")
    if not use_conda and not use_container:
        errors.append("no software environment available: install apptainer or conda "
                      "(mamba/micromamba), or pass --skip-checks if the tools are "
                      "already on PATH (e.g. inside the Docker image).")
    return errors


# Utility Function: Read CSV and Prepare Inputs
def _read_csv(path: Path):
    with open(path, newline="") as fh:
        return [[c.strip() for c in r] for r in csv.reader(fh) if any(c.strip() for c in r)]

# Utility Function: Prepare Inputs for genomes.csv and samples.csv
def prepare_inputs(genomes: Path, samples: Path, outdir: Path):
    """Validate both CSVs; write copies with absolute paths into <outdir>/inputs/."""
    errors, warnings = [], []

    def resolve(p, base):
        p = Path(p).expanduser()
        return (p if p.is_absolute() else base / p).resolve()

    g_out, s_out = [], []

    rows = _read_csv(genomes)
    hdr = rows[0] if rows else []
    if "id" not in hdr or "path" not in hdr:
        errors.append(f"{genomes}: header must contain 'id' and 'path'")
    else:
        i, j, seen = hdr.index("id"), hdr.index("path"), set()
        for n, r in enumerate(rows[1:], start=2):
            if len(r) <= max(i, j):
                errors.append(f"{genomes}:{n}: too few columns")
                continue
            gid, gp = r[i], resolve(r[j], genomes.parent)
            if gid in seen:
                errors.append(f"{genomes}:{n}: duplicate id '{gid}'")
            seen.add(gid)
            if not re.fullmatch(r"[A-Za-z0-9_-]+", gid):
                warnings.append(f"genome id '{gid}' has characters other than [A-Za-z0-9_-]; "
                                "this can confuse snakemake wildcards")
            if not gp.is_file():
                errors.append(f"{genomes}:{n}: genome not found: {gp}")
            g_out.append([gid, str(gp)])
        if not g_out:
            errors.append(f"{genomes}: no genomes listed")

    rows = _read_csv(samples)
    if len(rows) < 2 or len(rows[0]) < 2:
        errors.append(f"{samples}: need a header plus rows with 2 columns (sample_id,dir_path)")
    else:
        s_hdr, seen = rows[0][:2], set()
        fq_re = re.compile(r"(\.fq\.gz|\.fastq\.gz)$")
        for n, r in enumerate(rows[1:], start=2):
            if len(r) < 2:
                errors.append(f"{samples}:{n}: too few columns")
                continue
            sid, d = r[0], resolve(r[1], samples.parent)
            if sid in seen:
                errors.append(f"{samples}:{n}: duplicate sample '{sid}'")
            seen.add(sid)
            if d.is_file():
                if not fq_re.search(d.name):
                    errors.append(f"{samples}:{n}: file must end in .fq.gz or .fastq.gz: {d}")
            elif d.is_dir():
                if not any(fq_re.search(f.name) for f in d.rglob("*") if f.is_file()):
                    errors.append(f"{samples}:{n}: no *.fq.gz / *.fastq.gz under {d}")
            else:
                errors.append(f"{samples}:{n}: not found: {d}")
            s_out.append([sid, str(d)])

    if errors:
        return None, None, errors, warnings

    inputs = outdir / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    g_path, s_path = inputs / "genomes.csv", inputs / "samples.csv"
    for path, header, body in ((g_path, ["id", "path"], g_out), (s_path, s_hdr, s_out)):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(body)
    return g_path, s_path, errors, warnings


# Pipeline Function: Build Snakemake Command
def build_command(args, cfg_path: Path, outdir: Path, threads: int, sm_exe: str, sm_help: str):
    cmd = [sm_exe, "--snakefile", str(SNAKEFILE), "--directory", str(outdir),
           "--cores", str(threads), "--configfile", str(cfg_path),
           "--rerun-incomplete", "--printshellcmds"]

    methods = (["conda"] if args.use_conda else []) + (["apptainer"] if args.use_container else [])
    if methods:
        if "--software-deployment-method" in sm_help:       # newer snakemake
            cmd += ["--software-deployment-method", *methods]
        else:
            if args.use_conda:
                cmd.append("--use-conda")
            if args.use_container:
                cmd.append("--use-apptainer" if "--use-apptainer" in sm_help
                           else "--use-singularity")

    if args.unlock:
        cmd.append("--unlock")
    if args.dry_run:
        cmd.append("--dry-run")
    if args.force:
        cmd.append("--forceall")
    if args.keep_going:
        cmd.append("--keep-going")
    if args.retries:
        cmd += ["--retries" if "--retries" in sm_help else "--restart-times", str(args.retries)]
    if args.ignore_changes:
        cmd += ["--rerun-triggers", "mtime"]
    return cmd + args.snakemake_args


# Pipeline Function: Run Snakemake and Stream Output
def run_snakemake(cmd, quiet: bool):
    """Stream snakemake output to console (unless quiet) and the run log."""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                            bufsize=1, env={**os.environ, "PYTHONUNBUFFERED": "1"})
    locked = False
    try:
        for line in proc.stdout:
            line = line.rstrip("\n")
            sm_log.info(line)
            if not quiet:
                print(line, flush=True)
            if "LockException" in line or "Directory cannot be locked" in line:
                locked = True
        return proc.wait(), locked
    except KeyboardInterrupt:
        # Ctrl-C already reached snakemake (same process group); just wait for it to clean up.
        log.warning("Interrupted. Waiting for snakemake to stop running jobs ...")
        try:
            proc.wait(timeout=60)
        except subprocess.TimeoutExpired:
            proc.terminate()
            proc.wait()
        return 130, False


# Main Function: Parse Args, Validate Inputs, Run Snakemake
def main(argv=None) -> int:
    defaults = load_defaults()
    parser = build_parser(defaults)
    args = parser.parse_args(argv)

    # ---- merge config: packaged defaults < --configfile < CLI flags
    user_cfg = {}
    if args.configfile:
        if not args.configfile.is_file():
            parser.error(f"config file not found: {args.configfile}")
        user_cfg = yaml.safe_load(args.configfile.read_text()) or {}
    cfg = {**defaults, **user_cfg}
    for key in PARAMS:
        if getattr(args, key) is not None:
            cfg[key] = getattr(args, key)
    if args.threads is not None:
        cfg["threads"] = args.threads
    if args.memory is not None:
        cfg["max_mem"] = args.memory
    try:
        cfg["max_mem"] = parse_mem_gb(cfg["max_mem"])
    except argparse.ArgumentTypeError as e:
        parser.error(str(e))
    cfg["threads"] = int(cfg["threads"])
    if cfg["threads"] < 1:
        parser.error("threads must be >= 1")

    genomes = args.genomes or user_cfg.get("genomes")
    samples = args.samples or user_cfg.get("samples")
    if not genomes or not samples:
        parser.error("genomes (-g) and samples (-s) are required "
                     "(on the command line or in --configfile)")
    genomes, samples = Path(genomes).expanduser().resolve(), Path(samples).expanduser().resolve()
    for f in (genomes, samples):
        if not f.is_file():
            parser.error(f"file not found: {f}")

    # ---- logging
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    run_log = setup_logging(outdir / "logs" / "pang2p", args.verbose, args.quiet)
    log.info("pang2p %s | python %s | %s", __version__, sys.version.split()[0], sys.platform)
    log.debug("command line: %s", " ".join(sys.argv))
    log.info("run log: %s", run_log)

    unknown = set(user_cfg) - set(defaults) - {"genomes", "samples"}
    if unknown:
        log.warning("unknown keys in %s (ignored by the workflow): %s",
                    args.configfile, ", ".join(sorted(unknown)))

    # ---- choose a software environment if none was requested
    if not args.use_conda and not args.use_container and not args.skip_checks:
        if has_apptainer():
            args.use_container = True
        elif has_conda():
            args.use_conda = True
        if args.use_container or args.use_conda:
            log.info("no environment flag given; using %s",
                     "apptainer container" if args.use_container else "conda")

    # ---- dependency preflight
    sm_exe, sm_ver, sm_help = snakemake_info()
    if sm_exe:
        log.info("snakemake %s (%s)", sm_ver, sm_exe)
    if not args.skip_checks:
        errors = preflight(args.use_conda, args.use_container, bool(sm_exe))
        if errors:
            for e in errors:
                log.error(e)
            return 2
    elif not sm_exe:
        log.error("snakemake not found on PATH")
        return 2

    # ---- inputs
    g_csv, s_csv, errors, warnings = prepare_inputs(genomes, samples, outdir)
    for w in warnings:
        log.warning(w)
    if errors:
        for e in errors[:20]:
            log.error(e)
        if len(errors) > 20:
            log.error("... and %d more", len(errors) - 20)
        return 2
    cfg["genomes"], cfg["samples"] = str(g_csv), str(s_csv)

    cfg_path = outdir / "pang2p_config.yaml"
    if cfg_path.exists():
        old = yaml.safe_load(cfg_path.read_text()) or {}
        changed = sorted(k for k in set(old) | set(cfg) if old.get(k) != cfg.get(k)
                         and k not in ("genomes", "samples", "threads", "max_mem"))
        if changed:
            log.warning("settings differ from the previous run in this outdir: %s. "
                        "Affected steps will rerun (see --ignore-changes / --force).",
                        ", ".join(changed))
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    log.info("effective config: %s", cfg_path)

    # ---- run
    cmd = build_command(args, cfg_path, outdir, cfg["threads"], sm_exe, sm_help)
    log.info("running: %s", " ".join(cmd))
    t0 = time.monotonic()
    rc, locked = run_snakemake(cmd, args.quiet)
    elapsed = str(dt.timedelta(seconds=int(time.monotonic() - t0)))

    if rc == 0:
        log.info("finished in %s", elapsed)
    else:
        log.error("snakemake exited with code %d after %s", rc, elapsed)
        log.error("full log: %s ; per-rule logs: %s", run_log, outdir / "logs")
        if locked:
            log.error("the directory is locked by a previous run: "
                      "if none is running, use 'pang2p ... --unlock' and re-run.")
        elif rc != 130:
            log.error("fix the problem and re-run the same command to resume.")
    return rc

# Entry point for CLI
if __name__ == "__main__":
    sys.exit(main())
