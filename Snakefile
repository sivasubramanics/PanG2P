from snakemake.utils import validate, min_version
import pandas as pd
from os import path
from pathlib import Path

BASEDIR = Path(workflow.basedir)
SCRIPT_DIR = BASEDIR / "workflow" / "scripts"
ENV_DIR = BASEDIR / "workflow" / "envs"
RULES_DIR = BASEDIR / "workflow" / "rules"

configfile: "config.yaml"

genome_df = pd.read_csv(config["genomes"], sep=',').set_index("id").sort_index()
genome_map = genome_df["path"].to_dict()

include: RULES_DIR / "01.single_copy.smk"
include: RULES_DIR / "02.clustering.smk"
include: RULES_DIR / "03.pankmers.smk"
include: RULES_DIR / "04.genotyping.smk"


rule all:
    input:
        "work/01.single_copy/.done",
        "work/02.clustering/.done",
        "work/03.pankmers/.done",
        "work/04.genotyping/.done"
