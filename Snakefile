from snakemake.utils import validate, min_version
import pandas as pd
from os import path
from pathlib import Path

BASEDIR = Path(workflow.basedir)
SCRIPT_DIR = BASEDIR / "scripts"
ENV_DIR = BASEDIR / "envs"

configfile: "config.yaml"

genome_df = pd.read_csv(config["genomes"], sep=',').set_index("id").sort_index()
genome_map = genome_df["path"].to_dict()

include: "rules/01.single_copy.smk"
include: "rules/02.clustering.smk"
include: "rules/03.pankmers.smk"
include: "rules/04.genotyping.smk"


rule all:
    input:
        "work/01.single_copy/.done",
        "work/02.clustering/.done",
        "work/03.pankmers/.done",
        "work/04.genotyping/.done"
