from snakemake.utils import validate, min_version
import pandas as pd
from os import path

configfile: "config.yaml"

genome_df = pd.read_csv(config["genomes"], sep=',').set_index("id").sort_index()
genome_map = genome_df["path"].to_dict()
print(genome_map) # remove or comment-out in production

include: "rules/01.single_copy.smk"
include: "rules/02.clustering.smk"

rule all:
    input:
        "work/01.single_copy/.done",
        "work/02.clustering/.done"
