#!/usr/bin/env python3
# Parse benchmark files from snakemake output
# usage: parse_benchmarks.py <benchmark_dir>
# the rules are expected to be in this order benchmark_dir/rule_A/rule_a/sample.tsv

from pathlib import Path

import click
import pandas as pd


def get_time(benchmark_file: Path) -> float:
    """Return the benchmark time from a Snakemake benchmark file."""
    with benchmark_file.open() as f:
        next(f)  # skip header
        line = next(f).strip()
    return float(line.split()[0])


def get_rules() -> list[str]:
    # get list of rules in rules.order where the script is located
    script_dir = Path(__file__).parent
    rules_order = script_dir / "rules.order"
    if not rules_order.exists():
        print(f"Warning: {rules_order} does not exist")
        return []
    with rules_order.open() as f:
        return [line.strip() for line in f if line.strip()]


@click.command()
@click.option(
    "-d",
    "--bdir",
    default=Path("benchmarks"),
    type=click.Path(exists=True, path_type=Path),
    help="Directory containing benchmark files",
)
@click.option(
    "--outprefix",
    "-o",
    type=click.Path(path_type=Path),
    help="Output file prefix to write parsed benchmark results",
)
def main(bdir: Path, outprefix: Path | None):

    rows = []
    rule_order = {}

    for benchmark in sorted(bdir.rglob("*.tsv")):
        rel = benchmark.relative_to(bdir)
        parts = rel.parts

        # Must have at least layer_one/file.tsv
        if len(parts) < 2:
            continue

        layer_one = parts[0]

        if len(parts) == 2:
            # File directly under layer_one
            layer_two = benchmark.stem
            layer_three = ""
        else:
            layer_two = parts[1]

            remainder = list(parts[2:])
            remainder[-1] = Path(remainder[-1]).stem
            layer_three = "/".join(remainder)

        rows.append(
            {
                "layer_one": layer_one,
                "layer_two": layer_two,
                "layer_three": layer_three,
                "time": get_time(benchmark),
            }
        )

    rules = get_rules()
    rule_order = {r: i for i, r in enumerate(rules)}

    df = pd.DataFrame(rows)

    df["layer_two"] = pd.Categorical(df["layer_two"], categories=rules, ordered=True)

    # order layer_two by rules if available
    if rules:
        df["layer_two_order"] = df["layer_two"].map(rule_order)

    df = df.sort_values(
        by=["layer_one", "layer_two_order", "layer_three"], ascending=[True, True, True]
    )

    if outprefix:
        df.to_csv(outprefix.with_suffix(".tsv"), sep="\t", index=False)
    else:
        print(df.to_string(index=False))

    plot_df = df.groupby(["layer_one", "layer_two"], as_index=False, observed=False)[
        "time"
    ].sum()

    pivot = plot_df.pivot(index="layer_two", columns="layer_one", values="time").fillna(
        0
    )

    pivot = pivot.reindex(rules)

    import matplotlib.pyplot as plt

    ax = pivot.plot(kind="bar", figsize=(12, 6), width=0.85)

    # ax.set_xticks(np.arange(len(pivot.index)))
    ax.set_xticklabels(pivot.index, rotation=45, ha="right")

    ax.set_xlabel("Rule (layer_two)")
    ax.set_ylabel("Time (seconds)")
    ax.set_title("Benchmark grouped by pipeline stage")

    plt.tight_layout()
    plt.savefig(outprefix.with_suffix(".png"), dpi=300)


if __name__ == "__main__":
    main()
