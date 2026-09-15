from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

# file paths
ROOT = Path(__file__).resolve().parents[1]
GO2SUM_DIR = ROOT / "Data" / "raw" / "go2sum"
RESULTS_DIR = ROOT / "Results"
FIGURES_DIR = ROOT / "Figures"

INTEGRATED_FILE = RESULTS_DIR / "21_integrated_annotation_poor_predictions.csv"
PER_PROTEIN_FILE = RESULTS_DIR / "23_go2sum_integrated_analysis.csv"
TYPE_SUMMARY_FILE = RESULTS_DIR / "23_go2sum_output_type_summary.csv"
FLAG_SUMMARY_FILE = RESULTS_DIR / "23_go2sum_quality_flag_summary.csv"
FLAGGED_FILE = RESULTS_DIR / "23_go2sum_flagged_outputs.csv"
CORRELATION_FILE = RESULTS_DIR / "23_go2sum_correlation_summary.csv"
FIGURE_FILE = FIGURES_DIR / "23_go2sum_quality_flags.png"

OUTPUT_TYPES = ("function", "pathway", "subunit")

# lexical patterns used only to flag outputs for manual review
REVIEW_RULES = {
    "obsolete_placeholder": r"\bobsolete\b|obisoluteute",
    "photosynthesis_or_chloroplast": r"\bphotosynth\w*\b|\bchloroplast\w*\b|\bplastid\w*\b|\bC4 acid pathway\b",
    "metazoan_or_clinical_context": r"\b(?:human|mammal\w*|mucosal|allergic|allergy|embryo\w*|neuron\w*|brain|blood|thyroid|tumou?r|cancer)\b",
    "bacterial_or_viral_context": r"\b(?:bacteri\w*|gram-positive|gram-negative|viral|virus|phage|host immune|host cell)\b",
}


def read_go2sum(kind):
    path = GO2SUM_DIR / f"{kind}_merged.tab"
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["Protein", "Description"],
        dtype=str, keep_default_na=False,
    )

    if df["Protein"].duplicated().any():
        raise ValueError(f"Duplicate protein IDs in {path.name}")

    df["Description"] = (
        df["Description"]
        .str.replace(rf"^{kind.upper()}:\s*\*?\s*", "", regex=True, flags=re.I)
        .str.strip()
    )
    return df


def correlation_row(data, predictor, metric, kind):
    usable = data[[predictor, metric]].apply(pd.to_numeric, errors="coerce").dropna()
    rho, p_value = spearmanr(usable[predictor], usable[metric])

    return {
        "Output_Type": kind,
        "ProtNote_Metric": predictor,
        "Text_Metric": metric,
        "N_used": len(usable),
        "Spearman_rho": rho,
        "P_value": p_value,
    }


def main():
    if not INTEGRATED_FILE.exists():
        raise FileNotFoundError(INTEGRATED_FILE)

    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)

    integrated = pd.read_csv(INTEGRATED_FILE)
    if len(integrated) != 521 or integrated["Protein"].nunique() != 521:
        raise ValueError("Integrated input must contain 521 unique proteins")

    expected_ids = set(integrated["Protein"])
    result = integrated.copy()
    long_tables = []
    type_rows = []

    for kind in OUTPUT_TYPES:
        raw = read_go2sum(kind)

        if set(raw["Protein"]) != expected_ids:
            raise ValueError(f"{kind} GO2SUM IDs do not match the integrated dataset")

        # measure length and exact text repetition
        raw["Word_Count"] = raw["Description"].str.split().str.len()
        normalised = raw["Description"].str.casefold().str.replace(r"\s+", " ", regex=True).str.strip()
        raw["Exact_Duplicate_Count"] = normalised.map(normalised.value_counts())

        # very long outputs are flagged using a conservative 3 x IQR rule
        q1, q3 = raw["Word_Count"].quantile([0.25, 0.75])
        length_threshold = q3 + 3 * (q3 - q1)
        extreme = raw["Word_Count"] > length_threshold

        reasons = []
        for text, is_extreme in zip(raw["Description"], extreme):
            matched = [name for name, pattern in REVIEW_RULES.items() if re.search(pattern, text, re.I)]
            if is_extreme:
                matched.append("extreme_length_outlier")
            reasons.append("; ".join(matched))

        raw["Review_Reasons"] = reasons
        raw["Review_Flag"] = raw["Review_Reasons"].ne("")
        raw["Output_Type"] = kind
        long_tables.append(raw.copy())

        # add each output type to the integrated protein table
        prefix = f"GO2SUM_{kind.title()}"
        wide = raw.rename(columns={
            "Description": f"{prefix}_Text",
            "Word_Count": f"{prefix}_Word_Count",
            "Exact_Duplicate_Count": f"{prefix}_Exact_Duplicate_Count",
            "Review_Reasons": f"{prefix}_Review_Reasons",
            "Review_Flag": f"{prefix}_Review_Flag",
        }).drop(columns="Output_Type")

        result = result.merge(wide, on="Protein", how="left", validate="one_to_one")

        type_rows.append({
            "Output_Type": kind,
            "N_missing": int(raw["Description"].eq("").sum()),
            "N_unique_exact_descriptions": int(normalised.nunique()),
            "Percent_unique_exact_descriptions": normalised.nunique() / len(raw) * 100,
            "Median_word_count": raw["Word_Count"].median(),
            "Maximum_word_count": raw["Word_Count"].max(),
            "N_extreme_length_outliers": int(extreme.sum()),
            "N_any_review_flag": int(raw["Review_Flag"].sum()),
            "Percent_any_review_flag": raw["Review_Flag"].mean() * 100,
        })

    long = pd.concat(long_tables, ignore_index=True)
    type_summary = pd.DataFrame(type_rows)

    flag_columns = [f"GO2SUM_{kind.title()}_Review_Flag" for kind in OUTPUT_TYPES]
    result["GO2SUM_Any_Review_Flag"] = result[flag_columns].any(axis=1)
    result["GO2SUM_Number_Of_Flagged_Output_Types"] = result[flag_columns].sum(axis=1)

    # count each review reason by GO2SUM output type
    flag_rows = []
    for kind in OUTPUT_TYPES:
        subset = long[long["Output_Type"] == kind]

        for flag_name in [*REVIEW_RULES, "extreme_length_outlier"]:
            matched = subset["Review_Reasons"].str.contains(
                rf"(?:^|; ){re.escape(flag_name)}(?:;|$)", regex=True
            )
            flag_rows.append({
                "Output_Type": kind,
                "Quality_Flag": flag_name,
                "N_flagged": int(matched.sum()),
                "Percent_of_521": matched.mean() * 100,
            })

    flag_summary = pd.DataFrame(flag_rows)

    # keep model context beside outputs selected for manual review
    metadata = [
        column for column in [
            "Protein", "ProtNLM_Top1", "ProtNLM_Top1_Score",
            "ProtNote_Top_GO_Description", "ProtNote_Top_Probability",
            "ProtNote_Specific_Count", "ProtNote_HighConf_Specific_Count",
            "Ortholog count", "Protein Length",
        ]
        if column in integrated.columns
    ]

    flagged = (
        long[long["Review_Flag"]]
        .merge(integrated[metadata], on="Protein", how="left")
        .sort_values(["Output_Type", "Review_Reasons", "Protein"])
    )

    # exploratory links between ProtNote input and GO2SUM wording
    correlation_rows = []
    predictors = [
        "ProtNote_Specific_Count",
        "ProtNote_HighConf_Specific_Count",
        "ProtNote_Top_Probability",
    ]

    for kind in OUTPUT_TYPES:
        metrics = [
            f"GO2SUM_{kind.title()}_Word_Count",
            f"GO2SUM_{kind.title()}_Exact_Duplicate_Count",
        ]
        for predictor in predictors:
            if predictor in result.columns:
                for metric in metrics:
                    correlation_rows.append(
                        correlation_row(result, predictor, metric, kind)
                    )

    correlations = pd.DataFrame(correlation_rows)

    result.to_csv(PER_PROTEIN_FILE, index=False)
    type_summary.to_csv(TYPE_SUMMARY_FILE, index=False)
    flag_summary.to_csv(FLAG_SUMMARY_FILE, index=False)
    flagged.to_csv(FLAGGED_FILE, index=False)
    correlations.to_csv(CORRELATION_FILE, index=False)

    # plot review-rule frequencies
    plot_data = (
        flag_summary
        .pivot(index="Quality_Flag", columns="Output_Type", values="Percent_of_521")
        .reindex(columns=OUTPUT_TYPES)
    )

    ax = plot_data.plot.bar(
        figsize=(9, 5.4),
        color=["#286D8E", "#D18F3B", "#6A8F5B"],
    )
    ax.set_ylabel("Proteins flagged (%)")
    ax.set_xlabel("")
    ax.set_xticklabels(
        [label.replace("_", " ") for label in plot_data.index],
        rotation=25, ha="right",
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(title="GO2SUM output")

    plt.tight_layout()
    plt.savefig(FIGURE_FILE, dpi=300, bbox_inches="tight")
    plt.savefig(FIGURE_FILE.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    print(f"proteins analysed: {len(result)}")
    print(f"proteins with review flags: {int(result['GO2SUM_Any_Review_Flag'].sum())}")
    print(f"flagged output rows: {len(flagged)} / {len(long)}")
    print("\noutput type summary:")
    print(type_summary.to_string(index=False))


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used extensively for lines 91-135 because i had not previously carried
# out quality checks on generated biological text. i asked how i could measure
# repeated descriptions, identify unusually long outputs and apply several
# transparent lexical review rules without treating the flags as proof of an
# incorrect annotation. GAI gave examples including:
#
# counts = cleaned_text.value_counts()
# df["duplicate_count"] = cleaned_text.map(counts)
#
# q1, q3 = df["word_count"].quantile([0.25, 0.75])
# threshold = q3 + 3 * (q3 - q1)
# outlier = df["word_count"] > threshold
#
# matched = [
#     name for name, pattern in rules.items()
#     if re.search(pattern, text, re.I)
# ]
#
# i adapted these approaches for GO2SUM function, pathway and subunit descriptions.
# the resulting rules are only used to identify outputs for manual review.

# lines 145-160 were also developed with GAI support. i asked how to count each
# individual review reason when several reasons could be stored together in one
# semicolon-separated field. GAI suggested looping through the known rule names
# and using str.contains() with a boundary-aware regular expression. i adapted
# this to create a separate count and percentage for every review rule and every
# GO2SUM output type.

# GAI was also used for lines 186-197 when setting up the exploratory correlation
# analysis. i had three ProtNote variables, three GO2SUM output types and two text
# properties, and asked how to avoid writing every Spearman test separately.
# ChatGPT suggested nested loops over predictors and output metrics, calling the
# same correlation function each time. i used that structure here to compare
# ProtNote prediction properties with GO2SUM word count and description
# repetition. these correlations describe generated-text behaviour only and are
# not treated as measurements of biological correctness.