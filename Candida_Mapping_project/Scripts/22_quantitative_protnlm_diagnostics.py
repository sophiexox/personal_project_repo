from pathlib import Path
import re
import unicodedata
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

# file paths
ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "Results"
FIGURES_DIR = ROOT / "Figures"

INTEGRATED_FILE = RESULTS_DIR / "21_integrated_annotation_poor_predictions.csv"
RAW_FILE = ROOT / "Data" / "raw" / "protnlm" / "ProtNLM_annotation_poor_521_predictions.csv"

DIAGNOSTIC_FILE = RESULTS_DIR / "22_protnlm_quantitative_diagnostics.csv"
CORRELATION_FILE = RESULTS_DIR / "22_protnlm_correlation_summary.csv"
BINNED_FILE = RESULTS_DIR / "22_protnlm_binned_summary.csv"
GROUP_FILE = RESULTS_DIR / "22_protnlm_generic_comparison.csv"
EXTREMES_FILE = RESULTS_DIR / "22_protnlm_manual_review_extremes.csv"
SUMMARY_FILE = RESULTS_DIR / "22_protnlm_quantitative_summary.md"

EXPECTED_PROTEINS = 521
EXPECTED_BEAMS = 10
REVIEW_N = 5


def normalise_prediction(text):
    # conservative cleanup for exact lexical comparison
    value = unicodedata.normalize("NFKC", str(text)).casefold().strip()
    value = re.sub(r"[_\-‐-‒–—]+", " ", value)
    value = re.sub(r"^[\s\.,;:|/\\]+|[\s\.,;:|/\\]+$", "", value)
    return re.sub(r"\s+", " ", value).strip()


def as_bool(series):
    return series.fillna(False).astype(str).str.strip().str.casefold().eq("true")


def validate_inputs(integrated, raw):
    integrated_needed = {
        "Protein", "Protein Length", "Ortholog count", "Paralog count",
        "ProtNLM_Top1", "ProtNLM_Top1_Score", "ProtNLM_Top1_Generic",
        "ProtNLM_Top10",
    }
    raw_needed = {"Protein", "Rank", "Prediction", "Score"}

    for name, frame, needed in [
        ("integrated", integrated, integrated_needed),
        ("raw ProtNLM", raw, raw_needed),
    ]:
        missing = needed - set(frame.columns)
        if missing:
            raise KeyError(f"Missing {name} columns: {', '.join(sorted(missing))}")

    if integrated["Protein"].duplicated().any():
        raise ValueError("Integrated file must contain one row per protein")

    if integrated["Protein"].nunique() != EXPECTED_PROTEINS:
        raise ValueError(f"Expected {EXPECTED_PROTEINS} proteins")

    if (raw.groupby("Protein").size() != EXPECTED_BEAMS).any():
        raise ValueError("Every protein must have 10 ProtNLM predictions")

    if raw.duplicated(["Protein", "Rank"]).any():
        raise ValueError("Duplicate Protein/Rank combinations found")

    if set(raw["Protein"]) != set(integrated["Protein"]):
        raise ValueError("Protein IDs differ between integrated and raw files")


def build_beam_diagnostics(raw):
    # summarise repeated wording across each protein's 10 predictions
    beams = raw.sort_values(["Protein", "Rank"]).copy()
    beams["Normalised_Prediction"] = beams["Prediction"].map(normalise_prediction)

    rows = []

    for protein, group in beams.groupby("Protein", sort=False):
        labels = group["Normalised_Prediction"].tolist()
        counts = Counter(labels)
        top_count = max(counts.values())

        repeated = sorted(
            label
            for label, count in counts.items()
            if count == top_count
        )

        rows.append({
            "Protein": protein,
            "Beam_Count": len(labels),
            "Beam_Unique_Normalised_Labels": len(counts),
            "Beam_Top_Label_Repetition_Count": top_count,
            "Beam_Exact_Lexical_Coherence": top_count / len(labels),
            "Beam_Most_Repeated_Normalised_Label": " | ".join(repeated),
            "Beam_Normalised_Labels": " | ".join(labels),
        })

    return pd.DataFrame(rows)


def correlation_row(data, predictor):
    usable = data[[predictor, "ProtNLM_Top1_Score"]].dropna()

    rho, p_value = spearmanr(
        usable[predictor],
        usable["ProtNLM_Top1_Score"],
    )

    return {
        "Outcome": "ProtNLM_Top1_Score",
        "Predictor": predictor,
        "Method": "Spearman rank correlation",
        "N_used": len(usable),
        "N_total": len(data),
        "Missing_predictor": int(data[predictor].isna().sum()),
        "Missing_outcome": int(data["ProtNLM_Top1_Score"].isna().sum()),
        "Spearman_rho": rho,
        "P_value": p_value,
    }


def add_quantile_bin(data, column, output):
    # quartiles come from the observed data rather than fixed cut-offs
    valid = data[column].dropna()
    bins = pd.qcut(valid, q=4, duplicates="drop").astype(str)
    data.loc[valid.index, output] = bins.to_numpy()


def binned_summary(data, bin_column, variable):
    rows = []

    for label, group in data.dropna(subset=[bin_column]).groupby(
        bin_column,
        sort=False,
    ):
        scores = group["ProtNLM_Top1_Score"].dropna()
        q1, q3 = scores.quantile([0.25, 0.75])

        rows.append({
            "Binning_Variable": variable,
            "Bin": str(label),
            "Count": len(scores),
            "Observed_Min": group[variable].min(),
            "Observed_Max": group[variable].max(),
            "Median_ProtNLM_Score": scores.median(),
            "Mean_ProtNLM_Score": scores.mean(),
            "Q1_ProtNLM_Score": q1,
            "Q3_ProtNLM_Score": q3,
            "IQR_ProtNLM_Score": q3 - q1,
        })

    return pd.DataFrame(rows)


def suspicious_reason(label):
    # narrow lexical flags for manual review only
    text = str(label).strip()
    reasons = []

    if re.match(r"^\([^)]{2,80}\)\s+", text):
        reasons.append("parenthetical organism-like prefix")

    if re.search(r"\bWGS\s+project\b|\bcontig\b|\bscaffold\b", text, re.I):
        reasons.append("WGS/contig/scaffold wording")

    return "; ".join(reasons)


def select_manual_review(data):
    score_hi = data["ProtNLM_Top1_Score"].quantile(0.75)
    score_lo = data["ProtNLM_Top1_Score"].quantile(0.25)
    coherence_hi = data["Beam_Exact_Lexical_Coherence"].quantile(0.75)
    coherence_lo = data["Beam_Exact_Lexical_Coherence"].quantile(0.25)

    definitions = [
        (
            "high score + high lexical coherence",
            (data["ProtNLM_Top1_Score"] >= score_hi)
            & (data["Beam_Exact_Lexical_Coherence"] >= coherence_hi),
            ["ProtNLM_Top1_Score", "Beam_Exact_Lexical_Coherence"],
            [False, False],
        ),
        (
            "high score + low lexical coherence",
            (data["ProtNLM_Top1_Score"] >= score_hi)
            & (data["Beam_Exact_Lexical_Coherence"] <= coherence_lo),
            ["ProtNLM_Top1_Score", "Beam_Exact_Lexical_Coherence"],
            [False, True],
        ),
        (
            "low score + high lexical coherence",
            (data["ProtNLM_Top1_Score"] <= score_lo)
            & (data["Beam_Exact_Lexical_Coherence"] >= coherence_hi),
            ["Beam_Exact_Lexical_Coherence", "ProtNLM_Top1_Score"],
            [False, True],
        ),
    ]

    selections = []

    for category, mask, sort_columns, ascending in definitions:
        chosen = (
            data[mask]
            .sort_values(sort_columns, ascending=ascending)
            .head(REVIEW_N)
            .copy()
        )

        chosen["Review_Category"] = category
        selections.append(chosen)

    suspicious = data[data["Suspicious_Output_Rule_Match"]].copy()
    suspicious = suspicious.sort_values(
        "ProtNLM_Top1_Score",
        ascending=False,
    )
    suspicious["Review_Category"] = "taxonomic/database-like lexical output"
    selections.append(suspicious)

    columns = [
        "Review_Category",
        "Protein",
        "ProtNLM_Top1",
        "ProtNLM_Top1_Score",
        "Beam_Exact_Lexical_Coherence",
        "Beam_Unique_Normalised_Labels",
        "Beam_Top_Label_Repetition_Count",
        "Ortholog count",
        "Paralog count",
        "Protein Length",
        "ProtNLM_Top1_Generic",
        "Suspicious_Output_Reason",
        "ProtNLM_Top10",
    ]

    return pd.concat(selections, ignore_index=True)[columns]


def generic_comparison(data):
    rows = []

    for generic, group in data.groupby("ProtNLM_Top1_Generic"):
        scores = group["ProtNLM_Top1_Score"].dropna()
        q1, q3 = scores.quantile([0.25, 0.75])

        rows.append({
            "Comparison": "generic vs non-generic top-1 label",
            "Group": "Generic" if generic else "Non-generic",
            "N": len(scores),
            "Mean_ProtNLM_Score": scores.mean(),
            "Median_ProtNLM_Score": scores.median(),
            "Q1_ProtNLM_Score": q1,
            "Q3_ProtNLM_Score": q3,
            "IQR_ProtNLM_Score": q3 - q1,
        })

    generic_scores = data.loc[
        data["ProtNLM_Top1_Generic"],
        "ProtNLM_Top1_Score",
    ].dropna()

    nongeneric_scores = data.loc[
        ~data["ProtNLM_Top1_Generic"],
        "ProtNLM_Top1_Score",
    ].dropna()

    u_stat, p_value = mannwhitneyu(
        generic_scores,
        nongeneric_scores,
        alternative="two-sided",
    )

    rank_biserial = (
        2 * u_stat
        / (len(generic_scores) * len(nongeneric_scores))
        - 1
    )

    result = pd.DataFrame(rows)
    result["Mann_Whitney_U"] = u_stat
    result["P_value"] = p_value
    result["Rank_Biserial_Generic_vs_NonGeneric"] = rank_biserial

    return result


def make_figures(data):
    FIGURES_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(22)

    specs = [
        (
            "Ortholog count",
            "Ortholog count (log1p scale)",
            "22_protnlm_score_vs_ortholog_count.png",
            True,
        ),
        (
            "Protein Length",
            "Protein length (amino acids)",
            "22_protnlm_score_vs_protein_length.png",
            False,
        ),
        (
            "Beam_Exact_Lexical_Coherence",
            "Exact lexical beam coherence (max repeated label / 10)",
            "22_protnlm_score_vs_lexical_beam_coherence.png",
            False,
        ),
    ]

    for column, label, filename, log_transform in specs:
        fig, ax = plt.subplots(figsize=(6.4, 4.5))

        x = data[column].astype(float)

        if log_transform:
            x = np.log1p(x)

        # fixed jitter only helps separate repeated x values visually
        if data[column].nunique() < 80:
            span = max(float(x.max() - x.min()), 1.0)
            x = x + rng.normal(0, span * 0.004, len(x))

        ax.scatter(
            x,
            data["ProtNLM_Top1_Score"],
            s=24,
            alpha=0.58,
            color="#286D8E",
            edgecolors="none",
        )

        ax.set_xlabel(label)
        ax.set_ylabel("ProtNLM top-1 score")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if log_transform:
            ticks = np.array([0, 1, 5, 10, 25, 100, 500, 2000])
            ticks = ticks[ticks <= data[column].max()]

            ax.set_xticks(
                np.log1p(ticks),
                labels=[str(value) for value in ticks],
            )

        path = FIGURES_DIR / filename

        fig.tight_layout()
        fig.savefig(path, dpi=300, bbox_inches="tight")
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


def format_p(value):
    return f"{value:.3g}" if value >= 0.001 else f"{value:.2e}"


def write_summary(diagnostic, correlations, group, review):
    correlation_lines = [
        f"- {row.Predictor}: rho = {row.Spearman_rho:.3f}, "
        f"p = {format_p(row.P_value)}, N = {row.N_used}"
        for row in correlations.itertuples(index=False)
    ]

    generic = group[group["Group"] == "Generic"].iloc[0]
    nongeneric = group[group["Group"] == "Non-generic"].iloc[0]

    summary = f"""# ProtNLM quantitative diagnostics

This exploratory analysis covers {len(diagnostic)} annotation-poor proteins.
Associations do not demonstrate model mechanism, training-data leakage or biological
correctness. Lexical beam coherence measures repeated wording only.

## Spearman correlations with ProtNLM top-1 score

{chr(10).join(correlation_lines)}

## Generic versus non-generic labels

- Generic: N = {int(generic['N'])}, median = {generic['Median_ProtNLM_Score']:.4f}.
- Non-generic: N = {int(nongeneric['N'])}, median = {nongeneric['Median_ProtNLM_Score']:.4f}.
- Mann-Whitney U = {generic['Mann_Whitney_U']:.1f}, p = {format_p(generic['P_value'])},
  rank-biserial = {generic['Rank_Biserial_Generic_vs_NonGeneric']:.3f}.

## Beam review

- Median unique normalised labels: {diagnostic['Beam_Unique_Normalised_Labels'].median():.1f}.
- Median exact lexical coherence: {diagnostic['Beam_Exact_Lexical_Coherence'].median():.2f}.
- Suspicious lexical rule matches: {int(diagnostic['Suspicious_Output_Rule_Match'].sum())}.
- Manual-review rows: {len(review)}.
"""

    SUMMARY_FILE.write_text(summary, encoding="utf-8")


def main():
    for path in [INTEGRATED_FILE, RAW_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    RESULTS_DIR.mkdir(exist_ok=True)

    integrated = pd.read_csv(INTEGRATED_FILE)
    raw = pd.read_csv(RAW_FILE)

    validate_inputs(integrated, raw)

    # restore numeric and boolean types after reading CSV files
    for column in [
        "Protein Length",
        "Ortholog count",
        "Paralog count",
        "ProtNLM_Top1_Score",
    ]:
        integrated[column] = pd.to_numeric(
            integrated[column],
            errors="coerce",
        )

    integrated["ProtNLM_Top1_Generic"] = as_bool(
        integrated["ProtNLM_Top1_Generic"]
    )

    beam = build_beam_diagnostics(raw)

    diagnostic = integrated.merge(
        beam,
        on="Protein",
        how="left",
        validate="one_to_one",
    )

    diagnostic["Suspicious_Output_Reason"] = (
        diagnostic["ProtNLM_Top1"]
        .map(suspicious_reason)
    )

    diagnostic["Suspicious_Output_Rule_Match"] = (
        diagnostic["Suspicious_Output_Reason"].ne("")
    )

    predictors = [
        "Ortholog count",
        "Paralog count",
        "Protein Length",
        "Beam_Exact_Lexical_Coherence",
        "Beam_Unique_Normalised_Labels",
    ]

    correlations = pd.DataFrame(
        [
            correlation_row(diagnostic, predictor)
            for predictor in predictors
        ]
    )

    add_quantile_bin(
        diagnostic,
        "Ortholog count",
        "Ortholog_Count_Quartile_Bin",
    )

    add_quantile_bin(
        diagnostic,
        "Protein Length",
        "Protein_Length_Quartile_Bin",
    )

    binned = pd.concat(
        [
            binned_summary(
                diagnostic,
                "Ortholog_Count_Quartile_Bin",
                "Ortholog count",
            ),
            binned_summary(
                diagnostic,
                "Protein_Length_Quartile_Bin",
                "Protein Length",
            ),
        ],
        ignore_index=True,
    )

    group = generic_comparison(diagnostic)
    review = select_manual_review(diagnostic)

    diagnostic.to_csv(DIAGNOSTIC_FILE, index=False)
    correlations.to_csv(CORRELATION_FILE, index=False)
    binned.to_csv(BINNED_FILE, index=False)
    group.to_csv(GROUP_FILE, index=False)
    review.to_csv(EXTREMES_FILE, index=False)

    make_figures(diagnostic)
    write_summary(diagnostic, correlations, group, review)

    print(f"proteins analysed: {len(diagnostic)}")

    print(
        "median lexical beam coherence: "
        f"{diagnostic['Beam_Exact_Lexical_Coherence'].median():.2f}"
    )

    print(
        "suspicious lexical outputs: "
        f"{int(diagnostic['Suspicious_Output_Rule_Match'].sum())}"
    )

    print("\ncorrelations:")

    print(
        correlations[
            ["Predictor", "Spearman_rho", "P_value", "N_used"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used extensively during the development of this script since the
# quantitative diagnostics went beyond the types of analysis i have previously
# written independently. i used GAI to help work out suitable pandas/scipy
# approaches, understand what the statistical outputs meant, and turn the
# analysis plan into working Python. the biological interpretation and decisions
# about which comparisons were appropriate for this project were then checked
# separately and kept deliberately cautious.

# for lines 31-35, i used GAI to help construct a conservative way of normalising
# ProtNLM text predictions before comparing the 10 beam outputs. i wanted labels
# such as "F-box", "F_box" and labels with different spacing or punctuation to be
# treated consistently, without doing semantic matching or assuming that similar
# biological phrases meant the same thing. GAI suggested using Unicode
# normalisation together with regular expressions, for example:
#
# value = unicodedata.normalize("NFKC", str(your_text)).lower().strip()
# value = re.sub(r"[_\-]+", " ", value)
# value = re.sub(r"\s+", " ", value)
#
# i adapted this into normalise_prediction(). the function only standardises
# formatting, so the resulting beam coherence measurement remains an exact
# lexical comparison rather than a biological similarity score.

# GAI was used heavily for lines 43-72 when adding input validation. i asked how
# to check that two related dataframes contained all required columns, that one
# table had exactly one row per protein, and that every protein in the raw
# ProtNLM output contained the expected number of ranked predictions. example
# code supplied by GAI included patterns such as:
#
# missing = required_columns - set(df.columns)
# if missing:
#     raise ValueError(...)
#
# counts = raw.groupby("your_id").size()
# if (counts != expected_count).any():
#     raise ValueError(...)
#
# if set(df1["your_id"]) != set(df2["your_id"]):
#     raise ValueError(...)
#
# these ideas were adapted into validate_inputs() so the quantitative analysis
# would stop rather than silently continue if the 521-protein dataset or the
# expected 10 ProtNLM beams per protein had been altered.

# lines 75-103 were one of the main sections for which i relied on GAI. i had not
# previously calculated within-protein consistency across multiple generated
# predictions and asked how the 10 beam labels could be summarised numerically.
# GAI suggested grouping predictions by protein, using Counter() to count the
# normalised labels, identifying the largest repetition count and dividing this
# by the total number of predictions. a simplified example was:
#
# counts = Counter(labels)
# most_common_count = max(counts.values())
# coherence = most_common_count / len(labels)
# unique_labels = len(counts)
#
# i adapted this to create Beam_Unique_Normalised_Labels,
# Beam_Top_Label_Repetition_Count and Beam_Exact_Lexical_Coherence. i also kept
# the normalised labels themselves so unusual beam behaviour could be inspected
# manually. GAI explanations were useful here for understanding that this is a
# measure of repeated wording only, not evidence that a prediction is correct.

# GAI also supported the statistical analysis in lines 106-158. i asked how to
# test for monotonic associations between ProtNLM score and variables such as
# ortholog count, paralog count, protein length and beam coherence. GAI explained
# the use of scipy.stats.spearmanr for a rank-based correlation and showed a
# general structure such as:
#
# usable = df[[x, y]].dropna()
# rho, p_value = spearmanr(usable[x], usable[y])
#
# i adapted this into correlation_row() and also recorded the number of usable
# observations and missing values. i separately asked how to divide continuous
# variables into data-derived quartiles and summarise the score distribution in
# each group. GAI suggested pd.qcut(..., q=4) and calculating median, mean,
# quartiles and IQR for each resulting bin. this became add_quantile_bin() and
# binned_summary(). these analyses were treated as exploratory associations and
# not as evidence of model mechanism or training-data leakage.

# lines 161-241 were developed with substantial GAI help for selecting examples
# for manual review. i wanted a reproducible way to find proteins showing
# contrasting combinations of ProtNLM score and beam coherence rather than
# manually choosing whichever examples looked interesting. GAI suggested using
# the 25th and 75th percentiles of each variable to define data-driven extremes,
# for example:
#
# score_high = df["score"].quantile(0.75)
# score_low = df["score"].quantile(0.25)
# coherence_high = df["coherence"].quantile(0.75)
#
# mask = (
#     (df["score"] >= score_high)
#     & (df["coherence"] <= coherence_low)
# )
#
# selected = df[mask].sort_values(...).head(n)
#
# i adapted this idea into the high-score/high-coherence,
# high-score/low-coherence and low-score/high-coherence review categories.
# GAI also helped with the regular-expression checks used by suspicious_reason()
# to flag narrow database-like or organism-like wording. these flags were used
# only to identify outputs worth inspecting and were not treated as evidence that
# the corresponding annotation was biologically wrong.

# GAI was used for lines 244-289 because i had not previously implemented a
# Mann-Whitney U comparison or rank-biserial effect size. i asked how two
# independent score distributions could be compared without assuming they were
# normally distributed. GAI provided an example using:
#
# u_stat, p_value = mannwhitneyu(
#     group_1,
#     group_2,
#     alternative="two-sided"
# )
#
# rank_biserial = (
#     2 * u_stat / (len(group_1) * len(group_2)) - 1
# )
#
# this was adapted to compare ProtNLM top-1 scores for generic and non-generic
# predictions. GAI was also used to help interpret the direction of the effect
# size and to make sure the output was described as a difference in score
# distributions rather than a validation of either prediction group.

# lines 292-360 were written with GAI assistance for producing the three
# diagnostic scatter plots from one reusable plotting routine. i asked how to
# loop over several x variables, optionally log-transform one variable, add a
# small reproducible jitter where points overlap heavily, and save both PNG and
# PDF versions. GAI supplied general examples using np.log1p(),
# np.random.default_rng(), ax.scatter() and a list of plotting specifications.
# these were adapted to ortholog count, protein length and lexical beam
# coherence. no regression line was added because the figures are intended to
# show the observed data rather than imply a fitted predictive relationship.

# finally, GAI helped with lines 367-402 and the corresponding calls in main()
# when assembling the results into a short markdown report. i asked how values
# stored in the correlation and group-comparison dataframes could be inserted
# automatically into a text summary using an f-string and then written to a
# .md file. an example structure was:
#
# report = f"""
# # Analysis summary
# - rho = {rho:.3f}
# - p = {p_value:.3g}
# """
#
# Path("your_summary.md").write_text(report, encoding="utf-8")
#
# i adapted this so the markdown file records the correlation results,
# generic/non-generic comparison and beam-review summaries generated by the
# script. the wording was deliberately written to emphasise that the analyses
# are exploratory and that statistical association or lexical consistency does
# not establish biological correctness.