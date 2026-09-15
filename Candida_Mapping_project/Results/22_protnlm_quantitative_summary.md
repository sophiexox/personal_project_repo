# ProtNLM quantitative diagnostics

## Scope and interpretation

This exploratory analysis covers all 521 annotation-poor proteins. It tests associations, not model mechanism or training-set memorisation. A relationship with ortholog count could be consistent with family representation in known sequence space, but cannot demonstrate training-data leakage. Exact lexical beam coherence is internal wording stability only; it is not biological correctness or semantic agreement.

## Spearman correlations with ProtNLM top-1 score

- Ortholog count: rho = 0.050, p = 0.258, N = 521; missing predictor/outcome = 0/0.
- Paralog count: rho = 0.052, p = 0.232, N = 521; missing predictor/outcome = 0/0.
- Protein Length: rho = 0.122, p = 0.00519, N = 521; missing predictor/outcome = 0/0.
- Beam_Exact_Lexical_Coherence: rho = 0.131, p = 0.00264, N = 521; missing predictor/outcome = 0/0.
- Beam_Unique_Normalised_Labels: rho = -0.134, p = 0.00223, N = 521; missing predictor/outcome = 0/0.

Effect sizes should be considered alongside p-values. A small p-value does not make a weak association biologically important.

## Generic versus non-generic top-1 labels

- Generic: N = 85, median = 0.0587, mean = 0.0898.
- Non-generic: N = 436, median = 0.0867, mean = 0.1418.
- Mann-Whitney U = 13776.0, p = 1.81e-04, rank-biserial effect (generic minus non-generic orientation) = -0.257.

This comparison uses the pre-existing generic flag. It describes score distributions and does not validate either class of prediction.

## Beam and review outputs

- Median unique normalised labels among 10 beams: 10.0.
- Median exact lexical beam coherence: 0.10.
- Narrow suspicious-output rule matches: 13. These are lexical prompts for manual review, not biological classifications.
- Manual-review table rows: 28. Extremes use observed score/coherence quartiles, with up to 5 rows for each corner; all narrow suspicious matches are retained.

## Files

- `22_protnlm_quantitative_diagnostics.csv`: one row per protein with beam diagnostics and review flags.
- `22_protnlm_correlation_summary.csv`: effect sizes, p-values, N and missingness.
- `22_protnlm_binned_summary.csv`: data-derived quartile summaries with median, mean and IQR.
- `22_protnlm_generic_comparison.csv`: generic/non-generic descriptive and rank-based comparison.
- `22_protnlm_manual_review_extremes.csv`: illustrative extremes for manual review.
- Figures are saved as PNG and PDF without fitted regression lines.
