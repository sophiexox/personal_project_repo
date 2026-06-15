# C. auris AI-assisted protein annotation prototyping

# scripts for testing how protein language model outputs could support the existing C. auris functional annotation pipeline.

## Current focus::::

- Check HPC access and environment requirements.
- Review ProtNote / ProtNLM-style protein annotation tools.
- Prepare inputs for GO2Sum from predicted GO terms.
- Compare AI-derived summaries against existing annotated C. auris genes.
- Identify where this could fit alongside InterProScan, FoldSeek/CATH and Pfam-based naming.

## Planned workflow

1. Start with C. auris proteins currently labelled hypothetical / uncharacterised.
2. Collect existing evidence:
   - InterProScan domains
   - GO terms
   - CATH/FoldSeek structural hits
   - Pfam/domain labels
   - UniProt/CGD annotations where available
3. Run or import predictions from ProtNote / ProtNLM-style tools.
4. Convert GO predictions to summaries using GO2Sum.
5. Compare generated summaries against curated or semi-curated C. auris annotations.
6. Flag high-confidence cases where multiple methods agree.

## Notes

This repo is currently a scaffold rather than a finished production pipeline.