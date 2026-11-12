# Ancestry-component allele frequencies of variants previously proposed as markers of lead susceptibility in the Mexico City Prospective Study

Analysis code and supplementary material for the manuscript.

The work describes the allele frequencies of the 27 loci in 10 genes proposed
by Ding et al. 2016 as markers of lead susceptibility, separated by Indigenous
Mexican and European ancestry component through local ancestry inference in the
Mexico City Prospective Study, and evaluates whether the differentiation
observed is exceptional against an empirical background of variants of
comparable frequency in 219 control genes.

## Checking that the figures follow from the data

    python3 scripts/28_reproduccion_independiente.py

That script reads only `data/processed/mcps_frecuencias.tsv` and
`data/raw/regiones.bed`, reimplements from scratch what the Methods section
describes without importing any other script in the pipeline, and compares each
figure it obtains with the one published in the supplementary tables. It
reproduces seventeen figures with no differences, including the frequencies
recomputed from AC and AN, Hudson's FST estimator, the three variant counts, the
outlier blocks of the three sets, the rate ratio, the binomial test and the
percentiles of the targeted analysis.

## Supplementary tables

`Table_1_body` is Table 1 as it appears in the article: ten columns with the
gene, the rsID, the GRCh38 coordinate, the reference and alternate alleles, the
name under which the literature cites the locus, the named allele the alternate
allele corresponds to where the literature names one, the position of the
alternate allele relative to the minor allele of the source study, the support
for that correspondence, and the frequencies in the Indigenous Mexican and
European components.

`Table_1_allelic_orientation` is the complete version, with 24 columns. It keeps
apart two things worth not conflating. The first is whether the alternate allele
corresponds to an allele the literature names, which happens at nine of the 24
loci. The second is whether it corresponds to the allele whose effect the source
study reported, graded on three levels: confirmed at 19 loci, provisional at 4
and undetermined at 1. Of the 23 resolved, at 19 the alternate allele is the
minor allele of that study and at 4 the major one. The direction of effect is
reported with its full context, and the ratios are hazard ratios from a Cox
model expressing the association between bone lead and coronary heart disease
per doubling of the biomarker within each genotype stratum, not the effect of
carrying the allele.

`Table_S1_results_per_locus` holds the main result per locus: AC allele counts
and AN chromosome numbers for both components, frequencies, the matching
frequency under both definitions, an indicator of major allele inversion,
Hudson's FST, the number of neighbours and blocks, the matching distance
reached, the percentile with its resampling interval, and the percentile under
the alternative definition. Counts are fractional because they come from local
ancestry inference, and frequencies below 1e-4 are given in scientific notation:
HFE rs1800562 has an estimated count of 0.2 out of 170,634 chromosomes in the
Indigenous component, that is 1.17e-6, which to four decimal places would print
as zero.

`Table_S2_sensitivity` gives the sensitivity analyses in long format, 240 rows
across four scenarios. `Table_S2b_control_subsets` summarizes by control gene
subset. `Table_S3_scan_grid` is the complete grid of 54 comparisons; in the
restricted sets the percentiles are recomputed within the set rather than
inheriting the cutoff obtained with all variants. `Table_S4_block_extent` gives
the physical extent of the blocks, counted within each stratum just as in the
grid.

`Table_S5_matched_null` holds the null distribution of outlier blocks in control
regions matched to the candidates on length, density of evaluable common
variants and whole-genome fraction. It replaces the denominator of the binomial
test, which distributed the expected blocks by region length and assumed that a
candidate megabase and a control megabase offer the same opportunity to produce
them. It equalizes opportunity, not architecture. The table holds two runs of
each variant set, one with the nine candidate genes and one with the worst
matched of them removed from both the candidate set and the matching,
distinguished by the `excluded_genes` column; the nominal value of the complete
set does not survive that removal. `Table_S6_region_matching` reports the
quality of that matching gene by gene.

`Table_S7_fst_simulation` holds the simulation of Hudson's estimator over
fractional counts. It answers two things: the sample size correction remains
practically unbiased when the denominator is a sum of probabilities, with a
maximum absolute bias of 0.0008; and assignment error attenuates FST toward zero
rather than inflating it under the schemes simulated, which are symmetric
between components and independent across variants. `Table_S8_tail_calibration`
extends that to the upper tail, which is what the percentile threshold uses: the
99th percentile of the null FST distribution under assignment error does not
exceed the one obtained with perfect assignment at either of the two sample
sizes evaluated. Neither simulation covers error correlated along a haplotype or
asymmetric between components, so the values are not presented as bounds.

## Structure

    LEEME.txt                        full index, table by table (Spanish)
    Tabla_*.tsv                      tables as generated (Spanish headers)
    en/                              the same tables with English headers
    verificacion/                    outlier blocks, control genes,
                                     regiones.bed and provenance with SHA-256
    verificacion/codigo/             the complete pipeline, with its manifest

`verificacion/codigo/MANIFIESTO_CODIGO.txt` carries the library versions, the
execution order, the seeds and the intermediate files needed to regenerate the
supplementary material. Script comments are in Spanish.

## Environment

Python 3.11 or later, with numpy 1.24, pandas 2.0, scipy 1.10, requests 2.28 and
matplotlib 3.7. Variant extraction uses htslib 1.17, and converting the
manuscript to Word uses pandoc 3.0 and python-docx 1.1.

    pip install "numpy>=1.24" "pandas>=2.0" "scipy>=1.10" requests \
        "matplotlib>=3.7" "python-docx>=1.1"

## Data not redistributed

Mexico City Prospective Study frequencies are publicly available at
https://rgc-mcps.regeneron.com/ and are not included here. Coordinates and the
resolution of rs identifiers come from gnomAD v4.1 through its GraphQL API.
Query dates and the SHA-256 checksums of each download are in
`verificacion/procedencia.txt`.

## License

The code is distributed under the MIT license. The result tables are
distributed under CC BY 4.0.
