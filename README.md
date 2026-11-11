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

`Table_1_body` is Table 1 as it appears in the article. `Table_1_allelic_orientation`
is the complete version, with 24 columns, which keeps apart whether the
alternate allele corresponds to an allele the literature names and whether it
corresponds to the allele whose effect the source study reported: confirmed at
19 loci, provisional at 4 and undetermined at 1.

`Table_S1_results_per_locus` holds the main result per locus. Counts are
fractional because they come from local ancestry inference: HFE rs1800562 has an
estimated count of 0.2 out of 170,634 chromosomes in the Indigenous component.

`Table_S2_sensitivity` and `Table_S2b_control_subsets` hold the sensitivity
analyses. `Table_S3_scan_grid` is the grid of 54 comparisons, with percentiles
recomputed within each restricted set. `Table_S4_block_extent` gives the
physical extent of the blocks.

`Table_S5_matched_null` holds the null distribution of outlier blocks in control
regions matched on length, density of evaluable common variants and
whole-genome fraction. It replaces the denominator of the binomial test, which
assumed that a candidate megabase and a control megabase offer the same
opportunity. It equalizes opportunity, not architecture.
`Table_S6_region_matching` reports the quality of that matching gene by gene.

`Table_S7_fst_simulation` holds the simulation of Hudson's estimator over
fractional counts: the sample size correction remains practically unbiased when
the denominator is a sum of probabilities, with a maximum absolute bias of
0.0008, and assignment error attenuates FST toward zero rather than inflating
it, so the values reported are lower bounds.

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
matplotlib 3.7. Variant extraction uses htslib 1.17.

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
