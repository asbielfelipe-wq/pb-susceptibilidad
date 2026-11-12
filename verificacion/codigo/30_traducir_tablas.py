#!/usr/bin/env python3
"""
Traduce al ingles los encabezados y los valores categoricos de las tablas.

No toca ninguna cifra: solo sustituye nombres de columna y cadenas de texto
por su equivalente en ingles, segun un diccionario explicito. Lo que no
aparece en el diccionario se deja tal cual, y al final el script avisa de
toda cadena que siga pareciendo espanola, para que no se quede nada sin
traducir sin que se note.

Correr: python3 scripts/30_traducir_tablas.py
Entradas: results/suplementario/*.tsv
Salidas:  results/suplementario/en/*.tsv
"""

import csv
import re
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
SUP = RAIZ / "results" / "suplementario"
SAL = SUP / "en"

# Nombre del archivo traducido. Se renombra para que el paquete en ingles no
# mezcle idiomas en los nombres de archivo.
ARCHIVOS = {
    "Tabla_1_cuerpo.tsv": "Table_1_body.tsv",
    "Tabla_1_orientacion_alelica.tsv": "Table_1_allelic_orientation.tsv",
    "Tabla_S1_resultados_por_locus.tsv": "Table_S1_results_per_locus.tsv",
    "Tabla_S2_sensibilidades.tsv": "Table_S2_sensitivity.tsv",
    "Tabla_S2b_subconjuntos_control.tsv": "Table_S2b_control_subsets.tsv",
    "Tabla_S3_rejilla_barrido.tsv": "Table_S3_scan_grid.tsv",
    "Tabla_S4_extension_bloques.tsv": "Table_S4_block_extent.tsv",
    "Tabla_S5_nula_apareada.tsv": "Table_S5_matched_null.tsv",
    "Tabla_S6_apareamiento_regiones.tsv": "Table_S6_region_matching.tsv",
    "Tabla_S7_simulacion_fst.tsv": "Table_S7_fst_simulation.tsv",
    "Tabla_S8_calibracion_cola.tsv": "Table_S8_tail_calibration.tsv",
    "bloques_atipicos.tsv": "outlier_blocks.tsv",
}

COLUMNAS = {
    "Gen": "Gene", "gen": "gene", "gen_candidato": "candidate_gene",
    "Nombre": "Name", "Posicion": "Position", "Certeza": "Certainty",
    "Alelo nombrado": "Named allele", "ALT en el estudio": "ALT in source study",
    "ALT_frente_al_alelo_menor_de_Ding": "ALT_vs_Ding_minor_allele",
    "alelo_de_literatura_que_corresponde_a_ALT": "literature_allele_matching_ALT",
    "base_del_alelo_nombrado": "basis_of_named_allele",
    "base_de_la_orientacion_frente_a_Ding": "basis_of_orientation_vs_Ding",
    "certeza_de_la_orientacion": "orientation_certainty",
    "biomarcador_y_desenlace": "biomarker_and_outcome",
    "modelo_genetico": "genetic_model",
    "direccion_en_el_estudio_de_origen": "direction_in_source_study",
    "significacion_tras_FDR": "significance_after_FDR",
    "respaldo_del_gen": "gene_support",
    "respaldo_de_la_variante": "variant_support",
    "transcrito_de_la_anotacion": "annotation_transcript",
    "coordenada_GRCh38": "GRCh38_coordinate", "consecuencia": "consequence",
    "nomenclatura": "nomenclature", "nota": "note", "fuente": "source_of_field",
    "clase": "class", "grupo": "group", "estrato": "stratum",
    "conjunto": "set", "escenario": "scenario", "condicion": "condition",
    "subconjunto": "subset", "locus": "locus",
    "maf_apareamiento_principal": "matching_maf_principal",
    "maf_apareamiento_alternativa": "matching_maf_alternative",
    "maf_media": "mean_maf", "dif_abs": "abs_difference",
    "invierte_mayoritario": "major_allele_inverted",
    "n_vecinos": "n_neighbours", "n_bloques": "n_blocks",
    "distancia_max_al_fondo": "max_distance_to_background",
    "distancia_max": "max_distance",
    "percentil": "percentile", "pct_def_alternativa": "pct_alternative_definition",
    "cambio_entre_definiciones": "change_between_definitions",
    "IC95_inferior": "CI95_lower", "IC95_superior": "CI95_upper",
    "IC_cand_Mb": "CI_cand_per_Mb",
    "bloques_cand": "blocks_cand", "bloques_ctrl": "blocks_ctrl",
    "var_atipicas_cand": "outlier_var_cand", "var_atipicas_ctrl": "outlier_var_ctrl",
    "tasa_cand_Mb": "rate_cand_per_Mb", "tasa_ctrl_Mb": "rate_ctrl_per_Mb",
    "razon": "ratio", "p_binomial_unilateral": "p_one_sided_binomial",
    "umbral_percentil": "percentile_threshold", "ventana_kb": "window_kb",
    "extension_kb_mediana": "extent_kb_median", "extension_kb_p90": "extent_kb_p90",
    "extension_kb_maxima": "extent_kb_max",
    "variantes_por_bloque_mediana": "variants_per_block_median",
    "variantes_por_bloque_maxima": "variants_per_block_max",
    "bloques_candidatos_observados": "observed_candidate_blocks",
    "nula_media": "null_mean", "nula_mediana": "null_median",
    "nula_IC95": "null_CI95", "nula_maxima": "null_maximum",
    "p_empirica": "p_empirical", "replicas": "replicates",
    "genes_candidatos": "candidate_genes",
    "genes_sin_comparadores_suficientes": "genes_without_enough_comparators",
    "relajacion_maxima": "max_relaxation", "relajacion": "relaxation",
    "replicas_con_reutilizacion": "replicates_with_reuse",
    "largo_kb": "length_kb", "variantes_evaluables": "evaluable_variants",
    "densidad_por_kb": "density_per_kb", "fraccion_wgs": "wgs_fraction",
    "comparadores": "comparators", "tolerancia_largo": "length_tolerance",
    "tolerancia_densidad": "density_tolerance",
    "tolerancia_fraccion_wgs": "wgs_fraction_tolerance",
    "mediana_largo_kb_comparadores": "median_length_kb_comparators",
    "mediana_densidad_comparadores": "median_density_comparators",
    "p_indigena": "p_Indigenous", "p_europea": "p_European",
    "FST_verdadero": "true_FST", "FST_estimado_medio": "mean_estimated_FST",
    "cromosomas": "chromosomes",
    "exactitud_asignacion": "mean_posterior_to_true_ancestry",
    "variantes": "variants", "FST_mediana": "FST_median",
    "FST_p95": "FST_p95", "FST_p99": "FST_p99", "FST_p999": "FST_p999",
    "FST_maxima": "FST_max",
    "proporcion_sobre_p99_sin_error": "fraction_above_error_free_p99",
    "sesgo": "bias", "sesgo_relativo": "relative_bias",
    "desviacion_estandar": "standard_deviation",
    "AN_indigena_medio": "mean_AN_Indigenous",
    "loci_evaluables": "evaluable_loci",
    "loci_percentil_95_o_mas": "loci_at_or_above_p95",
    "mediana": "median", "maximo": "maximum", "n_genes": "n_genes",
    "evaluable": "evaluable", "bloque": "block", "pos": "pos",
}

VALORES = {
    "si": "yes", "no": "no", "alta": "high", "wgs": "wgs", "WGS": "WGS",
    "candidato": "candidate", "control": "control", "completo": "complete",
    "restringido_asimetrico": "restricted_asymmetric",
    "restringido_simetrico": "restricted_symmetric",
    "menor": "minor", "mayor": "major", "indet.": "undet.",
    "confirmada": "confirmed", "provisional": "provisional",
    "indeterminada": "undetermined", "no determinable": "not determinable",
    "promotor": "promoter", "alelo F": "F allele",
    "ALT es el alelo menor de Ding": "ALT is Ding's minor allele",
    "ALT es el alelo mayor de Ding": "ALT is Ding's major allele",
    "par de alelos no ambiguo": "unambiguous allele pair",
    "par de alelos no ambiguo, concordante con la anotacion proteica":
        "unambiguous allele pair, concordant with the protein annotation",
    "anotacion proteica, que resuelve la transversion ambigua":
        "protein annotation, which resolves the ambiguous transversion",
    "concordancia de frecuencia en transversion ambigua":
        "frequency concordance in an ambiguous transversion",
    "transversion ambigua con lecturas discordantes entre si":
        "ambiguous transversion with mutually discordant readings",
    "identidad del identificador del estudio de origen inferida del par de "
    "alelos y de la concordancia de frecuencia, sin fusion registrada en dbSNP":
        "identity of the source study identifier inferred from the allele pair "
        "and from frequency concordance, with no merge recorded in dbSNP",
    "confirmada por alelos o por anotacion molecular":
        "confirmed by alleles or by molecular annotation",
    "provisional: concordancia de frecuencia, sin verificacion de hebra":
        "provisional: frequency concordance, without strand verification",
    "provisional: identidad con el marcador del estudio de origen inferida, "
    "equivalencia de identificadores no documentada":
        "provisional: identity with the source study marker inferred, "
        "identifier equivalence not documented",
    "anotacion proteica del transcrito canonico":
        "protein annotation of the canonical transcript",
    "no aplicable: la literatura no nombra un alelo":
        "not applicable: the literature does not name an allele",
    "no aplicable: de promotor, sin alelo nombrado en la literatura":
        "not applicable: promoter, no allele named in the literature",
    "no aplicable: intronica, sin alelo nombrado en la literatura":
        "not applicable: intronic, no allele named in the literature",
    "no aplicable: intronica, la literatura la nombra por el sitio BsmI y no "
    "por la base":
        "not applicable: intronic, the literature names it by the BsmI site "
        "and not by the base",
    "no aplicable: intronica, la literatura la nombra por el sitio Tru9I":
        "not applicable: intronic, the literature names it by the Tru9I site",
    "no aplicable: intronica, la literatura la nombra por el sitio ApaI":
        "not applicable: intronic, the literature names it by the ApaI site",
    "no aplicable: sinonima, la literatura la nombra por el sitio TaqI":
        "not applicable: synonymous, the literature names it by the TaqI site",
    "Thr268, que la literatura llama T235 o M268T segun la numeracion":
        "Thr268, which the literature calls T235 or M268T depending on the "
        "numbering",
    "ALAD2 (Asn59; C en la hebra del gen, G en la hebra de referencia)":
        "ALAD2 (Asn59; C on the gene strand, G on the reference strand)",
    "el alelo Arg130 que, junto con el alelo C de rs7412, define e4; no es e4 "
    "por si solo":
        "the Arg130 allele which, together with the C allele of rs7412, defines "
        "e4; it is not e4 on its own",
    "el alelo Cys176 que, junto con el alelo T de rs429358, define e2; no es "
    "e2 por si solo":
        "the Cys176 allele which, together with the T allele of rs429358, "
        "defines e2; it is not e2 on its own",
    "alelo F, el que carece del primer ATG y da la proteina de 424 residuos; "
    "ALT produce p.Met1? y su frecuencia europea de 0.65 coincide con la "
    "descrita para F":
        "F allele, the one lacking the first ATG and yielding the 424-residue "
        "protein; ALT produces p.Met1? and its European frequency of 0.65 "
        "matches the one described for F",
    "plomo en tibia, cardiopatia coronaria incidente; el HR es por duplicacion "
    "del plomo oseo dentro de cada estrato genotipico y no el efecto de portar "
    "el alelo":
        "tibia lead, incident coronary heart disease; the HR is per doubling of "
        "bone lead within each genotype stratum and not the effect of carrying "
        "the allele",
    "plomo en rotula, cardiopatia coronaria incidente; el HR es por duplicacion "
    "del plomo oseo dentro de cada estrato genotipico y no el efecto de portar "
    "el alelo":
        "patella lead, incident coronary heart disease; the HR is per doubling "
        "of bone lead within each genotype stratum and not the effect of "
        "carrying the allele",
    "aditivo, genotipo como variable continua":
        "additive, genotype as a continuous variable",
    "dominante, portadores del alelo menor frente a no portadores":
        "dominant, minor allele carriers versus non-carriers",
    "no aplicable: la interaccion no fue significativa antes de la correccion":
        "not applicable: the interaction was not significant before correction",
    "no; FDR ajustada = 0.06, sin alcanzar el umbral de 0.05":
        "no; adjusted FDR = 0.06, not reaching the 0.05 threshold",
    "establecido en toxicocinetica del plomo":
        "established in lead toxicokinetics",
    "mecanistico, evidencia fragmentaria": "mechanistic, fragmentary evidence",
    "evaluada en el estudio de origen sin modificacion significativa":
        "evaluated in the source study with no significant modification",
    "interaccion limitrofe en el estudio de origen: FDR ajustada = 0.06, que no "
    "alcanza el umbral convencional de 0.05":
        "borderline interaction in the source study: adjusted FDR = 0.06, which "
        "does not reach the conventional 0.05 threshold",
    "no evaluada contra un desenlace de plomo fuera del estudio de origen":
        "not evaluated against a lead outcome outside the source study",
    "evaluada en el estudio de origen sin modificacion significativa; "
    "metaanalisis posterior de 27 estudios sin asociacion con plombemia":
        "evaluated in the source study with no significant modification; a "
        "later meta-analysis of 27 studies found no association with blood lead",
    "sin diferenciacion": "no differentiation",
    "diferenciacion leve": "mild differentiation",
    "diferenciacion fuerte": "strong differentiation",
    "alelo poco frecuente": "low-frequency allele",
    "alelo mayoritario invertido": "major allele inverted",
    "version de frecuencia": "frequency version",
    "plataformas separadas": "platforms separated",
    "subconjunto de genes control": "control gene subset",
    "definicion alternativa de la frecuencia de apareamiento":
        "alternative definition of the matching frequency",
    "AGTR1; sin estudios de asociacion con plomo fuera de Ding":
        "AGTR1; no lead association studies outside Ding",
    "Alelo ancestral; estructura poblacional documentada":
        "Ancestral allele; documented population structure",
    "HFE reconocido; Ding no hallo modificacion de efecto":
        "HFE recognized; Ding found no effect modification",
    "r2 alto con rs1544410 segun Ding": "high r2 with rs1544410 per Ding",
    "Intronica en el transcrito MANE Select; la consecuencia proteica de un "
    "transcrito alternativo no describe esta variante y se omite":
        "Intronic in the MANE Select transcript; the protein consequence of an "
        "alternative transcript does not describe this variant and is omitted",
    "Marcador reconocido; metaanalisis 2025 sin asociacion global con "
    "plombemia; La Tabla 3 de la fuente imprime este identificador como "
    "rs1833435, mientras que su seccion de genotipificacion y su Tabla 1 lo "
    "dan como rs1800435; se anota para localizar la fila correspondiente":
        "Recognized marker; a 2025 meta-analysis found no overall association "
        "with blood lead; Table 3 of the source prints this identifier as "
        "rs1833435, while its genotyping section and its Table 1 give it as "
        "rs1800435; noted here so the corresponding row can be located",
    "Ding reporta rs1073581; nomenclatura vigente rs2228570; FokI figura en el "
    "estudio de origen y en publicaciones previas con el identificador "
    "rs1073581. En dbSNP build 157 FokI corresponde a rs2228570, en "
    "chr12:47,879,112 (GRCh38), identificador en el que se fusiono rs10735810, "
    "y es el que se analizo aqui. Los alelos A y G de rs2228570 en la hebra de "
    "referencia coinciden con los genotipos G/G, G/A y A/A que reporta ese "
    "estudio, y la frecuencia europea del alelo alterno, 0.652, concuerda con "
    "la designacion de A como alelo menor, de modo que ALT corresponde al "
    "alelo mayor de ese estudio y al alelo F de la nomenclatura clasica. Esa "
    "correspondencia descansa en la concordancia de alelos y frecuencias, no "
    "en una equivalencia de identificadores registrada en dbSNP, y por eso la "
    "orientacion se reporta como provisional":
        "Ding reports rs1073581; the current nomenclature is rs2228570; FokI "
        "appears in the source study and in earlier publications under the "
        "identifier rs1073581. In dbSNP build 157 FokI corresponds to "
        "rs2228570, at chr12:47,879,112 (GRCh38), the identifier into which "
        "rs10735810 was merged, and that is the one analyzed here. The A and G "
        "alleles of rs2228570 on the reference strand match the G/G, G/A and "
        "A/A genotypes that study reports, and the European alternate allele "
        "frequency, 0.652, agrees with the designation of A as the minor "
        "allele, so that ALT corresponds to the major allele of that study and "
        "to the F allele of the classical nomenclature. That correspondence "
        "rests on the concordance of alleles and frequencies, not on an "
        "identifier equivalence recorded in dbSNP, and for that reason the "
        "orientation is reported as provisional",
}

# Fragmentos formulaicos. Se sustituyen por patron y no cadena a cadena,
# para que las cifras que llevan dentro pasen intactas y no haya que
# enumerar las veinticuatro variantes de la misma frase.
FRASES = [
    (r"\bsin modificacion de efecto detectada\b",
     "no effect modification detected"),
    (r"; no hubo homocigotos ALT suficientes para estimar ese estrato",
     "; there were not enough ALT homozygotes to estimate that stratum"),
    (r"\bla asociacion se atenua con la dosis de ALT\b",
     "the association attenuates with ALT dosage"),
    (r"\bla asociacion se refuerza con la dosis de ALT\b",
     "the association strengthens with ALT dosage"),
    (r"\brelacion no monotona\b", "non-monotonic relationship"),
    (r"\bcon atenuacion en portadores de uno o dos alelos ALT\b",
     "with attenuation in carriers of one or two ALT alleles"),
    (r"\bentre portadores de al menos un alelo ALT\b",
     "among carriers of at least one ALT allele"),
    (r"\bentre no portadores del alelo menor de Ding\b",
     "among non-carriers of Ding's minor allele"),
    (r"\bentre no portadores de ALT\b", "among non-carriers of ALT"),
    (r"\bcomo la correspondencia con ALT no es determinable, la direccion no "
     r"puede asignarse a ALT ni a REF\b",
     "since the correspondence with ALT is not determinable, the direction "
     "cannot be assigned to ALT or to REF"),
    (r"\ben rotula no hubo modificacion\b",
     "in patella there was no modification"),
    (r"\bel alelo mas diferenciado aqui es el que menos reforzaba la "
     r"asociacion alli\b",
     "the most differentiated allele here is the one that least strengthened "
     "the association there"),
    (r"\ben homocigotos REF\b", "in REF homozygotes"),
    (r"\ben homocigotos ALT\b", "in ALT homozygotes"),
    (r"\ben heterocigotos\b", "in heterozygotes"),
    (r"\bsin ALT\b", "without ALT"),
    (r"\bIC95 ", "95% CI "),
    (r", que es ", ", which is "),
    (r"\bintronica\b", "intronic"),
    (r"El cambio p\.Asn14Lys corresponde a (\S+), que no es MANE",
     r"The p.Asn14Lys change corresponds to \1, which is not MANE"),
    (r"\btranscrito canonico\b", "canonical transcript"),
    (r"\bTabla (\d)\b", r"Table \1"),
    (r"\by dbSNP\b", "and dbSNP"),
    (r"\by Ding et al\.", "and Ding et al."),
]


def aplicar_frases(v):
    for pat, rep in FRASES:
        v = re.sub(pat, rep, v)
    return v


# Patrones que delatan una cadena espanola sin traducir.
# "no" y "si" se dejan fuera a proposito: aparecen en ingles y disparaban
# falsos positivos sobre cadenas ya traducidas.
SOSPECHA = re.compile(
    r"(?<![A-Za-z])(de|la|el|los|las|que|con|sin|por|para|segun|entre|"
    r"una|un|del|al|son|su|sus|mas|este|esta|hubo|fue|como|aqui|alli)"
    r"(?![A-Za-z])", re.I)
EXENTAS = re.compile(r"^(no|si|de novo)$", re.I)


def traducir(v):
    s = v.strip()
    if s in VALORES:
        return VALORES[s]
    return aplicar_frases(v)


def main():
    if not SUP.exists():
        sys.exit(f"Falta {SUP}")
    SAL.mkdir(parents=True, exist_ok=True)
    pendientes, n_arch = [], 0
    for f in sorted(SUP.glob("*.tsv")):
        with open(f, encoding="utf-8", newline="") as fh:
            r = list(csv.reader(fh, delimiter="\t"))
        if not r:
            continue
        enc = [COLUMNAS.get(c, c) for c in r[0]]
        cuerpo = [[traducir(c) for c in fila] for fila in r[1:]]
        destino = SAL / ARCHIVOS.get(f.name, f.name)
        with open(destino, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(enc)
            w.writerows(cuerpo)
        n_arch += 1
        for c in r[0]:
            if c not in COLUMNAS and SOSPECHA.search(c):
                pendientes.append((f.name, "columna", c))
        for fila in cuerpo:
            for c in fila:
                s = c.strip()
                limpio = re.sub(r"\bet al\.|doi:\S+|https?://\S+", "", s)
                if (len(s) > 3 and not EXENTAS.match(s)
                        and not re.fullmatch(r"[\d.,eE+-]+", s)
                        and len(SOSPECHA.findall(limpio)) >= 2):
                    pendientes.append((f.name, "valor", s[:90]))

    print(f"{n_arch} tablas traducidas en {SAL}")
    vistos, unicos = set(), []
    for a, t, s in pendientes:
        if (t, s) not in vistos:
            vistos.add((t, s)); unicos.append((a, t, s))
    if unicos:
        print(f"\n{len(unicos)} cadena(s) que siguen pareciendo espanolas:\n")
        for a, t, s in unicos:
            print(f"  [{a}] {t}: {s}")
        print("\nAnadirlas al diccionario y volver a correr.")
    else:
        print("\nNo queda ninguna cadena espanola sin traducir.")
    sys.exit(1 if unicos else 0)


if __name__ == "__main__":
    main()
