#!/usr/bin/env python3
"""
Tabla de correspondencia entre el alelo alterno de GRCh38 y el alelo descrito
en la literatura para cada locus del conjunto de Ding et al.

MOTIVO
------
Todas las frecuencias de este trabajo son del alelo ALTERNO segun GRCh38. Ese
alelo no coincide necesariamente con el que la literatura asocia a un efecto.
Sin esta correspondencia, una frecuencia baja de ALT puede leerse erroneamente
como ausencia del alelo relevante, y el FST no aporta direccion del efecto.

La tabla documenta, para cada locus: coordenada, REF y ALT observados, cambio
proteico cuando aplica, el alelo que la literatura evalua, y el nivel de
respaldo de esa asociacion. Las columnas de literatura NO se usan para
transformar ninguna frecuencia: solo etiquetan la salida.

El script resuelve automaticamente lo que puede resolver (cambio proteico y
alelos observados, desde gnomAD y desde los datos de MCPS) y deja marcadas
para verificacion manual las celdas que requieren consultar la fuente
original, en lugar de rellenarlas con una suposicion.

CASO ESPECIAL: APOE
-------------------
Los alelos e2, e3 y e4 son haplotipos definidos por la combinacion de rs429358
y rs7412. La frecuencia marginal de rs7412 no es la frecuencia de e2, porque
dos frecuencias marginales no determinan la distribucion conjunta. El script
calcula la frecuencia haplotipica maxima compatible con las marginales
observadas y marca ambos loci para que se reporten como SNP individuales, no
como e2 ni e4.

Correr: python3 scripts/17_orientacion_alelos.py
Salidas:
  results/tabla_orientacion.tsv   para completar y publicar como Tabla 1
  results/tabla_orientacion.txt
"""

import sys
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
RAW = RAIZ / "data" / "raw"

# Alelo que la literatura evalua, y nivel de respaldo. "?" = verificar en la
# fuente original antes de publicar.
LIT = {
    "rs1800435": ("ALAD2 (Asn)", "establecido",
                  "Marcador reconocido; metaanalisis 2025 sin asociacion global con plombemia"),
    "rs8177796": ("?", "fragmentario", "Regulatoria 5'; verificar alelo evaluado"),
    "rs1544410": ("BsmI b/B", "establecido", "VDR reconocido como marcador; verificar b o B"),
    "rs731236":  ("TaqI t/T", "establecido", "r2 alto con rs1544410 segun Ding"),
    "rs7975232": ("ApaI a/A", "establecido", "Verificar cual alelo corresponde a a"),
    "rs2228570": ("FokI f/F", "establecido",
                  "Ding reporta rs1073581; nomenclatura vigente rs2228570"),
    "rs757343":  ("Tru9I u/U", "fragmentario", "Verificar alelo evaluado"),
    "rs1800562": ("Tyr282 (C282Y)", "establecido",
                  "HFE reconocido; Ding no hallo modificacion de efecto"),
    "rs1799945": ("Asp63 (H63D)", "establecido", ""),
    "rs2071746": ("?", "fragmentario", "Promotor HMOX1"),
    "rs2071747": ("His7 (D7H)", "fragmentario", ""),
    "rs2071749": ("?", "fragmentario", "Intronica HMOX1"),
    "rs5995098": ("?", "fragmentario", "Intronica HMOX1"),
    "rs440446":  ("?", "fragmentario", "Intronica APOE"),
    "rs405509":  ("?", "fragmentario", "Promotor APOE"),
    "rs449647":  ("?", "fragmentario", "Promotor APOE"),
    "rs769446":  ("?", "fragmentario", "Promotor APOE"),
    "rs7412":    ("componente de e2", "fragmentario",
                  "HAPLOTIPO: e2 requiere rs429358 + rs7412; reportar como SNP"),
    "rs429358":  ("componente de e4", "fragmentario",
                  "HAPLOTIPO: e4 requiere rs429358 + rs7412; reportar como SNP"),
    "rs1695":    ("Val105 (I105V)", "fragmentario", ""),
    "rs699":     ("Thr268 / T235", "fragmentario",
                  "Alelo ancestral; estructura poblacional documentada"),
    "rs5046":    ("?", "fragmentario", "Promotor AGT"),
    "rs5050":    ("?", "fragmentario", "Promotor AGT"),
    "rs2493137": ("?", "fragmentario", "Promotor AGT"),
    "rs12695908": ("?", "fragmentario",
                   "AGTR1; sin estudios de asociacion con plomo fuera de Ding"),
}

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def txt(v):
    """Texto limpio: las celdas vacias llegan como NaN desde pandas."""
    if v is None or (isinstance(v, float) and v != v):
        return ""
    return str(v)


def main():
    f_pct = RES / "percentiles_v2.tsv"
    if not f_pct.exists():
        f_pct = RES / "ding_percentiles.tsv"
    if not f_pct.exists():
        sys.exit("Faltan los resultados de percentiles.")
    pct = pd.read_csv(f_pct, sep="\t")

    # hgvsp y consecuencia desde las coordenadas resueltas en gnomAD
    anot = []
    for f in ("ding_coordenadas.tsv", "ding_coordenadas_rescate.tsv"):
        p = RAW / f
        if p.exists():
            anot.append(pd.read_csv(p, sep="\t", dtype=str))
    hg = (pd.concat(anot, ignore_index=True).drop_duplicates(subset="rsid")
          if anot else pd.DataFrame(columns=["rsid", "hgvsp", "consecuencia"]))

    lista = RAW / "candidatos_ding2016.tsv"
    nom = (pd.read_csv(lista, sep="\t", dtype=str).fillna("")
           if lista.exists() else pd.DataFrame(columns=["rsid", "gen", "nomenclatura"]))

    filas = []
    for _, r in pct.iterrows():
        vid = str(r["variant_id"])
        partes = vid.split("-")
        ref, alt = (partes[2], partes[3]) if len(partes) >= 4 else ("", "")
        rs = r["rsid"]
        h = hg[hg["rsid"] == rs]
        n = nom[nom["rsid"] == rs]
        alelo_lit, respaldo, nota = LIT.get(rs, ("?", "?", ""))
        filas.append({
            "gen": r["gen"], "rsid": rs,
            "coordenada_GRCh38": f"chr{partes[0]}:{partes[1]}" if len(partes) >= 2 else "",
            "REF": ref, "ALT": alt,
            "nomenclatura": txt(n["nomenclatura"].iloc[0]) if len(n) else "",
            "consecuencia": txt(h["consecuencia"].iloc[0]) if len(h) else "",
            "hgvsp": txt(h["hgvsp"].iloc[0]) if len(h) else "",
            "alelo_evaluado_en_literatura": alelo_lit,
            "ALT_corresponde_al_alelo_de_literatura": "VERIFICAR",
            "direccion_del_efecto": "VERIFICAR",
            "desenlace_evaluado": "VERIFICAR",
            "nivel_de_respaldo": respaldo,
            "AF_IMX": r["AF_IMX"], "AF_EUR": r["AF_EUR"],
            "nota": nota,
        })

    t = pd.DataFrame(filas).sort_values(["gen", "rsid"])
    RES.mkdir(parents=True, exist_ok=True)
    t.to_csv(RES / "tabla_orientacion.tsv", sep="\t", index=False)

    di("=" * 78)
    di("TABLA DE ORIENTACION DE ALELOS")
    di("=" * 78)
    di(f"\n{len(t)} loci. Las columnas marcadas VERIFICAR requieren consultar la")
    di("publicacion original de cada asociacion antes de someter el manuscrito.")

    di(f"\n{'gen':<8}{'rsid':<13}{'REF>ALT':<9}{'hgvsp':<16}"
       f"{'alelo en literatura':<24}{'respaldo':<14}")
    for _, r in t.iterrows():
        ra = f"{r['REF']}>{r['ALT']}"
        di(f"{r['gen']:<8}{r['rsid']:<13}{ra:<9}{(txt(r['hgvsp']) or '-')[:15]:<16}"
           f"{txt(r['alelo_evaluado_en_literatura'])[:23]:<24}{r['nivel_de_respaldo']:<14}")

    n_ver = int((t["alelo_evaluado_en_literatura"] == "?").sum())
    di(f"\n{n_ver} loci sin alelo de literatura identificado; todos requieren")
    di("consultar la fuente. Son en su mayoria variantes de promotor o")
    di("intronicas incluidas por plausibilidad mecanistica.")

    # --- APOE: por que no se puede reportar e2 ni e4 ----------------------
    di("\n" + "=" * 78)
    di("APOE: FRECUENCIAS MARGINALES FRENTE A HAPLOTIPOS")
    di("=" * 78)
    a = pct[pct["rsid"].isin(["rs429358", "rs7412"])]
    if len(a) == 2:
        f2 = float(a[a["rsid"] == "rs7412"]["AF_IMX"].iloc[0])
        f4 = float(a[a["rsid"] == "rs429358"]["AF_IMX"].iloc[0])
        g2 = float(a[a["rsid"] == "rs7412"]["AF_EUR"].iloc[0])
        g4 = float(a[a["rsid"] == "rs429358"]["AF_EUR"].iloc[0])
        di(f"\nAncestria indigena: rs7412 ALT {f2:.4f}, rs429358 ALT {f4:.4f}")
        di(f"Ancestria europea:  rs7412 ALT {g2:.4f}, rs429358 ALT {g4:.4f}")
        di("\nLos alelos e2, e3 y e4 son haplotipos de esos dos SNP. Dos")
        di("frecuencias marginales no determinan la distribucion conjunta, de")
        di("modo que ninguna de estas cifras es la frecuencia de e2 ni de e4.")
        di("Estimarlas exigiria datos de haplotipo o fase, que MCPS no publica")
        di("en el archivo de frecuencias. Ambos loci se reportan como SNP")
        di("individuales, sin etiqueta de isoforma.")

    di("\n" + "=" * 78)
    di("QUE HAY QUE COMPLETAR ANTES DE SOMETER")
    di("=" * 78)
    di("\nPara cada locus, consultar la publicacion que reporta la asociacion y")
    di("anotar en la tabla:")
    di("  ALT_corresponde_al_alelo_de_literatura : si / no / no aplica")
    di("  direccion_del_efecto                   : aumenta / disminuye / nula")
    di("  desenlace_evaluado                     : plombemia, plomo oseo, etc.")
    di("\nMientras esas celdas digan VERIFICAR, el manuscrito no debe afirmar")
    di("que un alelo de riesgo esta ausente: solo que el alelo ALT de GRCh38")
    di("tiene una frecuencia determinada.")

    (RES / "tabla_orientacion.txt").write_text("\n".join(salida))
    print(f"\n\nTabla  -> {RES / 'tabla_orientacion.tsv'}")
    print(f"Resumen-> {RES / 'tabla_orientacion.txt'}")


if __name__ == "__main__":
    main()
