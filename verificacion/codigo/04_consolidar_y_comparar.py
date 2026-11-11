#!/usr/bin/env python3
"""
Consolida las tres salidas de gnomAD y compara la frecuencia del ALELO ALTERNO
entre AMR y los demas grupos poblacionales.

Convencion: se reportan los alelos tal como los define gnomAD sobre GRCh38.
REF es el alelo de referencia del ensamble y ALT el alterno, SIN reasignarlos
segun cual sea mas frecuente ni segun cual se considere de riesgo en la
literatura. Todas las frecuencias son del alelo ALT.

La correspondencia con la nomenclatura clasica (C677T, C282Y, FokI, etc.) se
documenta en data/raw/anotacion_variantes.tsv y se usa solo para etiquetar las
tablas, nunca para transformar las frecuencias.

Entradas:
  data/processed/gnomad_frecuencias.tsv
  data/processed/gnomad_frecuencias_rescate.tsv
  data/processed/gnomad_frecuencias_ventanas.tsv
  data/raw/anotacion_variantes.tsv  (opcional, solo etiquetas)

Salidas:
  results/tabla1_frecuencias.tsv       AF del alelo ALT por poblacion
  results/tabla2_comparacion_amr.tsv   AMR vs comparadores, OR, Fisher, q BH

Correr: python3 scripts/04_consolidar_y_comparar.py
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

RAIZ = Path(__file__).resolve().parent.parent
PROC = RAIZ / "data" / "processed"
RES = RAIZ / "results"
ANOT = RAIZ / "data" / "raw" / "anotacion_variantes.tsv"

ENTRADAS = [
    PROC / "gnomad_frecuencias.tsv",
    PROC / "gnomad_frecuencias_rescate.tsv",
    PROC / "gnomad_frecuencias_ventanas.tsv",
]

REFERENCIA = "amr"
COMPARADORES = ["nfe", "afr", "eas", "sas"]
TODAS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining"]

# AC global minimo para conservar un alelo alterno en un sitio multialelico.
# Los alelos con conteos residuales son ruido de secuenciacion o alelos reales
# pero irrelevantes; en cualquier caso no son la variante de la literatura.
AC_MINIMO = 100


def alelos(variant_id):
    """'12-51005267-T-G' -> ('T', 'G')"""
    partes = str(variant_id).split("-")
    return (partes[2], partes[3]) if len(partes) >= 4 else ("", "")


def cargar():
    partes = []
    for ruta in ENTRADAS:
        if ruta.exists():
            partes.append(pd.read_csv(ruta, sep="\t"))
        else:
            print(f"aviso: no existe {ruta.name}", file=sys.stderr)
    if not partes:
        sys.exit("No hay archivos de frecuencias que consolidar.")

    df = pd.concat(partes, ignore_index=True)
    df = df.drop_duplicates(subset="variant_id").reset_index(drop=True)

    cols_ac = [f"AC_{p}" for p in TODAS if f"AC_{p}" in df.columns]
    df["AC_global"] = df[cols_ac].sum(axis=1)

    # --- sitios multialelicos ---------------------------------------------
    # Un rsID puede devolver varias filas, una por alelo alterno. Se conserva
    # la de mayor AC global y se registran las descartadas.
    df = df.sort_values("AC_global", ascending=False)
    descartadas = df[df.duplicated(subset="rsid", keep="first")].copy()
    df = df.drop_duplicates(subset="rsid", keep="first")

    # Alelos residuales aun sin duplicado.
    residuales = df[df["AC_global"] < AC_MINIMO].copy()
    df = df[df["AC_global"] >= AC_MINIMO]

    df[["alelo_ref", "alelo_alt"]] = df["variant_id"].apply(
        lambda v: pd.Series(alelos(v)))

    return df.reset_index(drop=True), descartadas, residuales


def benjamini_hochberg(p):
    s = pd.Series(p).dropna()
    n = len(s)
    orden = s.sort_values()
    q = (orden.values * n) / (pd.RangeIndex(1, n + 1).values)
    q = pd.Series(q, index=orden.index)[::-1].cummin()[::-1]
    return q.reindex(pd.Series(p).index)


def main():
    df, descartadas, residuales = cargar()

    print(f"{len(df)} variantes conservadas (una por rsID)")
    if len(descartadas):
        print(f"\n{len(descartadas)} alelo(s) alterno(s) descartado(s) "
              f"por sitio multialelico:")
        for _, r in descartadas.iterrows():
            print(f"  {r['gen']:<9} {r['rsid']:<12} {r['variant_id']:<18} "
                  f"AC global {int(r['AC_global'])}")
    if len(residuales):
        print(f"\n{len(residuales)} variante(s) con AC global < {AC_MINIMO}, "
              f"excluida(s):")
        for _, r in residuales.iterrows():
            print(f"  {r['gen']:<9} {r['rsid']:<12} AC {int(r['AC_global'])}")

    # Etiquetas de la literatura, si existen. Nunca tocan las frecuencias.
    if ANOT.exists():
        anot = pd.read_csv(ANOT, sep="\t", dtype=str).fillna("")
        cols = [c for c in ["rsid", "nomenclatura_literatura", "nota"]
                if c in anot.columns]
        df = df.merge(anot[cols], on="rsid", how="left")
    else:
        df["nomenclatura_literatura"] = ""
        print(f"\naviso: no se encontro {ANOT.name}; "
              f"las tablas saldran sin etiquetas de literatura")

    RES.mkdir(parents=True, exist_ok=True)

    # ---- Tabla 1 ---------------------------------------------------------
    cols_t1 = (["gen", "rsid", "variant_id", "alelo_ref", "alelo_alt",
                "nomenclatura_literatura", "consecuencia", "hgvsp"]
               + [f"AF_{p}" for p in TODAS if f"AF_{p}" in df.columns]
               + [f"AN_{p}" for p in TODAS if f"AN_{p}" in df.columns])
    t1 = df[[c for c in cols_t1 if c in df.columns]].sort_values(["gen", "rsid"])
    t1.to_csv(RES / "tabla1_frecuencias.tsv", sep="\t", index=False)

    # ---- Tabla 2 ---------------------------------------------------------
    filas = []
    for _, r in df.iterrows():
        ac_r, an_r = r.get(f"AC_{REFERENCIA}"), r.get(f"AN_{REFERENCIA}")
        if not an_r:
            continue
        for comp in COMPARADORES:
            ac_c, an_c = r.get(f"AC_{comp}"), r.get(f"AN_{comp}")
            if not an_c or not ac_c:
                continue
            tabla = [[int(ac_r), int(an_r - ac_r)], [int(ac_c), int(an_c - ac_c)]]
            orr, p = fisher_exact(tabla)
            filas.append({
                "gen": r["gen"], "rsid": r["rsid"],
                "alelo_alt": r["alelo_alt"],
                "nomenclatura_literatura": r.get("nomenclatura_literatura", ""),
                "comparador": comp,
                f"AF_{REFERENCIA}": round(ac_r / an_r, 5),
                f"AF_{comp}": round(ac_c / an_c, 5),
                "razon_AF": round((ac_r / an_r) / (ac_c / an_c), 3),
                "OR": round(orr, 3),
                "p": p,
            })

    t2 = pd.DataFrame(filas)
    t2["q_BH"] = benjamini_hochberg(t2["p"])
    t2 = t2.sort_values(["gen", "rsid", "comparador"])
    t2.to_csv(RES / "tabla2_comparacion_amr.tsv", sep="\t", index=False)

    print(f"\nTabla 1 -> {RES / 'tabla1_frecuencias.tsv'}")
    print(f"Tabla 2 -> {RES / 'tabla2_comparacion_amr.tsv'}")

    # ---- Resumen AMR vs NFE ----------------------------------------------
    vs_nfe = t2[t2["comparador"] == "nfe"].sort_values("razon_AF")
    print(f"\nAMR vs NFE, frecuencia del alelo ALT ({len(vs_nfe)} variantes):")
    print(f"{'gen':<10}{'rsid':<12}{'ALT':<5}{'AF_amr':>9}{'AF_nfe':>9}"
          f"{'razon':>8}  {'q<0.05':<7}")
    for _, r in vs_nfe.iterrows():
        marca = "si" if r["q_BH"] < 0.05 else "no"
        print(f"{r['gen']:<10}{r['rsid']:<12}{r['alelo_alt']:<5}"
              f"{r['AF_amr']:>9.3f}{r['AF_nfe']:>9.3f}{r['razon_AF']:>8.2f}  {marca:<7}")

    print(f"\nCon AN de este tamano casi todo alcanza q<0.05. "
          f"Interpretar por razon_AF, no por q.")


if __name__ == "__main__":
    main()
