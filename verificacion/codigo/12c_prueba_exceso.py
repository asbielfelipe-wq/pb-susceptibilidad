#!/usr/bin/env python3
"""
Prueba formal del exceso de bloques atipicos en genes candidatos, con control
del artefacto de alelos casi privados de una ancestria.

DOS PROBLEMAS QUE RESUELVE
--------------------------
1. POTENCIA. El script 12b reporta razones crudas de bloques (5 frente a 9).
   Con conteos tan pequenos una razon de 5 es compatible con valores muy
   distintos. Aqui se calculan intervalos de confianza exactos de Poisson para
   la tasa de bloques por Mb, y una prueba binomial condicional sobre el total
   de bloques, que es la prueba apropiada para comparar dos tasas de conteo.

2. ALELOS CASI PRIVADOS. Varias senales corresponden a variantes presentes en
   una ancestria y practicamente ausentes en la otra (p. ej. 0.62 frente a
   0.008). Ese perfil es el esperado bajo el cuello de botella fundador
   americano, no bajo seleccion local, y produce FST alto por construccion. Se
   repite todo el analisis excluyendo variantes con AF < UMBRAL_PRIVADO en
   cualquiera de las dos ancestrias.

El analisis se repite ademas con varios umbrales de percentil y ventanas de
agrupamiento, para verificar que la conclusion no dependa de esas elecciones.

Correr: python3 scripts/12c_prueba_exceso.py
Salidas:
  results/exceso_pruebas.tsv
  results/exceso_resumen.txt
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, chi2

RAIZ = Path(__file__).resolve().parent.parent
TODAS = RAIZ / "results" / "barrido_todas.tsv"
BED = RAIZ / "data" / "raw" / "regiones.bed"
RES = RAIZ / "results"

UMBRALES_PCT = [99.0, 97.5, 95.0]
VENTANAS = [25000, 50000, 100000]
UMBRAL_PRIVADO = 0.02

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def ic_poisson(k, alpha=0.05):
    """Intervalo exacto de Poisson para un conteo k."""
    lo = 0.0 if k == 0 else chi2.ppf(alpha / 2, 2 * k) / 2
    hi = chi2.ppf(1 - alpha / 2, 2 * (k + 1)) / 2
    return lo, hi


def n_bloques(df, ventana):
    if not len(df):
        return 0, []
    d = df.sort_values(["chrom", "pos"])
    tam, actual, pc, pp, n = [], 0, None, None, 0
    for c, p in zip(d["chrom"].astype(str), d["pos"]):
        if pc != c or (p - pp) > ventana:
            actual += 1
            tam.append(0)
        tam[-1] += 1
        pc, pp = c, p
    return actual, tam


def mb_por_clase():
    tot = {}
    if not BED.exists():
        return tot
    with open(BED, encoding="utf-8") as fh:
        for linea in fh:
            if not linea.strip():
                continue
            c, i, f, nombre = linea.split("\t")[:4]
            _, _, clase = nombre.strip().rpartition("_")
            tot[clase] = tot.get(clase, 0) + (int(f) - int(i))
    return {k: v / 1e6 for k, v in tot.items()}


def analizar(df, mb, etiqueta, pct, ventana):
    """Devuelve una fila con los conteos y la prueba para un estrato."""
    out = []
    for est in sorted(df["estrato"].unique()):
        s = df[df["estrato"] == est]
        c = s[s["grupo"] == "candidato"]
        k = s[s["grupo"] == "control"]
        if not len(c) or not len(k):
            continue
        bc, tam_c = n_bloques(c[c["pct"] >= pct], ventana)
        bk, _ = n_bloques(k[k["pct"] >= pct], ventana)
        mb_c, mb_k = mb.get("candidato", np.nan), mb.get("control", np.nan)

        tasa_c = bc / mb_c if mb_c == mb_c and mb_c else np.nan
        tasa_k = bk / mb_k if mb_k == mb_k and mb_k else np.nan
        lo_c, hi_c = ic_poisson(bc)
        lo_k, hi_k = ic_poisson(bk)

        # Prueba binomial condicional: dado el total de bloques, la proporcion
        # que cae en candidatos se compara con la esperada por longitud.
        p_esperada = mb_c / (mb_c + mb_k) if (mb_c + mb_k) else np.nan
        if bc + bk > 0 and p_esperada == p_esperada:
            p_val = binomtest(bc, bc + bk, p_esperada, alternative="greater").pvalue
        else:
            p_val = np.nan

        out.append({
            "analisis": etiqueta, "pct": pct, "ventana_kb": ventana // 1000,
            "estrato": est,
            "n_var_cand": len(c), "n_var_ctrl": len(k),
            "bloques_cand": bc, "bloques_ctrl": bk,
            "tasa_cand_Mb": round(tasa_c, 3) if tasa_c == tasa_c else "",
            "IC_cand_Mb": f"{lo_c/mb_c:.2f}-{hi_c/mb_c:.2f}" if mb_c else "",
            "tasa_ctrl_Mb": round(tasa_k, 3) if tasa_k == tasa_k else "",
            "IC_ctrl_Mb": f"{lo_k/mb_k:.2f}-{hi_k/mb_k:.2f}" if mb_k else "",
            "razon": round(tasa_c / tasa_k, 2) if tasa_k else "",
            "p_binomial": p_val,
            "tam_bloque_max": max(tam_c) if tam_c else 0,
        })
    return out


def main():
    if not TODAS.exists():
        sys.exit(f"Falta {TODAS}. Correr antes el script 12.")

    df = pd.read_csv(TODAS, sep="\t", low_memory=False)
    col_a = [c for c in df.columns if c.startswith("AF_IMX")][0]
    col_b = [c for c in df.columns if c.startswith("AF_EUR")][0]
    mb = mb_por_clase()

    di("=" * 78)
    di("PRUEBA FORMAL DEL EXCESO DE BLOQUES ATIPICOS")
    di("=" * 78)
    di(f"\nRegion: candidato {mb.get('candidato', float('nan')):.2f} Mb, "
       f"control {mb.get('control', float('nan')):.2f} Mb")
    di(f"Proporcion esperada de bloques en candidatos por longitud: "
       f"{100*mb.get('candidato',0)/(mb.get('candidato',0)+mb.get('control',1)):.1f}%")

    # --- conjunto completo y conjunto sin alelos casi privados ------------
    privadas = ((df[col_a] < UMBRAL_PRIVADO) | (df[col_b] < UMBRAL_PRIVADO))
    di(f"\nVariantes con AF < {UMBRAL_PRIVADO} en alguna ancestria: "
       f"{privadas.sum():,} de {len(df):,} ({100*privadas.mean():.1f}%)")
    di("Esas variantes producen FST alto por el cuello de botella fundador,")
    di("no necesariamente por diferenciacion local.")

    conjuntos = [("completo", df), ("sin_casi_privadas", df[~privadas])]

    filas = []
    for etq, d in conjuntos:
        for pct in UMBRALES_PCT:
            for v in VENTANAS:
                filas.extend(analizar(d, mb, etq, pct, v))

    t = pd.DataFrame(filas)
    RES.mkdir(parents=True, exist_ok=True)
    t.to_csv(RES / "exceso_pruebas.tsv", sep="\t", index=False)

    # --- resultado principal ----------------------------------------------
    for etq, _ in conjuntos:
        di("\n" + "=" * 78)
        di(f"CONJUNTO: {etq}")
        di("=" * 78)
        s = t[(t["analisis"] == etq)]
        for est in sorted(s["estrato"].unique()):
            di(f"\nEstrato {est}:")
            di(f"{'pct':>6}{'vent':>6}{'blq cand':>10}{'blq ctrl':>10}"
               f"{'tasa cand/Mb':>14}{'IC95':>14}{'tasa ctrl/Mb':>14}"
               f"{'razon':>8}{'p':>10}")
            for _, r in s[s["estrato"] == est].iterrows():
                p = r["p_binomial"]
                di(f"{r['pct']:>6.1f}{r['ventana_kb']:>6}"
                   f"{r['bloques_cand']:>10}{r['bloques_ctrl']:>10}"
                   f"{str(r['tasa_cand_Mb']):>14}{r['IC_cand_Mb']:>14}"
                   f"{str(r['tasa_ctrl_Mb']):>14}{str(r['razon']):>8}"
                   f"{p:>10.3f}" if p == p else "")

    # --- lectura -----------------------------------------------------------
    di("\n" + "=" * 78)
    di("LECTURA")
    di("=" * 78)
    princ = t[(t["pct"] == 99.0) & (t["ventana_kb"] == 50)]
    for _, r in princ.iterrows():
        sig = "si" if r["p_binomial"] < 0.05 else "no"
        di(f"\n{r['analisis']:<20} estrato {r['estrato']:<5} "
           f"{r['bloques_cand']} vs {r['bloques_ctrl']} bloques, "
           f"razon {r['razon']}, p = {r['p_binomial']:.3f}  -> significativo: {sig}")

    di("\nSi el exceso desaparece al excluir las variantes casi privadas de una")
    di("ancestria, se explica por historia demografica (cuello de botella")
    di("fundador) y no por diferenciacion particular de estos genes.")
    di("\nLos intervalos de Poisson muestran la incertidumbre real: con pocos")
    di("bloques, razones aparentemente grandes son compatibles con ausencia")
    di("de exceso.")

    (RES / "exceso_resumen.txt").write_text("\n".join(salida))
    print(f"\n\nPruebas -> {RES / 'exceso_pruebas.tsv'}")
    print(f"Resumen -> {RES / 'exceso_resumen.txt'}")


if __name__ == "__main__":
    main()
