#!/usr/bin/env python3
"""
Reproduccion independiente de las cifras publicadas.

QUE HACE
--------
Recalcula desde los datos crudos, sin importar ninguno de los scripts del
canal, las cifras principales del manuscrito y las compara con las tablas
suplementarias. Lee solo dos archivos de entrada, los mismos que leeria quien
descargara el deposito: las frecuencias procesadas de MCPS y el BED de
regiones. Todo lo demas lo reimplementa a partir de lo que dice la seccion de
Metodos.

Si una cifra publicada no se puede reproducir asi, hay dos posibilidades y
ambas importan antes de publicar: que los Metodos no basten para reproducirla,
o que la cifra este mal. El script no decide cual de las dos es; dice en que
difiere y cuanto.

QUE NO ES
---------
No es una replica ciega. Quien escribe este script ha visto el canal original,
de modo que no demuestra que los Metodos sean suficientes para un lector que
parta de cero. Lo que demuestra es algo mas estrecho y aun asi util: que las
cifras publicadas se siguen de los datos de entrada por un camino distinto del
que las produjo, escrito aparte y sin reutilizar una sola linea de aquel.

Tampoco verifica lo que no esta en los datos de entrada: la curaduria alelica
de la Tabla 1, que viene de literatura y de dbSNP, y la simulacion del
estimador, que no usa datos de MCPS.

Correr: python3 scripts/28_reproduccion_independiente.py
No escribe nada. Solo imprime el resultado de la comparacion.
"""

import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
PROC = RAIZ / "data" / "processed"
RAW = RAIZ / "data" / "raw"
SUP = RAIZ / "results" / "suplementario"

# Parametros declarados en Metodos
MAF_MIN_BARRIDO = 0.01
N_VECINOS_BARRIDO = 300
PCT_UMBRAL = 99.0
VENTANA_PB = 50_000
N_VECINOS_DIRIGIDO = 300
DIST_MAX_DIRIGIDO = 0.02
MIN_VECINOS_DIRIGIDO = 50

ok = fallo = 0


def chk(etiqueta, obtenido, publicado, tol=0.0):
    """Compara una cifra recalculada con la publicada."""
    global ok, fallo
    try:
        igual = abs(float(obtenido) - float(publicado)) <= tol
    except (TypeError, ValueError):
        igual = str(obtenido) == str(publicado)
    if igual:
        ok += 1
        print(f"  coincide  {etiqueta:<52} {obtenido}")
    else:
        fallo += 1
        print(f"  DIFIERE   {etiqueta:<52} recalculado {obtenido}  "
              f"publicado {publicado}")


def fst_hudson(p1, n1, p2, n2):
    """Hudson en la formulacion de Bhatia et al. 2013, sin truncar."""
    with np.errstate(invalid="ignore", divide="ignore"):
        num = ((p1 - p2) ** 2
               - p1 * (1 - p1) / np.maximum(n1 - 1, 1e-12)
               - p2 * (1 - p2) / np.maximum(n2 - 1, 1e-12))
        den = p1 * (1 - p2) + p2 * (1 - p1)
        return np.where(den > 0, num / den, np.nan)


def percentil_en_vecindario(valores_vecinos, valor):
    """Porcentaje de vecinos con FST estrictamente menor."""
    v = np.asarray(valores_vecinos)
    v = v[~np.isnan(v)]
    if not len(v):
        return np.nan
    return 100.0 * np.sum(v < valor) / len(v)


def vecinos_por_posicion(orden_fondo, pos, k):
    """Los k vecinos mas proximos en una lista ordenada, por posicion."""
    n = len(orden_fondo)
    i = int(np.searchsorted(orden_fondo, pos))
    lo, hi = max(0, i - k // 2), min(n, i + k // 2)
    if hi - lo < k:
        if lo == 0:
            hi = min(n, k)
        elif hi == n:
            lo = max(0, n - k)
    return lo, hi


def bloques(df, ventana=VENTANA_PB):
    """Encadena variantes consecutivas mientras el hueco no supere la ventana."""
    if not len(df):
        return 0
    d = df.sort_values(["chrom", "pos"])
    cort = (d["chrom"].astype(str).ne(d["chrom"].astype(str).shift())
            | (d["pos"].diff() > ventana))
    return int(cort.sum())


def binomial_unilateral(k, n, p):
    """P(X >= k) con X binomial(n, p), por suma directa."""
    from math import comb
    return float(sum(comb(n, i) * p ** i * (1 - p) ** (n - i)
                     for i in range(k, n + 1)))


def main():
    f_mc, f_bed = PROC / "mcps_frecuencias.tsv", RAW / "regiones.bed"
    for f in (f_mc, f_bed):
        if not f.exists():
            sys.exit(f"Falta {f}")
    mc = pd.read_csv(f_mc, sep="\t", low_memory=False)

    def tabla(n):
        f = SUP / n
        if not f.exists():
            sys.exit(f"Falta {f}")
        with open(f, encoding="utf-8") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    s1, s3 = tabla("Tabla_S1_resultados_por_locus.tsv"), tabla("Tabla_S3_rejilla_barrido.tsv")

    # El archivo de frecuencias no siempre trae rsid. Se empareja por la
    # primera clave disponible que exista en los dos lados.
    def elegir_clave(cols_mc, filas_s1):
        for c in ("rsid", "variant_id"):
            if c in cols_mc and filas_s1 and filas_s1[0].get(c):
                return c
        return None

    clave = elegir_clave(set(mc.columns), s1)
    if clave is None and {"chrom", "pos"} <= set(mc.columns) and s1[0].get("variant_id"):
        clave = "_pos"
        mc["_pos"] = (mc["chrom"].astype(str) + "-" + mc["pos"].astype(str))
        for r in s1:
            partes = str(r["variant_id"]).split("-")
            r["_pos"] = f"{partes[0]}-{partes[1]}" if len(partes) >= 2 else ""
    print(f"Clave de emparejamiento de loci: {clave or 'ninguna disponible'}")

    # --- columnas de la version principal de frecuencias ---
    ver = "PLANA"
    a, b = "IMX", "EUR"
    afa, ana = f"AF_{a}_{ver}", f"AN_{a}_{ver}"
    afb, anb = f"AF_{b}_{ver}", f"AN_{b}_{ver}"
    aca = f"AC_{a}_{ver}"
    falta = [c for c in (afa, ana, afb, anb) if c not in mc.columns]
    if falta:
        sys.exit(f"El archivo de frecuencias no trae {falta}. Columnas: "
                 f"{sorted(c for c in mc.columns if c.startswith(('AF_','AN_','AC_')))[:12]}")

    print("=" * 84)
    print("REPRODUCCION INDEPENDIENTE DESDE LOS DATOS CRUDOS")
    print("=" * 84)
    print(f"\nEntradas: {f_mc.name} ({len(mc):,} filas) y {f_bed.name}")

    # ---------- 1. frecuencias y FST de los 24 loci ----------
    print("\n--- Los 24 loci de la Tabla S1 ---")
    por_rs = {r[clave]: r for r in s1} if clave else {}
    sub = mc[mc[clave].isin(por_rs)] if clave else pd.DataFrame()
    if not len(sub):
        print("  (no se pudo emparejar los loci con el archivo de frecuencias; se omite)")
    else:
        dif_af = dif_fst = 0.0
        n_comp = 0
        for _, r in sub.iterrows():
            p = por_rs.get(str(r[clave]))
            if p is None:
                continue
            if aca in mc.columns and float(r[ana]) > 0:
                af_i = float(r[aca]) / float(r[ana])
                dif_af = max(dif_af, abs(af_i - float(p["AF_IMX"])))
            f = float(fst_hudson(np.array([float(r[afa])]), np.array([float(r[ana])]),
                                 np.array([float(r[afb])]), np.array([float(r[anb])]))[0])
            dif_fst = max(dif_fst, abs(f - float(p["FST"])))
            n_comp += 1
        chk(f"frecuencias recalculadas desde AC/AN ({n_comp} loci), dif maxima",
            f"{dif_af:.2e}", 0.0, tol=1e-5)
        chk(f"FST de Hudson recalculado ({n_comp} loci), dif maxima",
            f"{dif_fst:.2e}", 0.0, tol=1e-4)

    # ---------- 2. conjunto analitico del barrido ----------
    print("\n--- Conjunto analitico del barrido ---")
    df = mc.copy()
    for col, val in (("reparto_default", 0), ("es_indel", 0)):
        if col in df.columns:
            df = df[df[col] == val]
    df = df.dropna(subset=[afa, afb, ana, anb])
    df = df[(df[ana] > 0) & (df[anb] > 0)].copy()
    df["maf_media"] = ((np.minimum(df[afa], 1 - df[afa])
                        + np.minimum(df[afb], 1 - df[afb])) / 2)
    df = df[df["maf_media"] >= MAF_MIN_BARRIDO].copy()
    df["estrato"] = np.where(df["source"] == "WGS", "wgs", "alta")
    df["FST"] = fst_hudson(df[afa].values, df[ana].values,
                           df[afb].values, df[anb].values)
    df = df.dropna(subset=["FST"]).reset_index(drop=True)

    pub = {r["conjunto"]: r for r in s3
           if r["umbral_percentil"] == "99.0" and r["ventana_kb"] == "50"}
    tot_pub = {}
    for c in ("completo", "restringido_asimetrico", "restringido_simetrico"):
        fs = [r for r in s3 if r["conjunto"] == c
              and r["umbral_percentil"] == "99.0" and r["ventana_kb"] == "50"]
        tot_pub[c] = sum(int(r["n_var_cand"]) + int(r["n_var_ctrl"]) for r in fs)

    asim = df[(df[afa] >= 0.02) & (df[afb] >= 0.02)]
    sim = df[(np.minimum(df[afa], 1 - df[afa]) >= 0.02)
             & (np.minimum(df[afb], 1 - df[afb]) >= 0.02)]
    chk("variantes evaluadas, conjunto completo", len(df), tot_pub["completo"])
    chk("variantes, criterio asimetrico", len(asim), tot_pub["restringido_asimetrico"])
    chk("variantes, criterio simetrico", len(sim), tot_pub["restringido_simetrico"])
    chk("diferencia entre criterios", len(asim) - len(sim),
        tot_pub["restringido_asimetrico"] - tot_pub["restringido_simetrico"])

    # ---------- 3. bloques atipicos por conjunto ----------
    print("\n--- Bloques atipicos al percentil 99, ventana 50 kb ---")
    for nombre, conj in (("completo", df), ("restringido_asimetrico", asim),
                         ("restringido_simetrico", sim)):
        p = pub.get(nombre)
        if p is None:
            continue
        for est in ("wgs",):
            s = conj[conj["estrato"] == est]
            ctrl = s[s["clase"] == "control"]
            cand = s[s["clase"] == "candidato"]
            if len(ctrl) < N_VECINOS_BARRIDO or not len(cand):
                continue
            ordc = ctrl.sort_values("maf_media")
            maf_c = ordc["maf_media"].values
            fst_c = ordc["FST"].values

            pct_cand = []
            for m, f in zip(cand["maf_media"].values, cand["FST"].values):
                lo, hi = vecinos_por_posicion(maf_c, m, N_VECINOS_BARRIDO)
                pct_cand.append(percentil_en_vecindario(fst_c[lo:hi], f))
            at_cand = cand[np.array(pct_cand) >= PCT_UMBRAL]

            pct_ctrl = np.empty(len(ordc))
            genes = ordc["gen"].values
            for g in pd.unique(genes):
                m_g = genes == g
                otros_maf = maf_c[~m_g]
                otros_fst = fst_c[~m_g]
                for i in np.flatnonzero(m_g):
                    lo, hi = vecinos_por_posicion(otros_maf, maf_c[i],
                                                  N_VECINOS_BARRIDO)
                    pct_ctrl[i] = percentil_en_vecindario(otros_fst[lo:hi],
                                                          fst_c[i])
            at_ctrl = ordc[pct_ctrl >= PCT_UMBRAL]

            b_cand, b_ctrl = bloques(at_cand), bloques(at_ctrl)
            chk(f"{nombre}: bloques candidatos ({est})", b_cand,
                int(p["bloques_cand"]))
            chk(f"{nombre}: bloques control ({est})", b_ctrl,
                int(p["bloques_ctrl"]))
            if nombre == "completo" and b_cand + b_ctrl > 0:
                r = pd.read_csv(f_bed, sep="\t", header=None,
                                names=["chrom", "ini", "fin", "nombre"])
                r["mb"] = (r["fin"] - r["ini"]) / 1e6
                es_cand = r["nombre"].astype(str).str.contains("candidato",
                                                              case=False)
                mb_c, mb_k = r.loc[es_cand, "mb"].sum(), r.loc[~es_cand, "mb"].sum()
                esperada = mb_c / (mb_c + mb_k)
                chk("proporcion esperada por longitud", round(esperada, 4),
                    0.0343, tol=0.0006)
                pv = binomial_unilateral(b_cand, b_cand + b_ctrl, esperada)
                chk("p binomial unilateral, conjunto completo", round(pv, 3),
                    round(float(p["p_binomial_unilateral"]), 3), tol=0.002)
                razon = (b_cand / mb_c) / (b_ctrl / mb_k)
                chk("razon de tasas, conjunto completo", round(razon, 2),
                    round(float(p["razon"]), 2), tol=0.02)

    # ---------- 4. percentiles del analisis dirigido ----------
    print("\n--- Analisis dirigido: percentiles de los 24 loci ---")
    d = mc.copy()
    for col, val in (("reparto_default", 0), ("es_indel", 0)):
        if col in d.columns:
            d = d[d[col] == val]
    d = d.dropna(subset=[afa, afb, ana, anb])
    d = d[(d[ana] > 0) & (d[anb] > 0)].copy()
    d["pbar"] = (d[afa] + d[afb]) / 2
    d["maf_ap"] = np.minimum(d["pbar"], 1 - d["pbar"])
    d["estrato"] = np.where(d["source"] == "WGS", "wgs", "alta")
    d["FST"] = fst_hudson(d[afa].values, d[ana].values,
                          d[afb].values, d[anb].values)
    d = d.dropna(subset=["FST"])
    # Los Metodos lo dicen: cada locus se compara con las variantes de control
    # de frecuencia mas cercana DENTRO DE SU MISMO ESTRATO. Ignorarlo cambia
    # los percentiles hasta en 30 puntos.
    fondo_por_estrato = {e: g.sort_values("maf_ap")
                         for e, g in d[d["clase"] == "control"].groupby("estrato")}

    recalc, n_ok = {}, 0
    for r in s1:
        fila = (d[(d["clase"] == "candidato") & (d[clave] == r[clave])]
                if clave else pd.DataFrame())
        if not len(fila):
            continue
        est = fila.iloc[0]["estrato"]
        fondo = fondo_por_estrato.get(est)
        if fondo is None or not len(fondo):
            continue
        maf_k, fst_k = fondo["maf_ap"].values, fondo["FST"].values
        m = float(fila.iloc[0]["maf_ap"])
        f = float(fila.iloc[0]["FST"])
        cerca = np.abs(maf_k - m) <= DIST_MAX_DIRIGIDO
        if cerca.sum() < MIN_VECINOS_DIRIGIDO:
            continue
        idx = np.flatnonzero(cerca)
        if len(idx) > N_VECINOS_DIRIGIDO:
            orden = np.argsort(np.abs(maf_k[idx] - m))[:N_VECINOS_DIRIGIDO]
            idx = idx[orden]
        recalc[r[clave]] = percentil_en_vecindario(fst_k[idx], f)
        n_ok += 1

    if recalc:
        difs = [abs(recalc[k] - float(por_rs[k]["percentil"])) for k in recalc]
        chk(f"percentiles del analisis dirigido ({n_ok} loci), dif maxima",
            round(max(difs), 3), 0.0, tol=0.05)
        med = float(np.median(list(recalc.values())))
        chk("mediana de los percentiles recalculados", round(med, 1),
            round(float(np.median([float(r["percentil"]) for r in s1])), 1),
            tol=0.1)
    else:
        print("  (no se pudo emparejar ningun locus candidato)")

    print("\n" + "=" * 84)
    print(f"{ok} cifras reproducidas, {fallo} con diferencia")
    print("=" * 84)
    if fallo:
        print("\nUna diferencia no significa necesariamente que la cifra "
              "publicada este mal.\nPuede venir de una decision que los Metodos "
              "no fijan con suficiente detalle,\ny eso tambien conviene saberlo "
              "antes de depositar.")
    sys.exit(1 if fallo else 0)


if __name__ == "__main__":
    main()
