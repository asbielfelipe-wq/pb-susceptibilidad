#!/usr/bin/env python3
"""
Sensibilidades recalculadas bajo el esquema del analisis principal.

MOTIVO
------
Las sensibilidades previas se calcularon con la definicion de apareamiento
antigua, el promedio de las frecuencias menores de cada componente, mientras
que el analisis principal usa min(pbar, 1-pbar). Comparar unas con otro mezcla
dos esquemas y no verifica lo que pretende verificar. Aqui todo se recalcula
con la definicion principal, dejando fija cualquier otra decision salvo la que
se somete a prueba en cada caso.

Incluye ademas:
  - intervalo condicional de la razon de tasas del barrido, calculado por
    inversion de la prueba binomial sobre el total de bloques
  - comparacion del remuestreo por variante y por bloque bajo el mismo fondo
    y el mismo apareamiento
  - subconjuntos de control con semilla declarada y numero de loci con
    comparadores suficientes en cada uno

Correr: python3 scripts/20_sensibilidades_v2.py
Salidas: results/sens2_*.tsv y results/sensibilidades_v2.txt
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom, binomtest

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
PROC = RAIZ / "data" / "processed"
RAW = RAIZ / "data" / "raw"

MAF_MIN, D_MAX = 0.01, 0.02
N_OBJETIVO, N_MINIMO = 300, 50
BLOQUE_KB, N_BOOT = 50_000, 2000
SEMILLA = 20261006

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def fst_hudson(p1, n1, p2, n2):
    with np.errstate(invalid="ignore", divide="ignore"):
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)
        return np.where(den > 0, num / den, np.nan)


def bloques(df):
    d = df.sort_values(["chrom", "pos"]).copy()
    b, actual, pc, pp = [], 0, None, None
    for c, p in zip(d["chrom"].astype(str), d["pos"]):
        if pc != c or (p - pp) > BLOQUE_KB:
            actual += 1
        b.append(actual)
        pc, pp = c, p
    d["bloque"] = b
    return d


def preparar(mc, version, estratos):
    af_a, an_a = f"AF_IMX_{version}", f"AN_IMX_{version}"
    af_b, an_b = f"AF_EUR_{version}", f"AN_EUR_{version}"
    if af_a not in mc.columns:
        return None
    d = mc[(mc["reparto_default"] == 0) & (mc["es_indel"] == 0)].copy()
    d = d.dropna(subset=[af_a, af_b, an_a, an_b])
    d = d[(d[an_a] > 0) & (d[an_b] > 0)]
    p1, p2 = d[af_a].values, d[af_b].values
    pb = (p1 + p2) / 2
    d["maf"] = np.minimum(pb, 1 - pb)          # definicion principal
    d["FST"] = fst_hudson(p1, d[an_a].values, p2, d[an_b].values)
    d["estrato"] = (np.where(d["source"] == "WGS", "wgs", "alta")
                    if estratos == "agrupado" else d["source"].str.lower())
    return d[d["maf"] >= MAF_MIN].dropna(subset=["FST"])


def percentiles(d, cand, genes=None, boot=None, rng=None):
    fondo = d[d["clase"] == "control"]
    if genes is not None:
        fondo = fondo[fondo["gen"].isin(genes)]
    fondo = bloques(fondo)
    out = []
    for _, c in cand.iterrows():
        f = d[d["variant_id"] == c["variant_id"]]
        if not len(f):
            out.append({"rsid": c["rsid"], "gen": c["gen"], "pct": np.nan, "n_vec": 0})
            continue
        r = f.iloc[0]
        g = fondo[fondo["estrato"] == r["estrato"]]
        vec = g[(g["maf"] - r["maf"]).abs() <= D_MAX]
        if len(vec) > N_OBJETIVO:
            vec = vec.loc[(vec["maf"] - r["maf"]).abs().nsmallest(N_OBJETIVO).index]
        fila = {"rsid": c["rsid"], "gen": c["gen"], "estrato": r["estrato"],
                "n_vec": len(vec)}
        if len(vec) < N_MINIMO:
            fila["pct"] = np.nan
            out.append(fila)
            continue
        fila["pct"] = round((vec["FST"].values < r["FST"]).mean() * 100, 1)
        if boot and rng is not None:
            fst_v = vec["FST"].values
            if boot == "bloque":
                grupos = [x["FST"].values for _, x in vec.groupby("bloque")]
                nb = len(grupos)
                ps = np.array([(np.concatenate([grupos[j] for j in
                                rng.integers(0, nb, nb)]) < r["FST"]).mean() * 100
                               for _ in range(N_BOOT)])
            else:
                idx = rng.integers(0, len(fst_v), size=(N_BOOT, len(fst_v)))
                ps = (fst_v[idx] < r["FST"]).mean(axis=1) * 100
            lo, hi = np.percentile(ps, [2.5, 97.5])
            fila["IC_inf"], fila["IC_sup"] = round(lo, 1), round(hi, 1)
        out.append(fila)
    return pd.DataFrame(out)


def ic_razon(bc, bk, mb_c, mb_k, alpha=0.05):
    """Intervalo condicional de la razon de tasas por inversion de la prueba
    binomial sobre el total de bloques: dado n = bc + bk, se buscan los valores
    de la razon r para los que bc no se rechaza, con p = r*mb_c/(r*mb_c + mb_k)."""
    n = bc + bk
    if n == 0:
        return (np.nan, np.nan)

    def acepta(r):
        p = r * mb_c / (r * mb_c + mb_k)
        lo = binom.ppf(alpha / 2, n, p)
        hi = binom.ppf(1 - alpha / 2, n, p)
        return lo <= bc <= hi

    rs = np.exp(np.linspace(np.log(1e-3), np.log(1e3), 4000))
    ok = [r for r in rs if acepta(r)]
    return (min(ok), max(ok)) if ok else (np.nan, np.nan)


def main():
    mc = pd.read_csv(PROC / "mcps_frecuencias.tsv", sep="\t", low_memory=False)
    cand = pd.read_csv(RES / "principal_percentiles.tsv", sep="\t")
    rng = np.random.default_rng(SEMILLA)

    di("=" * 78)
    di("SENSIBILIDADES BAJO EL ESQUEMA DEL ANALISIS PRINCIPAL")
    di("=" * 78)
    di(f"\nApareamiento: min(pbar, 1-pbar) en todos los casos")
    di(f"Semilla de remuestreo y de seleccion de subconjuntos: {SEMILLA}")

    base = cand.set_index("rsid")["pct"]
    etiq = dict(zip(cand["rsid"], cand["gen"] + " " + cand["rsid"]))

    # ---- A: versiones de frecuencia, con la definicion principal ----
    di("\n" + "=" * 78)
    di("A. VERSIONES DE FRECUENCIA DE MCPS")
    di("=" * 78)
    cols = {}
    for v in ["PLANA", "UNREL", "IBD"]:
        d = preparar(mc, v, "agrupado")
        if d is None:
            di(f"\nversion {v}: columnas no disponibles")
            continue
        t = percentiles(d, cand)
        cols[v] = t.set_index("rsid")["pct"]
        ev = t.dropna(subset=["pct"])
        di(f"\n{v}: {len(ev)} evaluables, percentil >=95: "
           f"{int((ev['pct'] >= 95).sum())}, mediana {ev['pct'].median():.1f}")
    tv = pd.DataFrame(cols)
    tv["etiqueta"] = tv.index.map(etiq)
    tv["rango"] = tv[list(cols)].max(axis=1) - tv[list(cols)].min(axis=1)
    di(f"\n{'locus':<22}" + "".join(f"{k:>9}" for k in cols) + f"{'rango':>8}")
    for _, r in tv.sort_values("rango", ascending=False).iterrows():
        di(f"{r['etiqueta']:<22}" +
           "".join(f"{r[k]:>9.1f}" if r[k] == r[k] else f"{'-':>9}" for k in cols) +
           f"{r['rango']:>8.1f}")
    di(f"\nCambio maximo entre versiones: {tv['rango'].max():.1f} puntos")
    tv.to_csv(RES / "sens2_versiones.tsv", sep="\t")

    # ---- B: estratos separados ----
    di("\n" + "=" * 78)
    di("B. ARREGLO Y EXOMA SEPARADOS")
    di("=" * 78)
    dS = preparar(mc, "PLANA", "separado")
    tS = percentiles(dS, cand, boot="bloque", rng=rng)
    tS["etiqueta"] = tS["rsid"].map(etiq)
    tS["pct_principal"] = tS["rsid"].map(base)
    di(f"\n{'locus':<22}{'estrato':<8}{'vec':>5}{'pct sep':>9}{'pct ppal':>10}"
       f"{'cambio':>8}")
    for _, r in tS.sort_values("pct", ascending=False, na_position="last").iterrows():
        if r["pct"] != r["pct"]:
            di(f"{r['etiqueta']:<22}{'-':<8}{r['n_vec']:>5}{'no eval':>9}"
               f"{r['pct_principal']:>10.1f}")
            continue
        di(f"{r['etiqueta']:<22}{r['estrato']:<8}{r['n_vec']:>5}{r['pct']:>9.1f}"
           f"{r['pct_principal']:>10.1f}{r['pct']-r['pct_principal']:>+8.1f}")
    ev = tS.dropna(subset=["pct"])
    di(f"\nEvaluables {len(ev)} de {len(tS)}; percentil >=95: "
       f"{int((ev['pct'] >= 95).sum())}")
    di("Los loci del estrato de genoma completo no cambian, porque la division")
    di("afecta solo a arreglo y exoma.")
    tS.to_csv(RES / "sens2_plataforma.tsv", sep="\t", index=False)

    # ---- C: remuestreo por variante frente a por bloque, mismo fondo ----
    di("\n" + "=" * 78)
    di("C. REMUESTREO POR VARIANTE FRENTE A POR BLOQUE")
    di("=" * 78)
    d = preparar(mc, "PLANA", "agrupado")
    rng_v = np.random.default_rng(SEMILLA)
    tvar = percentiles(d, cand, boot="variante", rng=rng_v)
    rng_b = np.random.default_rng(SEMILLA)
    tblq = percentiles(d, cand, boot="bloque", rng=rng_b)
    amp_v = (tvar["IC_sup"] - tvar["IC_inf"]).mean()
    amp_b = (tblq["IC_sup"] - tblq["IC_inf"]).mean()
    di(f"\nAmplitud media del IC95 con el mismo fondo y apareamiento:")
    di(f"  remuestreo por variante: {amp_v:.1f} puntos")
    di(f"  remuestreo por bloque:   {amp_b:.1f} puntos")
    di(f"  razon: {amp_b/amp_v:.2f}")
    comp = tvar[["rsid"]].copy()
    comp["etiqueta"] = comp["rsid"].map(etiq)
    comp["amp_variante"] = tvar["IC_sup"] - tvar["IC_inf"]
    comp["amp_bloque"] = tblq["IC_sup"] - tblq["IC_inf"]
    comp.to_csv(RES / "sens2_remuestreo.tsv", sep="\t", index=False)

    # ---- D: subconjuntos de control, con semilla y elegibilidad ----
    di("\n" + "=" * 78)
    di("D. SUBCONJUNTOS DE GENES CONTROL")
    di("=" * 78)
    genes = sorted(d[d["clase"] == "control"]["gen"].dropna().unique())
    di(f"\n{len(genes)} genes control. Las mitades se obtienen por muestreo sin")
    di(f"reemplazo de {len(genes)//2} genes con semilla {SEMILLA}; cambian a la vez")
    di("el tamano y la composicion del fondo, de modo que el efecto observado no")
    di("es atribuible solo al numero de genes.")
    rng_s = np.random.default_rng(SEMILLA)
    sub = {"completo": genes}
    for i in range(3):
        sub[f"mitad_{i+1}"] = sorted(rng_s.choice(genes, len(genes)//2, replace=False))
    f74 = RAW / "genes_control_v1_74.tsv"
    if f74.exists():
        g74 = pd.read_csv(f74, sep="\t")["gen"].tolist()
        sub["v1_74genes"] = [g for g in g74 if g in genes]
    filas = []
    for nombre, gs in sub.items():
        t = percentiles(d, cand, genes=gs)
        ev = t.dropna(subset=["pct"])
        filas.append({"subconjunto": nombre, "n_genes": len(gs),
                      "loci_con_comparadores": len(ev),
                      "pct95": int((ev["pct"] >= 95).sum()),
                      "mediana": round(float(ev["pct"].median()), 1),
                      "maximo": round(float(ev["pct"].max()), 1),
                      "AGT_rs2493137": float(t[t["rsid"] == "rs2493137"]["pct"].iloc[0])
                      if len(t[t["rsid"] == "rs2493137"]) else None})
    ts = pd.DataFrame(filas)
    di(f"\n{'subconjunto':<14}{'genes':>7}{'con compar.':>13}{'pct>=95':>9}"
       f"{'mediana':>9}{'maximo':>8}{'rs2493137':>11}")
    for _, r in ts.iterrows():
        v = f"{r['AGT_rs2493137']:.1f}" if r["AGT_rs2493137"] == r["AGT_rs2493137"] else "-"
        di(f"{r['subconjunto']:<14}{r['n_genes']:>7}{r['loci_con_comparadores']:>13}"
           f"{r['pct95']:>9}{r['mediana']:>9}{r['maximo']:>8}{v:>11}")
    ts.to_csv(RES / "sens2_controles.tsv", sep="\t", index=False)

    # ---- E: intervalo condicional de la razon de tasas ----
    di("\n" + "=" * 78)
    di("E. INTERVALO CONDICIONAL DE LA RAZON DE TASAS DEL BARRIDO")
    di("=" * 78)
    di("\nSe obtiene por inversion de la prueba binomial: dado el total de")
    di("bloques, se conservan los valores de la razon para los que el conteo")
    di("candidato no se rechaza al 5%.")
    mb_c, mb_k = 0.494, 13.91
    di(f"\n{'comparacion':<24}{'razon':>8}{'IC95 condicional':>22}")
    for etq, bc, bk in [("4 frente a 51 bloques", 4, 51),
                        ("2 frente a 37 bloques", 2, 37),
                        ("4 frente a 49 bloques", 4, 49)]:
        lo, hi = ic_razon(bc, bk, mb_c, mb_k)
        razon = (bc / mb_c) / (bk / mb_k)
        di(f"{etq:<24}{razon:>8.2f}{f'{lo:.2f}-{hi:.2f}':>22}")

    di("\n" + "=" * 78)
    di("RESUMEN")
    di("=" * 78)
    di(f"\nPercentil de AGT rs2493137 bajo cada condicion:")
    di(f"  analisis principal                 {base.get('rs2493137', float('nan')):.1f}")
    for v in cols:
        val = cols[v].get("rs2493137")
        di(f"  version {v:<26}{val:.1f}" if val == val else f"  version {v:<26}-")
    sep_val = tS[tS["rsid"] == "rs2493137"]["pct"]
    if len(sep_val):
        v = sep_val.iloc[0]
        di(f"  estratos separados                 "
           f"{v:.1f}" if v == v else "  estratos separados                 no eval")
    for _, r in ts.iterrows():
        if r["AGT_rs2493137"] == r["AGT_rs2493137"]:
            di(f"  control {r['subconjunto']:<26}{r['AGT_rs2493137']:.1f}")

    (RES / "sensibilidades_v2.txt").write_text("\n".join(salida))
    print(f"\n\nResumen -> {RES / 'sensibilidades_v2.txt'}")


if __name__ == "__main__":
    main()
