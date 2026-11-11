#!/usr/bin/env python3
"""
Analisis de sensibilidad solicitados en revision. Comprueba si el resultado
negativo depende de decisiones concretas del diseno.

A. DEFINICION DE LA FRECUENCIA DE APAREAMIENTO
   Se compara la implementada, promedio de las frecuencias menores de cada
   componente, con la alternativa min(pbar, 1-pbar). Divergen cuando el alelo
   mayoritario cambia entre componentes, lo que ocurre en siete loci.

B. SEPARACION DE ARREGLO Y EXOMA
   Se agruparon por equivalencia de tamano muestral efectivo, pero sus
   mecanismos de seleccion de variantes difieren. Se repite el analisis con
   los dos estratos separados.

C. VERSIONES ALTERNATIVAS DE FRECUENCIA
   MCPS publica tres: con todas las muestras (usada), con el subconjunto de no
   emparentados, y corregida por segmentos identicos por descendencia. Se
   repiten los percentiles con las tres.

D. FILTRO SIMETRICO DE ALELOS POCO FRECUENTES
   El filtro original excluye por la frecuencia del alelo alterno, de modo que
   no atrapa los casos en que el escaso es el de referencia. Aqui se aplica el
   criterio simetrico: se excluye si min(p, 1-p) < umbral en algun componente.

E. SUBCONJUNTOS DE CONTROLES
   Se repiten los percentiles con mitades aleatorias del conjunto control y con
   los 74 genes originales, para comprobar que la conclusion no depende de unos
   pocos genes.

F. PROPIEDADES DE LOS CONJUNTOS
   Longitud, densidad de variantes y proporcion codificante en candidatos y
   controles.

G. DETALLE POR LOCUS
   Vecinos, bloques independientes y distancia maxima de frecuencia alcanzada.

Correr: python3 scripts/18_sensibilidades.py
Salidas: results/sensibilidad_*.tsv y results/sensibilidades.txt
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
PROC = RAIZ / "data" / "processed"
RAW = RAIZ / "data" / "raw"

D_MAX, N_OBJETIVO, N_MINIMO = 0.02, 300, 50
BLOQUE_KB, N_BOOT = 50_000, 2000
MAF_MIN, UMBRAL_PRIVADO = 0.01, 0.02

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def fst_hudson(p1, n1, p2, n2):
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


def preparar(mc, version, maf_def, estrato_modo):
    af_a, an_a = f"AF_IMX_{version}", f"AN_IMX_{version}"
    af_b, an_b = f"AF_EUR_{version}", f"AN_EUR_{version}"
    if af_a not in mc.columns:
        return None
    d = mc[(mc["reparto_default"] == 0) & (mc["es_indel"] == 0)].copy()
    d = d.dropna(subset=[af_a, af_b, an_a, an_b])
    d = d[(d[an_a] > 0) & (d[an_b] > 0)]
    p1, p2 = d[af_a].values, d[af_b].values
    if maf_def == "A":
        d["maf"] = (np.minimum(p1, 1 - p1) + np.minimum(p2, 1 - p2)) / 2
    else:
        pb = (p1 + p2) / 2
        d["maf"] = np.minimum(pb, 1 - pb)
    d["FST"] = fst_hudson(p1, d[an_a].values, p2, d[an_b].values)
    if estrato_modo == "agrupado":
        d["estrato"] = np.where(d["source"] == "WGS", "wgs", "alta")
    else:
        d["estrato"] = d["source"].str.lower()
    d = d[d["maf"] >= MAF_MIN].dropna(subset=["FST"])
    return d


def percentiles(d, cand_ids, boot=False, genes_fondo=None, rng=None):
    fondo = d[d["clase"] == "control"]
    if genes_fondo is not None:
        fondo = fondo[fondo["gen"].isin(genes_fondo)]
    fondo = bloques(fondo)
    out = []
    for vid in cand_ids:
        f = d[d["variant_id"] == vid]
        if not len(f):
            continue
        r = f.iloc[0]
        g = fondo[fondo["estrato"] == r["estrato"]]
        dd = (g["maf"] - r["maf"]).abs()
        vec = g[dd <= D_MAX]
        if len(vec) > N_OBJETIVO:
            vec = vec.loc[(vec["maf"] - r["maf"]).abs().nsmallest(N_OBJETIVO).index]
        if len(vec) < N_MINIMO:
            out.append({"variant_id": vid, "pct": None, "n_vec": len(vec)})
            continue
        pct = (vec["FST"].values < r["FST"]).mean() * 100
        fila = {"variant_id": vid, "pct": round(pct, 1), "n_vec": len(vec),
                "n_blq": int(vec["bloque"].nunique()),
                "d_max": round(float((vec["maf"] - r["maf"]).abs().max()), 5),
                "estrato": r["estrato"], "maf": round(float(r["maf"]), 4)}
        if boot and rng is not None:
            grupos = [x["FST"].values for _, x in vec.groupby("bloque")]
            nb = len(grupos)
            pcts = np.empty(N_BOOT)
            for i in range(N_BOOT):
                sel = rng.integers(0, nb, nb)
                v = np.concatenate([grupos[j] for j in sel])
                pcts[i] = (v < r["FST"]).mean() * 100
            lo, hi = np.percentile(pcts, [2.5, 97.5])
            fila["IC_inf"], fila["IC_sup"] = round(lo, 1), round(hi, 1)
        out.append(fila)
    return pd.DataFrame(out)


def main():
    f_mcps = PROC / "mcps_frecuencias.tsv"
    f_pct = RES / "percentiles_v2.tsv"
    if not f_mcps.exists() or not f_pct.exists():
        sys.exit("Faltan mcps_frecuencias.tsv o percentiles_v2.tsv")

    mc = pd.read_csv(f_mcps, sep="\t", low_memory=False)
    base = pd.read_csv(f_pct, sep="\t")
    cand_ids = base["variant_id"].tolist()
    etiq = dict(zip(base["variant_id"], base["gen"] + " " + base["rsid"]))
    rng = np.random.default_rng(20261006)

    di("=" * 78)
    di("ANALISIS DE SENSIBILIDAD")
    di("=" * 78)

    resumen = []

    # ---- A: definicion de MAF ----
    di("\n" + "=" * 78)
    di("A. DEFINICION DE LA FRECUENCIA DE APAREAMIENTO")
    di("=" * 78)
    dA = preparar(mc, "PLANA", "A", "agrupado")
    dB = preparar(mc, "PLANA", "B", "agrupado")
    pA = percentiles(dA, cand_ids).rename(columns={"pct": "pct_defA"})
    pB = percentiles(dB, cand_ids).rename(columns={"pct": "pct_defB"})
    comp = pA[["variant_id", "pct_defA", "maf"]].merge(
        pB[["variant_id", "pct_defB"]], on="variant_id")
    comp["etiqueta"] = comp["variant_id"].map(etiq)
    comp["cambio"] = comp["pct_defB"] - comp["pct_defA"]
    di(f"\n{'locus':<22}{'maf A':>8}{'pct A':>8}{'pct B':>8}{'cambio':>9}")
    for _, r in comp.sort_values("cambio", key=abs, ascending=False).iterrows():
        if r["pct_defA"] is None or r["pct_defB"] is None:
            continue
        di(f"{r['etiqueta']:<22}{r['maf']:>8.4f}{r['pct_defA']:>8.1f}"
           f"{r['pct_defB']:>8.1f}{r['cambio']:>+9.1f}")
    cam = comp["cambio"].abs()
    di(f"\nCambio medio {cam.mean():.1f} puntos, maximo {cam.max():.1f}")
    di(f"Loci que superarian el percentil 95 con la definicion B: "
       f"{int((comp['pct_defB'] >= 95).sum())}")
    comp.to_csv(RES / "sensibilidad_maf.tsv", sep="\t", index=False)
    resumen.append(("Definicion de MAF (B)", int((comp["pct_defB"] >= 95).sum())))

    # ---- B: arreglo y exoma separados ----
    di("\n" + "=" * 78)
    di("B. ARREGLO Y EXOMA SEPARADOS")
    di("=" * 78)
    dS = preparar(mc, "PLANA", "A", "separado")
    pS = percentiles(dS, cand_ids, boot=True, rng=rng)
    pS["etiqueta"] = pS["variant_id"].map(etiq)
    di(f"\n{'locus':<22}{'estrato':<8}{'vec':>5}{'blq':>5}{'pct':>7}{'IC95':>13}")
    for _, r in pS.sort_values("pct", ascending=False, na_position="last").iterrows():
        if r["pct"] is None or r["pct"] != r["pct"]:
            di(f"{r['etiqueta']:<22}{'-':<8}{r['n_vec']:>5}{'-':>5}{'no eval':>7}")
            continue
        ic = f"[{r['IC_inf']:.0f}-{r['IC_sup']:.0f}]"
        di(f"{r['etiqueta']:<22}{r['estrato']:<8}{r['n_vec']:>5}{r['n_blq']:>5}"
           f"{r['pct']:>7.1f}{ic:>13}")
    ev = pS.dropna(subset=["pct"])
    di(f"\nEvaluables {len(ev)} de {len(pS)}; percentil >=95: "
       f"{int((ev['pct'] >= 95).sum())}")
    pS.to_csv(RES / "sensibilidad_plataforma.tsv", sep="\t", index=False)
    resumen.append(("Arreglo y exoma separados", int((ev["pct"] >= 95).sum())))

    # ---- C: versiones de frecuencia ----
    di("\n" + "=" * 78)
    di("C. VERSIONES ALTERNATIVAS DE FRECUENCIA")
    di("=" * 78)
    tabs = {}
    for v, nombre in [("PLANA", "todas las muestras"),
                      ("UNREL", "no emparentados"),
                      ("IBD", "corregida por IBD")]:
        d = preparar(mc, v, "A", "agrupado")
        if d is None:
            di(f"\nversion {v}: columnas no disponibles")
            continue
        t = percentiles(d, cand_ids)
        tabs[v] = t.set_index("variant_id")["pct"]
        ev = t.dropna(subset=["pct"])
        di(f"\n{nombre}: {len(ev)} evaluables, "
           f"percentil >=95: {int((ev['pct'] >= 95).sum())}, "
           f"mediana {ev['pct'].median():.1f}")
        resumen.append((f"Version {v}", int((ev["pct"] >= 95).sum())))
    if len(tabs) > 1:
        tv = pd.DataFrame(tabs)
        tv["etiqueta"] = tv.index.map(etiq)
        di(f"\n{'locus':<22}" + "".join(f"{k:>10}" for k in tabs))
        for _, r in tv.sort_values(list(tabs)[0], ascending=False).iterrows():
            di(f"{r['etiqueta']:<22}" +
               "".join(f"{r[k]:>10.1f}" if r[k] == r[k] else f"{'-':>10}"
                       for k in tabs))
        tv.to_csv(RES / "sensibilidad_versiones.tsv", sep="\t")

    # ---- D: filtro simetrico ----
    di("\n" + "=" * 78)
    di("D. FILTRO SIMETRICO DE ALELOS POCO FRECUENTES")
    di("=" * 78)
    d = preparar(mc, "PLANA", "A", "agrupado")
    af_a, af_b = "AF_IMX_PLANA", "AF_EUR_PLANA"
    asim = (d[af_a] < UMBRAL_PRIVADO) | (d[af_b] < UMBRAL_PRIVADO)
    sim = ((np.minimum(d[af_a], 1 - d[af_a]) < UMBRAL_PRIVADO) |
           (np.minimum(d[af_b], 1 - d[af_b]) < UMBRAL_PRIVADO))
    di(f"\nVariantes excluidas con el criterio asimetrico (AF_ALT < "
       f"{UMBRAL_PRIVADO}): {int(asim.sum()):,} de {len(d):,} "
       f"({100*asim.mean():.1f}%)")
    di(f"Variantes excluidas con el criterio simetrico (min(p,1-p) < "
       f"{UMBRAL_PRIVADO}): {int(sim.sum()):,} ({100*sim.mean():.1f}%)")
    solo_sim = int((sim & ~asim).sum())
    di(f"Capturadas solo por el simetrico: {solo_sim:,}")
    di("Son variantes donde el alelo escaso es el de referencia, que el filtro")
    di("original no atrapaba.")

    # ---- E: subconjuntos de controles ----
    di("\n" + "=" * 78)
    di("E. SUBCONJUNTOS DE GENES CONTROL")
    di("=" * 78)
    genes_ctrl = sorted(d[d["clase"] == "control"]["gen"].dropna().unique())
    di(f"\n{len(genes_ctrl)} genes control disponibles")
    subconjuntos = {"completo": genes_ctrl}
    for i in range(3):
        sel = list(rng.choice(genes_ctrl, size=len(genes_ctrl) // 2, replace=False))
        subconjuntos[f"mitad_{i+1}"] = sel
    f74 = RAW / "genes_control_v1_74.tsv"
    if f74.exists():
        g74 = pd.read_csv(f74, sep="\t")["gen"].tolist()
        subconjuntos["v1_74genes"] = [g for g in g74 if g in genes_ctrl]
    filas = []
    for nombre, genes in subconjuntos.items():
        t = percentiles(d, cand_ids, genes_fondo=genes)
        ev = t.dropna(subset=["pct"])
        filas.append({"subconjunto": nombre, "n_genes": len(genes),
                      "n_evaluables": len(ev),
                      "n_pct95": int((ev["pct"] >= 95).sum()),
                      "mediana_pct": round(float(ev["pct"].median()), 1) if len(ev) else None,
                      "max_pct": round(float(ev["pct"].max()), 1) if len(ev) else None})
        resumen.append((f"Control {nombre}", int((ev["pct"] >= 95).sum())))
    ts = pd.DataFrame(filas)
    di(f"\n{'subconjunto':<14}{'genes':>7}{'evaluables':>12}{'pct>=95':>9}"
       f"{'mediana':>9}{'maximo':>9}")
    for _, r in ts.iterrows():
        di(f"{r['subconjunto']:<14}{r['n_genes']:>7}{r['n_evaluables']:>12}"
           f"{r['n_pct95']:>9}{r['mediana_pct']:>9}{r['max_pct']:>9}")
    ts.to_csv(RES / "sensibilidad_controles.tsv", sep="\t", index=False)

    # ---- F: propiedades de los conjuntos ----
    di("\n" + "=" * 78)
    di("F. PROPIEDADES DE LOS CONJUNTOS DE REGIONES")
    di("=" * 78)
    bed = RAW / "regiones.bed"
    if bed.exists():
        reg = pd.read_csv(bed, sep="\t", header=None,
                          names=["chrom", "ini", "fin", "nombre"])
        reg["clase"] = reg["nombre"].str.rsplit("_", n=1).str[-1]
        reg["kb"] = (reg["fin"] - reg["ini"]) / 1000
        di(f"\n{'clase':<12}{'regiones':>10}{'Mb':>8}{'kb mediana':>12}"
           f"{'kb min':>9}{'kb max':>9}")
        for cl, g in reg.groupby("clase"):
            di(f"{cl:<12}{len(g):>10}{g['kb'].sum()/1000:>8.2f}"
               f"{g['kb'].median():>12.1f}{g['kb'].min():>9.1f}{g['kb'].max():>9.1f}")
        tot = reg.groupby("clase")["kb"].sum() / 1000
        di(f"\n{'clase':<12}{'variantes':>11}{'var/Mb':>10}{'comunes':>10}{'com/Mb':>9}")
        for cl in ["candidato", "control"]:
            tod = mc[mc["clase"] == cl]
            com = d[d["clase"] == cl]
            if cl in tot.index and tot[cl]:
                di(f"{cl:<12}{len(tod):>11,}{len(tod)/tot[cl]:>10,.0f}"
                   f"{len(com):>10,}{len(com)/tot[cl]:>9,.0f}")
        if "consecuencia" in mc.columns:
            di("")
            for cl in ["candidato", "control"]:
                sub = mc[mc["clase"] == cl]
                cod = sub["consecuencia"].fillna("").str.contains(
                    "missense|synonymous|stop|frameshift", case=False).mean()
                di(f"{cl}: {100*cod:.1f}% de variantes con consecuencia codificante")
        reg.to_csv(RES / "sensibilidad_regiones.tsv", sep="\t", index=False)

    # ---- G: detalle por locus ----
    di("\n" + "=" * 78)
    di("G. DETALLE DEL FONDO LOCAL POR LOCUS")
    di("=" * 78)
    det = percentiles(d, cand_ids, boot=True, rng=rng)
    det["etiqueta"] = det["variant_id"].map(etiq)
    di(f"\n{'locus':<22}{'estrato':<7}{'maf':>8}{'vecinos':>9}{'bloques':>9}"
       f"{'d_max':>8}{'pct':>7}{'IC95':>13}")
    for _, r in det.sort_values("pct", ascending=False, na_position="last").iterrows():
        ic = f"[{r['IC_inf']:.0f}-{r['IC_sup']:.0f}]" if "IC_inf" in r and r.get("IC_inf") == r.get("IC_inf") else "-"
        di(f"{r['etiqueta']:<22}{r['estrato']:<7}{r['maf']:>8.4f}{r['n_vec']:>9}"
           f"{r['n_blq']:>9}{r['d_max']:>8.4f}{r['pct']:>7.1f}{ic:>13}")
    det.to_csv(RES / "sensibilidad_detalle_locus.tsv", sep="\t", index=False)

    # ---- resumen ----
    di("\n" + "=" * 78)
    di("RESUMEN: LOCI EN PERCENTIL >= 95 BAJO CADA VARIANTE DEL ANALISIS")
    di("=" * 78)
    di("")
    for nombre, n in resumen:
        di(f"  {nombre:<34}{n}")
    di("\nSi todas las filas son cero, el resultado negativo del analisis")
    di("dirigido no depende de las decisiones evaluadas.")

    (RES / "sensibilidades.txt").write_text("\n".join(salida))
    print(f"\n\nResumen -> {RES / 'sensibilidades.txt'}")


if __name__ == "__main__":
    main()
