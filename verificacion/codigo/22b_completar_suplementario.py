#!/usr/bin/env python3
"""
Completa el material suplementario en los tres puntos que quedaron pendientes
tras la ultima revision.

QUE CORRIGE
-----------
1. Etiqueta del escenario de sensibilidad. La tabla S2 traia un escenario
   llamado "fondo local, analisis principal" cuyos 24 percentiles coinciden
   exactamente con la columna pct_def_alternativa de S1. El nombre era
   incorrecto por dos razones: el fondo no era local sino el de los genes
   control, y la definicion de la frecuencia de apareamiento no era la
   principal sino la alternativa. Procedia del script 18, anterior al cambio
   de definicion. Aqui esas filas se recalculan desde cero bajo la definicion
   alternativa con el fondo de genes control, de modo que el numero de
   vecinos, bloques y la distancia alcanzada correspondan al escenario que la
   etiqueta declara.

2. Conteos alelicos en S1. Se agregan AC y AN de IMX y EUR para que el lector
   pueda reproducir cada frecuencia, verificar la correccion por tamano de
   muestra del estimador de Hudson y comprobar que las frecuencias no
   detectables corresponden a conteos estimados de cero y no a redondeo. AC y
   AN son fraccionarios porque provienen de inferencia de ancestria local:
   AC_k = suma de p_ijk * G_ij y AN_k = suma de p_ijk, de modo que AN es ya el
   conteo de cromosomas atribuidos al componente k.

3. Resultados por locus de los subconjuntos de control. El script 20 calculaba
   los percentiles de cada locus para los cinco subconjuntos pero solo
   conservaba el resumen. Aqui se guardan las filas por locus, que es lo que
   permite ver cual locus se mueve y cual deja de tener comparadores.

Requiere que ya se hayan corrido 10, 19 y 20.

Correr: python3 scripts/22b_completar_suplementario.py
Salida: sobrescribe en results/suplementario/
  Tabla_S1_resultados_por_locus.tsv
  Tabla_S2_sensibilidades.tsv
  Tabla_S2b_subconjuntos_control.tsv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
PROC = RAIZ / "data" / "processed"
RAW = RAIZ / "data" / "raw"
SAL = RES / "suplementario"

VERSION = "PLANA"
MAF_MIN, D_MAX = 0.01, 0.02
N_OBJETIVO, N_MINIMO = 300, 50
BLOQUE_KB, N_BOOT = 50_000, 2000
SEMILLA = 20261006


def fst_hudson(p1, n1, p2, n2):
    with np.errstate(invalid="ignore", divide="ignore"):
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)
        return np.where(den > 0, num / den, np.nan)


def asignar_bloques(df):
    d = df.sort_values(["chrom", "pos"]).copy()
    b, actual, pc, pp = [], 0, None, None
    for c, p in zip(d["chrom"].astype(str), d["pos"]):
        if pc != c or (p - pp) > BLOQUE_KB:
            actual += 1
        b.append(actual)
        pc, pp = c, p
    d["bloque"] = b
    return d


def preparar(mc):
    """Marco analitico con ambas definiciones de frecuencia de apareamiento."""
    af_a, an_a = f"AF_IMX_{VERSION}", f"AN_IMX_{VERSION}"
    af_b, an_b = f"AF_EUR_{VERSION}", f"AN_EUR_{VERSION}"
    d = mc[(mc["reparto_default"] == 0) & (mc["es_indel"] == 0)].copy()
    d = d.dropna(subset=[af_a, af_b, an_a, an_b])
    d = d[(d[an_a] > 0) & (d[an_b] > 0)]
    p1, p2 = d[af_a].values, d[af_b].values
    pbar = (p1 + p2) / 2
    d["maf_principal"] = np.minimum(pbar, 1 - pbar)
    d["maf_alternativa"] = (np.minimum(p1, 1 - p1) + np.minimum(p2, 1 - p2)) / 2
    d["FST"] = fst_hudson(p1, d[an_a].values, p2, d[an_b].values)
    d["estrato"] = np.where(d["source"] == "WGS", "wgs", "alta")
    return d.dropna(subset=["FST"])


def percentiles(d, cand, col_maf, genes=None, boot=False, rng=None):
    """Percentil de cada candidata con la columna de apareamiento indicada."""
    dd = d[d[col_maf] >= MAF_MIN]
    fondo = dd[dd["clase"] == "control"]
    if genes is not None:
        fondo = fondo[fondo["gen"].isin(genes)]
    fondo = asignar_bloques(fondo)
    out = []
    for _, c in cand.iterrows():
        etiqueta = f"{c['gen']} {c['rsid']}"
        f = dd[dd["variant_id"] == c["variant_id"]]
        if not len(f):
            out.append({"rsid": c["rsid"], "gen": c["gen"], "locus": etiqueta,
                        "estrato": "", "percentil": np.nan, "n_vecinos": 0,
                        "n_bloques": np.nan, "distancia_max": np.nan})
            continue
        r = f.iloc[0]
        g = fondo[fondo["estrato"] == r["estrato"]]
        dist = (g[col_maf] - r[col_maf]).abs()
        vec = g[dist <= D_MAX]
        if len(vec) > N_OBJETIVO:
            vec = vec.loc[(vec[col_maf] - r[col_maf]).abs()
                          .nsmallest(N_OBJETIVO).index]
        fila = {"rsid": c["rsid"], "gen": c["gen"], "locus": etiqueta,
                "estrato": r["estrato"], "n_vecinos": len(vec)}
        if len(vec) < N_MINIMO:
            fila.update({"percentil": np.nan, "n_bloques": np.nan,
                         "distancia_max": np.nan})
            out.append(fila)
            continue
        # dos decimales: con uno, un percentil de 95.0 no permite saber si
        # cumple la regla de >= 95, que se aplica sobre el valor sin redondear
        fila["percentil"] = round(float((vec["FST"].values < r["FST"]).mean() * 100), 2)
        fila["n_bloques"] = int(vec["bloque"].nunique())
        fila["distancia_max"] = round(float((vec[col_maf] - r[col_maf]).abs().max()), 5)
        if boot and rng is not None:
            grupos = [x["FST"].values for _, x in vec.groupby("bloque")]
            nb = len(grupos)
            ps = np.array([(np.concatenate([grupos[j] for j in
                            rng.integers(0, nb, nb)]) < r["FST"]).mean() * 100
                           for _ in range(N_BOOT)])
            lo, hi = np.percentile(ps, [2.5, 97.5])
            fila["IC95_inferior"], fila["IC95_superior"] = round(lo, 1), round(hi, 1)
        out.append(fila)
    return pd.DataFrame(out)


def main():
    f_mcps = PROC / "mcps_frecuencias.tsv"
    f_pr = RES / "principal_percentiles.tsv"
    if not f_mcps.exists() or not f_pr.exists():
        sys.exit("Faltan mcps_frecuencias.tsv o principal_percentiles.tsv; "
                 "correr primero los scripts 10 y 19.")
    if not SAL.exists():
        sys.exit(f"No existe {SAL}; correr primero el script 21.")

    mc = pd.read_csv(f_mcps, sep="\t", low_memory=False)
    pr = pd.read_csv(f_pr, sep="\t")
    d = preparar(mc)
    cand = pr[["gen", "rsid", "variant_id"]].copy()

    print("=" * 76)
    print("COMPLETAR MATERIAL SUPLEMENTARIO")
    print("=" * 76)

    # ---------------- 1. S1 con conteos alelicos ----------------
    cols_cuenta = [f"AC_IMX_{VERSION}", f"AN_IMX_{VERSION}",
                   f"AC_EUR_{VERSION}", f"AN_EUR_{VERSION}"]
    faltan = [c for c in cols_cuenta if c not in mc.columns]
    if faltan:
        sys.exit(f"mcps_frecuencias.tsv no trae {faltan}; "
                 "volver a correr el script 10.")

    cuentas = mc[["variant_id"] + cols_cuenta].drop_duplicates("variant_id")
    cuentas = cuentas.rename(columns={
        f"AC_IMX_{VERSION}": "AC_IMX", f"AN_IMX_{VERSION}": "AN_IMX",
        f"AC_EUR_{VERSION}": "AC_EUR", f"AN_EUR_{VERSION}": "AN_EUR"})

    f_s1 = SAL / "Tabla_S1_resultados_por_locus.tsv"
    s1 = pd.read_csv(f_s1, sep="\t")
    for c in ("AC_IMX", "AN_IMX", "AC_EUR", "AN_EUR"):
        s1 = s1.drop(columns=[c], errors="ignore")
    s1 = s1.merge(cuentas, on="variant_id", how="left")

    orden = list(s1.columns)
    for c in ("AC_IMX", "AN_IMX", "AC_EUR", "AN_EUR"):
        orden.remove(c)
    i = orden.index("AF_EUR") + 1
    orden = orden[:i] + ["AC_IMX", "AN_IMX", "AC_EUR", "AN_EUR"] + orden[i:]
    s1 = s1[orden]
    s1["AC_IMX"] = s1["AC_IMX"].round(2)
    s1["AC_EUR"] = s1["AC_EUR"].round(2)
    s1["AN_IMX"] = s1["AN_IMX"].round(1)
    s1["AN_EUR"] = s1["AN_EUR"].round(1)

    # Las frecuencias venian redondeadas a cuatro decimales, de modo que
    # HFE rs1800562 se imprimia como 0.0 pese a tener un conteo estimado de
    # 0.2 cromosomas. Se recalculan desde AC y AN, que es lo que permite
    # verificarlas, y las muy pequenas van en notacion cientifica para no
    # perder la magnitud.
    def fmt_af(x):
        if pd.isna(x):
            return ""
        if x == 0:
            return "0"
        if x < 1e-4 or x > 1 - 1e-4:
            return f"{x:.3e}"
        return f"{x:.6f}".rstrip("0").rstrip(".")

    for comp in ("IMX", "EUR"):
        ac, an, af = f"AC_{comp}", f"AN_{comp}", f"AF_{comp}"
        crudo = pd.to_numeric(s1[ac], errors="coerce") / pd.to_numeric(
            s1[an], errors="coerce").replace(0, np.nan)
        s1[af] = crudo.map(fmt_af)
    chicas = s1[pd.to_numeric(s1["AF_IMX"], errors="coerce") < 1e-4]
    if len(chicas):
        print("\n  frecuencias IMX por debajo de 1e-4, recalculadas de AC/AN:")
        for _, q in chicas.iterrows():
            print(f"    {q['gen']} {q['rsid']:<12} AF={q['AF_IMX']:<10} "
                  f"AC={q['AC_IMX']} / AN={q['AN_IMX']}")
    s1.to_csv(f_s1, sep="\t", index=False)

    sin = s1[s1["AC_IMX"].fillna(1) == 0]
    print(f"\nS1: {len(s1)} loci con AC y AN de IMX y EUR")
    print(f"  loci con AC_IMX estimado en cero: {len(sin)}"
          + (f" ({', '.join(sin['gen'] + ' ' + sin['rsid'])})" if len(sin) else ""))
    print(f"  AN_IMX mediano {s1['AN_IMX'].median():,.0f} cromosomas, "
          f"AN_EUR mediano {s1['AN_EUR'].median():,.0f}")

    # ---------------- 2. escenario de definicion alternativa ----------------
    rng = np.random.default_rng(SEMILLA)
    alt = percentiles(d, cand, "maf_alternativa", boot=True, rng=rng)
    alt["escenario"] = "definicion alternativa de la frecuencia de apareamiento"
    alt["condicion"] = alt["estrato"]

    f_s2 = SAL / "Tabla_S2_sensibilidades.tsv"
    s2 = pd.read_csv(f_s2, sep="\t")
    # El script tiene que poder correrse mas de una vez sobre el mismo archivo
    # sin acumular filas. Se retiran tanto el escenario mal etiquetado como los
    # que este script genera, porque si solo se quitara el primero una segunda
    # ejecucion duplicaria los nuevos.
    A_RETIRAR = ("fondo local, analisis principal",
                 "definicion alternativa de la frecuencia de apareamiento",
                 "subconjunto de genes control")
    quitadas = {e: int((s2["escenario"] == e).sum()) for e in A_RETIRAR}
    s2 = s2[~s2["escenario"].isin(A_RETIRAR)]
    print("\nS2: filas retiradas antes de reescribir")
    for e, k in quitadas.items():
        if k:
            print(f"  {e:<52}{k:>4}")
    if not any(quitadas.values()):
        print("  ninguna: el archivo no traia estos escenarios")
    print(f"  se escriben {len(alt)} filas de la definicion alternativa")

    # ---------------- 3. subconjuntos de control por locus ----------------
    genes = sorted(d[d["clase"] == "control"]["gen"].dropna().unique())
    rng_s = np.random.default_rng(SEMILLA)
    sub = {"completo": genes}
    for i in range(3):
        sub[f"mitad_{i+1}"] = sorted(rng_s.choice(genes, len(genes) // 2,
                                                  replace=False))
    f74 = RAW / "genes_control_v1_74.tsv"
    if f74.exists():
        g74 = pd.read_csv(f74, sep="\t")["gen"].tolist()
        sub["v1_74genes"] = [g for g in g74 if g in genes]
    else:
        print(f"  aviso: no se encontro {f74.name}; el subconjunto de 74 genes "
              "no se incluye")

    bloques_sub, resumen, pct_completo = [], [], None
    for nombre, gs in sub.items():
        t = percentiles(d, cand, "maf_principal", genes=gs)
        if nombre == "completo":
            # El subconjunto "completo" es el analisis principal: los mismos
            # 219 genes, la misma frecuencia de apareamiento y los mismos
            # parametros. Se guarda para llevar el percentil de S1 a dos
            # decimales, de modo que S1, S2b y el texto den la misma mediana y
            # el caso limite de la regla del umbral sea verificable.
            pct_completo = {r: v for r, v in zip(t["rsid"], t["percentil"])
                            if pd.notna(v)}
        t["escenario"] = "subconjunto de genes control"
        t["condicion"] = nombre
        bloques_sub.append(t)
        ev = t.dropna(subset=["percentil"])
        resumen.append({"subconjunto": nombre, "n_genes": len(gs),
                        "loci_evaluables": len(ev),
                        "loci_percentil_95_o_mas": int((ev["percentil"] >= 95).sum()),
                        "mediana": round(float(ev["percentil"].median()), 2),
                        "maximo": round(float(ev["percentil"].max()), 2),
                        "AGT_rs2493137": float(
                            t.loc[t["rsid"] == "rs2493137", "percentil"].iloc[0])
                        if (t["rsid"] == "rs2493137").any() else None})
        lim = ev[(ev["percentil"] >= 94.5) & (ev["percentil"] < 95.5)]
        print(f"  control {nombre:<12} {len(gs):>4} genes  "
              f"{len(ev):>3} evaluables  mediana {ev['percentil'].median():.2f}  "
              f">=95: {int((ev['percentil'] >= 95).sum())}")
        for _, q in lim.iterrows():
            print(f"      en el limite: {q['locus']} = {q['percentil']:.2f} "
                  f"({'cuenta' if q['percentil'] >= 95 else 'no cuenta'})")

    pd.DataFrame(resumen).to_csv(SAL / "Tabla_S2b_subconjuntos_control.tsv",
                                 sep="\t", index=False)

    # ---------------- S1: percentil con dos decimales ----------------
    if pct_completo:
        s1 = pd.read_csv(f_s1, sep="\t")
        antes = s1["percentil"].copy()
        s1["percentil"] = s1["rsid"].map(pct_completo).combine_first(
            s1["percentil"])
        dif = (s1["percentil"] - antes).abs()
        s1.to_csv(f_s1, sep="\t", index=False)
        ev = s1[s1["evaluable"] == "si"]["percentil"].dropna()
        print(f"\nS1: percentil a dos decimales")
        print(f"  cambio maximo respecto del valor a un decimal: {dif.max():.2f}")
        print(f"  mediana de los {len(ev)} loci evaluables: {ev.median():.2f}")
        print(f"  rango: {ev.min():.2f} a {ev.max():.2f}")
        print(f"  loci con percentil >= 95: {int((ev >= 95).sum())}")
        print("  el intervalo de remuestreo conserva el decimal del analisis "
              "principal")

    # ---------------- ensamblar S2 ----------------
    cols = ["escenario", "condicion", "rsid", "locus", "percentil",
            "IC95_inferior", "IC95_superior", "n_vecinos", "n_bloques",
            "distancia_max"]
    nuevo = pd.concat([alt] + bloques_sub, ignore_index=True)
    nuevo = nuevo.reindex(columns=cols)
    s2 = s2.reindex(columns=cols)

    # completar rsid donde venia vacio, a partir de la etiqueta del locus
    def rs_de(locus):
        for parte in str(locus).split():
            if parte.startswith("rs"):
                return parte
        return ""

    for marco in (s2, nuevo):
        falta = marco["rsid"].isna() | marco["rsid"].astype(str).isin(("", "-", "nan"))
        marco.loc[falta, "rsid"] = marco.loc[falta, "locus"].map(rs_de)

    tS2 = pd.concat([nuevo, s2], ignore_index=True)
    orden_esc = {"definicion alternativa de la frecuencia de apareamiento": 0,
                 "version de frecuencia": 1, "plataformas separadas": 2,
                 "subconjunto de genes control": 3}
    tS2["_o"] = tS2["escenario"].map(orden_esc).fillna(9)
    tS2 = (tS2.sort_values(["_o", "condicion", "percentil"],
                           ascending=[True, True, False], na_position="last")
              .drop(columns="_o"))
    tS2.to_csv(f_s2, sep="\t", index=False)

    print(f"\nS2: {len(tS2)} filas en {tS2['escenario'].nunique()} escenarios")
    for e, g in tS2.groupby("escenario"):
        print(f"  {e:<52}{len(g):>4} filas")
    print(f"\nEscrito en {SAL}")


if __name__ == "__main__":
    main()
