#!/usr/bin/env python3
"""
Sensibilidad del barrido al criterio de exclusion de alelos poco frecuentes, y
extension fisica de los bloques.

POR QUE IMPORTA EL CRITERIO
---------------------------
El analisis restringido excluye variantes cuya frecuencia del alelo ALTERNO
queda por debajo de 0.02 en alguno de los dos componentes. Ese criterio es
asimetrico: una variante con frecuencia alterna de 0.99 sobrevive al filtro
aunque el alelo de referencia sea igual de escaso. El criterio simetrico,
min(p, 1-p) < 0.02, trata los dos alelos por igual. Hasta ahora solo se
cuantificaba cuantas variantes mas excluiria el simetrico; aqui se corre el
barrido completo bajo los tres conjuntos, que es lo que permite decir si
cambian los bloques, las razones de tasas y la interpretacion.

POR QUE ESTE SCRIPT IMPORTA EL 12 EN LUGAR DE REIMPLEMENTARLO
-------------------------------------------------------------
Una version anterior de este script reimplementaba el barrido y no reproducia
los resultados del 12: definia atipica a la variante que supera un percentil
global de FST dentro de su estrato, cuando el 12 le asigna a cada variante un
percentil contra sus vecinas en frecuencia dentro del fondo de control, y usaba
la frecuencia de apareamiento principal cuando el 12 usa la alternativa. Con
esas dos diferencias el numero de bloques de control cambiaba de 51 a 18 y
aparecian comparaciones significativas donde el analisis original no tiene
ninguna.

Para que eso no vuelva a pasar, este script carga 12_barrido.py como modulo y
reutiliza sus funciones y sus constantes. La unica diferencia entre los tres
conjuntos es el criterio de exclusion. Al final se comprueba que, en el umbral
y la ventana principales del 12, el conjunto completo reproduce exactamente los
bloques candidatos que el 12 dejo en barrido_atipicos.tsv; si no coincide, el
script lo dice en lugar de callarlo.

NOTA SOBRE LOS BLOQUES
----------------------
Los bloques se construyen agrupando variantes atipicas consecutivas mientras la
distancia con la anterior no supere la ventana, de modo que el parametro es la
distancia maxima entre variantes vecinas y no el tamano del bloque: la
extension de un bloque no esta acotada por la ventana. Se reporta la
distribucion real.

Correr: python3 scripts/24_sensibilidad_simetrica.py
Salidas:
  results/suplementario/Tabla_S3_rejilla_barrido.tsv
  results/suplementario/Tabla_S4_extension_bloques.tsv
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, chi2

RAIZ = Path(__file__).resolve().parent.parent
SCRIPTS = RAIZ / "scripts"
RES = RAIZ / "results"
PROC = RAIZ / "data" / "processed"
RAW = RAIZ / "data" / "raw"
SAL = RES / "suplementario"

UMBRALES = (95.0, 97.5, 99.0)
VENTANAS_KB = (25, 50, 100)
CASI_PRIVADO = 0.02


def cargar_12():
    """Carga 12_barrido.py como modulo.

    El nombre empieza con digito, de modo que no se puede importar por la via
    normal. main() esta protegido por __main__, asi que cargarlo no ejecuta el
    barrido.
    """
    f = SCRIPTS / "12_barrido.py"
    if not f.exists():
        sys.exit(f"Falta {f}: este script reutiliza sus funciones.")
    spec = importlib.util.spec_from_file_location("barrido12", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    faltan = [n for n in ("fst_hudson", "percentiles_contra_fondo",
                          "N_VECINOS", "VENTANA_LD", "PCT_UMBRAL", "MAF_MIN",
                          "VERSION", "A", "B")
              if not hasattr(mod, n)]
    if faltan:
        sys.exit(f"12_barrido.py no expone {faltan}; revisar sus nombres.")
    return mod


def agrupar(df_at, ventana_pb):
    """Agrupamiento por distancia maxima entre variantes vecinas.

    Misma logica que agrupar_bloques del script 12, con la ventana como
    argumento para poder recorrer la rejilla. Devuelve el marco con la columna
    de bloque, sin colapsar, para poder contar y medir extension.
    """
    if not len(df_at):
        return df_at.assign(bloque=pd.Series(dtype=int))
    d = df_at.sort_values(["chrom", "pos"]).copy()
    b, actual, pc, pp = [], 0, None, None
    for c, p in zip(d["chrom"].astype(str), d["pos"]):
        if pc != c or (p - pp) > ventana_pb:
            actual += 1
        b.append(actual)
        pc, pp = c, p
    d["bloque"] = b
    return d


def ic_poisson(k, expo, alpha=0.05):
    """Intervalo exacto de Poisson para una tasa, por su relacion con chi2."""
    lo = chi2.ppf(alpha / 2, 2 * k) / 2 if k > 0 else 0.0
    hi = chi2.ppf(1 - alpha / 2, 2 * k + 2) / 2
    return lo / expo, hi / expo


def mb_por_clase(bed):
    """Megabases evaluadas por clase, desde el BED de regiones."""
    r = pd.read_csv(bed, sep="\t", header=None,
                    names=["chrom", "ini", "fin", "nombre"])
    r["clase"] = np.where(r["nombre"].astype(str).str.contains("candidato"),
                          "candidato", "control")
    r["mb"] = (r["fin"] - r["ini"]) / 1e6
    return r.groupby("clase")["mb"].sum().to_dict()


def preparar(b12, mc):
    """Marco analitico del barrido, igual que en main() del script 12."""
    af_a, an_a = f"AF_{b12.A}_{b12.VERSION}", f"AN_{b12.A}_{b12.VERSION}"
    af_b, an_b = f"AF_{b12.B}_{b12.VERSION}", f"AN_{b12.B}_{b12.VERSION}"
    df = mc[(mc["reparto_default"] == 0) & (mc["es_indel"] == 0)]
    df = df.dropna(subset=[af_a, af_b, an_a, an_b])
    df = df[(df[an_a] > 0) & (df[an_b] > 0)].copy()
    # frecuencia de apareamiento del barrido: promedio de las frecuencias
    # menores de cada componente, que es la definicion alternativa
    df["maf_media"] = ((np.minimum(df[af_a], 1 - df[af_a])
                        + np.minimum(df[af_b], 1 - df[af_b])) / 2)
    df = df[df["maf_media"] >= b12.MAF_MIN]
    df["estrato"] = np.where(df["source"] == "WGS", "wgs", "alta")
    df["FST"] = b12.fst_hudson(df[af_a].values, df[an_a].values,
                               df[af_b].values, df[an_b].values)
    return df.dropna(subset=["FST"]).reset_index(drop=True), af_a, af_b


def percentilar(b12, df):
    """Percentil de cada variante, igual que en el script 12.

    Las candidatas se evaluan contra el fondo de control completo; las de
    control, contra el control de los demas genes, que es la tasa de fondo.
    """
    partes = []
    for est in ("alta", "wgs"):
        sub = df[df["estrato"] == est]
        ctrl = sub[sub["clase"] == "control"]
        cand = sub[sub["clase"] == "candidato"].copy()
        if len(ctrl) < b12.N_VECINOS or not len(cand):
            continue
        cand["pct"] = b12.percentiles_contra_fondo(cand, ctrl)
        cand["grupo"] = "candidato"
        ctrl_ev = ctrl.copy()
        pcts = []
        for gen, g in ctrl.groupby("gen"):
            otros = ctrl[ctrl["gen"] != gen]
            pcts.append(pd.Series(b12.percentiles_contra_fondo(g, otros),
                                  index=g.index))
        ctrl_ev["pct"] = pd.concat(pcts).reindex(ctrl_ev.index)
        ctrl_ev["grupo"] = "control"
        partes.append(pd.concat([cand, ctrl_ev]))
    return pd.concat(partes, ignore_index=True) if partes else pd.DataFrame()


def comparar(todo, mb, umbral, ventana_kb):
    """Bloques atipicos por grupo, tasas, intervalos y prueba binomial."""
    filas = []
    mb_c, mb_k = mb.get("candidato", np.nan), mb.get("control", np.nan)
    prop = mb_c / (mb_c + mb_k)
    for est in ("alta", "wgs"):
        s = todo[todo["estrato"] == est]
        if not len(s):
            continue
        at = s[s["pct"] >= umbral]
        bc = agrupar(at[at["grupo"] == "candidato"],
                     ventana_kb * 1000)["bloque"].nunique()
        bk = agrupar(at[at["grupo"] == "control"],
                     ventana_kb * 1000)["bloque"].nunique()
        tc, tk = bc / mb_c, bk / mb_k
        lo_c, hi_c = ic_poisson(bc, mb_c)
        n = bc + bk
        pval = (binomtest(bc, n, prop, alternative="greater").pvalue
                if n else np.nan)
        filas.append({
            "estrato": est, "umbral_percentil": umbral,
            "ventana_kb": ventana_kb,
            "n_var_cand": int((s["grupo"] == "candidato").sum()),
            "n_var_ctrl": int((s["grupo"] == "control").sum()),
            "var_atipicas_cand": int((at["grupo"] == "candidato").sum()),
            "var_atipicas_ctrl": int((at["grupo"] == "control").sum()),
            "bloques_cand": bc, "bloques_ctrl": bk,
            "tasa_cand_Mb": round(tc, 3),
            "IC_cand_Mb": f"{lo_c:.2f}-{hi_c:.2f}",
            "tasa_ctrl_Mb": round(tk, 3),
            "razon": round(tc / tk, 3) if tk else np.nan,
            "p_binomial_unilateral": pval})
    return filas


def main():
    b12 = cargar_12()
    f_mcps = PROC / "mcps_frecuencias.tsv"
    bed = RAW / "regiones.bed"
    for f in (f_mcps, bed):
        if not f.exists():
            sys.exit(f"Falta {f}; correr primero los scripts 08 y 10.")
    SAL.mkdir(parents=True, exist_ok=True)

    print("=" * 76)
    print("SENSIBILIDAD DEL BARRIDO AL CRITERIO DE EXCLUSION")
    print("=" * 76)
    print(f"\nParametros tomados de 12_barrido.py:")
    print(f"  MAF_MIN      {b12.MAF_MIN}")
    print(f"  N_VECINOS    {b12.N_VECINOS}")
    print(f"  PCT_UMBRAL   {b12.PCT_UMBRAL}")
    print(f"  VENTANA_LD   {b12.VENTANA_LD} pb")
    print(f"  VERSION      {b12.VERSION}, {b12.A} frente a {b12.B}")

    mc = pd.read_csv(f_mcps, sep="\t", low_memory=False)
    df, af_a, af_b = preparar(b12, mc)

    p1, p2 = df[af_a].values, df[af_b].values
    df["pasa_asimetrico"] = (p1 >= CASI_PRIVADO) & (p2 >= CASI_PRIVADO)
    df["pasa_simetrico"] = ((np.minimum(p1, 1 - p1) >= CASI_PRIVADO)
                            & (np.minimum(p2, 1 - p2) >= CASI_PRIVADO))

    n_tot = len(df)
    n_a = int(df["pasa_asimetrico"].sum())
    n_s = int(df["pasa_simetrico"].sum())
    print(f"\nVariantes comunes evaluadas            {n_tot:>8,}")
    print(f"  conserva el criterio asimetrico      {n_a:>8,}  "
          f"({100 * n_a / n_tot:.1f}%)")
    print(f"  conserva el criterio simetrico       {n_s:>8,}  "
          f"({100 * n_s / n_tot:.1f}%)")
    print(f"  excluidas de mas por el simetrico    {n_a - n_s:>8,}  "
          f"({100 * (n_a - n_s) / n_tot:.1f}% del total)")

    mb = mb_por_clase(bed)
    print(f"\nMegabases evaluadas: candidatos {mb.get('candidato', 0):.2f}, "
          f"control {mb.get('control', 0):.2f}")

    conjuntos = {"completo": df,
                 "restringido_asimetrico": df[df["pasa_asimetrico"]],
                 "restringido_simetrico": df[df["pasa_simetrico"]]}

    filas, extension, percentilados = [], [], {}
    for nombre, sub in conjuntos.items():
        todo = percentilar(b12, sub)
        if not len(todo):
            print(f"\n{nombre}: ningun estrato evaluable")
            continue
        percentilados[nombre] = todo
        for u in UMBRALES:
            for v in VENTANAS_KB:
                for f in comparar(todo, mb, u, v):
                    f["conjunto"] = nombre
                    filas.append(f)
        # Extension de los bloques de variantes atipicas al umbral principal.
        # Se recorre estrato por estrato, igual que comparar(), porque los
        # bloques se cuentan dentro de cada estrato: mezclarlos encadena
        # variantes de plataformas distintas y el numero de bloques deja de
        # coincidir con el de la rejilla.
        at = todo[todo["pct"] >= b12.PCT_UMBRAL]
        for est in ("alta", "wgs"):
            at_e = at[at["estrato"] == est]
            for v in VENTANAS_KB:
                for grp in ("candidato", "control"):
                    g = agrupar(at_e[at_e["grupo"] == grp], v * 1000)
                    if not len(g):
                        continue
                    tam = g.groupby(["chrom", "bloque"])["pos"].agg(
                        ["min", "max", "count"])
                    largo = (tam["max"] - tam["min"]) / 1000
                    extension.append({
                        "conjunto": nombre, "estrato": est, "grupo": grp,
                        "ventana_kb": v,
                        "umbral_percentil": b12.PCT_UMBRAL,
                        "n_bloques": len(tam),
                        "extension_kb_mediana": round(float(largo.median()), 1),
                        "extension_kb_p90": round(float(largo.quantile(0.9)), 1),
                        "extension_kb_maxima": round(float(largo.max()), 1),
                        "variantes_por_bloque_mediana":
                            float(tam["count"].median()),
                        "variantes_por_bloque_maxima":
                            int(tam["count"].max())})

    cols = ["conjunto", "umbral_percentil", "ventana_kb", "estrato",
            "n_var_cand", "n_var_ctrl", "var_atipicas_cand",
            "var_atipicas_ctrl", "bloques_cand", "bloques_ctrl",
            "tasa_cand_Mb", "IC_cand_Mb", "tasa_ctrl_Mb", "razon",
            "p_binomial_unilateral"]
    tb = (pd.DataFrame(filas)[cols]
          .sort_values(["conjunto", "estrato", "umbral_percentil",
                        "ventana_kb"]))
    tb.to_csv(SAL / "Tabla_S3_rejilla_barrido.tsv", sep="\t", index=False)

    ventana_principal = int(b12.VENTANA_LD / 1000)
    print(f"\n{len(tb)} comparaciones. Umbral {b12.PCT_UMBRAL:.0f}, ventana "
          f"{ventana_principal} kb:")
    pr = tb[(tb["umbral_percentil"] == b12.PCT_UMBRAL)
            & (tb["ventana_kb"] == ventana_principal)]
    print(f"\n  {'conjunto':<24}{'est':<5}{'cand':>6}{'ctrl':>6}{'razon':>8}"
          f"{'p':>8}")
    for _, r in pr.iterrows():
        raz = "-" if pd.isna(r["razon"]) else f"{r['razon']:.2f}"
        print(f"  {r['conjunto']:<24}{r['estrato']:<5}{r['bloques_cand']:>6}"
              f"{r['bloques_ctrl']:>6}{raz:>8}"
              f"{r['p_binomial_unilateral']:>8.3f}")

    sig = tb[tb["p_binomial_unilateral"] < 0.05]
    print(f"\nComparaciones con p < 0.05: {len(sig)} de {len(tb)}")
    print(f"Rango de razones: {tb['razon'].min():.2f} a {tb['razon'].max():.2f}")
    print(f"Rango de p: {tb['p_binomial_unilateral'].min():.3f} a "
          f"{tb['p_binomial_unilateral'].max():.3f}")

    # ---- comprobacion contra la salida del script 12 ----
    f_at = RES / "barrido_atipicos.tsv"
    print("\n" + "=" * 76)
    print("COMPROBACION CONTRA EL SCRIPT 12")
    print("=" * 76)
    if not f_at.exists():
        print(f"\nNo esta {f_at.name}; no se puede comprobar. Correr el 12.")
    else:
        orig = pd.read_csv(f_at, sep="\t")
        print(f"\nbarrido_atipicos.tsv trae {len(orig)} bloques candidatos")
        ok = True
        for est, g in orig.groupby("estrato"):
            aqui = pr[(pr["estrato"] == est)
                      & (pr["conjunto"] == "completo")]["bloques_cand"]
            v = int(aqui.iloc[0]) if len(aqui) else 0
            marca = "coincide" if v == len(g) else "NO COINCIDE"
            if v != len(g):
                ok = False
            print(f"  estrato {est:<5} script 12: {len(g):>3}   "
                  f"aqui: {v:>3}   {marca}")
        print("\nReproduce el barrido original." if ok else
              "\nNO reproduce el barrido original: no usar estas cifras en el "
              "manuscrito hasta resolver la diferencia.")

    te = pd.DataFrame(extension)[
        ["conjunto", "estrato", "grupo", "ventana_kb", "umbral_percentil",
         "n_bloques", "extension_kb_mediana", "extension_kb_p90",
         "extension_kb_maxima", "variantes_por_bloque_mediana",
         "variantes_por_bloque_maxima"]].sort_values(
        ["conjunto", "estrato", "grupo", "ventana_kb"])
    te.to_csv(SAL / "Tabla_S4_extension_bloques.tsv", sep="\t", index=False)

    # ---- las dos tablas deben dar el mismo numero de bloques ----
    # Se cuentan en dos lugares del script, con el mismo criterio pero
    # recorriendo la rejilla de distinta manera. Si alguna vez dejan de
    # coincidir, el paquete queda con dos conteos distintos del mismo objeto.
    desajustes = []
    for _, r in te.iterrows():
        col = "bloques_cand" if r["grupo"] == "candidato" else "bloques_ctrl"
        m = tb[(tb["conjunto"] == r["conjunto"])
               & (tb["estrato"] == r["estrato"])
               & (tb["ventana_kb"] == r["ventana_kb"])
               & (tb["umbral_percentil"] == r["umbral_percentil"])]
        if not len(m):
            desajustes.append((r["conjunto"], r["estrato"], r["grupo"],
                               r["ventana_kb"], "sin fila en la rejilla", ""))
        elif int(m.iloc[0][col]) != int(r["n_bloques"]):
            desajustes.append((r["conjunto"], r["estrato"], r["grupo"],
                               r["ventana_kb"], int(m.iloc[0][col]),
                               int(r["n_bloques"])))
    print("\n" + "=" * 76)
    print("COHERENCIA ENTRE LA REJILLA Y LA EXTENSION")
    print("=" * 76)
    if desajustes:
        print(f"\n{len(desajustes)} desajuste(s). No usar estas tablas:")
        for d in desajustes:
            print(f"  {d[0]:<24}{d[1]:<5}{d[2]:<11}{d[3]:>4} kb   "
                  f"rejilla {d[4]}   extension {d[5]}")
    else:
        print(f"\nLos {len(te)} conteos de bloques coinciden con la rejilla.")
    print("\n" + "=" * 76)
    print(f"EXTENSION FISICA DE LOS BLOQUES, ventana {ventana_principal} kb")
    print("=" * 76)
    print(f"\n{'conjunto':<24}{'est':<5}{'grupo':<11}{'blq':>5}{'mediana':>10}"
          f"{'p90':>9}{'maxima':>10}")
    for _, r in te[te["ventana_kb"] == ventana_principal].iterrows():
        print(f"{r['conjunto']:<24}{r['estrato']:<5}"
              f"{r['grupo']:<11}{r['n_bloques']:>5}"
              f"{r['extension_kb_mediana']:>9.1f}k"
              f"{r['extension_kb_p90']:>8.1f}k"
              f"{r['extension_kb_maxima']:>9.1f}k")
    print("\nLa ventana acota la distancia entre variantes vecinas, no la")
    print("extension del bloque, que puede ser varias veces mayor.")

    print(f"\nEscrito en {SAL}")


if __name__ == "__main__":
    main()
