#!/usr/bin/env python3
"""
Distribucion nula de bloques atipicos con regiones control apareadas.

POR QUE
-------
La prueba binomial del barrido usa la longitud fisica de las regiones como
denominador: reparte los bloques esperados en proporcion a las megabases de
cada clase. Ese denominador supone que una megabase candidata y una megabase
control ofrecen la misma oportunidad de generar un bloque atipico, y no es
asi: la densidad de variantes comunes evaluables es alrededor de un 51% mayor
en las regiones candidatas. Con mas variantes por kilobase hay mas ocasiones
de que alguna caiga en la cola alta de FST y mas ocasiones de que dos caigan
lo bastante cerca para encadenarse en un bloque.

Este script sustituye ese denominador por uno empirico. En lugar de suponer
que la longitud basta, construye la distribucion nula del numero de bloques
que se observaria en nueve regiones control elegidas para parecerse a las
nueve candidatas en las tres cosas que determinan la oportunidad: su longitud,
su densidad de variantes comunes evaluables y la plataforma de la que
proceden esas variantes. El valor de p sale de comparar el numero observado
con esa distribucion, sin suponer nada sobre la forma de la nula.

QUE NO RESUELVE
---------------
El apareamiento iguala oportunidad, no arquitectura. Dos regiones con la misma
longitud, la misma densidad y la misma plataforma pueden diferir en tasa de
recombinacion, estructura de haplotipos y regimen selectivo, y eso tambien
afecta al encadenamiento de variantes atipicas. La nula de este script es mas
exigente que la binomial por longitud, no es una nula perfecta.

ELECCIONES DECLARADAS
---------------------
El apareamiento se hace gen a gen, no por ventanas arbitrarias, porque la
unidad del analisis original es la region genica y los percentiles de control
se calcularon dejando fuera el propio gen. Cada replica toma un gen control
distinto por cada gen candidato, sin reemplazo dentro de la replica, para que
una misma region no cuente dos veces. Los percentiles de las variantes control
son los de dejar un gen fuera, que ya trae el barrido: son los que
corresponden, porque son los que miden la cola alta contra el resto del fondo.

Las tolerancias se relajan por pasos solo cuando un gen candidato no reune el
minimo de comparadores, y la relajacion aplicada se registra por gen en la
tabla de apareamiento, para que se vea en que casos el fondo no daba de si.

Correr: python3 scripts/26_nula_apareada.py
Salidas:
  results/suplementario/Tabla_S5_nula_apareada.tsv
  results/suplementario/Tabla_S6_apareamiento_regiones.tsv
"""

import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
SCRIPTS = RAIZ / "scripts"
RES = RAIZ / "results"
PROC = RAIZ / "data" / "processed"
RAW = RAIZ / "data" / "raw"
SAL = RES / "suplementario"

N_REPLICAS = 10000
SEMILLA = 20261110
VENTANA_KB = 50
UMBRAL_PCT = 99.0
MIN_COMPARADORES = 5

# Conjuntos de genes candidatos que se evaluan. El segundo retira AGT, que
# es el peor apareado y a la vez aporta uno de los bloques observados.
EXCLUSIONES = [(), ("AGT",)]

# Tolerancias, de la mas estricta a la mas laxa. Cada entrada es la razon
# maxima de longitud, la razon maxima de densidad y la diferencia maxima en
# la fraccion de variantes de genoma completo.
TOLERANCIAS = [(1.5, 1.25, 0.10),
               (2.0, 1.50, 0.20),
               (3.0, 2.00, 0.35),
               (5.0, 3.00, 1.00)]


def cargar_24():
    """Reutiliza 24b_sensibilidad_simetrica.py, que a su vez reutiliza el 12.

    Asi el marco analitico, los percentiles y el agrupamiento son exactamente
    los del barrido original y no una reimplementacion.
    """
    for nombre in ("24b_sensibilidad_simetrica.py", "24_sensibilidad_simetrica.py"):
        f = SCRIPTS / nombre
        if f.exists():
            break
    else:
        sys.exit("Falta 24b_sensibilidad_simetrica.py en scripts/.")
    spec = importlib.util.spec_from_file_location("sens24", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    faltan = [n for n in ("cargar_12", "preparar", "percentilar", "agrupar",
                          "CASI_PRIVADO")
              if not hasattr(mod, n)]
    if faltan:
        sys.exit(f"{f.name} no expone {faltan}; revisar sus nombres.")
    return mod


def longitudes_por_gen(bed):
    """Longitud total y clase de cada region, desde el BED.

    Un gen puede aparecer en varios intervalos; se suman. El nombre del BED
    trae el simbolo y la clase, con separadores que han variado entre
    versiones, de modo que se extrae el simbolo como el primer token
    alfanumerico y la clase por la presencia de la palabra candidato.
    """
    r = pd.read_csv(bed, sep="\t", header=None,
                    names=["chrom", "ini", "fin", "nombre"])
    nom = r["nombre"].astype(str)
    r["gen"] = nom.map(lambda s: re.split(r"[|;,_\s]", s)[0].strip())
    r["clase"] = np.where(nom.str.contains("candidato", case=False),
                          "candidato", "control")
    r["largo"] = r["fin"] - r["ini"]
    g = r.groupby(["gen", "clase"], as_index=False)["largo"].sum()
    dup = g["gen"].duplicated(keep=False)
    if dup.any():
        sys.exit("Hay genes con las dos clases en el BED: "
                 f"{sorted(g.loc[dup, 'gen'].unique())}")
    return g.set_index("gen")


def perfil_de_regiones(todo, largos):
    """Longitud, densidad y composicion de plataforma de cada region.

    La densidad se mide sobre las variantes que el barrido evalua de verdad,
    no sobre todas las del VCF: son esas las que pueden llegar a ser atipicas.
    """
    n = todo.groupby("gen").size().rename("n_variantes")
    wgs = (todo.assign(e=(todo["estrato"] == "wgs").astype(float))
           .groupby("gen")["e"].mean().rename("fraccion_wgs"))
    p = pd.concat([n, wgs], axis=1).join(largos, how="inner")
    p["densidad_por_kb"] = p["n_variantes"] / (p["largo"] / 1000)
    return p


def comparadores(perf, gen, tol):
    """Genes control que se parecen al candidato dado, con una tolerancia."""
    r_largo, r_dens, d_wgs = tol
    c = perf.loc[gen]
    ctrl = perf[perf["clase"] == "control"]
    razon_l = ctrl["largo"] / c["largo"]
    razon_d = ctrl["densidad_por_kb"] / c["densidad_por_kb"]
    ok = ((razon_l.between(1 / r_largo, r_largo))
          & (razon_d.between(1 / r_dens, r_dens))
          & ((ctrl["fraccion_wgs"] - c["fraccion_wgs"]).abs() <= d_wgs))
    return list(ctrl.index[ok])


def emparejar(perf):
    """Para cada gen candidato, su grupo de comparadores y la tolerancia usada."""
    cands = list(perf.index[perf["clase"] == "candidato"])
    grupos, filas = {}, []
    for gen in cands:
        for i, tol in enumerate(TOLERANCIAS):
            pool = comparadores(perf, gen, tol)
            if len(pool) >= MIN_COMPARADORES:
                break
        grupos[gen] = pool
        c = perf.loc[gen]
        filas.append({
            "gen_candidato": gen,
            "largo_kb": round(c["largo"] / 1000, 1),
            "variantes_evaluables": int(c["n_variantes"]),
            "densidad_por_kb": round(c["densidad_por_kb"], 2),
            "fraccion_wgs": round(c["fraccion_wgs"], 3),
            "comparadores": len(pool),
            "tolerancia_largo": tol[0],
            "tolerancia_densidad": tol[1],
            "tolerancia_fraccion_wgs": tol[2],
            "relajacion": i,
            "mediana_largo_kb_comparadores":
                round(float(perf.loc[pool, "largo"].median() / 1000), 1) if pool else np.nan,
            "mediana_densidad_comparadores":
                round(float(perf.loc[pool, "densidad_por_kb"].median()), 2) if pool else np.nan,
        })
    return grupos, pd.DataFrame(filas)


def bloques_de(sub, ventana_pb, agrupar):
    """Numero de bloques atipicos en un subconjunto de variantes."""
    if not len(sub):
        return 0
    g = agrupar(sub, ventana_pb)
    return int(g.groupby(["chrom", "bloque"]).ngroups)


def nula(atipicas_por_gen, grupos, rng, ventana_pb, agrupar):
    """Distribucion nula del numero de bloques en regiones apareadas."""
    genes = list(grupos)
    cuentas = np.empty(N_REPLICAS, dtype=int)
    sin_hueco = 0
    for r in range(N_REPLICAS):
        tomados, piezas = set(), []
        for gen in genes:
            pool = [x for x in grupos[gen] if x not in tomados]
            if not pool:
                pool = grupos[gen]
                sin_hueco += 1
            elegido = pool[rng.integers(len(pool))]
            tomados.add(elegido)
            piezas.append(atipicas_por_gen.get(elegido))
        piezas = [p for p in piezas if p is not None and len(p)]
        cuentas[r] = (bloques_de(pd.concat(piezas), ventana_pb, agrupar)
                      if piezas else 0)
    return cuentas, sin_hueco


def main():
    s24 = cargar_24()
    b12 = s24.cargar_12()

    f_mc = PROC / "mcps_frecuencias.tsv"
    f_bed = RAW / "regiones.bed"
    for f in (f_mc, f_bed):
        if not f.exists():
            sys.exit(f"Falta {f}")
    mc = pd.read_csv(f_mc, sep="\t", low_memory=False)
    largos = longitudes_por_gen(f_bed)

    df, af_a, af_b = s24.preparar(b12, mc)
    rng = np.random.default_rng(SEMILLA)
    ventana_pb = VENTANA_KB * 1000

    conjuntos = {
        "completo": df,
        "restringido_asimetrico": df[(df[af_a] >= s24.CASI_PRIVADO)
                                     & (df[af_b] >= s24.CASI_PRIVADO)],
        "restringido_simetrico": df[
            (np.minimum(df[af_a], 1 - df[af_a]) >= s24.CASI_PRIVADO)
            & (np.minimum(df[af_b], 1 - df[af_b]) >= s24.CASI_PRIVADO)],
    }

    print("=" * 76)
    print("NULA APAREADA POR LONGITUD, DENSIDAD Y PLATAFORMA")
    print("=" * 76)
    print(f"\n{N_REPLICAS:,} replicas, semilla {SEMILLA}, umbral {UMBRAL_PCT}, "
          f"ventana {VENTANA_KB} kb")

    resultados, apareamiento = [], None
    for nombre, sub in conjuntos.items():
        todo = s24.percentilar(b12, sub.copy())
        if not len(todo):
            print(f"\n{nombre}: sin variantes evaluables")
            continue
        falta = set(todo["gen"]) - set(largos.index)
        if falta:
            sys.exit(f"Genes sin region en el BED: {sorted(falta)[:5]}")

        perf_todo = perfil_de_regiones(todo, largos)
        at = todo[todo["pct"] >= UMBRAL_PCT]
        atipicas_por_gen = {g: s for g, s in at.groupby("gen")}

        # El apareamiento no sale igual de bueno en los nueve genes: AGT
        # necesita relajar las tolerancias y sus comparadores tienen menor
        # densidad que el propio gen, y ademas aporta uno de los bloques
        # observados. Se repite todo sin el para ver si el resultado depende
        # de ese unico emparejamiento deficiente.
        for excluidos in EXCLUSIONES:
            perf = perf_todo.drop(index=[g for g in excluidos
                                         if g in perf_todo.index])
            grupos, tabla = emparejar(perf)
            if nombre == "completo" and not excluidos:
                apareamiento = tabla

            cand = at[(at["clase"] == "candidato")
                      & (~at["gen"].isin(excluidos))]
            obs = bloques_de(cand, ventana_pb, s24.agrupar)

            cuentas, sin_hueco = nula(atipicas_por_gen, grupos, rng,
                                      ventana_pb, s24.agrupar)
            p = (1 + int((cuentas >= obs).sum())) / (N_REPLICAS + 1)
            lo, hi = np.percentile(cuentas, [2.5, 97.5])
            pocos = int((tabla["comparadores"] < MIN_COMPARADORES).sum())

            resultados.append({
                "conjunto": nombre,
                "genes_excluidos": ",".join(excluidos) if excluidos else "ninguno",
                "umbral_percentil": UMBRAL_PCT,
                "ventana_kb": VENTANA_KB,
                "bloques_candidatos_observados": obs,
                "nula_media": round(float(cuentas.mean()), 2),
                "nula_mediana": float(np.median(cuentas)),
                "nula_IC95": f"{lo:.0f}-{hi:.0f}",
                "nula_maxima": int(cuentas.max()),
                "p_empirica": round(p, 4),
                "replicas": N_REPLICAS,
                "genes_candidatos": len(grupos),
                "genes_sin_comparadores_suficientes": pocos,
                "relajacion_maxima": int(tabla["relajacion"].max()),
                "replicas_con_reutilizacion": sin_hueco,
            })
            etq = "sin " + ",".join(excluidos) if excluidos else "los nueve"
            print(f"\n{nombre}  ({etq})")
            print(f"  observado           {obs}")
            print(f"  nula                media {cuentas.mean():.2f}, mediana "
                  f"{np.median(cuentas):.0f}, IC95 {lo:.0f}-{hi:.0f}, "
                  f"maxima {cuentas.max()}")
            print(f"  p empirica          {p:.4f}")

    if not resultados:
        sys.exit("No se obtuvo ningun resultado.")

    SAL.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(resultados).to_csv(SAL / "Tabla_S5_nula_apareada.tsv",
                                    sep="\t", index=False)
    apareamiento.to_csv(SAL / "Tabla_S6_apareamiento_regiones.tsv",
                        sep="\t", index=False)

    print("\n" + "=" * 76)
    print("CALIDAD DEL APAREAMIENTO, CONJUNTO COMPLETO")
    print("=" * 76)
    print(f"\n{'gen':8}{'largo kb':>10}{'dens/kb':>9}{'wgs':>7}"
          f"{'comp':>6}{'relaj':>7}{'kb comp':>9}{'dens comp':>11}")
    for _, r in apareamiento.iterrows():
        print(f"{r['gen_candidato']:8}{r['largo_kb']:>10}"
              f"{r['densidad_por_kb']:>9}{r['fraccion_wgs']:>7}"
              f"{r['comparadores']:>6}{r['relajacion']:>7}"
              f"{r['mediana_largo_kb_comparadores']:>9}"
              f"{r['mediana_densidad_comparadores']:>11}")
    pocos = apareamiento[apareamiento["comparadores"] < MIN_COMPARADORES]
    if len(pocos):
        print(f"\nAviso: {len(pocos)} gen(es) no reunen {MIN_COMPARADORES} "
              f"comparadores ni con la tolerancia mas laxa: "
              f"{', '.join(pocos['gen_candidato'])}")
    print(f"\nEscrito en {SAL}")


if __name__ == "__main__":
    main()
