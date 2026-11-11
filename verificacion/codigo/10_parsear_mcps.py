#!/usr/bin/env python3
"""
Convierte el VCF extraido de MCPS en una tabla plana lista para analizar.

Decisiones, todas tomadas a partir de los Metodos de Ziyatdinov et al. 2023
(Nature 622:784-793, doi 10.1038/s41586-023-06595-3):

1. CAMPOS PLANOS, no _UNREL ni _IBD.
   Los autores estiman las frecuencias con todas las muestras, porque
   restringir a no emparentados reduce en 8.6% el numero de variantes
   polimorficas, y reportan poca diferencia entre las tres versiones. La
   version plana es la principal; _UNREL e _IBD quedan como sensibilidad.
   Este script guarda las tres para poder repetir el analisis con cualquiera.

2. AN ESTA EN CROMOSOMAS.
   Segun los Metodos, AC_k = suma de p_ijk * G_ij, AN_k = suma de p_ijk, y el
   tamano de muestra efectivo es n_k = AN_k / 2. Es decir AN ya es el conteo
   efectivo de cromosomas, que es la escala que pide el estimador de Hudson.

3. SE MARCAN LOS ALELOS REPARTIDOS POR DEFECTO.
   El 46% de los singletons de exoma cayeron en tramos de ancestria
   heterocigota; a esos los autores les dieron peso igual a dos ancestrias.
   Eso produce AC exactamente 0.5 en dos ancestrias. No es informacion: es un
   reparto ante la imposibilidad de asignar. Se marcan con reparto_default y
   se excluyen del analisis.

4. SOURCE ES VARIABLE DE ESTRATIFICACION OBLIGATORIA.
   Tamanos de muestra efectivos medianos para ancestria indigena mexicana:
   91,856 en variantes de exoma frente a 6,549 en variantes de genoma. Un
   fondo que mezcle ambas no es homogeneo. Ademas, cuando una variante
   aparecia en arreglo y en secuenciacion, se usaron los datos del arreglo.

5. AF SE RECALCULA COMO AC/AN.
   El AF publicado se trunca a cero cuando AC < 1.

Correr: python3 scripts/10_parsear_mcps.py
Salida: data/processed/mcps_frecuencias.tsv
"""

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
ENTRADA = RAIZ / "data" / "raw" / "mcps" / "mcps_regiones.vcf"
REGIONES = RAIZ / "data" / "raw" / "regiones.bed"
SALIDA = RAIZ / "data" / "processed" / "mcps_frecuencias.tsv"

ANCESTRIAS = ["IMX", "EUR", "AFR"]
VERSIONES = {"": "plana", "_UNREL": "no emparentados", "_IBD": "corregida IBD"}

TOL = 1e-9


def cargar_regiones():
    regiones = []
    if not REGIONES.exists():
        print(f"aviso: falta {REGIONES.name}; las variantes saldran sin gen",
              file=sys.stderr)
        return regiones
    with open(REGIONES, encoding="utf-8") as fh:
        for linea in fh:
            if not linea.strip():
                continue
            c, i, f, nombre = linea.split("\t")[:4]
            gen, _, clase = nombre.strip().rpartition("_")
            regiones.append((c.replace("chr", ""), int(i), int(f), gen, clase))
    return regiones


def gen_de(chrom, pos, regiones):
    for c, i, f, gen, clase in regiones:
        if c == chrom and i < pos <= f:
            return gen, clase
    return "", ""


def info_a_dict(campo):
    d = {}
    for par in campo.split(";"):
        if "=" in par:
            k, v = par.split("=", 1)
            d[k] = v
    return d


def num(info, clave):
    try:
        return float(info[clave])
    except (KeyError, ValueError):
        return float("nan")


def main():
    if not ENTRADA.exists():
        sys.exit(f"Falta {ENTRADA}. Correr antes: bash scripts/09_extraer_mcps.sh sin")

    regiones = cargar_regiones()

    cols = ["variant_id", "chrom", "pos", "ref", "alt", "gen", "clase", "source",
            "es_indel", "reparto_default", "n_efectivo_IMX"]
    for suf in VERSIONES:
        etq = suf.lstrip("_") or "PLANA"
        for a in ANCESTRIAS:
            cols += [f"AC_{a}_{etq}", f"AN_{a}_{etq}", f"AF_{a}_{etq}"]

    filas = []
    n_lineas = 0

    with open(ENTRADA, encoding="utf-8") as fh:
        for linea in fh:
            if linea.startswith("#") or not linea.strip():
                continue
            n_lineas += 1
            campos = linea.rstrip("\n").split("\t")
            if len(campos) < 8:
                continue

            chrom, pos, ref, alt = campos[0], campos[1], campos[3], campos[4]
            info = info_a_dict(campos[7])
            gen, clase = gen_de(chrom, int(pos), regiones)

            fila = {
                "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
                "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                "gen": gen, "clase": clase,
                "source": info.get("SOURCE", ""),
                "es_indel": int(len(ref) != 1 or len(alt) != 1),
            }

            # Reparto por defecto: el singleton no asignable que los autores
            # dividen en partes iguales entre dos ancestrias (AC = 0.5 cada una).
            ac_planos = [num(info, f"AC_{a}") for a in ANCESTRIAS]
            reparto = int(sum(1 for v in ac_planos if abs(v - 0.5) < TOL) >= 2)
            fila["reparto_default"] = reparto

            for suf in VERSIONES:
                etq = suf.lstrip("_") or "PLANA"
                for a in ANCESTRIAS:
                    ac = num(info, f"AC_{a}{suf}")
                    an = num(info, f"AN_{a}{suf}")
                    af = ac / an if an == an and an > 0 and ac == ac else float("nan")
                    fila[f"AC_{a}_{etq}"] = round(ac, 4) if ac == ac else ""
                    fila[f"AN_{a}_{etq}"] = round(an, 1) if an == an else ""
                    fila[f"AF_{a}_{etq}"] = round(af, 8) if af == af else ""

            an_imx = num(info, "AN_IMX")
            fila["n_efectivo_IMX"] = round(an_imx / 2, 0) if an_imx == an_imx else ""

            filas.append(fila)

    if not filas:
        sys.exit("No se parseo ninguna variante.")

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    with open(SALIDA, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for f in filas:
            fh.write("\t".join(str(f.get(c, "")) for c in cols) + "\n")

    # ------------------------------------------------------------------
    print(f"{n_lineas} lineas leidas, {len(filas)} variantes parseadas\n")

    # Tamano de muestra efectivo por fuente: la razon de estratificar
    print("Tamano de muestra efectivo para ancestria indigena (n = AN/2), por fuente:")
    print(f"{'fuente':<10}{'n variantes':>13}{'n_ef mediano':>15}{'min':>10}{'max':>10}")
    por_source = {}
    for f in filas:
        por_source.setdefault(f["source"], []).append(f)
    for src, g in sorted(por_source.items(), key=lambda x: -len(x[1])):
        ns = sorted(x["n_efectivo_IMX"] for x in g
                    if isinstance(x["n_efectivo_IMX"], float))
        if ns:
            print(f"{src or '(sin)':<10}{len(g):>13}{ns[len(ns)//2]:>15,.0f}"
                  f"{ns[0]:>10,.0f}{ns[-1]:>10,.0f}")
        else:
            print(f"{src or '(sin)':<10}{len(g):>13}{'-':>15}")

    n_rep = sum(f["reparto_default"] for f in filas)
    n_ind = sum(f["es_indel"] for f in filas)
    print(f"\nIndels: {n_ind}")
    print(f"Con reparto por defecto (AC=0.5 en dos ancestrias): {n_rep} "
          f"({100*n_rep/len(filas):.1f}%)  -> excluir del analisis")

    por_clase = {}
    for f in filas:
        por_clase[f["clase"] or "(fuera)"] = por_clase.get(f["clase"] or "(fuera)", 0) + 1
    print("\nPor clase de gen:")
    for k, v in sorted(por_clase.items(), key=lambda x: -x[1]):
        print(f"  {k:<14} {v:>8}")

    # Variantes analizables: comunes, sin reparto, no indel, por fuente
    print("\nVariantes analizables (AF_IMX entre 0.01 y 0.99, sin reparto, SNV):")
    print(f"{'fuente':<10}{'candidato':>11}{'control':>10}{'total':>8}")
    for src, g in sorted(por_source.items(), key=lambda x: -len(x[1])):
        util = [f for f in g
                if not f["reparto_default"] and not f["es_indel"]
                and isinstance(f["AF_IMX_PLANA"], float)
                and 0.01 <= f["AF_IMX_PLANA"] <= 0.99]
        nc = sum(1 for f in util if f["clase"] == "candidato")
        nk = sum(1 for f in util if f["clase"] == "control")
        print(f"{src or '(sin)':<10}{nc:>11}{nk:>10}{len(util):>8}")

    print(f"\nSalida -> {SALIDA}")
    print("\nEl analisis debe estratificarse por 'source': el tamano de muestra")
    print("efectivo difiere en mas de un orden de magnitud entre WGS y WES/Array.")


if __name__ == "__main__":
    main()
