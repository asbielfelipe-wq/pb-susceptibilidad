#!/usr/bin/env python3
"""
Rescate de las variantes que la consulta por gen no recupero.

Dos causas distintas:
  1. Variantes de promotor/regulatorias que caen rio arriba del gen y por eso
     no aparecen en una consulta acotada a los limites genicos.
  2. Genes que fallaron por limite de peticiones (HTTP 429).

Estrategia: obtener las coordenadas del gen y consultar la REGION con padding
rio arriba y rio abajo. Guarda cada gen en cache para poder reanudar sin
repetir consultas ya hechas.

Correr: python3 scripts/03b_rescate_regiones.py
Salida: data/processed/gnomad_frecuencias_rescate.tsv
"""

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
REFERENCE = "GRCh38"

PADDING_5P = 20000   # pb rio arriba: cubre promotores tipo -617, -308, -262
PADDING_3P = 5000    # pb rio abajo

PAUSA = 8            # segundos entre consultas; el 429 aparece con pausas cortas
REINTENTOS = 4

RAIZ = Path(__file__).resolve().parent.parent
CANDIDATOS = RAIZ / "data" / "raw" / "candidatos_literatura.tsv"
CACHE = RAIZ / "data" / "raw" / "gnomad_cache"
SALIDA = RAIZ / "data" / "processed" / "gnomad_frecuencias_rescate.tsv"

POBLACIONES = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining"]

# Genes pendientes: los que dieron cero por ser regulatorios y el que fallo por 429.
PENDIENTES = ["ALAD", "SLC30A10", "NFE2L2", "CAT", "SOD2", "TNF"]

Q_COORDENADAS = """
query Coords($simbolo: String!, $referencia: ReferenceGenomeId!) {
  gene(gene_symbol: $simbolo, reference_genome: $referencia) {
    gene_id
    symbol
    chrom
    start
    stop
  }
}
"""

Q_REGION = """
query Region($chrom: String!, $start: Int!, $stop: Int!,
             $referencia: ReferenceGenomeId!, $dataset: DatasetId!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: $referencia) {
    variants(dataset: $dataset) {
      variant_id
      rsids
      consequence
      hgvsp
      genome { ac an populations { id ac an } }
      exome  { ac an populations { id ac an } }
    }
  }
}
"""


def peticion(consulta, variables, etiqueta):
    carga = {"query": consulta, "variables": variables}
    for intento in range(1, REINTENTOS + 1):
        try:
            r = requests.post(API, json=carga, timeout=300)
            if r.status_code == 429:
                espera = 30 * intento
                print(f"  429 en {etiqueta}; esperando {espera}s", file=sys.stderr)
                time.sleep(espera)
                continue
            r.raise_for_status()
            datos = r.json()
            if "errors" in datos:
                print(f"  error de API en {etiqueta}: {datos['errors']}", file=sys.stderr)
                return None
            return datos["data"]
        except Exception as exc:
            espera = 10 * intento
            print(f"  intento {intento} fallo en {etiqueta} ({exc}); espera {espera}s",
                  file=sys.stderr)
            time.sleep(espera)
    return None


def rsids_objetivo(ruta):
    objetivo = {}
    with open(ruta, encoding="utf-8") as fh:
        enc = fh.readline().rstrip("\n").split("\t")
        i_gen, i_rs, i_tipo = (enc.index(c) for c in ("gen", "rsid", "tipo_variante"))
        for linea in fh:
            c = linea.rstrip("\n").split("\t")
            if c[i_tipo] in ("CNV", "STR") or not c[i_rs] or c[i_rs] == "NA":
                continue
            objetivo.setdefault(c[i_gen].split(" ")[0], set()).add(c[i_rs])
    return objetivo


def frecuencias(bloque):
    if not bloque:
        return {}
    return {p["id"]: (p["ac"], p["an"])
            for p in (bloque.get("populations") or [])
            if p["id"] in POBLACIONES}


def combinar(a, b):
    salida = {}
    for p in set(a) | set(b):
        salida[p] = (a.get(p, (0, 0))[0] + b.get(p, (0, 0))[0],
                     a.get(p, (0, 0))[1] + b.get(p, (0, 0))[1])
    return salida


def main():
    objetivo = rsids_objetivo(CANDIDATOS)
    CACHE.mkdir(parents=True, exist_ok=True)
    filas = []

    for simbolo in PENDIENTES:
        buscados = objetivo.get(simbolo, set())
        if not buscados:
            print(f"{simbolo}: sin rsID objetivo, se omite")
            continue

        archivo_cache = CACHE / f"{simbolo}_region.json"

        if archivo_cache.exists():
            print(f"{simbolo}: usando cache")
            region = json.loads(archivo_cache.read_text())["region"]
        else:
            print(f"{simbolo}: coordenadas ...")
            d = peticion(Q_COORDENADAS,
                         {"simbolo": simbolo, "referencia": REFERENCE},
                         f"{simbolo} coords")
            if not d or not d.get("gene"):
                print(f"  no se obtuvieron coordenadas para {simbolo}", file=sys.stderr)
                continue
            g = d["gene"]
            inicio = max(1, g["start"] - PADDING_5P)
            fin = g["stop"] + PADDING_3P
            print(f"  chr{g['chrom']}:{inicio}-{fin} "
                  f"({fin - inicio} pb, gen {g['start']}-{g['stop']})")
            time.sleep(PAUSA)

            print(f"{simbolo}: variantes de la region ...")
            d = peticion(Q_REGION,
                         {"chrom": g["chrom"], "start": inicio, "stop": fin,
                          "referencia": REFERENCE, "dataset": DATASET},
                         f"{simbolo} region")
            if not d or not d.get("region"):
                print(f"  sin variantes para {simbolo}", file=sys.stderr)
                continue
            region = d["region"]
            archivo_cache.write_text(json.dumps({
                "fecha_consulta_utc": datetime.now(timezone.utc).isoformat(),
                "simbolo": simbolo,
                "gene_id": g["gene_id"],
                "chrom": g["chrom"],
                "start_consulta": inicio,
                "stop_consulta": fin,
                "region": region,
            }))
            time.sleep(PAUSA)

        encontrados = 0
        for v in region.get("variants") or []:
            coincide = set(v.get("rsids") or []) & buscados
            if not coincide:
                continue
            encontrados += 1
            comb = combinar(frecuencias(v.get("genome")), frecuencias(v.get("exome")))
            fila = {
                "gen": simbolo,
                "rsid": ";".join(sorted(coincide)),
                "variant_id": v["variant_id"],
                "consecuencia": v.get("consequence") or "",
                "hgvsp": v.get("hgvsp") or "",
            }
            for p in POBLACIONES:
                ac, an = comb.get(p, (0, 0))
                fila[f"AC_{p}"] = ac
                fila[f"AN_{p}"] = an
                fila[f"AF_{p}"] = round(ac / an, 6) if an else ""
            filas.append(fila)

        print(f"  {encontrados} fila(s) para {len(buscados)} rsID objetivo")
        faltan = buscados - {r for f in filas if f["gen"] == simbolo
                             for r in f["rsid"].split(";")}
        if faltan:
            print(f"  sin recuperar: {', '.join(sorted(faltan))}")

    if not filas:
        print("\nNo se recupero nada. Revisar el esquema en https://gnomad.broadinstitute.org/api",
              file=sys.stderr)
        sys.exit(1)

    columnas = list(filas[0].keys())
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    with open(SALIDA, "w", encoding="utf-8") as fh:
        fh.write("\t".join(columnas) + "\n")
        for f in filas:
            fh.write("\t".join(str(f[c]) for c in columnas) + "\n")

    print(f"\n{len(filas)} filas escritas en {SALIDA}")
    print(f"Cache por gen en {CACHE} (borrar un archivo fuerza su reconsulta)")


if __name__ == "__main__":
    main()
