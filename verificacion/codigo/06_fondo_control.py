#!/usr/bin/env python3
"""
Descarga las variantes de un conjunto de genes control no relacionados con
metales, para construir un fondo independiente de los genes candidatos.

Motivacion: el fondo del script 05 sale de los mismos 13 genes candidatos. Si
esos genes comparten historia selectiva, el resultado nulo podria ser circular:
las candidatas no destacarian porque todo su entorno esta igualmente desplazado.
Un fondo de genes sin relacion con el fenotipo cierra esa objecion.

Criterios del conjunto control (data/raw/genes_control.tsv): genes autosomicos
de expresion constitutiva, sin funcion conocida en transporte de metales,
metabolismo de xenobioticos ni respuesta inmune, lejos de la region HLA y sin
senales de seleccion reportadas.

Consulta por ventanas de 30 kb con cache por ventana, igual que el script 03c,
porque las regiones grandes devuelven HTTP 500.

Correr: python3 scripts/06_fondo_control.py
Salida: data/raw/gnomad_cache_control/ (una ventana por archivo)
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

VENTANA = 30000
PADDING = 5000
PAUSA = 8
REINTENTOS = 4

RAIZ = Path(__file__).resolve().parent.parent
LISTA = RAIZ / "data" / "raw" / "genes_control.tsv"
CACHE = RAIZ / "data" / "raw" / "gnomad_cache_control"

Q_COORDENADAS = """
query Coords($simbolo: String!, $referencia: ReferenceGenomeId!) {
  gene(gene_symbol: $simbolo, reference_genome: $referencia) {
    gene_id symbol chrom start stop
  }
}
"""

Q_REGION = """
query Region($chrom: String!, $start: Int!, $stop: Int!,
             $referencia: ReferenceGenomeId!, $dataset: DatasetId!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: $referencia) {
    variants(dataset: $dataset) {
      variant_id rsids consequence
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
            if r.status_code in (429, 500, 502, 503):
                espera = 20 * intento
                print(f"    HTTP {r.status_code} en {etiqueta}; espera {espera}s",
                      file=sys.stderr)
                time.sleep(espera)
                continue
            r.raise_for_status()
            d = r.json()
            if "errors" in d:
                print(f"    error de API en {etiqueta}: {d['errors']}", file=sys.stderr)
                return None
            return d["data"]
        except Exception as exc:
            espera = 10 * intento
            print(f"    intento {intento} fallo en {etiqueta} ({exc}); espera {espera}s",
                  file=sys.stderr)
            time.sleep(espera)
    return None


def main():
    if not LISTA.exists():
        sys.exit(f"Falta {LISTA}")

    genes = []
    with open(LISTA, encoding="utf-8") as fh:
        fh.readline()
        for linea in fh:
            if linea.strip():
                genes.append(linea.split("\t")[0].strip())

    CACHE.mkdir(parents=True, exist_ok=True)
    print(f"{len(genes)} genes control\n")

    total_v = 0
    for i, simbolo in enumerate(genes, 1):
        hechas = sorted(CACHE.glob(f"{simbolo}_v*.json"))
        marca = CACHE / f"{simbolo}_COMPLETO.txt"
        if marca.exists():
            n = sum(len(json.loads(f.read_text())["variants"]) for f in hechas)
            print(f"[{i}/{len(genes)}] {simbolo}: completo en cache ({n} variantes)")
            total_v += n
            continue

        print(f"[{i}/{len(genes)}] {simbolo}: coordenadas ...")
        d = peticion(Q_COORDENADAS, {"simbolo": simbolo, "referencia": REFERENCE},
                     f"{simbolo} coords")
        if not d or not d.get("gene"):
            print(f"  sin coordenadas para {simbolo}", file=sys.stderr)
            continue
        g = d["gene"]
        inicio = max(1, g["start"] - PADDING)
        fin = g["stop"] + PADDING
        n_vent = (fin - inicio + VENTANA - 1) // VENTANA
        print(f"  chr{g['chrom']}:{inicio}-{fin} ({fin-inicio} pb, {n_vent} ventanas)")
        time.sleep(PAUSA)

        n_gen, completo = 0, True
        for k in range(n_vent):
            v_ini = inicio + k * VENTANA
            v_fin = min(fin, v_ini + VENTANA - 1)
            cache_v = CACHE / f"{simbolo}_v{k:02d}.json"

            if cache_v.exists():
                n_gen += len(json.loads(cache_v.read_text())["variants"])
                continue

            d = peticion(Q_REGION,
                         {"chrom": g["chrom"], "start": v_ini, "stop": v_fin,
                          "referencia": REFERENCE, "dataset": DATASET},
                         f"{simbolo} v{k}")
            if not d or not d.get("region"):
                print(f"    ventana {k+1}/{n_vent} sin resultado", file=sys.stderr)
                completo = False
                time.sleep(PAUSA)
                continue
            vs = d["region"].get("variants") or []
            cache_v.write_text(json.dumps({
                "fecha_consulta_utc": datetime.now(timezone.utc).isoformat(),
                "simbolo": simbolo, "gene_id": g["gene_id"], "chrom": g["chrom"],
                "start": v_ini, "stop": v_fin, "variants": vs,
            }))
            n_gen += len(vs)
            print(f"    ventana {k+1}/{n_vent}: {len(vs)} variantes")
            time.sleep(PAUSA)

        if completo:
            marca.write_text(datetime.now(timezone.utc).isoformat())
        print(f"  {simbolo}: {n_gen} variantes"
              f"{'' if completo else '  (INCOMPLETO, volver a correr)'}")
        total_v += n_gen

    print(f"\nTotal: {total_v} variantes en {CACHE}")
    print("Si algun gen quedo incompleto, volver a correr: retoma donde se quedo.")


if __name__ == "__main__":
    main()
