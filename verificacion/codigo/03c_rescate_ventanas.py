#!/usr/bin/env python3
"""
Tercer intento para los genes cuya region es demasiado grande para una sola
consulta: gnomAD devuelve HTTP 500 cuando la region contiene demasiadas
variantes. Se trocea en ventanas de 30 kb y se unen los resultados.

Cada ventana se guarda en cache por separado, asi que una interrupcion no
obliga a repetir lo ya consultado.

Correr: python3 scripts/03c_rescate_ventanas.py
Salida: data/processed/gnomad_frecuencias_ventanas.tsv
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

VENTANA = 30000      # pb por consulta; por encima de ~60 kb el servidor da 500
PADDING_5P = 20000
PADDING_3P = 5000
PAUSA = 8
REINTENTOS = 4

RAIZ = Path(__file__).resolve().parent.parent
CANDIDATOS = RAIZ / "data" / "raw" / "candidatos_literatura.tsv"
CACHE = RAIZ / "data" / "raw" / "gnomad_cache"
SALIDA = RAIZ / "data" / "processed" / "gnomad_frecuencias_ventanas.tsv"

POBLACIONES = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining"]

PENDIENTES = ["SLC30A10", "NFE2L2", "SOD2"]

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
      variant_id rsids consequence hgvsp
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
            datos = r.json()
            if "errors" in datos:
                print(f"    error de API en {etiqueta}: {datos['errors']}", file=sys.stderr)
                return None
            return datos["data"]
        except Exception as exc:
            espera = 10 * intento
            print(f"    intento {intento} fallo en {etiqueta} ({exc}); espera {espera}s",
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
    return {p: (a.get(p, (0, 0))[0] + b.get(p, (0, 0))[0],
                a.get(p, (0, 0))[1] + b.get(p, (0, 0))[1])
            for p in set(a) | set(b)}


def main():
    objetivo = rsids_objetivo(CANDIDATOS)
    CACHE.mkdir(parents=True, exist_ok=True)
    filas = []

    for simbolo in PENDIENTES:
        buscados = objetivo.get(simbolo, set())
        if not buscados:
            continue

        print(f"\n{simbolo}: coordenadas ...")
        d = peticion(Q_COORDENADAS, {"simbolo": simbolo, "referencia": REFERENCE},
                     f"{simbolo} coords")
        if not d or not d.get("gene"):
            print(f"  no se obtuvieron coordenadas", file=sys.stderr)
            continue
        g = d["gene"]
        inicio = max(1, g["start"] - PADDING_5P)
        fin = g["stop"] + PADDING_3P
        total = fin - inicio
        n_vent = (total + VENTANA - 1) // VENTANA
        print(f"  chr{g['chrom']}:{inicio}-{fin} ({total} pb en {n_vent} ventanas)")
        time.sleep(PAUSA)

        variantes = []
        for i in range(n_vent):
            v_ini = inicio + i * VENTANA
            v_fin = min(fin, v_ini + VENTANA - 1)
            cache_v = CACHE / f"{simbolo}_v{i:02d}.json"

            if cache_v.exists():
                print(f"  ventana {i+1}/{n_vent}: cache")
                variantes.extend(json.loads(cache_v.read_text())["variants"])
                continue

            print(f"  ventana {i+1}/{n_vent}: {v_ini}-{v_fin} ...")
            d = peticion(Q_REGION,
                         {"chrom": g["chrom"], "start": v_ini, "stop": v_fin,
                          "referencia": REFERENCE, "dataset": DATASET},
                         f"{simbolo} v{i}")
            if not d or not d.get("region"):
                print(f"    ventana {i+1} sin resultado", file=sys.stderr)
                time.sleep(PAUSA)
                continue
            vs = d["region"].get("variants") or []
            cache_v.write_text(json.dumps({
                "fecha_consulta_utc": datetime.now(timezone.utc).isoformat(),
                "simbolo": simbolo, "chrom": g["chrom"],
                "start": v_ini, "stop": v_fin, "variants": vs,
            }))
            print(f"    {len(vs)} variantes")
            variantes.extend(vs)
            time.sleep(PAUSA)

        encontrados = 0
        vistos = set()
        for v in variantes:
            if v["variant_id"] in vistos:
                continue
            vistos.add(v["variant_id"])
            coincide = set(v.get("rsids") or []) & buscados
            if not coincide:
                continue
            encontrados += 1
            comb = combinar(frecuencias(v.get("genome")), frecuencias(v.get("exome")))
            fila = {"gen": simbolo, "rsid": ";".join(sorted(coincide)),
                    "variant_id": v["variant_id"],
                    "consecuencia": v.get("consequence") or "",
                    "hgvsp": v.get("hgvsp") or ""}
            for p in POBLACIONES:
                ac, an = comb.get(p, (0, 0))
                fila[f"AC_{p}"] = ac
                fila[f"AN_{p}"] = an
                fila[f"AF_{p}"] = round(ac / an, 6) if an else ""
            filas.append(fila)

        print(f"  {encontrados} fila(s) de {len(buscados)} rsID objetivo "
              f"({len(vistos)} variantes revisadas)")

    if not filas:
        print("\nNo se recupero nada.", file=sys.stderr)
        sys.exit(1)

    columnas = list(filas[0].keys())
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    with open(SALIDA, "w", encoding="utf-8") as fh:
        fh.write("\t".join(columnas) + "\n")
        for f in filas:
            fh.write("\t".join(str(f[c]) for c in columnas) + "\n")
    print(f"\n{len(filas)} filas escritas en {SALIDA}")


if __name__ == "__main__":
    main()
