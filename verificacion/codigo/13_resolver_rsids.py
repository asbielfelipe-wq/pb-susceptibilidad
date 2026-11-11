#!/usr/bin/env python3
"""
Resuelve las coordenadas GRCh38 de los rsID del conjunto de Ding et al. 2016 y
los busca en la tabla de frecuencias de MCPS.

El VCF de MCPS no incluye rsID: identifica las variantes por posicion y alelos.
Para cruzarlo con una lista de rsID hace falta resolver primero cada rsID a sus
coordenadas. Se usa gnomAD, que si indexa rsID, consultando por gen y filtrando
localmente (la busqueda directa por rsID en esa API es menos estable).

Correr: python3 scripts/13_resolver_rsids.py
Salidas:
  data/raw/ding_coordenadas.tsv      rsID -> variant_id, con alelos
  results/ding_en_mcps.tsv           las que se encontraron en MCPS, con AF
  results/ding_cobertura.txt
"""

import sys
import time
from pathlib import Path

import pandas as pd
import requests

API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
REFERENCE = "GRCh38"
PAUSA = 6
REINTENTOS = 4

RAIZ = Path(__file__).resolve().parent.parent
LISTA = RAIZ / "data" / "raw" / "candidatos_ding2016.tsv"
MCPS = RAIZ / "data" / "processed" / "mcps_frecuencias.tsv"
COORD = RAIZ / "data" / "raw" / "ding_coordenadas.tsv"
RES = RAIZ / "results"

VERSION = "PLANA"

Q = """
query G($s: String!, $r: ReferenceGenomeId!, $d: DatasetId!) {
  gene(gene_symbol: $s, reference_genome: $r) {
    symbol
    variants(dataset: $d) {
      variant_id rsids consequence hgvsp
      genome { ac an }
      exome  { ac an }
    }
  }
}
"""

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def consultar(simbolo):
    carga = {"query": Q, "variables": {"s": simbolo, "r": REFERENCE, "d": DATASET}}
    for intento in range(1, REINTENTOS + 1):
        try:
            r = requests.post(API, json=carga, timeout=300)
            if r.status_code in (429, 500, 502, 503):
                espera = 20 * intento
                print(f"    HTTP {r.status_code}; espera {espera}s", file=sys.stderr)
                time.sleep(espera)
                continue
            r.raise_for_status()
            d = r.json()
            if "errors" in d:
                print(f"    error API: {d['errors']}", file=sys.stderr)
                return None
            return d["data"]["gene"]
        except Exception as exc:
            espera = 10 * intento
            print(f"    intento {intento} fallo ({exc}); espera {espera}s",
                  file=sys.stderr)
            time.sleep(espera)
    return None


def main():
    if not LISTA.exists():
        sys.exit(f"Falta {LISTA}")

    dj = pd.read_csv(LISTA, sep="\t", dtype=str).fillna("")
    snv = dj[dj["tipo"] == "SNV"].copy()
    objetivo = {}
    for _, r in snv.iterrows():
        objetivo.setdefault(r["gen"], set()).add(r["rsid"])

    di("=" * 72)
    di("RESOLUCION DE COORDENADAS DE LOS LOCI DE DING ET AL. 2016")
    di("=" * 72)
    di(f"\n{len(snv)} loci SNV en {len(objetivo)} genes")
    excl = dj[dj["tipo"] != "SNV"]
    if len(excl):
        di(f"{len(excl)} excluidos por tipo (CNV/STR): "
           f"{', '.join(excl['gen'] + ' ' + excl['nomenclatura'])}")

    # --- cache -----------------------------------------------------------
    if COORD.exists():
        prev = pd.read_csv(COORD, sep="\t", dtype=str)
        di(f"\nUsando cache: {COORD.name} ({len(prev)} filas)")
        coords = prev
    else:
        filas = []
        for simbolo, rsids in objetivo.items():
            print(f"consultando {simbolo} ({len(rsids)} rsID) ...")
            g = consultar(simbolo)
            if not g:
                print(f"  sin resultado para {simbolo}", file=sys.stderr)
                continue
            hallados = set()
            for v in g.get("variants") or []:
                coincide = set(v.get("rsids") or []) & rsids
                if not coincide:
                    continue
                ac = ((v.get("genome") or {}).get("ac") or 0) + \
                     ((v.get("exome") or {}).get("ac") or 0)
                for rs in coincide:
                    filas.append({
                        "gen": simbolo, "rsid": rs,
                        "variant_id": v["variant_id"],
                        "consecuencia": v.get("consequence") or "",
                        "hgvsp": v.get("hgvsp") or "",
                        "AC_global_gnomad": ac,
                    })
                    hallados.add(rs)
            faltan = rsids - hallados
            print(f"  {len(hallados)} de {len(rsids)}"
                  + (f"; sin resolver: {', '.join(sorted(faltan))}" if faltan else ""))
            time.sleep(PAUSA)

        if not filas:
            sys.exit("No se resolvio ningun rsID.")
        coords = pd.DataFrame(filas)
        # En sitios multialelicos se conserva el alelo con mayor AC global
        coords["AC_global_gnomad"] = pd.to_numeric(coords["AC_global_gnomad"])
        coords = (coords.sort_values("AC_global_gnomad", ascending=False)
                        .drop_duplicates(subset="rsid", keep="first"))
        coords.to_csv(COORD, sep="\t", index=False)
        di(f"\nCoordenadas -> {COORD}")

    di(f"\n{len(coords)} de {len(snv)} rsID resueltos a coordenadas GRCh38")
    sin_coord = set(snv["rsid"]) - set(coords["rsid"])
    if sin_coord:
        di(f"Sin coordenadas: {', '.join(sorted(sin_coord))}")

    # --- cruce con MCPS --------------------------------------------------
    if not MCPS.exists():
        sys.exit(f"\nFalta {MCPS}. Correr antes el script 10.")

    mc = pd.read_csv(MCPS, sep="\t", low_memory=False)
    cols_af = [f"AF_{a}_{VERSION}" for a in ("IMX", "EUR", "AFR")]
    j = coords.merge(
        mc[["variant_id", "gen", "source", "clase", "reparto_default",
            "n_efectivo_IMX"] + cols_af],
        on="variant_id", how="left", suffixes=("", "_mcps"))

    hallados = j[j["source"].notna()]
    perdidos = j[j["source"].isna()]

    di("\n" + "=" * 72)
    di("COBERTURA EN MCPS")
    di("=" * 72)
    di(f"\n{len(hallados)} de {len(coords)} presentes en MCPS")

    if len(hallados):
        di(f"\n{'gen':<9}{'rsid':<13}{'variante':<20}{'fuente':<7}"
           f"{'AF_IMX':>8}{'AF_EUR':>8}{'rep':>5}")
        for _, r in hallados.sort_values(["gen", "rsid"]).iterrows():
            di(f"{r['gen']:<9}{r['rsid']:<13}{r['variant_id']:<20}"
               f"{r['source']:<7}{r[cols_af[0]]:>8.4f}{r[cols_af[1]]:>8.4f}"
               f"{int(r['reparto_default']):>5}")

    if len(perdidos):
        di(f"\n{len(perdidos)} no estan en las regiones extraidas de MCPS:")
        for _, r in perdidos.iterrows():
            di(f"  {r['gen']:<9}{r['rsid']:<13}{r['variant_id']}")
        di("\nSi caen fuera de las regiones del BED, ampliar el padding del")
        di("script 08 y repetir la extraccion.")

    por_fuente = hallados["source"].value_counts().to_dict() if len(hallados) else {}
    di(f"\nPor fuente: " + ", ".join(f"{k}={v}" for k, v in por_fuente.items()))

    RES.mkdir(parents=True, exist_ok=True)
    j.to_csv(RES / "ding_en_mcps.tsv", sep="\t", index=False)
    (RES / "ding_cobertura.txt").write_text("\n".join(salida))
    print(f"\n\nTabla  -> {RES / 'ding_en_mcps.tsv'}")
    print(f"Resumen-> {RES / 'ding_cobertura.txt'}")


if __name__ == "__main__":
    main()
