#!/usr/bin/env python3
"""
Consulta frecuencias alelicas por poblacion en gnomAD v4.1 para los genes
candidatos de susceptibilidad al plomo.

Estrategia: se consulta por GEN (no por rsID) porque la busqueda por rsID en la
API es menos estable entre versiones. Despues se filtra localmente por los rsID
de la tabla de candidatos.

Correr: python3 scripts/03_gnomad_frecuencias.py
Salida: data/processed/gnomad_frecuencias.tsv

ADVERTENCIA: el esquema GraphQL de gnomAD cambia entre versiones. Si la consulta
falla, abrir https://gnomad.broadinstitute.org/api y revisar el esquema antes de
asumir que el problema es de red. Dejar registrada la fecha de consulta.
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

RAIZ = Path(__file__).resolve().parent.parent
ENTRADA = RAIZ / "data" / "raw" / "candidatos_literatura.tsv"
SALIDA = RAIZ / "data" / "processed" / "gnomad_frecuencias.tsv"
CRUDO = RAIZ / "data" / "raw" / "gnomad_respuestas.json"

# Poblaciones de gnomAD v4. AMR es la de interes; las demas son comparadores.
POBLACIONES = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining"]

CONSULTA = """
query VariantesDeGen($simbolo: String!, $referencia: ReferenceGenomeId!, $dataset: DatasetId!) {
  gene(gene_symbol: $simbolo, reference_genome: $referencia) {
    gene_id
    symbol
    variants(dataset: $dataset) {
      variant_id
      rsids
      consequence
      hgvsp
      genome {
        ac
        an
        populations { id ac an }
      }
      exome {
        ac
        an
        populations { id ac an }
      }
    }
  }
}
"""


def leer_candidatos(ruta):
    genes, rsids = [], set()
    with open(ruta, encoding="utf-8") as fh:
        encabezado = fh.readline().rstrip("\n").split("\t")
        i_gen, i_rs, i_tipo = (encabezado.index(c) for c in ("gen", "rsid", "tipo_variante"))
        for linea in fh:
            campos = linea.rstrip("\n").split("\t")
            if campos[i_tipo] in ("CNV", "STR"):
                continue  # no consultables como SNV
            gen = campos[i_gen].split(" ")[0]  # "SLC11A2 (DMT1)" -> "SLC11A2"
            if gen not in genes:
                genes.append(gen)
            if campos[i_rs] and campos[i_rs] != "NA":
                rsids.add(campos[i_rs])
    return genes, rsids


def consultar_gen(simbolo, intentos=3):
    carga = {
        "query": CONSULTA,
        "variables": {"simbolo": simbolo, "referencia": REFERENCE, "dataset": DATASET},
    }
    for intento in range(1, intentos + 1):
        try:
            r = requests.post(API, json=carga, timeout=180)
            r.raise_for_status()
            datos = r.json()
            if "errors" in datos:
                print(f"  error de la API en {simbolo}: {datos['errors']}", file=sys.stderr)
                return None
            return datos["data"]["gene"]
        except Exception as exc:
            espera = 5 * intento
            print(f"  intento {intento} fallo ({exc}); reintentando en {espera}s", file=sys.stderr)
            time.sleep(espera)
    return None


def frecuencias(bloque):
    """Devuelve dict poblacion -> (ac, an) desde un bloque genome o exome."""
    if not bloque:
        return {}
    salida = {}
    for p in bloque.get("populations") or []:
        pid = p["id"]
        # gnomAD devuelve tambien estratos por sexo (ej. "amr_XX"); se omiten.
        if pid in POBLACIONES:
            salida[pid] = (p["ac"], p["an"])
    return salida


def combinar(gen_fr, exo_fr):
    """Suma conteos de genomas y exomas por poblacion."""
    todas = set(gen_fr) | set(exo_fr)
    salida = {}
    for p in todas:
        ac = gen_fr.get(p, (0, 0))[0] + exo_fr.get(p, (0, 0))[0]
        an = gen_fr.get(p, (0, 0))[1] + exo_fr.get(p, (0, 0))[1]
        salida[p] = (ac, an)
    return salida


def main():
    genes, rsids = leer_candidatos(ENTRADA)
    print(f"{len(genes)} genes, {len(rsids)} rsID objetivo")

    respuestas, filas = {}, []

    for simbolo in genes:
        print(f"consultando {simbolo} ...")
        gen = consultar_gen(simbolo)
        if gen is None:
            print(f"  sin resultado para {simbolo}", file=sys.stderr)
            continue
        respuestas[simbolo] = gen

        encontrados = 0
        for v in gen.get("variants") or []:
            coincide = set(v.get("rsids") or []) & rsids
            if not coincide:
                continue
            encontrados += 1
            comb = combinar(frecuencias(v.get("genome")), frecuencias(v.get("exome")))
            fila = {
                "gen": gen["symbol"],
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
        print(f"  {encontrados} variante(s) candidata(s)")
        time.sleep(1)  # cortesia con la API

    CRUDO.parent.mkdir(parents=True, exist_ok=True)
    with open(CRUDO, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "fecha_consulta_utc": datetime.now(timezone.utc).isoformat(),
                "dataset": DATASET,
                "reference_genome": REFERENCE,
                "respuestas": respuestas,
            },
            fh,
        )

    if not filas:
        print("No se recupero ninguna variante. Revisar el esquema de la API.", file=sys.stderr)
        sys.exit(1)

    columnas = list(filas[0].keys())
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    with open(SALIDA, "w", encoding="utf-8") as fh:
        fh.write("\t".join(columnas) + "\n")
        for fila in filas:
            fh.write("\t".join(str(fila[c]) for c in columnas) + "\n")

    print(f"\n{len(filas)} filas escritas en {SALIDA}")
    print(f"Respuestas crudas en {CRUDO} (incluye fecha de consulta para Metodos)")


if __name__ == "__main__":
    main()
