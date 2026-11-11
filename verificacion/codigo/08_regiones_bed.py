#!/usr/bin/env python3
"""
Obtiene las coordenadas GRCh38 de los genes candidatos y control y escribe un
archivo BED para extraer esas regiones de los VCF de MCPS con tabix.

Los VCF de MCPS estan indexados, asi que no hace falta descargar los 4.4 GB:
tabix baja solo los bloques de las regiones pedidas.

Correr: python3 scripts/08_regiones_bed.py
Salida: data/raw/regiones.bed  y  data/raw/regiones.txt (formato tabix)
"""

import sys
import time
from pathlib import Path

import requests

API = "https://gnomad.broadinstitute.org/api"
REFERENCE = "GRCh38"
PADDING_5P = 20000
PADDING_3P = 5000
PAUSA = 3

RAIZ = Path(__file__).resolve().parent.parent
RAW = RAIZ / "data" / "raw"

# Genes candidatos: se leen de data/raw/genes_candidatos.tsv si existe
# (conjunto derivado de Ding et al. 2016, PMID 27584680); si no, se usa la
# lista previa, conservada aqui como analisis exploratorio.
CANDIDATOS_PREVIOS = ["ALAD", "VDR", "HFE", "SLC11A2", "SLC30A10", "MTHFR",
                      "GSTP1", "APOE", "NFE2L2", "CAT", "SOD2", "TNF", "ABCB1"]


def cargar_candidatos():
    f = RAW / "genes_candidatos.tsv"
    if not f.exists():
        print("aviso: no esta genes_candidatos.tsv; se usa la lista previa",
              file=sys.stderr)
        return CANDIDATOS_PREVIOS
    genes = []
    with open(f, encoding="utf-8") as fh:
        fh.readline()
        for linea in fh:
            if linea.strip():
                genes.append(linea.split("\t")[0].strip())
    return genes

Q = """
query C($s: String!, $r: ReferenceGenomeId!) {
  gene(gene_symbol: $s, reference_genome: $r) { gene_id symbol chrom start stop }
}
"""


def coords(simbolo):
    for intento in range(1, 4):
        try:
            r = requests.post(API, json={"query": Q,
                                         "variables": {"s": simbolo, "r": REFERENCE}},
                              timeout=120)
            if r.status_code in (429, 500, 502, 503):
                time.sleep(15 * intento)
                continue
            r.raise_for_status()
            d = r.json()
            if "errors" in d:
                return None
            return d["data"]["gene"]
        except Exception as exc:
            print(f"  intento {intento} fallo ({exc})", file=sys.stderr)
            time.sleep(10 * intento)
    return None


def main():
    ctrl_f = RAW / "genes_control.tsv"
    controles = []
    if ctrl_f.exists():
        with open(ctrl_f, encoding="utf-8") as fh:
            fh.readline()
            controles = [l.split("\t")[0].strip() for l in fh if l.strip()]
    else:
        print(f"aviso: no esta {ctrl_f.name}; solo se haran los candidatos",
              file=sys.stderr)

    candidatos = cargar_candidatos()
    todos = [(g, "candidato") for g in candidatos] + \
            [(g, "control") for g in controles]
    print(f"{len(todos)} genes ({len(candidatos)} candidatos, "
          f"{len(controles)} control)\n")

    filas = []
    for simbolo, clase in todos:
        g = coords(simbolo)
        if not g:
            print(f"  SIN COORDENADAS: {simbolo}", file=sys.stderr)
            continue
        ini = max(1, g["start"] - PADDING_5P)
        fin = g["stop"] + PADDING_3P
        filas.append((g["chrom"], ini, fin, simbolo, clase, fin - ini))
        print(f"  {simbolo:<10} chr{g['chrom']}:{ini}-{fin}  "
              f"{fin-ini:>9} pb  [{clase}]")
        time.sleep(PAUSA)

    if not filas:
        sys.exit("No se obtuvo ninguna coordenada.")

    def orden(f):
        c = f[0]
        return (int(c) if c.isdigit() else 99, f[1])
    filas.sort(key=orden)

    with open(RAW / "regiones.bed", "w") as fh:
        for c, i, f, s, cl, n in filas:
            fh.write(f"chr{c}\t{i-1}\t{f}\t{s}_{cl}\n")

    with open(RAW / "regiones.txt", "w") as fh:
        for c, i, f, s, cl, n in filas:
            fh.write(f"chr{c}:{i}-{f}\n")

    total = sum(f[5] for f in filas)
    print(f"\n{len(filas)} regiones, {total/1e6:.1f} Mb en total")
    print(f"BED    -> {RAW / 'regiones.bed'}")
    print(f"tabix  -> {RAW / 'regiones.txt'}")
    print("\nOJO: revisar si los VCF de MCPS usan 'chr1' o '1' como nombre de")
    print("contig (tabix -H ... | grep contig) y ajustar el prefijo si hace falta.")


if __name__ == "__main__":
    main()
