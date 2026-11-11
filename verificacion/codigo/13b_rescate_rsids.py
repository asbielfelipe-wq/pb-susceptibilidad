#!/usr/bin/env python3
"""
Recupera los rsID del conjunto de Ding que la consulta por gen no resolvio.

Causa: la consulta por simbolo devuelve variantes dentro de los limites del
gen, y la mayoria de los loci pendientes son de promotor o regiones
reguladoras rio arriba. Se consulta por REGION con padding amplio y en
ventanas de 30 kb, porque las regiones grandes devuelven HTTP 500.

Los que sigan sin aparecer pueden ser rsID retirados o fusionados en dbSNP;
el script lo senala para verificacion manual.

Correr: python3 scripts/13b_rescate_rsids.py
Salidas:
  data/raw/ding_coordenadas_rescate.tsv
  results/ding_cobertura_final.txt
"""

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
REFERENCE = "GRCh38"

PADDING_5P = 50000     # amplio: cubre promotores lejanos
PADDING_3P = 20000
VENTANA = 30000
PAUSA = 8
REINTENTOS = 4

RAIZ = Path(__file__).resolve().parent.parent
LISTA = RAIZ / "data" / "raw" / "candidatos_ding2016.tsv"
PREV = RAIZ / "data" / "raw" / "ding_coordenadas.tsv"
SAL = RAIZ / "data" / "raw" / "ding_coordenadas_rescate.tsv"
MCPS = RAIZ / "data" / "processed" / "mcps_frecuencias.tsv"
CACHE = RAIZ / "data" / "raw" / "gnomad_cache_ding"
RES = RAIZ / "results"

VERSION = "PLANA"

Q_COORD = """
query C($s: String!, $r: ReferenceGenomeId!) {
  gene(gene_symbol: $s, reference_genome: $r) { symbol chrom start stop }
}
"""

Q_REG = """
query R($c: String!, $i: Int!, $f: Int!, $r: ReferenceGenomeId!, $d: DatasetId!) {
  region(chrom: $c, start: $i, stop: $f, reference_genome: $r) {
    variants(dataset: $d) {
      variant_id rsids consequence hgvsp
      genome { ac an } exome { ac an }
    }
  }
}
"""

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def peticion(q, v, etq):
    for intento in range(1, REINTENTOS + 1):
        try:
            r = requests.post(API, json={"query": q, "variables": v}, timeout=300)
            if r.status_code in (429, 500, 502, 503):
                esp = 20 * intento
                print(f"    HTTP {r.status_code} en {etq}; espera {esp}s",
                      file=sys.stderr)
                time.sleep(esp)
                continue
            r.raise_for_status()
            d = r.json()
            if "errors" in d:
                print(f"    error API en {etq}: {d['errors']}", file=sys.stderr)
                return None
            return d["data"]
        except Exception as exc:
            esp = 10 * intento
            print(f"    intento {intento} fallo en {etq} ({exc}); espera {esp}s",
                  file=sys.stderr)
            time.sleep(esp)
    return None


def main():
    dj = pd.read_csv(LISTA, sep="\t", dtype=str).fillna("")
    snv = dj[dj["tipo"] == "SNV"]

    resueltos = set()
    if PREV.exists():
        resueltos = set(pd.read_csv(PREV, sep="\t", dtype=str)["rsid"])

    pend = {}
    for _, r in snv.iterrows():
        if r["rsid"] not in resueltos:
            pend.setdefault(r["gen"], set()).add(r["rsid"])

    di("=" * 72)
    di("RESCATE DE LOCI DE DING NO RESUELTOS POR CONSULTA POR GEN")
    di("=" * 72)
    if not pend:
        di("\nNo quedan pendientes.")
        return
    di(f"\n{sum(len(v) for v in pend.values())} rsID pendientes en "
       f"{len(pend)} genes")
    for g, rs in pend.items():
        di(f"  {g:<9} {', '.join(sorted(rs))}")

    CACHE.mkdir(parents=True, exist_ok=True)
    filas = []

    for simbolo, rsids in pend.items():
        di(f"\n{simbolo}: coordenadas ...")
        d = peticion(Q_COORD, {"s": simbolo, "r": REFERENCE}, f"{simbolo} coords")
        if not d or not d.get("gene"):
            di(f"  sin coordenadas para {simbolo}")
            continue
        g = d["gene"]
        ini = max(1, g["start"] - PADDING_5P)
        fin = g["stop"] + PADDING_3P
        n_v = (fin - ini + VENTANA - 1) // VENTANA
        di(f"  chr{g['chrom']}:{ini}-{fin} ({fin-ini} pb, {n_v} ventanas)")
        time.sleep(PAUSA)

        variantes, vistos = [], set()
        for k in range(n_v):
            v_i = ini + k * VENTANA
            v_f = min(fin, v_i + VENTANA - 1)
            cf = CACHE / f"{simbolo}_v{k:02d}.json"
            if cf.exists():
                vs = json.loads(cf.read_text())["variants"]
            else:
                d = peticion(Q_REG, {"c": g["chrom"], "i": v_i, "f": v_f,
                                     "r": REFERENCE, "d": DATASET},
                             f"{simbolo} v{k}")
                if not d or not d.get("region"):
                    print(f"    ventana {k+1}/{n_v} sin resultado", file=sys.stderr)
                    time.sleep(PAUSA)
                    continue
                vs = d["region"].get("variants") or []
                cf.write_text(json.dumps({
                    "fecha_consulta_utc": datetime.now(timezone.utc).isoformat(),
                    "simbolo": simbolo, "chrom": g["chrom"],
                    "start": v_i, "stop": v_f, "variants": vs}))
                time.sleep(PAUSA)
            for v in vs:
                if v["variant_id"] not in vistos:
                    vistos.add(v["variant_id"])
                    variantes.append(v)

        hallados = set()
        for v in variantes:
            coincide = set(v.get("rsids") or []) & rsids
            if not coincide:
                continue
            ac = ((v.get("genome") or {}).get("ac") or 0) + \
                 ((v.get("exome") or {}).get("ac") or 0)
            for rs in coincide:
                filas.append({"gen": simbolo, "rsid": rs,
                              "variant_id": v["variant_id"],
                              "consecuencia": v.get("consequence") or "",
                              "hgvsp": v.get("hgvsp") or "",
                              "AC_global_gnomad": ac})
                hallados.add(rs)
        faltan = rsids - hallados
        di(f"  {len(hallados)} de {len(rsids)} recuperados"
           + (f"; sin resolver: {', '.join(sorted(faltan))}" if faltan else ""))

    if filas:
        res = pd.DataFrame(filas)
        res["AC_global_gnomad"] = pd.to_numeric(res["AC_global_gnomad"])
        res = (res.sort_values("AC_global_gnomad", ascending=False)
                  .drop_duplicates(subset="rsid", keep="first"))
        res.to_csv(SAL, sep="\t", index=False)
        di(f"\n{len(res)} coordenadas nuevas -> {SAL}")
    else:
        res = pd.DataFrame(columns=["gen", "rsid", "variant_id", "consecuencia",
                                    "hgvsp", "AC_global_gnomad"])
        di("\nNo se recupero ninguna coordenada nueva.")

    # --- cobertura final --------------------------------------------------
    todas = pd.concat([pd.read_csv(PREV, sep="\t", dtype=str), res.astype(str)]
                      if PREV.exists() else [res.astype(str)], ignore_index=True)
    todas = todas.drop_duplicates(subset="rsid")

    di("\n" + "=" * 72)
    di("COBERTURA FINAL")
    di("=" * 72)
    di(f"\n{len(todas)} de {len(snv)} loci SNV con coordenadas")
    sin = set(snv["rsid"]) - set(todas["rsid"])
    if sin:
        di(f"\nSin resolver tras el rescate ({len(sin)}): {', '.join(sorted(sin))}")
        di("Verificar en dbSNP si estan retirados o fusionados con otro rsID.")

    if MCPS.exists():
        mc = pd.read_csv(MCPS, sep="\t", low_memory=False)
        cols_af = [f"AF_{a}_{VERSION}" for a in ("IMX", "EUR", "AFR")]
        j = todas.merge(mc[["variant_id", "source", "clase", "reparto_default",
                            "n_efectivo_IMX"] + cols_af],
                        on="variant_id", how="left")
        hall = j[j["source"].notna()]
        di(f"\n{len(hall)} de {len(todas)} presentes en MCPS")
        if len(hall):
            di(f"\n{'gen':<9}{'rsid':<13}{'variante':<22}{'fuente':<7}"
               f"{'AF_IMX':>8}{'AF_EUR':>8}{'dif':>8}")
            h = hall.copy()
            h["dif"] = (h[cols_af[0]] - h[cols_af[1]]).abs()
            for _, r in h.sort_values("dif", ascending=False).iterrows():
                di(f"{r['gen']:<9}{r['rsid']:<13}{r['variant_id']:<22}"
                   f"{r['source']:<7}{r[cols_af[0]]:>8.4f}"
                   f"{r[cols_af[1]]:>8.4f}{r['dif']:>8.4f}")
        perd = j[j["source"].isna()]
        if len(perd):
            di(f"\n{len(perd)} con coordenadas pero fuera de las regiones de MCPS:")
            for _, r in perd.iterrows():
                di(f"  {r['gen']:<9}{r['rsid']:<13}{r['variant_id']}")
            di("Ampliar el padding del script 08 y repetir la extraccion.")
        j.to_csv(RES / "ding_en_mcps_final.tsv", sep="\t", index=False)

    RES.mkdir(parents=True, exist_ok=True)
    (RES / "ding_cobertura_final.txt").write_text("\n".join(salida))
    print(f"\n\nResumen -> {RES / 'ding_cobertura_final.txt'}")


if __name__ == "__main__":
    main()
