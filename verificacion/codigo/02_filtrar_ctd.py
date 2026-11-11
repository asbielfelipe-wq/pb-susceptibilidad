#!/usr/bin/env python3
"""
Deriva de la Comparative Toxicogenomics Database una lista de genes que
interactuan con plomo, con criterios declarados y reproducibles.

POR QUE
-------
Los 13 genes analizados hasta ahora se eligieron por conocimiento de la
literatura, sin busqueda sistematica ni criterio declarado. Eso no es
defendible en un manuscrito. Este script sustituye esa eleccion por una
derivacion a partir de evidencia curada.

QUE CORRIGE RESPECTO AL USO HABITUAL DE CTD
-------------------------------------------
1. FILTRO DE ORGANISMO. CTD mezcla interacciones de humano, rata, raton, pez
   cebra y hasta plantas. Sin filtrar, la lista de genes es una quimera
   interespecie que no puede interpretarse como mecanismo humano.

2. EVIDENCIA POR PUBLICACIONES, NO POR FILAS. Una sola publicacion genera
   varias filas de interaccion (una por accion reportada). Contar filas
   sobreestima la evidencia.

3. ESPECIFICIDAD. El numero de interacciones de un gen refleja cuanto se ha
   estudiado ese gen, no su relevancia. TNF encabeza la lista de casi
   cualquier quimico porque TNF se mide en todo. Se reporta tambien el
   cociente entre interacciones con plomo e interacciones totales del gen.

4. DIRECCION DEL EFECTO. El campo InteractionActions distingue
   increases^expression de decreases^activity. Agruparlas borra la biologia.

5. DESCENDIENTES MeSH. Lead (D007854) tiene descendientes en la jerarquia
   (acetato, nitrato, cromato, tetraetilo). Incluirlos o no cambia el conteo
   por completo. Aqui se decide explicitamente y se reporta.

Correr: python3 scripts/02_filtrar_ctd.py
Salidas:
  results/ctd_genes.tsv     genes con sus metricas de evidencia
  results/ctd_resumen.txt
"""

import gzip
import sys
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
CTD = RAIZ / "data" / "raw" / "ctd"
RES = RAIZ / "results"

PLOMO = "D007854"
TAXON_HUMANO = 9606
INCLUIR_DESCENDIENTES = True   # decidido explicitamente; se reporta cuales
MIN_PMIDS = 3                  # publicaciones independientes para entrar

# Los 13 genes usados hasta ahora, para contrastar
USADOS = ["ALAD", "VDR", "HFE", "SLC11A2", "SLC30A10", "MTHFR", "GSTP1",
          "APOE", "NFE2L2", "CAT", "SOD2", "TNF", "ABCB1"]

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def leer_ctd(nombre, usecols=None):
    """Los archivos de CTD llevan la cabecera en una linea comentada.
    Se detecta la ultima linea '#' que contiene nombres de campo."""
    ruta = CTD / f"{nombre}.tsv.gz"
    if not ruta.exists():
        sys.exit(f"Falta {ruta}. Correr antes: bash scripts/01_descarga_ctd.sh")

    cabecera, n_coment = None, 0
    with gzip.open(ruta, "rt", encoding="utf-8", errors="replace") as fh:
        for linea in fh:
            if not linea.startswith("#"):
                break
            n_coment += 1
            campos = linea.lstrip("# ").rstrip("\n").split("\t")
            if len(campos) > 2:
                cabecera = [c.strip() for c in campos]
    if cabecera is None:
        sys.exit(f"No se pudo detectar la cabecera de {nombre}")

    return pd.read_csv(ruta, sep="\t", comment="#", names=cabecera,
                       usecols=usecols, low_memory=False, dtype=str)


def descendientes_de(chem_id):
    """Resuelve la jerarquia MeSH por prefijo de TreeNumbers."""
    q = leer_ctd("CTD_chemicals")
    col_id = next(c for c in q.columns if c.lower() == "chemicalid")
    col_tree = next((c for c in q.columns if "treenumber" in c.lower()
                     and "parent" not in c.lower()), None)
    col_nom = next(c for c in q.columns if c.lower() == "chemicalname")

    q[col_id] = q[col_id].str.replace("MESH:", "", regex=False)
    fila = q[q[col_id] == chem_id]
    if not len(fila) or col_tree is None:
        return {chem_id}, pd.DataFrame()

    arboles = [t for t in str(fila.iloc[0][col_tree]).split("|") if t and t != "nan"]
    if not arboles:
        return {chem_id}, pd.DataFrame()

    mask = q[col_tree].fillna("").apply(
        lambda s: any(t.startswith(a) for t in s.split("|") for a in arboles))
    desc = q[mask]
    return set(desc[col_id]) | {chem_id}, desc[[col_id, col_nom]]


def main():
    di("=" * 76)
    di("GENES QUE INTERACTUAN CON PLOMO SEGUN CTD")
    di("=" * 76)

    proc = CTD / "PROCEDENCIA.txt"
    if proc.exists():
        di("\nProcedencia de la descarga:")
        for linea in proc.read_text().splitlines()[:3]:
            di(f"  {linea}")

    # --- quimicos a incluir ----------------------------------------------
    if INCLUIR_DESCENDIENTES:
        ids, tabla = descendientes_de(PLOMO)
        di(f"\nQuimicos incluidos: {PLOMO} y {len(ids)-1} descendientes MeSH")
        if len(tabla):
            for _, r in tabla.head(15).iterrows():
                di(f"  {r.iloc[0]:<10} {r.iloc[1]}")
            if len(tabla) > 15:
                di(f"  ... y {len(tabla)-15} mas")
    else:
        ids = {PLOMO}
        di(f"\nQuimico incluido: solo {PLOMO}, sin descendientes")

    # --- interacciones ----------------------------------------------------
    di("\nLeyendo CTD_chem_gene_ixns ...")
    ix = leer_ctd("CTD_chem_gene_ixns")
    cols = {c.lower(): c for c in ix.columns}
    c_chem = cols.get("chemicalid")
    c_gen = cols.get("genesymbol")
    c_org = cols.get("organismid")
    c_act = cols.get("interactionactions")
    c_pmid = cols.get("pubmedids")
    c_forms = cols.get("geneforms")
    faltan = [n for n, c in [("ChemicalID", c_chem), ("GeneSymbol", c_gen),
                             ("OrganismID", c_org), ("PubMedIDs", c_pmid)]
              if c is None]
    if faltan:
        sys.exit(f"No se encontraron columnas: {', '.join(faltan)}. "
                 f"Columnas detectadas: {list(ix.columns)}")

    ix[c_chem] = ix[c_chem].str.replace("MESH:", "", regex=False)
    n_total = len(ix)

    pb = ix[ix[c_chem].isin(ids)].copy()
    di(f"{n_total:,} interacciones en total -> {len(pb):,} con plomo")

    pb_h = pb[pb[c_org] == str(TAXON_HUMANO)].copy()
    di(f"{len(pb_h):,} en Homo sapiens "
       f"({100*len(pb_h)/max(len(pb),1):.1f}%)")

    if c_forms:
        antes = len(pb_h)
        pb_h = pb_h[pb_h[c_forms].fillna("").str.contains("protein|mRNA", case=False)]
        di(f"{len(pb_h):,} con GeneForms de proteina o mRNA "
           f"({antes-len(pb_h):,} excluidas)")

    # --- evidencia por gen ------------------------------------------------
    def pmids(s):
        return set(str(s).split("|")) - {"nan", ""}

    agg = []
    for gen, g in pb_h.groupby(c_gen):
        p = set()
        for v in g[c_pmid]:
            p |= pmids(v)
        acciones = "|".join(g[c_act].fillna("")) if c_act else ""
        sube = acciones.count("increases")
        baja = acciones.count("decreases")
        agg.append({"gen": gen, "n_pmids": len(p), "n_filas": len(g),
                    "n_increases": sube, "n_decreases": baja})
    genes = pd.DataFrame(agg)

    # Especificidad: interacciones con plomo sobre interacciones totales
    tot = (ix[ix[c_org] == str(TAXON_HUMANO)]
           .groupby(c_gen).size().rename("n_total_quimicos"))
    genes = genes.join(tot, on="gen")
    genes["especificidad"] = (genes["n_filas"] /
                              genes["n_total_quimicos"]).round(4)
    genes = genes.sort_values("n_pmids", ascending=False)

    RES.mkdir(parents=True, exist_ok=True)
    genes.to_csv(RES / "ctd_genes.tsv", sep="\t", index=False)

    di(f"\n{len(genes):,} genes humanos con al menos una interaccion curada")
    robustos = genes[genes["n_pmids"] >= MIN_PMIDS]
    di(f"{len(robustos):,} con evidencia en >= {MIN_PMIDS} publicaciones "
       f"independientes")

    di("\n" + "=" * 76)
    di(f"20 GENES CON MAS PUBLICACIONES INDEPENDIENTES")
    di("=" * 76)
    di(f"\n{'gen':<12}{'PMIDs':>7}{'filas':>7}{'incr':>6}{'decr':>6}"
       f"{'tot quim':>10}{'especif':>9}")
    for _, r in genes.head(20).iterrows():
        di(f"{r['gen']:<12}{r['n_pmids']:>7}{r['n_filas']:>7}"
           f"{r['n_increases']:>6}{r['n_decreases']:>6}"
           f"{r['n_total_quimicos']:>10.0f}{r['especificidad']:>9.4f}")

    di("\nLa columna especificidad es la proporcion de las interacciones de ese")
    di("gen que corresponden a plomo. Un gen muy estudiado puede encabezar el")
    di("conteo bruto y tener especificidad baja.")

    # --- contraste con los 13 usados --------------------------------------
    di("\n" + "=" * 76)
    di("LOS 13 GENES ANALIZADOS, SEGUN ESTA DERIVACION")
    di("=" * 76)
    di(f"\n{'gen':<12}{'PMIDs':>7}{'rango':>8}{'percentil':>11}{'especif':>9}")
    g_idx = genes.reset_index(drop=True)
    for u in USADOS:
        f = g_idx[g_idx["gen"] == u]
        if not len(f):
            di(f"{u:<12}{'-':>7}{'ausente':>8}")
            continue
        i = f.index[0]
        pct = 100 * (1 - i / len(g_idx))
        di(f"{u:<12}{f.iloc[0]['n_pmids']:>7}{i+1:>8}{pct:>10.1f}%"
           f"{f.iloc[0]['especificidad']:>9.4f}")

    en_rob = [u for u in USADOS if u in set(robustos["gen"])]
    di(f"\n{len(en_rob)} de los 13 superan el umbral de {MIN_PMIDS} publicaciones: "
       f"{', '.join(en_rob)}")
    fuera = [u for u in USADOS if u not in en_rob]
    if fuera:
        di(f"No lo superan: {', '.join(fuera)}")

    di("\n" + "=" * 76)
    di("PARA METODOS")
    di("=" * 76)
    di(f"\nLas interacciones plomo-gen se obtuvieron de la Comparative")
    di(f"Toxicogenomics Database (ver PROCEDENCIA.txt para version y fecha)")
    di(f"para Lead (MeSH:{PLOMO}), "
       f"{'incluyendo' if INCLUIR_DESCENDIENTES else 'excluyendo'} compuestos")
    di(f"descendientes. Se restringio a Homo sapiens (NCBI taxon {TAXON_HUMANO})")
    di(f"y a interacciones con evidencia en al menos {MIN_PMIDS} publicaciones")
    di(f"independientes, resultando en {len(robustos):,} genes.")

    (RES / "ctd_resumen.txt").write_text("\n".join(salida))
    print(f"\n\nGenes   -> {RES / 'ctd_genes.tsv'}")
    print(f"Resumen -> {RES / 'ctd_resumen.txt'}")


if __name__ == "__main__":
    main()
