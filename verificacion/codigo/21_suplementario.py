#!/usr/bin/env python3
"""
Ensambla el material suplementario en tres tablas numeradas y un conjunto de
archivos de verificacion.

ESTRUCTURA
----------
En el cuerpo del articulo:
  Tabla 1    correspondencia entre el alelo alterno de GRCh38 y el alelo
             evaluado en la literatura, con el nivel de respaldo

Suplementario numerado:
  Tabla S1   resultado principal por locus y percentiles bajo ambas
             definiciones de la frecuencia de apareamiento
  Tabla S2   sensibilidades por locus, en formato largo: un renglon por locus
             y escenario
  Tabla S3   rejilla completa del barrido, 36 comparaciones

Archivos de verificacion, sin numeracion de tabla:
  bloques_atipicos.tsv, genes_control.tsv, regiones.bed, procedencia.txt
  y la carpeta codigo/

La version anterior de este script producia siete tablas numeradas, varias de
ellas materiales de verificacion mas que tablas de lectura, y consolidaba las
sensibilidades concatenando marcos de columnas distintas, lo que generaba una
tabla ancha y casi vacia. Aqui las sensibilidades se reorganizan en formato
largo y los materiales de verificacion salen de la numeracion.

Correr: python3 scripts/21_suplementario.py
Salida: results/suplementario/
"""

import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
RAW = RAIZ / "data" / "raw"
SAL = RES / "suplementario"

notas = []


def di(t=""):
    print(t)
    notas.append(t)


def leer(nombre, carpeta=RES):
    f = carpeta / nombre
    return pd.read_csv(f, sep="\t") if f.exists() else None


def main():
    SAL.mkdir(parents=True, exist_ok=True)
    ver = SAL / "verificacion"
    ver.mkdir(exist_ok=True)

    di("=" * 76)
    di("MATERIAL SUPLEMENTARIO")
    di("=" * 76)
    di(f"\nGenerado {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC\n")

    # ---------------- Tabla 1 ----------------
    t1 = leer("tabla_orientacion.tsv")
    if t1 is not None:
        orden = ["gen", "rsid", "coordenada_GRCh38", "REF", "ALT", "nomenclatura",
                 "consecuencia", "hgvsp", "alelo_evaluado_en_literatura",
                 "ALT_corresponde_al_alelo_de_literatura", "direccion_del_efecto",
                 "desenlace_evaluado", "nivel_de_respaldo", "AF_IMX", "AF_EUR", "nota"]
        t1 = t1[[c for c in orden if c in t1.columns]]
        t1.to_csv(SAL / "Tabla_1_orientacion_alelica.tsv", sep="\t", index=False)
        di(f"  Tabla 1    {len(t1):>4} loci   correspondencia alelica")
    else:
        di("  Tabla 1    FALTA (correr el script 17)")

    # ---------------- Tabla S1 ----------------
    pr = leer("principal_percentiles.tsv")
    s2 = leer("sensibilidad_maf.tsv")
    if pr is not None:
        s1 = pr.copy()
        if s2 is not None and "variant_id" in s2.columns:
            alt = s2.set_index("variant_id")["pct_defA"]
            s1["pct_def_alternativa"] = s1["variant_id"].map(alt)
            s1["cambio_entre_definiciones"] = (
                s1["pct"] - s1["pct_def_alternativa"]).round(1)
        s1 = s1.rename(columns={
            "maf_conjunta": "maf_apareamiento_principal",
            "maf_alternativa": "maf_apareamiento_alternativa",
            "pct": "percentil", "IC95_inf": "IC95_inferior",
            "IC95_sup": "IC95_superior", "d_max": "distancia_max_al_fondo"})
        for c in s1.select_dtypes("float").columns:
            s1[c] = s1[c].round(4)
        s1.to_csv(SAL / "Tabla_S1_resultados_por_locus.tsv", sep="\t", index=False)
        di(f"  Tabla S1   {len(s1):>4} loci   resultado principal y ambas definiciones")
    else:
        di("  Tabla S1   FALTA (correr el script 19)")

    # ---------------- Tabla S2: sensibilidades en formato largo ----------------
    largo = []

    v = leer("sens2_versiones.tsv")
    if v is not None:
        col_rs = "rsid" if "rsid" in v.columns else v.columns[0]
        for ver_ in [c for c in ("PLANA", "UNREL", "IBD") if c in v.columns]:
            for _, r in v.iterrows():
                largo.append({"escenario": "version de frecuencia",
                              "condicion": ver_, "rsid": r[col_rs],
                              "locus": r.get("etiqueta", ""),
                              "percentil": r[ver_]})

    pl = leer("sens2_plataforma.tsv")
    if pl is not None:
        for _, r in pl.iterrows():
            largo.append({"escenario": "plataformas separadas",
                          "condicion": r.get("estrato", ""), "rsid": r["rsid"],
                          "locus": r.get("etiqueta", ""), "percentil": r.get("pct"),
                          "n_vecinos": r.get("n_vec"),
                          "IC95_inferior": r.get("IC_inf"),
                          "IC95_superior": r.get("IC_sup")})

    det = leer("sensibilidad_detalle_locus.tsv")
    if det is not None:
        for _, r in det.iterrows():
            largo.append({"escenario": "fondo local, analisis principal",
                          "condicion": r.get("estrato", ""), "rsid": r.get("rsid", ""),
                          "locus": r.get("etiqueta", ""), "percentil": r.get("pct"),
                          "n_vecinos": r.get("n_vec"), "n_bloques": r.get("n_blq"),
                          "distancia_max": r.get("d_max")})

    if largo:
        tS2 = pd.DataFrame(largo)
        cols = ["escenario", "condicion", "rsid", "locus", "percentil",
                "IC95_inferior", "IC95_superior", "n_vecinos", "n_bloques",
                "distancia_max"]
        tS2 = tS2[[c for c in cols if c in tS2.columns]]
        tS2.to_csv(SAL / "Tabla_S2_sensibilidades.tsv", sep="\t", index=False)
        di(f"  Tabla S2   {len(tS2):>4} filas  sensibilidades por locus y escenario")

    ctr = leer("sens2_controles.tsv")
    if ctr is not None:
        ctr = ctr.rename(columns={"loci_con_comparadores": "loci_evaluables",
                                  "pct95": "loci_percentil_95_o_mas"})
        ctr.to_csv(SAL / "Tabla_S2b_subconjuntos_control.tsv", sep="\t", index=False)
        di(f"  Tabla S2b  {len(ctr):>4} filas  resumen por subconjunto de control")

    # ---------------- Tabla S3 ----------------
    rj = leer("exceso_pruebas.tsv")
    if rj is not None:
        rj = rj.rename(columns={"analisis": "conjunto", "pct": "umbral_percentil",
                                "p_binomial": "p_binomial_unilateral"})
        rj["conjunto"] = rj["conjunto"].replace(
            {"completo": "completo", "sin_casi_privadas": "restringido"})
        rj = rj.sort_values(["estrato", "conjunto", "umbral_percentil", "ventana_kb"])
        rj.to_csv(SAL / "Tabla_S3_rejilla_barrido.tsv", sep="\t", index=False)
        di(f"  Tabla S3   {len(rj):>4} filas  rejilla del barrido")

    # ---------------- verificacion ----------------
    di("")
    at = leer("barrido_atipicos.tsv")
    if at is not None:
        at = at[at.get("clase", at.get("grupo")) == "candidato"] \
            if "clase" in at.columns or "grupo" in at.columns else at
        for c in at.select_dtypes("float").columns:
            at[c] = at[c].round(5)
        at.to_csv(ver / "bloques_atipicos.tsv", sep="\t", index=False)
        di(f"  verificacion/bloques_atipicos.tsv        {len(at):>5} variantes")
    for origen, destino in [((RAW / "genes_control.tsv"), "genes_control.tsv"),
                            ((RAW / "regiones.bed"), "regiones.bed")]:
        if origen.exists():
            shutil.copy(origen, ver / destino)
            n = sum(1 for _ in open(origen, encoding="utf-8", errors="ignore"))
            di(f"  verificacion/{destino:<28}{n:>5} filas")

    proc = RAW / "ctd" / "PROCEDENCIA.txt"
    lineas = ["PROCEDENCIA DE LOS DATOS", "=" * 40, "",
              "MCPS: archivos de frecuencias por cromosoma, version del 30 de",
              "enero de 2023, de https://rgc-mcps.regeneron.com/", "",
              "gnomAD v4.1: API GraphQL, conjunto gnomad_r4, GRCh38", ""]
    lineas.append(proc.read_text(encoding="utf-8", errors="ignore")
                  if proc.exists() else "CTD: https://ctdbase.org/reports")
    (ver / "procedencia.txt").write_text("\n".join(lineas))
    di(f"  verificacion/procedencia.txt")

    cod = ver / "codigo"
    cod.mkdir(exist_ok=True)
    n = 0
    for f in sorted((RAIZ / "scripts").glob("*")):
        if f.suffix in (".py", ".sh", ".R"):
            shutil.copy(f, cod / f.name)
            n += 1
    di(f"  verificacion/codigo/                     {n:>5} scripts")

    # ---------------- indice ----------------
    (SAL / "LEEME.txt").write_text("""MATERIAL SUPLEMENTARIO

En el cuerpo del articulo

  Tabla 1   Correspondencia entre el alelo alterno de GRCh38 y el alelo
            evaluado en la literatura para cada locus, con el nivel de
            respaldo de la asociacion. Las celdas sin resolver indican, segun
            el caso, orientacion alelica no determinada en la fuente, ausencia
            de una asociacion especifica para ese locus, o asociacion sin
            direccion de efecto reportada.

Suplementario

  Tabla S1  Resultado principal por locus: frecuencias por componente de
            ancestria, frecuencia de apareamiento bajo ambas definiciones,
            FST, numero de vecinos, bloques fisicos de 50 kb, distancia maxima
            al fondo, percentil con intervalo por remuestreo de bloques, y
            percentil bajo la definicion alternativa de apareamiento.

  Tabla S2  Sensibilidades por locus y escenario, en formato largo: versiones
            de frecuencia de MCPS, separacion de plataformas y detalle del
            fondo local.

  Tabla S2b Resumen por subconjunto de genes control: numero de genes, loci
            con comparadores suficientes, loci en el percentil 95 o superior,
            mediana y maximo.

  Tabla S3  Rejilla completa del barrido: tres umbrales de percentil, tres
            ventanas de agrupamiento, dos conjuntos y dos estratos.

Archivos de verificacion

  verificacion/bloques_atipicos.tsv   coordenadas y frecuencias de las
                                      variantes atipicas de genes candidatos
  verificacion/genes_control.tsv      genes control con justificacion funcional
  verificacion/regiones.bed           intervalos genomicos evaluados
  verificacion/procedencia.txt        fechas, versiones y sumas de verificacion
  verificacion/codigo/                scripts, numerados en orden de ejecucion
""")
    di(f"\n  LEEME.txt                                indice del paquete")

    di("\n" + "=" * 76)
    di("PENDIENTE ANTES DE DEPOSITAR")
    di("=" * 76)
    di("\nCompletar en la Tabla 1 las columnas sin resolver. Para los loci con")
    di("respaldo establecido, consultar la publicacion original. Para el resto,")
    di("sustituir VERIFICAR por la etiqueta que corresponda:")
    di("  orientacion no resuelta   : la fuente no indica que alelo evaluo")
    di("  sin asociacion especifica : entro al panel por plausibilidad mecanistica")
    di("  direccion no reportada    : hay asociacion pero sin direccion de efecto")

    (SAL / "resumen_generacion.txt").write_text("\n".join(notas))
    print(f"\n\nPaquete en {SAL}")


if __name__ == "__main__":
    main()
