#!/usr/bin/env python3
"""
Version compacta de la Tabla 1 para el cuerpo del articulo.

La Tabla 1 completa tiene 24 columnas y es util como material consultable,
pero no es publicable en el cuerpo de un articulo: no cabe en una pagina y
obliga al lector a cruzar campos que solo necesita si va a verificar algo.
Este script produce la version de nueve columnas que responde lo que el lector
necesita saber mientras lee Resultados, y deja el resto en el suplemento.

Que conserva y por que:
  gen, rsid, REF>ALT      identifican el locus y el alelo que se reporta
  nomenclatura            el nombre con el que la literatura lo cita
  alelo_de_literatura     que alelo nombrado es ALT, cuando la literatura lo
                          nombra; es lo que decide si una frecuencia baja
                          describe el alelo relevante
  ALT_frente_a_Ding       si ALT es el alelo que el estudio de origen conto
  certeza                 con que respaldo se establecio esa correspondencia
  AF_IMX, AF_EUR          las dos frecuencias que discute el texto

Que sale al suplemento: biomarcador, desenlace, modelo genetico, direccion del
efecto, resultado de la correccion por FDR, respaldo del gen y de la variante,
transcrito, base de cada correspondencia, fuentes y notas.

Las etiquetas largas se abrevian para que la tabla entre en una pagina; la
leyenda al pie declara cada abreviatura.

Correr: python3 scripts/25b_tabla1_cuerpo.py
Salidas:
  results/suplementario/Tabla_1_cuerpo.tsv
  results/suplementario/Tabla_1_cuerpo.md
"""

import csv
import re
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
SAL = RAIZ / "results" / "suplementario"
ENTRADA = SAL / "Tabla_1_orientacion_alelica.tsv"

# Las dos columnas que pueden quedar sin resolver usan la misma abreviatura,
# "indet.", porque en ambas significa lo mismo: el dato no se pudo determinar.
VS = {"ALT es el alelo menor de Ding": "menor",
      "ALT es el alelo mayor de Ding": "mayor",
      "no determinable": "indet."}

CERTEZA = {"confirmada por alelos o por anotacion molecular": "confirmada",
           "provisional: concordancia de frecuencia, sin verificacion de hebra":
               "provisional",
           "provisional: identidad con el marcador del estudio de origen "
           "inferida, equivalencia de identificadores no documentada":
               "provisional",
           "indeterminada": "indet."}

# Nombre corto del alelo de la literatura. La columna completa explica el
# razonamiento; aqui basta el nombre.
CORTO = {
    "rs1800435": "ALAD2 (Asn59)",
    "rs699": "Thr268 (T235)",
    "rs1799945": "Asp63 (H63D)",
    "rs1800562": "Tyr282 (C282Y)",
    "rs2071747": "His7 (D7H)",
    "rs1695": "Val105 (I105V)",
    "rs429358": "Arg130",
    "rs7412": "Cys176",
    "rs2228570": "alelo F",
}

ENC = ["Gen", "rsID", "Posicion", "REF>ALT", "Nombre",
       "Alelo nombrado", "ALT en el estudio", "Certeza", "AF IMX", "AF EUR"]

PIE = (
    "Posicion, coordenada en GRCh38. REF>ALT, alelos de referencia y "
    "alterno en la hebra de referencia; todas las frecuencias corresponden al "
    "alelo alterno. Nombre, designacion con la que la literatura cita el "
    "locus. Alelo nombrado, alelo al que corresponde el alelo alterno cuando "
    "la literatura nombra uno; un guion indica que el "
    "locus se cita por sitio de restriccion o por posicion en el promotor y no "
    "por la base. ALT en el estudio, si el alelo alterno es el que el "
    "estudio de origen designa como menor o como mayor. Certeza, respaldo de "
    "esa correspondencia: confirmada, por el par de alelos o por la anotacion "
    "proteica; provisional, cuando se apoya en concordancia de frecuencia o en "
    "una identidad de identificador inferida. En las dos columnas, indet. "
    "senala que el dato no se pudo determinar, porque las lecturas posibles "
    "se contradicen. AF IMX y AF EUR, frecuencia del alelo "
    "alterno en los componentes indigena mexicano y europeo. Los dos sitios de "
    "APOE, rs429358 y rs7412, se nombran por su cambio de aminoacido: ninguno "
    "define por si solo una isoforma, porque epsilon-2, epsilon-3 y epsilon-4 "
    "son haplotipos que resultan de la combinacion de ambos. La version "
    "completa, con biomarcador, desenlace, modelo genetico, direccion del "
    "efecto, resultado de la correccion por tasa de descubrimientos falsos, "
    "respaldo del gen y de la variante, transcrito de la anotacion y fuente de "
    "cada campo, se incluye como material suplementario."
)


def frec(v):
    """Frecuencia con precision suficiente para no imprimir cero un valor que
    no lo es. Las muy pequenas van en notacion cientifica."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return ""
    if x == 0:
        return "0"
    if x < 1e-3:
        return f"{x:.2e}".replace("e-0", "e-")
    return f"{x:.4f}".rstrip("0").rstrip(".")


def main():
    if not ENTRADA.exists():
        sys.exit(f"Falta {ENTRADA}; correr antes el script 23.")
    with open(ENTRADA, encoding="utf-8") as fh:
        filas = list(csv.DictReader(fh, delimiter="\t"))

    # Las frecuencias de la Tabla 1 vienen redondeadas a cuatro decimales, de
    # modo que una del orden de 1e-6 se imprime como cero. S1 las trae
    # recalculadas desde AC y AN, asi que se prefieren esas cuando existen.
    s1f = SAL / "Tabla_S1_resultados_por_locus.tsv"
    precisas = {}
    if s1f.exists():
        with open(s1f, encoding="utf-8") as fh:
            for q in csv.DictReader(fh, delimiter="\t"):
                precisas[q["rsid"]] = (q.get("AF_IMX", ""), q.get("AF_EUR", ""))
        print(f"Frecuencias tomadas de {s1f.name}")
    else:
        print(f"Aviso: no esta {s1f.name}; se usan las frecuencias de la "
              "Tabla 1, redondeadas a cuatro decimales")

    out = []
    sin_mapear = []
    for r in filas:
        cert = CERTEZA.get(r["certeza_de_la_orientacion"])
        vs = VS.get(r["ALT_frente_al_alelo_menor_de_Ding"])
        if cert is None:
            sin_mapear.append(f"certeza: {r['certeza_de_la_orientacion'][:50]}")
        if vs is None:
            sin_mapear.append(
                f"orientacion: {r['ALT_frente_al_alelo_menor_de_Ding'][:50]}")
        out.append({
            "Gen": r["gen"],
            "rsID": r["rsid"],
            "Posicion": r["coordenada_GRCh38"].replace("chr", ""),
            "REF>ALT": f"{r['REF']}>{r['ALT']}",
            "Nombre": re.sub(r"\s*\((nomenclatura habitual|K59N)\)", "",
                             r["nomenclatura"]) or "-",
            "Alelo nombrado": CORTO.get(r["rsid"], "-"),
            "ALT en el estudio": vs or "?",
            "Certeza": cert or "?",
            "AF IMX": frec(precisas.get(r["rsid"], (None, None))[0]
                           or r["AF_IMX"]),
            "AF EUR": frec(precisas.get(r["rsid"], (None, None))[1]
                           or r["AF_EUR"]),
        })

    if sin_mapear:
        print("Valores sin abreviatura declarada; revisar el script:")
        for x in sorted(set(sin_mapear)):
            print(f"  {x}")
        sys.exit(1)

    out.sort(key=lambda f: (f["Gen"], f["rsID"]))

    with open(SAL / "Tabla_1_cuerpo.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=ENC, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    ali = ["---", "---", "---", "---", "---", "---", "---", "---",
           "---:", "---:"]
    md = ["| " + " | ".join(ENC) + " |", "| " + " | ".join(ali) + " |"]
    md += ["| " + " | ".join(f[c] for c in ENC) + " |" for f in out]
    md.append("")
    md.append(PIE)
    (SAL / "Tabla_1_cuerpo.md").write_text("\n".join(md) + "\n",
                                           encoding="utf-8")

    print("=" * 76)
    print("TABLA 1 PARA EL CUERPO DEL ARTICULO")
    print("=" * 76)
    print(f"\n{len(out)} loci, {len(ENC)} columnas "
          f"(la version completa tiene {len(filas[0])})")
    anchos = {c: max([len(c)] + [len(f[c]) for f in out]) for c in ENC}
    print("\nAncho maximo por columna, en caracteres:")
    for c in ENC:
        print(f"  {c:<28}{anchos[c]:>3}")
    print(f"\n  total                       {sum(anchos.values()):>3}")
    # Los recuentos se hacen sobre los valores que de verdad quedan escritos
    # en la tabla, que son las abreviaturas, y sobre el nombre real de la
    # columna, que es el de ENC.
    print("\nResumen de la columna de certeza:")
    for v in ("confirmada", "provisional", "indet."):
        print(f"  {v:<16}{sum(1 for f in out if f['Certeza'] == v):>3}")
    print("\nResumen de la orientacion (columna 'ALT en el estudio'):")
    for v in ("menor", "mayor", "indet."):
        print(f"  {v:<16}"
              f"{sum(1 for f in out if f['ALT en el estudio'] == v):>3}")
    print(f"\nEscrito en {SAL / 'Tabla_1_cuerpo.tsv'}")
    print(f"           {SAL / 'Tabla_1_cuerpo.md'}")


if __name__ == "__main__":
    main()
