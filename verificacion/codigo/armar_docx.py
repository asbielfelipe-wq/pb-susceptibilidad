#!/usr/bin/env python3
"""
Convierte manuscrito.md a Word y corrige el formato: todo el texto en negro,
sin lineas decorativas, y tablas que no se parten entre renglones ni entre
paginas.

El problema que resuelve en las tablas: pandoc entrega las tablas con ancho
automatico, de modo que Word reparte las columnas segun el contenido y las
celdas numericas estrechas acaban cortando los valores a mitad ("0.16" y "4"
en renglones distintos). Aqui se fija el ancho de cada columna en funcion de
la longitud real de su contenido, se desactiva el autoajuste, se repite el
renglon de encabezado en cada pagina y se impide que un renglon se divida.

Correr: python3 armar_docx.py
"""

import copy
import re
import subprocess
import sys
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Pt, RGBColor, Cm

AQUI = Path(__file__).resolve().parent
MD = AQUI / "manuscrito.md"
DOCX = AQUI / "Manuscrito_plomo_IMX.docx"

ANCHO_UTIL_CM = 16.0      # pagina carta con margenes de 2.5 cm
CM_POR_CARACTER = 0.168   # ancho medio de caracter a 8.5 pt
CM_MINIMO_CARACTER = 0.19 # holgura para la palabra que no se puede partir
RELLENO_CM = 0.22         # margen interno de celda

# Simbolos genicos que van en cursiva segun el estilo HGNC. La lista es
# explicita y no un patron de mayusculas, para no poner en cursiva siglas que
# no son genes (HLA designa una region, no un gen, y no se italiza; CAT y TF
# coincidirian con palabras corrientes y no aparecen aqui como genes).
GENES = ["AGTR1", "SLC11A2", "NFE2L2", "ABCB1", "HMOX1", "GSTP1", "GSTT1",
         "GSTM1", "MTHFR", "CYP1A1", "PEPT2", "ALAD", "APOE", "SOD2", "AGT",
         "VDR", "HFE", "REN", "TNF"]
RE_GENES = re.compile(r"(?<![A-Za-z0-9_-])(" + "|".join(GENES)
                      + r")(?![A-Za-z0-9_-])")


def convertir():
    subprocess.run(["pandoc", str(MD), "-o", str(DOCX)], check=True)


def _parrafos_de(contenedor):
    """Todos los parrafos, incluidos los de dentro de las tablas."""
    yield from contenedor.paragraphs
    for t in contenedor.tables:
        for fila in t.rows:
            for celda in fila.cells:
                yield from _parrafos_de(celda)


def cursiva_en_genes(doc):
    """Pone en cursiva los simbolos genicos, segun el estilo HGNC.

    Un run de Word lleva un solo formato, de modo que para italizar una palabra
    dentro de un run hay que partirlo. Se clona el run original para conservar
    su formato (negrita, tamano, color) y solo se cambia la cursiva en los
    fragmentos que son simbolos genicos.
    """
    n = 0
    for p in _parrafos_de(doc):
        for run in list(p.runs):
            texto = run.text
            if not texto or not RE_GENES.search(texto):
                continue
            partes = [x for x in RE_GENES.split(texto) if x != ""]
            if len(partes) < 2:
                # El run entero es un simbolo genico y no hay nada que partir.
                # Es el caso de las celdas de la columna "Gen" de la Tabla 1,
                # que quedaban en redonda porque el bucle las descartaba aqui.
                if RE_GENES.fullmatch(texto) and not run.italic:
                    run.italic = True
                    n += 1
                continue
            run.text = partes[0]
            run.italic = bool(RE_GENES.fullmatch(partes[0]))
            if run.italic:
                n += 1
            anterior = run._r
            for parte in partes[1:]:
                nuevo = copy.deepcopy(run._r)
                for hijo in nuevo.findall(qn("w:t")):
                    nuevo.remove(hijo)
                el = OxmlElement("w:t")
                el.set(qn("xml:space"), "preserve")
                el.text = parte
                nuevo.append(el)
                anterior.addnext(nuevo)
                anterior = nuevo
            # reasignar cursiva en los runs recien creados
            for run2 in p.runs:
                if RE_GENES.fullmatch(run2.text or ""):
                    if not run2.italic:
                        run2.italic = True
                        n += 1
                elif run2.italic and not RE_GENES.search(run2.text or ""):
                    run2.italic = False
    return n


def negro(doc):
    """Todo el texto en negro, incluidos encabezados y celdas."""
    def runs_de(contenedor):
        for p in contenedor.paragraphs:
            yield from p.runs
        for t in contenedor.tables:
            for fila in t.rows:
                for celda in fila.cells:
                    yield from runs_de(celda)

    n = 0
    for r in runs_de(doc):
        r.font.color.rgb = RGBColor(0, 0, 0)
        n += 1
    for s in doc.sections:
        for parte in (s.header, s.footer):
            for p in parte.paragraphs:
                for r in p.runs:
                    r.font.color.rgb = RGBColor(0, 0, 0)
    return n


def sin_bordes_decorativos(doc):
    """Quita bordes de parrafo, que es lo que pandoc usa como linea separadora."""
    n = 0
    for p in doc.paragraphs:
        pPr = p._p.get_or_add_pPr()
        for b in pPr.findall(qn("w:pBdr")):
            pPr.remove(b)
            n += 1
    return n


def _ancho_celda(celda, cm):
    celda.width = Cm(cm)
    for p in celda.paragraphs:
        p.paragraph_format.space_before = Pt(1)
        p.paragraph_format.space_after = Pt(1)
        for r in p.runs:
            r.font.size = Pt(8.5)


def _no_partir(fila, encabezado=False):
    trPr = fila._tr.get_or_add_trPr()
    cant = OxmlElement("w:cantSplit")
    trPr.append(cant)
    if encabezado:
        th = OxmlElement("w:tblHeader")
        trPr.append(th)


def ajustar_tablas(doc):
    """Ancho por columna proporcional al contenido, sin autoajuste."""
    for t in doc.tables:
        n_col = len(t.columns)
        if not n_col:
            continue
        # Dos medidas por columna: la palabra mas larga, que no se puede
        # partir y fija el ancho minimo, y la celda mas larga, que indica
        # cuanto espacio conviene darle del sobrante.
        minimos, deseados = [], []
        for j in range(n_col):
            textos = [t.cell(i, j).text.strip() for i in range(len(t.rows))]
            palabra = max((len(w) for txt in textos for w in txt.split()),
                          default=4)
            celda = max((len(txt) for txt in textos), default=4)
            minimos.append(palabra * CM_MINIMO_CARACTER + RELLENO_CM)
            deseados.append(celda * CM_POR_CARACTER + RELLENO_CM)

        if sum(minimos) >= ANCHO_UTIL_CM:
            f = ANCHO_UTIL_CM / sum(minimos)
            anchos = [m * f for m in minimos]
        else:
            sobrante = ANCHO_UTIL_CM - sum(minimos)
            extra = [max(0.0, d - m) for d, m in zip(deseados, minimos)]
            te = sum(extra)
            anchos = [m + (sobrante * e / te if te else sobrante / n_col)
                      for m, e in zip(minimos, extra)]

        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        t.autofit = False
        tblPr = t._tbl.tblPr
        for el in tblPr.findall(qn("w:tblLayout")):
            tblPr.remove(el)
        layout = OxmlElement("w:tblLayout")
        layout.set(qn("w:type"), "fixed")
        tblPr.append(layout)

        for i, fila in enumerate(t.rows):
            _no_partir(fila, encabezado=(i == 0))
            for j, celda in enumerate(fila.cells):
                if j < len(anchos):
                    _ancho_celda(celda, anchos[j])
        for j, col in enumerate(t.columns):
            if j < len(anchos):
                col.width = Cm(anchos[j])
    return len(doc.tables)


def main():
    if not MD.exists():
        sys.exit(f"No existe {MD}")
    convertir()
    doc = Document(DOCX)
    n_genes = cursiva_en_genes(doc)
    n_runs = negro(doc)
    n_bordes = sin_bordes_decorativos(doc)
    n_tablas = ajustar_tablas(doc)
    doc.save(DOCX)
    print(f"{DOCX.name}")
    print(f"  {n_genes} simbolos genicos en cursiva")
    print(f"  {n_runs} fragmentos de texto en negro")
    print(f"  {n_bordes} bordes de parrafo retirados")
    print(f"  {n_tablas} tabla(s) con ancho fijo, encabezado repetido y "
          f"renglones indivisibles")


if __name__ == "__main__":
    main()
