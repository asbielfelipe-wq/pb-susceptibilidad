#!/usr/bin/env python3
"""
Parche del panel 3B del script 14.

QUE CORRIGE
-----------
figura3 leia exceso_pruebas.tsv, la salida del script 12c, que calcula el
analisis restringido heredando el corte de percentil del conjunto completo. Al
recalcular los percentiles dentro de cada conjunto restringido los conteos
cambian, de modo que las barras "with control" del panel ya no corresponden a
lo que reporta el manuscrito. Este parche hace que el panel lea la rejilla
corregida.

Ademas el panel comparaba dos conjuntos y ahora son tres, porque se anadio el
criterio simetrico. Pasa a seis barras por umbral y la leyenda nombra los
conjuntos en lugar de decir "no control" y "with control", que con tres series
deja de distinguir.

Correr una sola vez: python3 scripts/parche_14_figura3B.py
Despues: python3 scripts/14_figuras.py
"""

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
F = RAIZ / "scripts" / "14_figuras.py"

VIEJO_LECTURA = '    exceso = leer("exceso_pruebas.tsv")'

NUEVO_LECTURA = '''    # La rejilla corregida vive en results/suplementario y trae los tres
    # conjuntos con los percentiles recalculados dentro de cada uno. Se
    # prefiere al exceso_pruebas.tsv del script 12c, que hereda el corte del
    # conjunto completo. Se aceptan los dos nombres posibles segun si ya se
    # renumeraron las tablas.
    exceso = None
    for nombre in ("Tabla_S3_rejilla_barrido.tsv",
                   "Tabla_S4_criterio_simetrico.tsv"):
        f = RAIZ / "results" / "suplementario" / nombre
        if not f.exists():
            continue
        t = pd.read_csv(f, sep="\\t")
        if "conjunto" not in t.columns or len(t) < 54:
            continue
        exceso = t.rename(columns={"conjunto": "analisis",
                                   "umbral_percentil": "pct"})
        print(f"  panel 3B: {nombre}, {len(exceso)} comparaciones")
        break
    if exceso is None:
        sys.exit("No se encontro la rejilla corregida de 54 comparaciones en "
                 "results/suplementario. Correr antes el script 24.")'''

VIEJO_PANEL = '''    series = [
        ("completo", "cand", C_CAND, "Candidate, no control"),
        ("completo", "ctrl", "#bdbdbd", "Control, no control"),
        ("sin_casi_privadas", "cand", C_ACENTO, "Candidate, with control"),
        ("sin_casi_privadas", "ctrl", "#e6b8b7", "Control, with control"),
    ]'''

NUEVO_PANEL = '''    series = [
        ("completo", "cand", C_CAND, "Candidate, all variants"),
        ("completo", "ctrl", "#bdbdbd", "Control, all variants"),
        ("restringido_asimetrico", "cand", C_ACENTO,
         "Candidate, asymmetric filter"),
        ("restringido_asimetrico", "ctrl", "#e6b8b7",
         "Control, asymmetric filter"),
        ("restringido_simetrico", "cand", "#4f7d6a",
         "Candidate, symmetric filter"),
        ("restringido_simetrico", "ctrl", "#b7d2c6",
         "Control, symmetric filter"),
    ]'''

CAMBIOS = [
    # seis series en lugar de cuatro: barra mas angosta y centrado nuevo
    ('    x = np.arange(len(orden)); w = 0.2',
     '    x = np.arange(len(orden)); w = 0.135'),
    ('        ax_b.bar(x + (i - 1.5) * w, vals, w, color=cl, label=lab)',
     '        ax_b.bar(x + (i - 2.5) * w, vals, w, color=cl, label=lab)'),
    ('    ax_b.legend(frameon=False, fontsize=6.3, loc="upper left")',
     '    ax_b.legend(frameon=False, fontsize=5.8, loc="upper left", ncol=2,\n'
     '                columnspacing=1.0, handlelength=1.1)'),
    ('    ax_b.set_title("Effect of the near-private allele control", loc="left")',
     '    ax_b.set_title("Effect of the near-private allele filter", loc="left")\n'
     '    ax_b.text(.5, -0.20, "whole-genome stratum, 50 kb grouping window",\n'
     '              fontsize=6.5, color="#616161", ha="center",\n'
     '              transform=ax_b.transAxes)'),
]


def main():
    if not F.exists():
        sys.exit(f"Falta {F}")
    t = F.read_text(encoding="utf-8")

    if "restringido_asimetrico" in t:
        print("El script 14 ya esta parchado. No se toca.")
        return

    faltan = [a for a, _ in CAMBIOS if a not in t]
    if VIEJO_LECTURA not in t:
        faltan.append(VIEJO_LECTURA)
    if VIEJO_PANEL not in t:
        faltan.append("el bloque de series del panel B")
    if faltan:
        print("No coincide el texto esperado; no se modifica nada:")
        for x in faltan:
            print(f"  {x[:70]}")
        sys.exit(1)

    t = t.replace(VIEJO_LECTURA, NUEVO_LECTURA)
    t = t.replace(VIEJO_PANEL, NUEVO_PANEL)
    for a, b in CAMBIOS:
        t = t.replace(a, b)

    copia = F.with_suffix(".py.antes_del_parche")
    if not copia.exists():
        copia.write_text(F.read_text(encoding="utf-8"), encoding="utf-8")
        print(f"Copia de seguridad en {copia.name}")
    F.write_text(t, encoding="utf-8")
    print(f"Parchado {F.name}")
    print("Ahora correr: python3 scripts/14_figuras.py")


if __name__ == "__main__":
    main()
