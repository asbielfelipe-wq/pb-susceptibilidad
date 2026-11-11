#!/usr/bin/env python3
"""
Construye la version final de la Tabla 1: correspondencia entre el alelo
alterno de GRCh38 y el alelo evaluado en el estudio de origen, con la base de
esa correspondencia, la direccion del efecto en su contexto completo y la
fuente de cada campo.

COMO SE RESUELVE LA CORRESPONDENCIA
-----------------------------------
Ding et al. 2016 reportan en su Tabla 1 el alelo de menor frecuencia de cada
SNP y en su Tabla 3 los genotipos nucleotidicos con la razon de riesgo de cada
uno. Eso permite cerrar la correspondencia sin recurrir a las letras de las
enzimas de restriccion, que nombran la presencia o ausencia del sitio y no la
base. Se usan tres vias, en este orden:

  1. Anotacion proteica. Si la variante es de cambio de aminoacido y la
     literatura nombra el alelo por el residuo, el residuo identifica el alelo
     sin ambiguedad: ALT de rs1800562 produce p.Cys282Tyr, de modo que ALT es
     el alelo Tyr282, que es C282Y. Esta via resuelve tambien las
     transversiones C/G de rs1799945 y rs2071747, donde la hebra por si sola
     no bastaria.

  2. Par de alelos no ambiguo. Si el par REF/ALT no coincide con su propio
     complemento (es decir, no es A/T ni C/G), el par de Ding identifica la
     hebra de forma unica y la correspondencia queda determinada.

  3. Transversiones A/T y C/G no codificantes. El par coincide con su
     complemento y la hebra no es verificable a partir de las fuentes. Cuando
     la frecuencia del alelo alterno en el componente europeo concuerda con la
     designacion de alelo menor de Ding, se declara concordancia por
     frecuencia, que es un argumento de consistencia y no una verificacion.
     Cuando no concuerda, la correspondencia se declara no determinable.

La direccion del efecto se reporta con su contexto completo, porque fuera de el
no significa nada: el biomarcador de exposicion (plomo en rotula o en tibia),
el desenlace (cardiopatia coronaria incidente), el modelo genetico y si la
interaccion alcanzo el umbral de la correccion por tasa de
descubrimientos falsos. Ninguna lo alcanza: cuatro quedaron en 0.06.

El respaldo se separa en dos columnas, porque el respaldo de un gen no es el
respaldo de una variante concreta: un gen puede tener mecanismo establecido en
la toxicocinetica del plomo mientras la variante especifica no se ha evaluado
nunca contra un desenlace de plomo fuera del panel de origen.

Correr: python3 scripts/23b_tabla1_orientacion.py
Salida: results/suplementario/Tabla_1_orientacion_alelica.tsv
"""

import csv
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
SAL = RES / "suplementario"
ENTRADA = RES / "tabla_orientacion.tsv"
SALIDA = SAL / "Tabla_1_orientacion_alelica.tsv"

DOI_DING = "10.1371/journal.pone.0161472"
F_DING_T1 = f"Ding et al. 2016, Tabla 1 (doi:{DOI_DING})"
F_DING_T3 = f"Ding et al. 2016, Tabla 3 (doi:{DOI_DING})"
F_DING_TXT = f"Ding et al. 2016, Resultados (doi:{DOI_DING})"
F_GNOMAD = "gnomAD v4.1, transcrito canonico"

VIA_PROT = "anotacion proteica del transcrito canonico"
VIA_PAR = "par de alelos no ambiguo frente a los genotipos de Ding"
VIA_FREC = "transversion A/T o C/G: hebra no verificable, concordante por frecuencia"
VIA_NO = "transversion A/T: hebra no verificable y lecturas discordantes"
VIA_INFER = ("identidad del identificador del estudio de origen inferida, no "
             "documentada")

MENOR = "ALT es el alelo menor de Ding"
MAYOR = "ALT es el alelo mayor de Ding"
INDET = "no determinable"

GEN_EST = "establecido en toxicocinetica del plomo"
GEN_FRA = "mecanistico, evidencia fragmentaria"

VAR_FDR = ("interaccion limitrofe en el estudio de origen: FDR ajustada = "
           "0.06, que no alcanza el umbral convencional de 0.05")
VAR_SIN = "evaluada en el estudio de origen sin modificacion significativa"
VAR_NO = "no evaluada contra un desenlace de plomo fuera del estudio de origen"

ROTULA = ("plomo en rotula, cardiopatia coronaria incidente; el HR es por "
          "duplicacion del plomo oseo dentro de cada estrato genotipico y no "
          "el efecto de portar el alelo")
TIBIA = ("plomo en tibia, cardiopatia coronaria incidente; el HR es por "
         "duplicacion del plomo oseo dentro de cada estrato genotipico y no "
         "el efecto de portar el alelo")
ADITIVO = "aditivo, genotipo como variable continua"
DOMINANTE = "dominante, portadores del alelo menor frente a no portadores"

# Curaduria por variante. Cada entrada declara de donde sale cada campo.
# alelo      : que alelo de la literatura es ALT, o por que no aplica
# via        : como se establecio
# vs_ding    : ALT frente al alelo menor de Ding
# biomarc    : biomarcador y desenlace del que proviene la direccion
# modelo     : modelo genetico
# direccion  : la direccion con sus cifras
# fdr        : resultado tras la correccion, con el valor ajustado
# gen_resp   : respaldo del gen
# var_resp   : respaldo de la variante
# fuente     : fuentes separadas por "; "
C = {
 "rs1800435": dict(
    alelo="ALAD2 (Asn59; C en la hebra del gen, G en la hebra de referencia)",
    via=VIA_PROT, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="sin modificacion de efecto detectada (p = 0.84)",
    fdr="no aplicable: la interaccion no fue significativa antes de la correccion",
    gen_resp=GEN_EST,
    var_resp="evaluada en el estudio de origen sin modificacion significativa; "
             "metaanalisis posterior de 27 estudios sin asociacion con plombemia",
    fuente=f"{F_GNOMAD}; {F_DING_T1}; {F_DING_T3}; "
           "Jiang y Xie 2025 (doi:10.1038/s41598-025-14850-y)",
    nota_propia=(
        "La Tabla 3 de la fuente imprime este identificador como rs1833435, "
        "mientras que su seccion de genotipificacion y su Tabla 1 lo dan como "
        "rs1800435; se anota para localizar la fila correspondiente"),
    ),

 "rs1544410": dict(
    alelo="no aplicable: intronica, la literatura la nombra por el sitio BsmI "
          "y no por la base",
    via=VIA_PAR, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="la asociacion se refuerza con la dosis de ALT: HR 0.97 sin ALT, "
              "1.59 en heterocigotos, 1.77 en homocigotos ALT (p = 0.004)",
    fdr="no; FDR ajustada = 0.06, sin alcanzar el umbral de 0.05", gen_resp=GEN_EST, var_resp=VAR_FDR,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs731236": dict(
    alelo="no aplicable: sinonima, la literatura la nombra por el sitio TaqI",
    via=VIA_PAR, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="la asociacion se refuerza con la dosis de ALT: HR 1.00 sin ALT, "
              "1.56 en heterocigotos, 1.70 en homocigotos ALT (p = 0.009)",
    fdr="no; FDR ajustada = 0.06, sin alcanzar el umbral de 0.05", gen_resp=GEN_EST, var_resp=VAR_FDR,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs7975232": dict(
    alelo="no aplicable: intronica, la literatura la nombra por el sitio ApaI",
    via=VIA_PAR, vs_ding=MAYOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="relacion no monotona: HR 1.48 en homocigotos ALT, 1.71 en "
              "heterocigotos y 0.84 en homocigotos REF (p = 0.007)",
    fdr="no; FDR ajustada = 0.06, sin alcanzar el umbral de 0.05", gen_resp=GEN_EST, var_resp=VAR_FDR,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs2228570": dict(
    alelo="alelo F, el que carece del primer ATG y da la proteina de 424 "
          "residuos; ALT produce p.Met1? y su frecuencia europea de 0.65 "
          "coincide con la descrita para F",
    via=VIA_INFER, vs_ding=MAYOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="sin modificacion de efecto detectada (p = 0.09)",
    fdr="no aplicable: la interaccion no fue significativa antes de la "
        "correccion",
    gen_resp=GEN_EST, var_resp=VAR_SIN,
    fuente=f"{F_GNOMAD}; dbSNP build 157; {F_DING_T1}; {F_DING_T3}",
    nota_propia=(
        "FokI figura en el estudio de origen y en publicaciones previas con el "
        "identificador rs1073581. En dbSNP build 157 FokI corresponde a "
        "rs2228570, en chr12:47,879,112 (GRCh38), identificador en el que se "
        "fusiono rs10735810, y es el que se analizo aqui. Los alelos A y G de "
        "rs2228570 en la hebra de referencia coinciden con los genotipos G/G, "
        "G/A y A/A que reporta ese estudio, y la frecuencia europea del alelo "
        "alterno, 0.652, concuerda con la designacion de A como alelo menor, "
        "de modo que ALT corresponde al alelo mayor de ese estudio y al alelo "
        "F de la nomenclatura clasica. Esa correspondencia descansa en la "
        "concordancia de alelos y frecuencias, no en una equivalencia de "
        "identificadores registrada en dbSNP, y por eso la orientacion se "
        "reporta como provisional"),
    ),

 "rs757343": dict(
    alelo="no aplicable: intronica, la literatura la nombra por el sitio Tru9I",
    via=VIA_PAR, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="la asociacion se refuerza con la dosis de ALT: HR 1.17 sin ALT, "
              "1.30 en heterocigotos, 2.09 en homocigotos ALT (p = 0.07)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs1799945": dict(
    alelo="Asp63, que es H63D", via=VIA_PROT, vs_ding=MENOR,
    biomarc=ROTULA, modelo=ADITIVO,
    direccion="sin modificacion de efecto detectada (p = 0.25)",
    fdr="no aplicable: la interaccion no fue significativa antes de la correccion",
    gen_resp=GEN_EST, var_resp=VAR_SIN,
    fuente=f"{F_GNOMAD}; {F_DING_T1}; {F_DING_T3}"),

 "rs1800562": dict(
    alelo="Tyr282, que es C282Y", via=VIA_PROT, vs_ding=MENOR,
    biomarc=ROTULA, modelo=ADITIVO,
    direccion="sin modificacion de efecto detectada (p = 0.22); no hubo "
              "homocigotos ALT suficientes para estimar ese estrato",
    fdr="no aplicable: la interaccion no fue significativa antes de la correccion",
    gen_resp=GEN_EST, var_resp=VAR_SIN,
    fuente=f"{F_GNOMAD}; {F_DING_T1}; {F_DING_T3}"),

 "rs2071746": dict(
    alelo="no aplicable: de promotor, sin alelo nombrado en la literatura",
    via=VIA_FREC, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="la asociacion se atenua con la dosis de ALT: HR 1.51 sin ALT, "
              "1.44 en heterocigotos, 1.02 en homocigotos ALT (p = 0.09)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs2071747": dict(
    alelo="His7, que es D7H", via=VIA_PROT, vs_ding=MENOR,
    biomarc=ROTULA, modelo=ADITIVO,
    direccion="sin modificacion de efecto detectada (p = 0.82); no hubo "
              "homocigotos ALT suficientes para estimar ese estrato",
    fdr="no aplicable: la interaccion no fue significativa antes de la correccion",
    gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_GNOMAD}; {F_DING_T1}; {F_DING_T3}"),

 "rs2071749": dict(
    alelo="no aplicable: intronica, sin alelo nombrado en la literatura",
    via=VIA_PAR, vs_ding=MAYOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="relacion no monotona: HR 1.02 en homocigotos ALT, 1.54 en "
              "heterocigotos y 1.34 en homocigotos REF (p = 0.10)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs5995098": dict(
    alelo="no aplicable: intronica, sin alelo nombrado en la literatura",
    via=VIA_FREC, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="relacion no monotona: HR 1.62 en homocigotos REF, 1.18 en "
              "heterocigotos y 1.40 en homocigotos ALT (p = 0.15)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs429358": dict(
    alelo="el alelo Arg130 que, junto con el alelo C de rs7412, define e4; no "
          "es e4 por si solo",
    via=VIA_PROT, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="la asociacion se atenua con la dosis de ALT: HR 1.43 sin ALT, "
              "1.05 en heterocigotos, 0.87 en homocigotos ALT (p = 0.10)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_GNOMAD}; {F_DING_T1}; {F_DING_T3}"),

 "rs7412": dict(
    alelo="el alelo Cys176 que, junto con el alelo T de rs429358, define e2; "
          "no es e2 por si solo",
    via=VIA_PROT, vs_ding=MENOR, biomarc=TIBIA, modelo=DOMINANTE,
    direccion="HR 2.24 (IC95 1.45-3.44) entre portadores de al menos un alelo "
              "ALT; en rotula no hubo modificacion (p = 0.30)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_GNOMAD}; {F_DING_T1}; {F_DING_T3}; {F_DING_TXT}"),

 "rs440446": dict(
    alelo="no aplicable: intronica, sin alelo nombrado en la literatura",
    via=VIA_FREC, vs_ding=MAYOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="sin modificacion de efecto detectada (p = 0.64)",
    fdr="no aplicable: la interaccion no fue significativa antes de la correccion",
    gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs405509": dict(
    alelo="no aplicable: de promotor, sin alelo nombrado en la literatura",
    via=VIA_PAR, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="la asociacion se refuerza con la dosis de ALT: HR 1.22 sin ALT, "
              "1.26 en heterocigotos, 1.58 en homocigotos ALT (p = 0.33)",
    fdr="no aplicable: la interaccion no fue significativa antes de la correccion",
    gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs449647": dict(
    alelo="no aplicable: de promotor, sin alelo nombrado en la literatura",
    via=VIA_NO, vs_ding=INDET, biomarc=TIBIA, modelo=DOMINANTE,
    direccion="HR 1.70 (IC95 1.36-2.12) entre no portadores del alelo menor de "
              "Ding; como la correspondencia con ALT no es determinable, la "
              "direccion no puede asignarse a ALT ni a REF",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}; {F_DING_TXT}"),

 "rs769446": dict(
    alelo="no aplicable: de promotor, sin alelo nombrado en la literatura",
    via=VIA_PAR, vs_ding=MENOR, biomarc=TIBIA, modelo=DOMINANTE,
    direccion="HR 2.81 (IC95 1.51-5.23) entre portadores de al menos un alelo "
              "ALT; en rotula no hubo modificacion (p = 0.46)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}; {F_DING_TXT}"),

 "rs1695": dict(
    alelo="Val105, que es I105V", via=VIA_PROT, vs_ding=MENOR,
    biomarc=TIBIA, modelo=DOMINANTE,
    direccion="HR 2.04 (IC95 1.47-2.84) entre portadores de al menos un alelo "
              "ALT; en rotula no hubo modificacion (p = 0.31)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_GNOMAD}; {F_DING_T1}; {F_DING_T3}; {F_DING_TXT}"),

 "rs12695908": dict(
    alelo="no aplicable: intronica, sin alelo nombrado en la literatura",
    via=VIA_PAR, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="sin modificacion de efecto detectada (p = 0.72); no hubo "
              "homocigotos ALT suficientes para estimar ese estrato",
    fdr="no aplicable: la interaccion no fue significativa antes de la correccion",
    gen_resp=GEN_FRA, var_resp=VAR_NO,
    fuente=f"{F_DING_T3}"),

 "rs699": dict(
    alelo="Thr268, que la literatura llama T235 o M268T segun la numeracion",
    via=VIA_PROT, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="la asociacion se atenua con la dosis de ALT: HR 2.17 en "
              "homocigotos REF, 1.29 en heterocigotos y 1.18 en homocigotos "
              "ALT (p = 0.02); el alelo mas diferenciado aqui es el que menos "
              "reforzaba la asociacion alli",
    fdr="no; FDR ajustada = 0.06, sin alcanzar el umbral de 0.05", gen_resp=GEN_FRA, var_resp=VAR_FDR,
    fuente=f"{F_GNOMAD}; {F_DING_T1}; {F_DING_T3}"),

 "rs5046": dict(
    alelo="no aplicable: de promotor, sin alelo nombrado en la literatura",
    via=VIA_PAR, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="la asociacion se atenua con la dosis de ALT: HR 1.57 sin ALT, "
              "1.13 en heterocigotos, 1.04 en homocigotos ALT (p = 0.05)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs5050": dict(
    alelo="no aplicable: de promotor, sin alelo nombrado en la literatura",
    via=VIA_PAR, vs_ding=MENOR, biomarc=ROTULA, modelo=ADITIVO,
    direccion="sin modificacion de efecto detectada (p = 0.75); no hubo "
              "homocigotos ALT suficientes para estimar ese estrato",
    fdr="no aplicable: la interaccion no fue significativa antes de la correccion",
    gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}"),

 "rs2493137": dict(
    alelo="no aplicable: de promotor, sin alelo nombrado en la literatura",
    via=VIA_PAR, vs_ding=MENOR, biomarc=TIBIA, modelo=DOMINANTE,
    direccion="HR 2.06 (IC95 1.48-2.87) entre no portadores de ALT, con "
              "atenuacion en portadores de uno o dos alelos ALT; en rotula no "
              "hubo modificacion (p = 0.25)",
    fdr="no", gen_resp=GEN_FRA, var_resp=VAR_SIN,
    fuente=f"{F_DING_T1}; {F_DING_T3}; {F_DING_TXT}"),
}

# Certeza de la orientacion. Separar estas tres categorias es necesario
# porque la concordancia de frecuencia es evidencia indirecta: sugiere que la
# hebra coincide, no lo demuestra.
CONF = "confirmada por alelos o por anotacion molecular"
# Hay dos motivos distintos para no dar una orientacion por confirmada y
# conviene no confundirlos: en un caso la hebra no es verificable, en el otro
# la hebra no esta en duda y lo que se infiere es a que polimorfismo designa
# el identificador del estudio de origen.
PROV_HEBRA = "provisional: concordancia de frecuencia, sin verificacion de hebra"
PROV_IDENT = ("provisional: identidad con el marcador del estudio de origen "
              "inferida, equivalencia de identificadores no documentada")
INDE = "indeterminada"

ORDEN = ["gen", "rsid", "coordenada_GRCh38", "REF", "ALT", "nomenclatura",
         "consecuencia", "hgvsp", "transcrito_de_la_anotacion",
         "alelo_de_literatura_que_corresponde_a_ALT", "base_del_alelo_nombrado",
         "ALT_frente_al_alelo_menor_de_Ding",
         "base_de_la_orientacion_frente_a_Ding", "certeza_de_la_orientacion",
         "biomarcador_y_desenlace", "modelo_genetico",
         "direccion_en_el_estudio_de_origen", "significacion_tras_FDR",
         "respaldo_del_gen", "respaldo_de_la_variante", "fuente",
         "AF_IMX", "AF_EUR", "nota"]

AMBIGUOS = ({"A", "T"}, {"C", "G"})


NOMENCLATURA_APOE = {"rs429358": "C130R", "rs7412": "R176C"}


def derivar_orientacion(ref, alt, via):
    """Base y certeza de la orientacion frente al estudio de origen.

    Un par A/T o C/G coincide con su propio complemento, de modo que la hebra
    no se deduce del par. Si la variante es de cambio de aminoacido, el residuo
    rompe esa ambiguedad; si no, solo queda la concordancia de frecuencia, que
    es indirecta.
    """
    ambiguo = {ref, alt} in AMBIGUOS
    if via == VIA_NO:
        return ("transversion ambigua con lecturas discordantes entre si", INDE)
    if via == VIA_FREC:
        return ("concordancia de frecuencia en transversion ambigua",
                PROV_HEBRA)
    if via == VIA_INFER:
        return ("identidad del identificador del estudio de origen inferida del "
                "par de alelos y de la concordancia de frecuencia, sin fusion "
                "registrada en dbSNP", PROV_IDENT)
    if via == VIA_PROT and ambiguo:
        return ("anotacion proteica, que resuelve la transversion ambigua", CONF)
    if via == VIA_PROT:
        return ("par de alelos no ambiguo, concordante con la anotacion "
                "proteica", CONF)
    return ("par de alelos no ambiguo", CONF)


def main():
    if not ENTRADA.exists():
        sys.exit(f"Falta {ENTRADA}; correr primero el script 17.")
    SAL.mkdir(parents=True, exist_ok=True)

    with open(ENTRADA, encoding="utf-8") as fh:
        base = list(csv.DictReader(fh, delimiter="\t"))

    faltan = [r["rsid"] for r in base if r["rsid"] not in C]
    sobran = [k for k in C if k not in {r["rsid"] for r in base}]
    if faltan:
        sys.exit(f"Sin curaduria para: {', '.join(faltan)}")
    if sobran:
        sys.exit(f"Curaduria sin locus correspondiente: {', '.join(sobran)}")

    # Las frecuencias del script 17 vienen redondeadas a cuatro decimales, de
    # modo que una del orden de 1e-6 se imprime como cero. S1 las trae
    # recalculadas desde AC y AN; se prefieren esas para que los tres archivos
    # que publican la frecuencia de un mismo locus digan lo mismo.
    s1f = SAL / "Tabla_S1_resultados_por_locus.tsv"
    precisas = {}
    if s1f.exists():
        with open(s1f, encoding="utf-8") as fh:
            for q in csv.DictReader(fh, delimiter="\t"):
                precisas[q["rsid"]] = (q.get("AF_IMX", ""), q.get("AF_EUR", ""))
        print(f"Frecuencias tomadas de {s1f.name}")
    else:
        print(f"Aviso: no esta {s1f.name}; correr antes el 22b. Se usan las "
              "frecuencias del script 17, redondeadas a cuatro decimales, y "
              "las menores de 1e-4 apareceran como cero.")

    filas = []
    for r in base:
        c = C[r["rsid"]]
        fila = {k: r.get(k, "") for k in
                ("gen", "rsid", "coordenada_GRCh38", "REF", "ALT",
                 "nomenclatura", "consecuencia", "hgvsp", "AF_IMX", "AF_EUR")}
        # Ni rs429358 ni rs7412 definen por si solos una isoforma: epsilon-2,
        # epsilon-3 y epsilon-4 son haplotipos de los dos sitios juntos. La
        # tabla los nombra por el cambio de aminoacido y la leyenda explica la
        # combinacion.
        if r["rsid"] in NOMENCLATURA_APOE:
            fila["nomenclatura"] = NOMENCLATURA_APOE[r["rsid"]]
        if r["rsid"] in precisas:
            af_i, af_e = precisas[r["rsid"]]
            if af_i:
                fila["AF_IMX"] = af_i
            if af_e:
                fila["AF_EUR"] = af_e
        base_or, certeza = derivar_orientacion(r["REF"], r["ALT"], c["via"])
        fila.update({
            "transcrito_de_la_anotacion": (
                "transcrito canonico de gnomAD v4.1"
                if (r.get("hgvsp") or "").strip() else ""),
            "alelo_de_literatura_que_corresponde_a_ALT": c["alelo"],
            "base_del_alelo_nombrado": (
                VIA_PROT if c["via"] in (VIA_PROT, VIA_INFER)
                else "no aplicable: la literatura no nombra un alelo"),
            "ALT_frente_al_alelo_menor_de_Ding": c["vs_ding"],
            "base_de_la_orientacion_frente_a_Ding": base_or,
            "certeza_de_la_orientacion": certeza,
            "biomarcador_y_desenlace": c["biomarc"],
            "modelo_genetico": c["modelo"],
            "direccion_en_el_estudio_de_origen": c["direccion"],
            "significacion_tras_FDR": c["fdr"],
            "respaldo_del_gen": c["gen_resp"],
            "respaldo_de_la_variante": c["var_resp"],
            "fuente": c["fuente"],
            "nota": "; ".join(x for x in (r.get("nota", ""),
                                          c.get("nota_propia", "")) if x),
        })
        # la anotacion proteica de rs440446 proviene de un transcrito
        # solapante y no describe esta variante; se retira para no inducir error
        if r["rsid"] == "rs440446":
            # Intronica en el transcrito MANE Select (c.-24+69C>G). El cambio
            # p.Asn14Lys solo aparece en NM_001302688.2, que no es MANE, de
            # modo que no describe esta variante en el transcrito de
            # referencia y se omite de la columna de consecuencia proteica.
            fila["hgvsp"] = ""
            fila["consecuencia"] = "intron_variant (MANE Select)"
            fila["transcrito_de_la_anotacion"] = (
                "NM_000041.4 / ENST00000252486.9 (MANE Select): c.-24+69C>G, "
                "intronica. El cambio p.Asn14Lys corresponde a "
                "NM_001302688.2, que no es MANE")
            fila["nota"] = ("Intronica en el transcrito MANE Select; la "
                            "consecuencia proteica de un transcrito "
                            "alternativo no describe esta variante y se omite")
        filas.append(fila)

    filas.sort(key=lambda f: (f["gen"], f["rsid"]))
    with open(SALIDA, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=ORDEN, delimiter="\t")
        w.writeheader()
        w.writerows(filas)

    def cuenta(col, val):
        return sum(1 for f in filas if f[col] == val)

    print("=" * 76)
    print("TABLA 1: ORIENTACION ALELICA")
    print("=" * 76)
    print(f"\n{len(filas)} loci\n")
    print("Correspondencia de ALT con el alelo de la literatura:")
    nombrado = [f for f in filas
                if not f["alelo_de_literatura_que_corresponde_a_ALT"]
                .startswith("no aplicable")]
    print(f"  alelo nombrado e identificado      {len(nombrado):>3}")
    print(f"  sin alelo nombrado en la literatura {len(filas)-len(nombrado):>3}")
    print("\nBase del alelo nombrado:")
    print(f"  anotacion proteica                "
          f"{cuenta('base_del_alelo_nombrado', VIA_PROT):>3}"
          "   (8 de cambio de aminoacido y 1 de perdida del codon de inicio)")

    print("\nCerteza de la orientacion frente al estudio de origen:")
    for v, e in ((CONF, "confirmada"),
                 (PROV_HEBRA, "provisional: hebra no verificable"),
                 (PROV_IDENT, "provisional: identidad del marcador"),
                 (INDE, "indeterminada")):
        print(f"  {e:<34}{cuenta('certeza_de_la_orientacion', v):>3}")

    print("\nRespaldo del gen:")
    for v, e in ((GEN_EST, "establecido"), (GEN_FRA, "fragmentario")):
        print(f"  {e:<34}{cuenta('respaldo_del_gen', v):>3}")
    print("\nRespaldo de la variante:")

    def por_prefijo(pref):
        return sum(1 for f in filas
                   if f["respaldo_de_la_variante"].startswith(pref))

    n_fdr = por_prefijo("interaccion limitrofe")
    n_sin = por_prefijo("evaluada en el estudio de origen")
    n_no = por_prefijo("no evaluada")
    print(f"  interaccion limitrofe, FDR 0.06       {n_fdr:>3}")
    print(f"  evaluada sin modificacion             {n_sin:>3}")
    print(f"  no evaluada fuera del panel           {n_no:>3}")
    if n_fdr + n_sin + n_no != len(filas):
        sys.exit("Las categorias de respaldo no suman el total de loci.")
    print(f"\nEscrito en {SALIDA}")


if __name__ == "__main__":
    main()
