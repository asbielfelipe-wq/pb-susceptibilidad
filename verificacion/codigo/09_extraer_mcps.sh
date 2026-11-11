#!/usr/bin/env bash
# Extrae las regiones de interes de los VCF de frecuencias de MCPS.
#
# Los VCF estan indexados en el servidor, asi que tabix descarga solo los
# bloques correspondientes a cada region en lugar de los 4.4 GB completos.
#
# El indice publicado es .tbi, pero htslib reciente busca .csi primero y falla.
# Se le indica el indice de forma explicita con la sintaxis ##idx##.
#
# Uso: bash scripts/09_extraer_mcps.sh [chr|sin]
#   chr  -> regiones como chr1:1000-2000   (por defecto)
#   sin  -> regiones como 1:1000-2000
set -uo pipefail

BASE="https://rgc-mcps.regeneron.com/downloads/20230130"
REGIONES="data/raw/regiones.txt"
DEST="data/raw/mcps"
PREFIJO="${1:-chr}"

[[ -f "$REGIONES" ]] || { echo "Falta $REGIONES (correr el script 08)"; exit 1; }
mkdir -p "$DEST"

echo "Modo de contig: $PREFIJO"
echo ""

n_ok=0; n_vacio=0; n_error=0
: > "$DEST/mcps_regiones.vcf"

while read -r region; do
  [[ -z "$region" ]] && continue

  # region viene como chr9:113366312-113406290
  chrom="${region%%:*}"
  resto="${region#*:}"
  if [[ "$PREFIJO" == "sin" ]]; then
    consulta="${chrom#chr}:$resto"
  else
    consulta="$region"
  fi

  archivo="$BASE/${chrom}.freq.vcf.gz"
  idx="$archivo.tbi"

  printf "%-28s " "$consulta"
  salida=$(tabix "${archivo}##idx##${idx}" "$consulta" 2>"$DEST/.err")

  if [[ -s "$DEST/.err" ]] && [[ -z "$salida" ]]; then
    echo "ERROR: $(head -c 120 "$DEST/.err")"
    n_error=$((n_error+1))
  elif [[ -z "$salida" ]]; then
    echo "sin variantes"
    n_vacio=$((n_vacio+1))
  else
    echo "$salida" >> "$DEST/mcps_regiones.vcf"
    n=$(echo "$salida" | wc -l | tr -d ' ')
    echo "$n variantes"
    n_ok=$((n_ok+1))
  fi
done < "$REGIONES"

rm -f "$DEST/.err"

echo ""
echo "regiones con datos: $n_ok   vacias: $n_vacio   con error: $n_error"
echo "total de lineas: $(wc -l < "$DEST/mcps_regiones.vcf" | tr -d ' ')"
echo "salida: $DEST/mcps_regiones.vcf"

if [[ $n_error -gt 0 && $n_ok -eq 0 ]]; then
  echo ""
  echo "Ninguna region funciono. Probar el otro modo de contig:"
  [[ "$PREFIJO" == "chr" ]] && echo "  bash scripts/09_extraer_mcps.sh sin" \
                            || echo "  bash scripts/09_extraer_mcps.sh chr"
fi
