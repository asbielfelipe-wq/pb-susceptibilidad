#!/usr/bin/env bash
set -euo pipefail

BASE="https://ctdbase.org/reports"
DEST="data/raw/ctd"
mkdir -p "$DEST"

ARCHIVOS=(
  CTD_chem_gene_ixns
  CTD_chemicals_diseases
  CTD_genes_pathways
  CTD_chemicals
  CTD_genes
)

tamano() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1"; }

hash256() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

descomprimir() {
  if command -v gzcat >/dev/null 2>&1; then gzcat "$1"; else zcat "$1"; fi
}

for f in "${ARCHIVOS[@]}"; do
  if [[ -f "$DEST/$f.tsv.gz" ]]; then
    echo "ya existe: $f.tsv.gz"
  else
    echo "descargando: $f.tsv.gz"
    curl -fL --retry 3 -o "$DEST/$f.tsv.gz" "$BASE/$f.tsv.gz"
  fi
done

{
  echo "fecha_registro_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "fuente: $BASE"
  echo ""
  printf "%s\t%s\t%s\n" "archivo" "bytes" "sha256"
  for f in "${ARCHIVOS[@]}"; do
    ruta="$DEST/$f.tsv.gz"
    if [[ -f "$ruta" ]]; then
      printf "%s\t%s\t%s\n" "$f.tsv.gz" "$(tamano "$ruta")" "$(hash256 "$ruta")"
    else
      printf "%s\t%s\t%s\n" "$f.tsv.gz" "FALTA" "FALTA"
    fi
  done
} > "$DEST/PROCEDENCIA.txt"

echo ""
echo "--- encabezado de CTD_chem_gene_ixns (contiene la version) ---"
descomprimir "$DEST/CTD_chem_gene_ixns.tsv.gz" 2>/dev/null | head -n 12 || true

echo ""
cat "$DEST/PROCEDENCIA.txt"
