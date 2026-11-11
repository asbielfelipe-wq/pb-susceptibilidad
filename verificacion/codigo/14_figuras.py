#!/usr/bin/env python3
"""
Manuscript figures, grouped into three composite panels by the argument each
one carries. All labels in English for submission.

Figure 1 - Allele frequency distribution at lead susceptibility loci
  A  Alternate allele frequency in the three ancestry components, one locus per
     row, ordered by the IMX minus EUR difference.
  B  IMX against EUR frequency for every common variant in the regions, with
     the candidate loci overlaid.

Figure 2 - Differentiation against the empirical background
  A  FST against mean minor allele frequency: control-gene background in grey,
     candidates in colour, 95th percentile curve of the background.
  B  Percentile of each locus with bootstrap confidence interval, by stratum.

Figure 3 - Demographic history and robustness
  A  IMX against EUR frequency coloured by FST, with the bands where either
     frequency falls below 0.02 shaded.
  B  Outlier block rate per Mb in candidate and control genes across the three
     percentile thresholds, with and without the near-private allele control.
  C  FST profile across the AGT region and rs699 frequency in population
     context.

Run: python3 scripts/14_figuras.py
Output: results/figuras/fig1.png, fig2.png, fig3.png (+ .pdf)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
PROC = RAIZ / "data" / "processed"
SAL = RES / "figuras"

C_IMX, C_EUR, C_AFR = "#1f3864", "#c0504d", "#8c8c8c"
C_CAND, C_FONDO, C_ACENTO = "#1f3864", "#cfcfcf", "#c0504d"

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def panel(ax, letra, dx=-0.14):
    ax.text(dx, 1.04, letra, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom", ha="left")


def leer(nombre, carpeta=RES):
    f = carpeta / nombre
    if not f.exists():
        sys.exit(f"Missing {f}. Run the analysis scripts first.")
    return pd.read_csv(f, sep="\t", low_memory=False)


def col(df, prefijo):
    c = [x for x in df.columns if x.startswith(prefijo)]
    if not c:
        sys.exit(f"No column starting with {prefijo}")
    return c[0]


def separar_etiquetas(ys, minimo):
    """Desplaza verticalmente etiquetas demasiado proximas, conservando el orden."""
    orden = np.argsort(ys)
    y = np.array(ys, dtype=float)
    for i in range(1, len(orden)):
        a, b = orden[i - 1], orden[i]
        if y[b] - y[a] < minimo:
            y[b] = y[a] + minimo
    return y


# ---------------------------------------------------------------- figure 1
def figura1(pct, ding, barrido):
    fig = plt.figure(figsize=(10, 6.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=.34)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # A
    a_imx, a_eur = col(ding, "AF_IMX"), col(ding, "AF_EUR")
    a_afr = col(ding, "AF_AFR")
    d = ding.dropna(subset=[a_imx, a_eur]).copy()
    d["delta"] = d[a_imx] - d[a_eur]
    d = d.sort_values("delta")
    y = np.arange(len(d))

    for yi, (_, r) in zip(y, d.iterrows()):
        ax_a.plot([r[a_eur], r[a_imx]], [yi, yi], color="#bdbdbd", lw=1.1, zorder=1)
    ax_a.scatter(d[a_eur], y, s=26, color=C_EUR, zorder=3, label="European")
    ax_a.scatter(d[a_imx], y, s=26, color=C_IMX, zorder=3, label="Indigenous Mexican")
    if a_afr in d:
        ax_a.scatter(d[a_afr], y, s=14, color=C_AFR, marker="s", zorder=2,
                     alpha=.85, label="African")
    ax_a.set_yticks(y)
    ax_a.set_yticklabels([f"{r['gen']} {r['rsid']}" for _, r in d.iterrows()])
    ax_a.set_xlim(-0.04, 1.03)
    ax_a.set_ylim(-0.7, len(d) - 0.3)
    ax_a.set_xlabel("Alternate allele frequency")
    ax_a.axvline(0, color="#e0e0e0", lw=.8, zorder=0)
    ax_a.legend(loc="lower right", frameon=False)
    ax_a.set_title("Allele frequency by ancestry component", loc="left")
    ax_a.text(0.0, -0.085, "red marks: frequency below 0.002 in Indigenous Mexican ancestry",
              transform=ax_a.transAxes, fontsize=6.3, color="#616161", ha="left")
    panel(ax_a, "A", dx=-0.30)

    for yi, (_, r) in zip(y, d.iterrows()):
        if r[a_imx] < 0.002:
            ax_a.plot([-0.032], [yi], marker="<", ms=5, color=C_ACENTO,
                      clip_on=False, zorder=5)

    # B, sin aspecto fijo para que ocupe la misma altura que A
    b_imx, b_eur = col(barrido, "AF_IMX"), col(barrido, "AF_EUR")
    fondo = barrido[barrido["grupo"] == "control"]
    ax_b.scatter(fondo[b_eur], fondo[b_imx], s=3, color=C_FONDO, alpha=.45,
                 lw=0, rasterized=True, label="Control-gene variants")
    ax_b.scatter(pct["AF_EUR"], pct["AF_IMX"], s=32, color=C_CAND,
                 edgecolor="white", lw=.5, zorder=3, label="Candidate loci")
    ax_b.plot([0, 1], [0, 1], color="#9e9e9e", ls="--", lw=.8, zorder=1)
    ax_b.set_xlabel("Frequency in European ancestry")
    ax_b.set_ylabel("Frequency in Indigenous Mexican ancestry")
    ax_b.set_xlim(-0.02, 1.02)
    ax_b.set_ylim(-0.02, 1.02)
    ax_b.legend(loc="upper left", frameon=False, markerscale=1.8)
    ax_b.set_title("Candidate loci within overall variation", loc="left")
    panel(ax_b, "B", dx=-0.17)

    return fig


# ---------------------------------------------------------------- figure 2
def figura2(pct, barrido):
    fig = plt.figure(figsize=(10, 6.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.05], wspace=.3)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # A
    fondo = barrido[barrido["grupo"] == "control"].dropna(subset=["FST", "maf_media"]).copy()
    # el barrido guarda maf_media con la definicion antigua; se recalcula con la
    # frecuencia conjunta, que es la del analisis principal
    bi, be = col(barrido, "AF_IMX"), col(barrido, "AF_EUR")
    pb = (fondo[bi].values + fondo[be].values) / 2
    fondo["maf_media"] = np.minimum(pb, 1 - pb)
    ax_a.scatter(fondo["maf_media"], fondo["FST"], s=3, color=C_FONDO,
                 alpha=.5, lw=0, rasterized=True)

    bins = np.linspace(0.01, 0.5, 40)
    centros, p95 = [], []
    for i in range(len(bins) - 1):
        g = fondo[(fondo["maf_media"] >= bins[i]) & (fondo["maf_media"] < bins[i + 1])]
        if len(g) >= 25:
            centros.append((bins[i] + bins[i + 1]) / 2)
            p95.append(np.percentile(g["FST"], 95))
    if len(p95) >= 5:
        suave = pd.Series(p95).rolling(5, center=True, min_periods=1).median()
        ax_a.plot(centros, suave, color=C_ACENTO, lw=1.6, zorder=4,
                  label="95th percentile of background")

    ax_a.scatter(pct["maf_media"], pct["FST"], s=34, color=C_CAND,
                 edgecolor="white", lw=.6, zorder=5, label="Candidate loci")

    # Etiquetas en vertical sobre o bajo cada punto, alternando. Sacarlas en
    # horizontal las hace invadir el espacio de los puntos vecinos, que estan
    # cerca en FST pero lejos en frecuencia.
    altos = pct.sort_values("FST", ascending=False).head(4)
    if len(altos):
        rango = max(float(pct["FST"].max()), 0.05)
        altos = altos.sort_values("maf_media")
        for k, (_, r) in enumerate(altos.iterrows()):
            arriba = (k % 2 == 0)
            dy = rango * (0.13 if arriba else -0.13)
            ax_a.annotate(f"{r['gen']} {r['rsid']}", (r["maf_media"], r["FST"]),
                          xytext=(r["maf_media"], r["FST"] + dy), fontsize=6.5,
                          color=C_CAND, ha="center",
                          va="bottom" if arriba else "top", zorder=6,
                          arrowprops=dict(arrowstyle="-", lw=.5, color=C_CAND,
                                          shrinkA=1, shrinkB=4))
    ax_a.margins(y=0.12)
    ax_a.set_xlim(-0.01, 0.52)
    ax_a.set_xlabel("Matching frequency, min(p\u0304, 1\u2212p\u0304)")
    ax_a.set_ylabel("Hudson's $F_{ST}$ (IMX vs EUR)")
    ax_a.legend(loc="upper right", frameon=False, markerscale=1.4)
    ax_a.set_title("Differentiation against control-gene background", loc="left")
    panel(ax_a, "A", dx=-0.17)

    # B
    d = pct.sort_values(["estrato", "pct_FST"])
    y, etiquetas, sep, filas = [], [], [], []
    pos = 0
    for est in ["alta", "wgs"]:
        s = d[d["estrato"] == est]
        if not len(s):
            continue
        if pos:
            sep.append(pos - 0.5)
            pos += 0.8
        for _, r in s.iterrows():
            y.append(pos)
            etiquetas.append(f"{r['gen']} {r['rsid']}")
            filas.append(r)
            pos += 1
    s = pd.DataFrame(filas)

    ax_b.hlines(y, s["IC95_inf"], s["IC95_sup"], color="#9e9e9e", lw=1.4, zorder=2)
    colores = [C_ACENTO if r["IC95_inf"] >= 95 else C_CAND for _, r in s.iterrows()]
    ax_b.scatter(s["pct_FST"], y, s=28, color=colores, zorder=3,
                 edgecolor="white", lw=.5)
    ax_b.axvline(95, color=C_ACENTO, ls="--", lw=1, zorder=1)
    for yy in sep:
        ax_b.axhline(yy, color="#e0e0e0", lw=.8)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(etiquetas)
    ax_b.set_xlim(-2, 106)
    ax_b.set_ylim(min(y) - 1.4, max(y) + 0.8)
    ax_b.set_xlabel("$F_{ST}$ percentile in frequency-matched background")
    ax_b.invert_yaxis()
    ax_b.set_title("Percentile with block-resampling interval", loc="left")
    ax_b.text(2, min(y) - 1.0, "Array and exome", fontsize=6.5, color="#616161")
    if sep:
        ax_b.text(2, sep[0] + 0.55, "Whole genome", fontsize=6.5, color="#616161")
    panel(ax_b, "B", dx=-0.30)

    return fig


# ---------------------------------------------------------------- figure 3
def figura3(barrido, exceso, mcps, pct):
    fig = plt.figure(figsize=(10, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=.42, wspace=.34)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    sub = gs[1, :].subgridspec(1, 2, width_ratios=[2.3, 1], wspace=.30)
    ax_c1 = fig.add_subplot(sub[0, 0])
    ax_c2 = fig.add_subplot(sub[0, 1])

    b_imx, b_eur = col(barrido, "AF_IMX"), col(barrido, "AF_EUR")

    # A
    d = barrido.dropna(subset=["FST", b_imx, b_eur])
    sc = ax_a.scatter(d[b_eur], d[b_imx], c=d["FST"], s=4, cmap="YlOrRd",
                      vmin=0, vmax=max(0.3, float(d["FST"].quantile(.999))),
                      lw=0, rasterized=True)
    for v in (0, 0.98):
        ax_a.axvspan(v, v + 0.02, color="#9e9e9e", alpha=.18, lw=0)
        ax_a.axhspan(v, v + 0.02, color="#9e9e9e", alpha=.18, lw=0)
    cb = fig.colorbar(sc, ax=ax_a, fraction=.046, pad=.03)
    cb.set_label("$F_{ST}$", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax_a.set_xlabel("Frequency in European ancestry")
    ax_a.set_ylabel("Frequency in Indigenous Mexican ancestry")
    ax_a.set_xlim(0, 1); ax_a.set_ylim(0, 1)
    ax_a.set_title("Near-private alleles and extreme $F_{ST}$", loc="left")
    ax_a.text(.5, -0.20, "shaded bands: allele frequency < 0.02", fontsize=6.5,
              color="#616161", ha="center", transform=ax_a.transAxes)
    panel(ax_a, "A", dx=-0.19)

    # B
    e = exceso[(exceso["ventana_kb"] == 50) & (exceso["estrato"] == "wgs")].copy()
    orden = [99.0, 97.5, 95.0]
    x = np.arange(len(orden)); w = 0.2
    series = [
        ("completo", "cand", C_CAND, "Candidate, no control"),
        ("completo", "ctrl", "#bdbdbd", "Control, no control"),
        ("sin_casi_privadas", "cand", C_ACENTO, "Candidate, with control"),
        ("sin_casi_privadas", "ctrl", "#e6b8b7", "Control, with control"),
    ]
    for i, (an, grp, cl, lab) in enumerate(series):
        vals = []
        for p in orden:
            f = e[(e["analisis"] == an) & (e["pct"] == p)]
            v = f[f"tasa_{grp}_Mb"]
            try:
                vals.append(float(v.iloc[0]))
            except (IndexError, ValueError):
                vals.append(np.nan)
        ax_b.bar(x + (i - 1.5) * w, vals, w, color=cl, label=lab)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([f"$\\geq${p:g}" for p in orden])
    ax_b.set_xlabel("Percentile threshold")
    ax_b.set_ylabel("Outlier blocks per Mb")
    ax_b.legend(frameon=False, fontsize=6.3, loc="upper left")
    ax_b.set_title("Effect of the near-private allele control", loc="left")
    panel(ax_b, "B", dx=-0.17)

    # C1
    agt = mcps[mcps["gen"] == "AGT"].copy() if "gen" in mcps.columns else pd.DataFrame()
    if len(agt):
        m_imx, m_eur = col(agt, "AF_IMX"), col(agt, "AF_EUR")
        m_ani, m_ane = col(agt, "AN_IMX"), col(agt, "AN_EUR")
        agt = agt.dropna(subset=[m_imx, m_eur, m_ani, m_ane])
        agt = agt[(agt["reparto_default"] == 0) & (agt["es_indel"] == 0)]
        p1, p2 = agt[m_imx].values, agt[m_eur].values
        n1, n2 = agt[m_ani].values, agt[m_ane].values
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)
        agt["FST"] = np.where(den > 0, num / den, np.nan)
        maf = (np.minimum(p1, 1 - p1) + np.minimum(p2, 1 - p2)) / 2
        agt = agt[maf >= 0.01].dropna(subset=["FST"])
        ax_c1.scatter(agt["pos"] / 1e6, agt["FST"], s=7, color="#b0b0b0", lw=0)

        marcas = pct[pct["gen"] == "AGT"].copy()
        if len(marcas):
            marcas["pos"] = marcas["variant_id"].apply(lambda v: int(v.split("-")[1]))
            marcas = marcas.sort_values("pos")
            rango = float(agt["FST"].max() - agt["FST"].min()) or 1.0
            for k, (_, r) in enumerate(marcas.iterrows()):
                ax_c1.scatter([r["pos"] / 1e6], [r["FST"]], s=44, color=C_ACENTO,
                              zorder=4, edgecolor="white", lw=.6)
                arriba = (k % 2 == 0)
                dy = rango * (0.16 if arriba else -0.16)
                ax_c1.annotate(r["rsid"], (r["pos"] / 1e6, r["FST"]),
                               xytext=(r["pos"] / 1e6, r["FST"] + dy), fontsize=6.5,
                               color=C_ACENTO, ha="center",
                               va="bottom" if arriba else "top", zorder=6,
                               arrowprops=dict(arrowstyle="-", lw=.5, color=C_ACENTO,
                                               shrinkA=1, shrinkB=4))
            ax_c1.set_xlim(float(agt["pos"].min() / 1e6) - 0.002,
                           float(agt["pos"].max() / 1e6) + 0.002)
            ax_c1.margins(y=0.14)
        ax_c1.set_xlabel("Position on chromosome 1 (Mb)")
        ax_c1.set_ylabel("$F_{ST}$")
    ax_c1.set_title("Differentiation profile across the AGT region", loc="left")
    panel(ax_c1, "C", dx=-0.11)

    # C2
    r699 = pct[pct["rsid"] == "rs699"]
    f_imx = float(r699["AF_IMX"].iloc[0]) if len(r699) else np.nan
    f_eur = float(r699["AF_EUR"].iloc[0]) if len(r699) else np.nan
    nombres = ["Indigenous\nMexican", "European", "Japanese", "French"]
    vals = [f_imx, f_eur, 0.76, 0.38]
    colores = [C_IMX, C_EUR, "#7f7f7f", "#c4c4c4"]
    barras = ax_c2.barh(range(4)[::-1], vals, color=colores, height=.62)
    for b, v in zip(barras, vals):
        if v == v:
            ax_c2.text(v + 0.02, b.get_y() + b.get_height() / 2, f"{v:.2f}",
                       va="center", fontsize=6.5, color="#424242")
    ax_c2.set_yticks(range(4)[::-1])
    ax_c2.set_yticklabels(nombres, fontsize=7)
    ax_c2.set_xlim(0, 1.12)
    ax_c2.set_xlabel("T235 allele frequency")
    ax_c2.set_title("rs699 in population context", loc="left", fontsize=8)
    ax_c2.text(0.0, -0.22, "Japanese and French values from Jeunemaitre et al. 1997",
               transform=ax_c2.transAxes, fontsize=6, color="#616161", ha="left")

    return fig


def main():
    pct = leer("principal_percentiles.tsv")
    # el analisis principal usa maf_conjunta, pct e IC95_*; se renombran a los
    # nombres que usan las figuras y se descartan los loci no evaluables
    pct = pct.rename(columns={"maf_conjunta": "maf_media", "pct": "pct_FST"})
    pct = pct[pct["evaluable"] == "si"].copy()
    ding = leer("ding_en_mcps_final.tsv")
    barrido = leer("barrido_todas.tsv")
    exceso = leer("exceso_pruebas.tsv")
    mcps = pd.read_csv(PROC / "mcps_frecuencias.tsv", sep="\t", low_memory=False)

    SAL.mkdir(parents=True, exist_ok=True)
    for n, f in [("fig1", figura1(pct, ding, barrido)),
                 ("fig2", figura2(pct, barrido)),
                 ("fig3", figura3(barrido, exceso, mcps, pct))]:
        f.savefig(SAL / f"{n}.png")
        f.savefig(SAL / f"{n}.pdf")
        plt.close(f)
        print(f"  {n}.png / {n}.pdf")
    print(f"\nFigures in {SAL}")


if __name__ == "__main__":
    main()
