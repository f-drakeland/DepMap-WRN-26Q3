from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr, t


# Keep generated SVG element IDs deterministic across runs.
plt.rcParams["svg.hashsalt"] = "depmap-wrn-26q3"


RECON = Path(
    "results/current-wgs/reconstruction_diagnostics.csv"
)

WRN = Path(
    "results/current-wgs/wrn_association_diagnostics.csv"
)

MSI = Path(
    "results/current-wgs/beyond_msi_diagnostics.csv"
)

OUT = Path(
    "figures/current-wgs"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


def save_figure(fig, stem):
    png = OUT / f"{stem}.png"
    svg = OUT / f"{stem}.svg"

    fig.savefig(
        png,
        dpi=240,
        bbox_inches="tight",
    )

    fig.savefig(
        svg,
        bbox_inches="tight",
        metadata={
            "Date": None,
        },
    )

    # Matplotlib SVG output includes trailing spaces in path data.
    # Normalize them so generated artifacts pass git diff --check
    # and remain byte-for-byte reproducible.
    svg_text = svg.read_text()
    svg.write_text(
        "\n".join(
            line.rstrip()
            for line in svg_text.splitlines()
        )
        + "\n"
    )

    plt.close(fig)

    print(f"Wrote: {png}")
    print(f"Wrote: {svg}")


def regression_line(x, y):
    fit = linregress(x, y)

    xs = np.linspace(
        float(np.min(x)),
        float(np.max(x)),
        200,
    )

    ys = (
        fit.intercept
        + fit.slope * xs
    )

    return fit, xs, ys


# ============================================================
# Figure 1
# Published FPBM vs reconstructed current-WGS FPBM
# ============================================================

recon = pd.read_csv(RECON)

x = recon["PublishedFPBM"].to_numpy(float)
y = recon["CurrentBrokenTAFPBM"].to_numpy(float)

r, p = pearsonr(x, y)
fit, xs, ys = regression_line(x, y)

fig, ax = plt.subplots(
    figsize=(7.2, 6.0)
)

ax.scatter(
    x,
    y,
    s=55,
)

identity_min = min(
    float(x.min()),
    float(y.min()),
)

identity_max = max(
    float(x.max()),
    float(y.max()),
)

ax.plot(
    [identity_min, identity_max],
    [identity_min, identity_max],
    linestyle=":",
    linewidth=1.5,
    label="Identity",
)

ax.plot(
    xs,
    ys,
    linewidth=1.8,
    label="Linear fit",
)

jhuem = recon[
    recon["ModelID"] == "ACH-000928"
].iloc[0]

ax.annotate(
    "JHUEM-1",
    (
        jhuem["PublishedFPBM"],
        jhuem["CurrentBrokenTAFPBM"],
    ),
    xytext=(8, 10),
    textcoords="offset points",
)

ax.set_xlabel(
    "Published broken-TA FPBM"
)

ax.set_ylabel(
    "Reconstructed current-WGS FPBM"
)

ax.set_title(
    "Current public WGS reproduces published broken-TA FPBM"
)

ax.text(
    0.03,
    0.97,
    (
        f"n = {len(recon)}\n"
        f"Pearson r = {r:.3f}\n"
        f"R² = {fit.rvalue ** 2:.3f}\n"
        f"fit: y = {fit.intercept:.2f} + "
        f"{fit.slope:.2f}x"
    ),
    transform=ax.transAxes,
    va="top",
)

ax.legend()

save_figure(
    fig,
    "01_published_vs_current_fpbm",
)


# ============================================================
# Figure 2
# Reconstructed current-WGS FPBM vs current WRN dependency
# ============================================================

wrn = pd.read_csv(WRN)

x = wrn[
    "CurrentBrokenTAFPBM"
].to_numpy(float)

y = wrn[
    "WRNGeneEffect"
].to_numpy(float)

r, p = pearsonr(x, y)
fit, xs, ys = regression_line(x, y)

fig, ax = plt.subplots(
    figsize=(7.2, 6.0)
)

ax.scatter(
    x,
    y,
    s=55,
)

ax.plot(
    xs,
    ys,
    linewidth=1.8,
    label="Linear fit",
)

for model_id, label, offset in [
    (
        "ACH-000969",
        "KM12",
        (10, 10),
    ),
    (
        "ACH-000943",
        "RKO",
        (8, -18),
    ),
]:
    row = wrn[
        wrn["ModelID"] == model_id
    ].iloc[0]

    ax.annotate(
        label,
        (
            row["CurrentBrokenTAFPBM"],
            row["WRNGeneEffect"],
        ),
        xytext=offset,
        textcoords="offset points",
    )

ax.set_xlabel(
    "Reconstructed current-WGS broken-TA FPBM"
)

ax.set_ylabel(
    "Current WRN GeneEffect"
)

ax.set_title(
    "Reconstructed broken-TA FPBM retains WRN association"
)

ax.text(
    0.03,
    0.97,
    (
        f"n = {len(wrn)}\n"
        f"Pearson r = {r:.3f}\n"
        f"R² = {fit.rvalue ** 2:.3f}"
    ),
    transform=ax.transAxes,
    va="top",
)

save_figure(
    fig,
    "02_current_fpbm_vs_wrn",
)


# ============================================================
# Figure 3
# MSI-adjusted residual relationship
# ============================================================

msi = pd.read_csv(MSI)

x = msi[
    "CurrentFPBMResidualAfterMSI"
].to_numpy(float)

y = msi[
    "WRNResidualAfterMSI"
].to_numpy(float)

r, _ = pearsonr(x, y)

# Partial-correlation significance with one covariate:
# df = n - 3, matching the canonical analysis.
df = len(msi) - 3
t_stat = r * np.sqrt(
    df / (1.0 - r ** 2)
)
p = 2.0 * t.sf(
    abs(t_stat),
    df,
)

fit, xs, ys = regression_line(x, y)

fig, ax = plt.subplots(
    figsize=(7.2, 6.0)
)

ax.scatter(
    x,
    y,
    s=55,
)

ax.plot(
    xs,
    ys,
    linewidth=1.8,
    label="Residual linear fit",
)

ax.axhline(
    0,
    linestyle=":",
    linewidth=1.0,
)

ax.axvline(
    0,
    linestyle=":",
    linewidth=1.0,
)

for model_id, label, offset in [
    (
        "ACH-000969",
        "KM12",
        (10, 10),
    ),
    (
        "ACH-000928",
        "JHUEM-1",
        (8, 10),
    ),
]:
    row = msi[
        msi["ModelID"] == model_id
    ].iloc[0]

    ax.annotate(
        label,
        (
            row[
                "CurrentFPBMResidualAfterMSI"
            ],
            row[
                "WRNResidualAfterMSI"
            ],
        ),
        xytext=offset,
        textcoords="offset points",
    )

ax.set_xlabel(
    "FPBM residual after MSI adjustment"
)

ax.set_ylabel(
    "WRN GeneEffect residual after MSI adjustment"
)

ax.set_title(
    "Broken-TA FPBM retains a modest MSI-adjusted relationship"
)

ax.text(
    0.03,
    0.97,
    (
        f"n = {len(msi)}\n"
        f"Partial Pearson r = {r:.3f}\n"
        f"Parametric p = {p:.4f}\n"
        "Permutation p = 0.0430\n"
        "Bootstrap 95% CI = -0.145 to 0.720"
    ),
    transform=ax.transAxes,
    va="top",
)

save_figure(
    fig,
    "03_fpbm_beyond_msi",
)

print()
print("PASS: current-WGS figure set generated.")
