import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# ── 1. Load data ─────────────────────────────────────────────────────────────

print("Loading data...")

# Gene expression matrix (rows = samples, cols = gene IDs)
expr = pd.read_csv(
    "/mnt/user-data/uploads/filtered_tsv.gz",
    sep="\t",
    header=0,
    compression="gzip",
)
# Strip whitespace from column names (the header had leading spaces)
expr.columns = expr.columns.str.strip()

# Class labels  (105 rows, no header; 1 = ER+, 0 = ER-)
labels = pd.read_csv(
    "/mnt/user-data/uploads/class.tsv",
    sep="\t",
    header=None,
    names=["label"],
).squeeze()

print(f"Expression matrix : {expr.shape}  (samples × genes)")
print(f"Labels            : {labels.shape}  (ER+={labels.sum()}, ER-={(labels==0).sum()})")

# ── 2. Identify XBP1 (4404) and GATA3 (4359) columns ────────────────────────

XBP1_ID  = "4404"
GATA3_ID = "4359"

xbp1  = expr[XBP1_ID].values
gata3 = expr[GATA3_ID].values

# ── 3. Colour palette ─────────────────────────────────────────────────────────
# Paper uses red for ER+ and blue/green for ER-
ER_POS_COLOR = "#D62728"   # red
ER_NEG_COLOR = "#1F77B4"   # blue

colors = np.where(labels == 1, ER_POS_COLOR, ER_NEG_COLOR)

# ── 4. Figure 1a – XBP1 vs GATA3 scatter plot ───────────────────────────────

fig, ax = plt.subplots(figsize=(6, 5))

ax.scatter(
    gata3[labels == 0], xbp1[labels == 0],
    c=ER_NEG_COLOR, s=40, alpha=0.85, linewidths=0.4,
    edgecolors="white", label="ER− (0)", zorder=3,
)
ax.scatter(
    gata3[labels == 1], xbp1[labels == 1],
    c=ER_POS_COLOR, s=40, alpha=0.85, linewidths=0.4,
    edgecolors="white", label="ER+ (1)", zorder=3,
)

ax.set_xlabel("GATA3 expression", fontsize=12)
ax.set_ylabel("XBP1 expression",  fontsize=12)
ax.set_title("Figure 1a – XBP1 vs GATA3\n(GSE5325 breast cancer cohort)", fontsize=12)

ax.legend(
    handles=[
        mpatches.Patch(color=ER_POS_COLOR, label="ER+"),
        mpatches.Patch(color=ER_NEG_COLOR, label="ER−"),
    ],
    frameon=True, fontsize=10,
)

ax.axhline(0, color="grey", linewidth=0.5, linestyle="--", zorder=1)
ax.axvline(0, color="grey", linewidth=0.5, linestyle="--", zorder=1)
ax.grid(True, alpha=0.25, zorder=0)
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig("/mnt/user-data/outputs/figure_1a_scatter.png", dpi=180, bbox_inches="tight")
plt.close()
print("Saved: figure_1a_scatter.png")

# ── 5. PCA on full expression matrix ─────────────────────────────────────────

print("\nRunning PCA...")

# Standardise: zero-mean, unit-variance per gene (standard for PCA on microarray)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(expr.values)   # shape: (105, 16174)

pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X_scaled)            # shape: (105, 2)

pc1 = X_pca[:, 0]
pc2 = X_pca[:, 1]

var_pc1 = pca.explained_variance_ratio_[0] * 100
var_pc2 = pca.explained_variance_ratio_[1] * 100

print(f"PC1 explains {var_pc1:.1f}% of variance")
print(f"PC2 explains {var_pc2:.1f}% of variance")

# ── 6. Figure 1c – Projection onto PC1 (strip / dot plot) ───────────────────
# The paper shows samples projected on the first principal component as a
# 1-D strip plot: each sample is a dot along the PC1 axis, jittered slightly
# on the y-axis to avoid overplotting, coloured by ER status.

rng  = np.random.default_rng(seed=0)
jitter_scale = 0.10

fig, ax = plt.subplots(figsize=(8, 3.5))

for cls, color, label_str in [(0, ER_NEG_COLOR, "ER−"), (1, ER_POS_COLOR, "ER+")]:
    mask   = labels.values == cls
    pc1_cl = pc1[mask]
    jitter = rng.uniform(-jitter_scale, jitter_scale, size=mask.sum())
    ax.scatter(
        pc1_cl, jitter,
        c=color, s=40, alpha=0.85,
        edgecolors="white", linewidths=0.4,
        label=label_str, zorder=3,
    )

ax.axhline(0, color="grey", linewidth=0.8, zorder=1)
ax.set_xlabel(f"PC1  ({var_pc1:.1f}% variance explained)", fontsize=12)
ax.set_ylabel("")
ax.set_yticks([])
ax.set_title("Figure 1c – Projection onto PC1\n(GSE5325 breast cancer cohort)", fontsize=12)

ax.legend(
    handles=[
        mpatches.Patch(color=ER_POS_COLOR, label="ER+"),
        mpatches.Patch(color=ER_NEG_COLOR, label="ER−"),
    ],
    frameon=True, fontsize=10, loc="upper left",
)

ax.grid(axis="x", alpha=0.25, zorder=0)
ax.spines[["top", "right", "left"]].set_visible(False)

plt.tight_layout()
plt.savefig("/mnt/user-data/outputs/figure_1c_pca.png", dpi=180, bbox_inches="tight")
plt.close()
print("Saved: figure_1c_pca.png")
print("\nDone.")
