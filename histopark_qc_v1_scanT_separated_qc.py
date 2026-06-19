#!/usr/bin/env python3
"""
===============================================================================
HistoPark MRI QC Pipeline v1
===============================================================================

This script performs first-pass MRI quality control (QC) for the HistoPark
pilot study.

The pipeline is intentionally conservative and transparent:
    - no silent resampling of quantitative maps
    - no voxelwise 3T-vs-7T comparisons in v1
    - all file matches are logged
    - all metrics are exported to CSV
    - screenshots are generated for visual inspection

Main goals
-----------
1. Quantitative QC
    - summary statistics
    - histogram analysis
    - error-map analysis
    - B1 homogeneity analysis

2. Visual QC
    - coronal screenshots
    - mask overlays
    - segmentation overlays

3. Cohort QC
    - outlier detection
    - cohort heatmaps
    - STX vs PTX comparisons

Key design principles
---------------------
- Histograms use fixed plausible value ranges where possible.
- Tissue masks are resampled into image space, NOT vice versa.
- Quantitative maps remain untouched for metric computation.
- Chimap is treated as signed/zero-centered and therefore does not use
  classical coefficient of variation (CV).

Outputs
-------
CSV:
    qc_metrics.csv
    qc_error_metrics.csv
    qc_histograms.csv
    qc_b1_homogeneity.csv
    qc_file_manifest.csv
    qc_outlier_flags.csv

Figures:
    screenshots/
    segmentation_overlays/
    histograms/
    cohort_plots/
    b1_homogeneity/

===============================================================================
"""

from __future__ import annotations

import argparse
import json
import math
import os
import stat
import subprocess
import shlex
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import nibabel as nib
    from nibabel.processing import resample_from_to
except Exception as exc:
    raise RuntimeError("This script requires nibabel. Try: pip install nibabel scipy") from exc

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise RuntimeError("This script requires matplotlib. Try: pip install matplotlib") from exc

try:
    from scipy.stats import wasserstein_distance
except Exception:
    wasserstein_distance = None

try:
    from tqdm import tqdm
except Exception:
    tqdm = lambda x, **kwargs: x


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

MAP_PATTERNS: Dict[str, List[str]] = {
    # qMRI final maps
    "R1": ["*R1map.nii", "*R1map.nii.gz"],
    "R2star": ["*R2starmap.nii", "*R2starmap.nii.gz"],
    "PD": ["*PDmap.nii", "*PDmap.nii.gz"],
    "MTsat": ["*MTsat.nii", "*MTsat.nii.gz"],

    # QSMxT: Chimap only. singlepass ignored by find_map_file().
    "Chimap": ["*Chimap.nii", "*Chimap.nii.gz"],

    # relax_R2
    "T2": ["*T2map.nii", "*T2map.nii.gz"],
    "R2": ["*R2map.nii", "*R2map.nii.gz"],
    "TB1_relax": ["*TB1map.nii", "*TB1map.nii.gz"],

    # r2prime
    "R2prime": ["*R2primemap.nii", "*R2primemap.nii.gz"],

    # qMRI Supplementary B1
    "B1_STX": ["*acq-tr1stx*TB1AFI_B1map.nii", "*acq-tr1stx*TB1AFI_B1map.nii.gz"],
    "B1_PTX": ["*acq-tr1ptx*TB1AFI_B1map.nii", "*acq-tr1ptx*TB1AFI_B1map.nii.gz"],
}

ERROR_PATTERNS: Dict[str, List[str]] = {
    "param_error": ["*param_error.nii", "*param_error.nii.gz"],
    "errorESTATICS": ["*errorESTATICS.nii", "*errorESTATICS.nii.gz"],
    "mSNR": ["*mSNR.nii", "*mSNR.nii.gz"],
    "rel_error": ["*rel_err*.nii", "*rel_err*.nii.gz", "*relative*error*.nii", "*relative*error*.nii.gz"],
    "generic_error": ["*error*.nii", "*error*.nii.gz", "*err*.nii", "*err*.nii.gz"],
}


def infer_error_source_from_name(filename: str) -> str:
    """Best-effort association of an error map with the image/parameter it describes."""
    n = filename.lower()
    if "r1" in n:
        return "R1"
    if "r2star" in n or "r2s" in n:
        return "R2star"
    if "pd" in n:
        return "PD"
    if "mtsat" in n or "mt" in n:
        return "MTsat"
    if "b1" in n or "tb1" in n:
        return "B1"
    if "t2" in n:
        return "T2"
    if "r2prime" in n:
        return "R2prime"
    if "r2" in n:
        return "R2"
    if "chi" in n or "qsm" in n:
        return "Chimap"
    if "estatics" in n:
        return "ESTATICS"
    if "param_error" in n:
        return "qMRI_params"
    if "msnr" in n:
        return "mSNR"
    return "unknown"

WEIGHTED_PATTERNS: Dict[str, List[str]] = {
    "PDw_undistorted": ["*acq-PDw*part-mag*desc-LCPCAout_desc-undistortedJac_MPM.nii",
                        "*acq-PDw*part-mag*desc-LCPCAout_desc-undistortedJac_MPM.nii.gz"],
    "T1w_undistorted": ["*acq-T1w*part-mag*desc-LCPCAout_desc-undistortedJac_MPM.nii",
                        "*acq-T1w*part-mag*desc-LCPCAout_desc-undistortedJac_MPM.nii.gz"],
    "MTw_undistorted": ["*acq-MTw*part-mag*desc-LCPCAout_desc-undistortedJac_MPM.nii",
                        "*acq-MTw*part-mag*desc-LCPCAout_desc-undistortedJac_MPM.nii.gz"],
}

DEFAULT_SCREENSHOT_MAPS = [
    "R1", "R2star", "PD", "MTsat", "Chimap", "B1_STX", "B1_PTX", "T2", "R2"
]

# Fixed display / plausibility ranges for screenshots and range-based outlier fractions.
# Values outside these ranges are not automatically excluded from summary statistics;
# they are counted as implausible/out-of-range voxels.
MAP_VALUE_RANGES = {
    "R1": (0.0, 2.0),
    "R2star": (0.0, 100.0),
    "PD": (0.0, 120.0),
    "MTsat": (0.0, 2.0),
    "Chimap": (-0.5, 0.5),
    "B1_STX": (0.0, 200.0),
    "B1_PTX": (0.0, 200.0),
    "TB1_relax": (0.0, 200.0),
    "T2": (0.0, 1.0),
    "R2": (0.0, 100.0),
    "R2prime": (-50.0, 100.0),
}

# Mean-based CV is misleading for signed / near-zero maps such as Chimap.
SIGNED_OR_ZERO_CENTERED_MAPS = {"Chimap"}


@dataclass
class FoundFile:
    subject: str
    session: str
    scanT: str
    scanner: str
    category: str
    name: str
    path: str
    status: str
    n_matches: int
    all_matches: str


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

# =============================================================================
# Basic helper functions
# =============================================================================

def str_or_empty(x) -> str:
    if pd.isna(x):
        return ""
    return str(x)


def apply_path_replacements(path: str, replacements: List[Tuple[str, str]]) -> str:
    out = str_or_empty(path)
    for old, new in replacements:
        out = out.replace(old, new)
    return out


def existing_dir(path: str) -> Optional[Path]:
    if not path:
        return None
    p = Path(path)
    return p if p.exists() and p.is_dir() else None


def glob_many(base: Optional[Path], patterns: Sequence[str], recursive: bool = False) -> List[Path]:
    if base is None or not base.exists():
        return []
    matches: List[Path] = []
    for pat in patterns:
        if recursive:
            matches.extend(base.rglob(pat))
        else:
            matches.extend(base.glob(pat))
    # unique, stable
    unique = sorted({str(p): p for p in matches}.values(), key=lambda p: str(p))
    return [p for p in unique if p.is_file()]


def find_single(
    base: Optional[Path],
    patterns: Sequence[str],
    *,
    recursive: bool = False,
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    prefer: Optional[Sequence[str]] = None,
) -> Tuple[Optional[Path], List[Path], str]:
    matches = glob_many(base, patterns, recursive=recursive)

    if include:
        matches = [p for p in matches if all(s.lower() in p.name.lower() for s in include)]
    if exclude:
        matches = [p for p in matches if not any(s.lower() in p.name.lower() for s in exclude)]

    if not matches:
        return None, [], "missing"

    if prefer:
        preferred = matches
        for token in prefer:
            tmp = [p for p in preferred if token.lower() in p.name.lower()]
            if tmp:
                preferred = tmp
        matches_for_choice = preferred
    else:
        matches_for_choice = matches

    if len(matches_for_choice) == 1:
        return matches_for_choice[0], matches, "ok"

    # deterministic fallback, but mark ambiguous
    return matches_for_choice[0], matches, "ambiguous"


# =============================================================================
# Image I/O and geometry utilities
# =============================================================================

def load_img(path: Path) -> nib.Nifti1Image:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return nib.load(str(path))


def finite_data(img: nib.Nifti1Image) -> np.ndarray:
    data = np.asanyarray(img.dataobj).astype(np.float32)
    data = np.squeeze(data)
    data[~np.isfinite(data)] = np.nan
    return data


def same_grid(img_a: nib.Nifti1Image, img_b: nib.Nifti1Image, atol: float = 1e-3) -> bool:
    return (
        img_a.shape[:3] == img_b.shape[:3]
        and np.allclose(img_a.affine, img_b.affine, atol=atol)
    )


def resample_mask_to_img(mask_img: nib.Nifti1Image, ref_img: nib.Nifti1Image, order: int = 0) -> np.ndarray:
    if same_grid(mask_img, ref_img):
        data = finite_data(mask_img)
    else:
        res = resample_from_to(mask_img, ref_img, order=order)
        data = finite_data(res)
    return data


# =============================================================================
# Statistical helper functions
# =============================================================================

def robust_mad(x: np.ndarray) -> float:
    med = np.nanmedian(x)
    return float(np.nanmedian(np.abs(x - med)))


def safe_stats(
    values: np.ndarray,
    *,
    value_range: Optional[Tuple[float, float]] = None,
    signed_or_zero_centered: bool = False,
) -> Dict[str, float]:
    # Remove NaN/Inf voxels before computing statistics.
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "n_vox": 0, "mean": np.nan, "median": np.nan, "sd": np.nan, "mad": np.nan,
            "cv": np.nan, "robust_cv": np.nan, "relative_spread": np.nan,
            "min": np.nan, "max": np.nan, "p05": np.nan, "p95": np.nan,
            "outlier_fraction_robust": np.nan, "out_of_range_fraction": np.nan,
        }

    # Classical statistics.
    mean = float(np.nanmean(values))
    sd = float(np.nanstd(values))
    med = float(np.nanmedian(values))
    mad = robust_mad(values)
    robust_sigma = 1.4826 * mad

    # Classical coefficient of variation (CV = SD / mean) is unstable
    # for maps centered around zero or containing negative values.
    #
    # Example:
    #   Chimap/QSM often has values around [-0.5, 0.5]
    #   -> mean may be close to zero
    #   -> CV explodes or becomes meaningless
    #
    # Therefore:
    #   - signed maps get cv = NaN
    #   - robust/spread metrics are preferred instead
    if signed_or_zero_centered or mean <= 0 or not np.isfinite(mean):
        cv = np.nan
    else:
        cv = float(sd / mean)

    robust_cv = float(robust_sigma / abs(med)) if abs(med) > 1e-8 else np.nan

    p05, p95 = np.nanpercentile(values, [5, 95])
    robust_range = float(p95 - p05)
    relative_spread = float(robust_sigma / robust_range) if robust_range > 0 else np.nan

    if mad > 0:
        outlier_fraction = float(np.mean(np.abs(values - med) > 5 * mad))
    else:
        outlier_fraction = 0.0

    if value_range is not None:
        lo, hi = value_range
        denom = hi - lo
        out_of_range_fraction = float(np.mean((values < lo) | (values > hi)))
        if denom > 0:
            mean_norm = float((mean - lo) / denom)
            median_norm = float((med - lo) / denom)
            sd_norm = float(sd / denom)
            mad_norm = float(mad / denom)
            p05_norm = float((p05 - lo) / denom)
            p95_norm = float((p95 - lo) / denom)
        else:
            mean_norm = median_norm = sd_norm = mad_norm = p05_norm = p95_norm = np.nan
    else:
        out_of_range_fraction = np.nan
        mean_norm = median_norm = sd_norm = mad_norm = p05_norm = p95_norm = np.nan

    return {
        "n_vox": int(values.size),
        "mean": mean,
        "median": med,
        "sd": sd,
        "mad": mad,
        "cv": cv,
        "robust_cv": robust_cv,
        "relative_spread": relative_spread,
        "mean_norm_range": mean_norm,
        "median_norm_range": median_norm,
        "sd_norm_range": sd_norm,
        "mad_norm_range": mad_norm,
        "p05_norm_range": p05_norm,
        "p95_norm_range": p95_norm,
        "min": float(np.nanmin(values)),
        "max": float(np.nanmax(values)),
        "p05": float(p05),
        "p95": float(p95),
        "outlier_fraction_robust": outlier_fraction,
        "out_of_range_fraction": out_of_range_fraction,
    }


def normalize_for_display(
    data: np.ndarray,
    mask: Optional[np.ndarray] = None,
    p_low=1,
    p_high=99,
    fixed_range: Optional[Tuple[float, float]] = None,
) -> Tuple[np.ndarray, float, float]:
    arr = data.copy()
    if fixed_range is not None:
        return arr, float(fixed_range[0]), float(fixed_range[1])

    if mask is not None:
        vals = arr[(mask > 0) & np.isfinite(arr)]
    else:
        vals = arr[np.isfinite(arr)]
    if vals.size == 0:
        return arr, 0.0, 1.0
    lo, hi = np.nanpercentile(vals, [p_low, p_high])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(vals)), float(np.nanmax(vals))
        if hi <= lo:
            hi = lo + 1.0
    return arr, float(lo), float(hi)


def select_coronal_slices(mask: np.ndarray, n_slices: int = 5) -> List[int]:
    """Return coronal y-indices spanning the non-zero mask."""
    if mask is None or np.sum(mask > 0) == 0:
        return []
    coords = np.where(mask > 0)
    y_min, y_max = int(np.min(coords[1])), int(np.max(coords[1]))
    if y_max <= y_min:
        return [y_min]
    # Use conservative central slices to avoid edge/wraparound artifacts.
    # Original was 0.20-0.80; 0.25-0.70 avoids the most anterior/posterior mask edges.
    qs = np.linspace(0.25, 0.70, n_slices)
    return [int(round(y_min + q * (y_max - y_min))) for q in qs]


# =============================================================================
# Screenshot and visualization functions
# =============================================================================

def save_coronal_screenshot(
    img_path: Path,
    out_png: Path,
    title: str,
    mask_data: Optional[np.ndarray] = None,
    n_slices: int = 5,
    fixed_range: Optional[Tuple[float, float]] = None,
    colorbar_label: str = "value",
) -> None:
    img = load_img(img_path)
    data = finite_data(img)
    if data.ndim != 3:
        return

    if mask_data is not None and mask_data.shape != data.shape:
        mask_data = None

    slices = select_coronal_slices(mask_data, n_slices=n_slices) if mask_data is not None else []
    if not slices:
        slices = [int(round(data.shape[1] * q)) for q in np.linspace(0.25, 0.70, n_slices)]

    _, vmin, vmax = normalize_for_display(data, mask_data, fixed_range=fixed_range)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, len(slices), figsize=(3.1 * len(slices), 3.8), constrained_layout=True)
    if len(slices) == 1:
        axes = [axes]

    last_im = None
    for ax, y in zip(axes, slices):
        sl = np.rot90(data[:, y, :])
        last_im = ax.imshow(sl, cmap="gray", vmin=vmin, vmax=vmax)
        if mask_data is not None:
            msl = np.rot90(mask_data[:, y, :] > 0)
            try:
                ax.contour(msl.astype(float), levels=[0.5], linewidths=0.6)
            except Exception:
                pass
        ax.set_title(f"y={y}", fontsize=8)
        ax.axis("off")

    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes, fraction=0.025, pad=0.02)
        cbar.set_label(colorbar_label, fontsize=8)
        cbar.ax.tick_params(labelsize=7)

    fig.suptitle(title, fontsize=10)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)



def save_overlay_screenshot(
    bg_path: Path,
    overlay_mask: np.ndarray,
    out_png: Path,
    title: str,
    n_slices: int = 5,
    bg_fixed_range: Optional[Tuple[float, float]] = None,
    alpha: float = 0.5,
) -> None:
    """Save red 50%-opacity mask overlay on a background image."""
    img = load_img(bg_path)
    bg = finite_data(img)
    if bg.ndim != 3 or overlay_mask.shape != bg.shape:
        return

    slices = select_coronal_slices(overlay_mask, n_slices=n_slices)
    if not slices:
        slices = [int(round(bg.shape[1] * q)) for q in np.linspace(0.25, 0.70, n_slices)]

    _, vmin, vmax = normalize_for_display(bg, overlay_mask, fixed_range=bg_fixed_range)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, len(slices), figsize=(3.1 * len(slices), 3.8), constrained_layout=True)
    if len(slices) == 1:
        axes = [axes]

    last_im = None
    for ax, y in zip(axes, slices):
        bg_sl = np.rot90(bg[:, y, :])
        ov_sl = np.rot90(overlay_mask[:, y, :] > 0)

        last_im = ax.imshow(bg_sl, cmap="gray", vmin=vmin, vmax=vmax)

        rgba = np.zeros((*ov_sl.shape, 4), dtype=float)
        rgba[..., 0] = 1.0  # red
        rgba[..., 3] = ov_sl.astype(float) * alpha
        ax.imshow(rgba)

        ax.set_title(f"y={y}", fontsize=8)
        ax.axis("off")

    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes, fraction=0.025, pad=0.02)
        cbar.set_label("R1", fontsize=8)
        cbar.ax.tick_params(labelsize=7)

    fig.suptitle(title, fontsize=10)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def make_histogram(
    values: np.ndarray,
    bins: int = 80,
    value_range: Optional[Tuple[float, float]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    # Remove NaN/Inf voxels before computing statistics.
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.array([]), np.array([])

    if value_range is not None:
        lo, hi = value_range
        # For map QC, histogram only the plausible display range. This prevents
        # extreme noise/background values from dominating the plot.
        values = values[(values >= lo) & (values <= hi)]
        if values.size == 0:
            return np.array([]), np.array([])
    else:
        lo, hi = np.nanpercentile(values, [0.5, 99.5])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
            if hi <= lo:
                hi = lo + 1.0

    counts, edges = np.histogram(values, bins=bins, range=(lo, hi), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts


def scanT_label(value) -> str:
    """Normalize scanT labels for filenames and grouping.

    Examples:
        3, 3.0, "3T" -> "3T"
        7, 7.0, "7T" -> "7T"
    """
    s = str_or_empty(value).strip()
    if not s:
        return "unknownT"
    s = s.replace(".0", "")
    if s.lower().endswith("t"):
        return s
    return f"{s}T"


def safe_plot_filename(text: str) -> str:
    """Make a safe string for plot filenames."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def iter_scanT_groups(df: pd.DataFrame):
    """Yield (scanT_label, subset) pairs. Keeps unknowns separate."""
    if df.empty or "scanT" not in df.columns:
        yield "unknownT", df
        return
    tmp = df.copy()
    tmp["_scanT_label"] = tmp["scanT"].apply(scanT_label)
    for label in sorted(tmp["_scanT_label"].dropna().unique()):
        yield str(label), tmp[tmp["_scanT_label"] == label].drop(columns=["_scanT_label"])



def plot_histogram_png(hist_df: pd.DataFrame, out_png: Path, title: str, max_lines: int = 30) -> None:
    if hist_df.empty:
        return
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))

    plotted = 0
    for _, row in hist_df.iterrows():
        centers = np.array(json.loads(row["bin_centers"]))
        density = np.array(json.loads(row["density"]))
        if centers.size == 0:
            continue
        label = f'{row["subject"]}_{row["session"]}'
        ax.plot(centers, density, linewidth=0.9, alpha=0.7, label=label)
        plotted += 1
        if plotted >= max_lines:
            break

    ax.set_title(title)
    ax.set_xlabel("value")
    ax.set_ylabel("density")
    if plotted <= 12:
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def plot_metric_boxplots(metrics: pd.DataFrame, outdir: Path) -> None:
    """Create metric boxplots separately for each field strength.

    This avoids mixing 3T and 7T distributions, which would make CV/error
    values look like QC outliers simply because the physics differs.
    """
    if metrics.empty:
        return
    outdir.mkdir(parents=True, exist_ok=True)

    for scanT_lab, scan_df in iter_scanT_groups(metrics):
        scan_dir = outdir / scanT_lab
        scan_dir.mkdir(parents=True, exist_ok=True)

        for metric in ["median", "cv", "robust_cv", "sd_norm_range", "out_of_range_fraction", "relative_spread"]:
            if metric not in scan_df.columns:
                continue

            for region in ["brain", "WM", "GM", "CSF"]:
                sub = scan_df[(scan_df["region"] == region) & np.isfinite(scan_df[metric])]
                if sub.empty:
                    continue

                maps = sorted(sub["map"].dropna().unique())
                if not maps:
                    continue

                fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(maps)), 4))
                data = [sub.loc[sub["map"] == m, metric].dropna().values for m in maps]
                if all(len(d) == 0 for d in data):
                    plt.close(fig)
                    continue

                ax.boxplot(data, tick_labels=maps, showfliers=True)
                ax.set_title(f"{metric} by map, region={region}, {scanT_lab}")
                ax.set_ylabel(metric)
                ax.tick_params(axis="x", rotation=45)
                fig.tight_layout()
                fig.savefig(scan_dir / f"boxplot_{metric}_{region}_{scanT_lab}.png", dpi=140)
                plt.close(fig)



def plot_heatmap(metrics: pd.DataFrame, out_png: Path, region: str = "brain", metric: str = "cv") -> None:
    """Create scanT-separated heatmaps.

    A single pooled heatmap can be misleading because 3T and 7T have different
    expected distributions. This function therefore writes:
        <stem>_3T.png
        <stem>_7T.png
    instead of one mixed heatmap.
    """
    if metrics.empty or metric not in metrics.columns:
        return

    for scanT_lab, scan_df in iter_scanT_groups(metrics):
        sub = scan_df[(scan_df["region"] == region) & np.isfinite(scan_df[metric])].copy()
        if sub.empty:
            continue

        sub["sub_ses"] = sub["subject"].astype(str) + "_" + sub["session"].astype(str)
        pivot = sub.pivot_table(index="sub_ses", columns="map", values=metric, aggfunc="median")
        if pivot.empty:
            continue

        # z-score per map column within this scanT group only.
        z = pivot.copy()
        for col in z.columns:
            vals = z[col].astype(float)
            sd = vals.std(skipna=True)
            if sd and np.isfinite(sd):
                z[col] = (vals - vals.mean(skipna=True)) / sd
            else:
                z[col] = np.nan

        out_png.parent.mkdir(parents=True, exist_ok=True)
        scan_out = out_png.with_name(f"{out_png.stem}_{scanT_lab}{out_png.suffix}")

        fig, ax = plt.subplots(figsize=(max(7, 0.7 * len(z.columns)), max(5, 0.25 * len(z.index))))
        im = ax.imshow(z.values, aspect="auto", interpolation="nearest")
        ax.set_xticks(range(len(z.columns)))
        ax.set_xticklabels(z.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(z.index)))
        ax.set_yticklabels(z.index, fontsize=6)
        ax.set_title(f"z-scored {metric}, region={region}, {scanT_lab}")
        fig.colorbar(im, ax=ax, label=f"{metric} z-score within {scanT_lab}")
        fig.tight_layout()
        fig.savefig(scan_out, dpi=150)
        plt.close(fig)




# ---------------------------------------------------------------------
# Main QC logic
# ---------------------------------------------------------------------

# =============================================================================
# Subject/session parsing and file discovery
# =============================================================================

def get_subject(row: pd.Series) -> str:
    return str_or_empty(row.get("bids_subjID")) or str_or_empty(row.get("subject"))


def get_session(row: pd.Series) -> str:
    return str_or_empty(row.get("bids_sesID")) or str_or_empty(row.get("session"))


def get_dirs(row: pd.Series, replacements: List[Tuple[str, str]]) -> Dict[str, Optional[Path]]:
    keys = ["qMRI", "QSMxT", "LCPCA_distCorr", "relax_R2", "r2prime", "nighres"]
    out: Dict[str, Optional[Path]] = {}
    for key in keys:
        out[key] = existing_dir(apply_path_replacements(str_or_empty(row.get(key)), replacements))
    return out


def find_maps_for_row(row: pd.Series, replacements: List[Tuple[str, str]]) -> Tuple[List[FoundFile], Dict[str, Path]]:
    subject = get_subject(row)
    session = get_session(row)
    scanT = str_or_empty(row.get("scanT"))
    scanner = str_or_empty(row.get("scanner"))
    dirs = get_dirs(row, replacements)

    found: List[FoundFile] = []
    chosen: Dict[str, Path] = {}

    def add(category: str, name: str, path: Optional[Path], matches: List[Path], status: str):
        if path is not None:
            chosen[f"{category}:{name}"] = path
        found.append(
            FoundFile(
                subject=subject,
                session=session,
                scanT=scanT,
                scanner=scanner,
                category=category,
                name=name,
                path=str(path) if path else "",
                status=status,
                n_matches=len(matches),
                all_matches=";".join(str(p) for p in matches),
            )
        )

    # qMRI maps
    for name in ["R1", "R2star", "PD", "MTsat"]:
        path, matches, status = find_single(dirs["qMRI"], MAP_PATTERNS[name])
        add("map", name, path, matches, status)

    # B1 maps in qMRI/Supplementary
    qsup = dirs["qMRI"] / "Supplementary" if dirs["qMRI"] else None
    for name in ["B1_STX", "B1_PTX"]:
        path, matches, status = find_single(qsup, MAP_PATTERNS[name])
        add("map", name, path, matches, status)

    # QSM Chimap: prefer pdf_phaseMask/subject/session/anat/transform_to_orig if present.
    qsm_dir = dirs["QSMxT"]
    qsm_candidates = []
    if qsm_dir:
        # If CSV points at QSMxT/sub/ses/anat, also check sibling pdf_phaseMask layout.
        p = qsm_dir
        # Direct transform_to_orig.
        qsm_candidates.append(p / "transform_to_orig")
        qsm_candidates.append(p)

        # Attempt to build pdf_phaseMask path from QSMxT root.
        parts = list(p.parts)
        if "QSMxT" in parts:
            idx = parts.index("QSMxT")
            # Find sub and ses in path.
            sub = subject
            ses = session
            root = Path(*parts[:idx + 1])
            qsm_candidates.insert(0, root / "pdf_phaseMask" / sub / ses / "anat" / "transform_to_orig")
            qsm_candidates.insert(1, root / "pdf_phaseMask" / sub / ses / "anat")

    qsm_path = None
    qsm_matches: List[Path] = []
    qsm_status = "missing"
    for cand in qsm_candidates:
        path, matches, status = find_single(
            cand,
            MAP_PATTERNS["Chimap"],
            exclude=["singlepass", "minIP", "swi", "R2star"],
            prefer=["Chimap"],
        )
        if path is not None:
            qsm_path, qsm_matches, qsm_status = path, matches, status
            break
    add("map", "Chimap", qsm_path, qsm_matches, qsm_status)

    # relax_R2 maps
    for name in ["T2", "R2", "TB1_relax"]:
        path, matches, status = find_single(dirs["relax_R2"], MAP_PATTERNS[name])
        add("map", name, path, matches, status)

    # r2prime
    path, matches, status = find_single(dirs["r2prime"], MAP_PATTERNS["R2prime"])
    add("map", "R2prime", path, matches, status)

    # Error maps: qMRI supplementary plus relax_R2.
    for err_name, pats in ERROR_PATTERNS.items():
        # Avoid generic_error duplicating param_error/errorESTATICS too much: still useful for manifest.
        path, matches, status = find_single(qsup, pats)
        if path:
            src = infer_error_source_from_name(path.name)
            add("error", f"qMRI_{src}_{err_name}", path, matches, status)

    if dirs["relax_R2"]:
        for err_name, pats in ERROR_PATTERNS.items():
            path, matches, status = find_single(dirs["relax_R2"], pats)
            if path:
                src = infer_error_source_from_name(path.name)
                add("error", f"relax_{src}_{err_name}", path, matches, status)

    # Weighted images for screenshots/mask creation reference; not used in stats by default.
    for name, pats in WEIGHTED_PATTERNS.items():
        path, matches, status = find_single(dirs["LCPCA_distCorr"], pats)
        add("weighted", name, path, matches, status)

    # Brain mask: prefer user-created/synthstrip-like masks, fallback to qMRI bmask.
    mask_search_dirs = [
        dirs["LCPCA_distCorr"],
        dirs["qMRI"],
        qsup,
        dirs["relax_R2"],
    ]
    brain_patterns = [
        "*brainmask*.nii", "*brainmask*.nii.gz",
        "*brain_mask*.nii", "*brain_mask*.nii.gz",
        "*desc-brain_mask*.nii", "*desc-brain_mask*.nii.gz",
        "*mask*.nii", "*mask*.nii.gz",
        "bmask*.nii", "bmask*.nii.gz",
    ]
    brain_path = None
    brain_matches: List[Path] = []
    brain_status = "missing"
    for md in mask_search_dirs:
        path, matches, status = find_single(
            md,
            brain_patterns,
            exclude=["seg", ".mat", "c1", "c2", "c3", "c4", "c5"],
            prefer=["brain", "mask"],
        )
        if path:
            brain_path, brain_matches, brain_status = path, matches, status
            break
    add("mask", "brain", brain_path, brain_matches, brain_status)

    # SPM tissue classes
    mpmcalc = qsup / "MPMCalc" if qsup else None
    tissue_specs = {
        "GM": ["c1*.nii", "c1*.nii.gz"],
        "WM": ["c2*.nii", "c2*.nii.gz"],
        "CSF": ["c3*.nii", "c3*.nii.gz"],
    }
    for tissue, pats in tissue_specs.items():
        path, matches, status = find_single(mpmcalc, pats)
        add("mask", tissue, path, matches, status)

    return found, chosen


def load_regions_for_image(
    img: nib.Nifti1Image,
    manifest_for_row: pd.DataFrame,
    tissue_threshold: float,
) -> Dict[str, np.ndarray]:
    regions: Dict[str, np.ndarray] = {}

    # Helper: retrieve the path for a specific mask from the manifest.
    def get_path(cat: str, name: str) -> Optional[Path]:
        sub = manifest_for_row[
            (manifest_for_row["category"] == cat)
            & (manifest_for_row["name"] == name)
            & (manifest_for_row["path"].astype(str) != "")
        ]
        if sub.empty:
            return None
        return Path(sub.iloc[0]["path"])

    # Brain mask:
    # Prefer externally generated masks (e.g. SynthStrip),
    # because they are often more robust than simple thresholding.
    brain_path = get_path("mask", "brain")
    if brain_path and brain_path.exists():
        brain = resample_mask_to_img(load_img(brain_path), img, order=0)
        regions["brain"] = brain > 0
    else:
        data = finite_data(img)
        regions["brain"] = np.isfinite(data) & (data != 0)

    for region in ["GM", "WM", "CSF"]:
        p = get_path("mask", region)
        if p and p.exists():
            # Tissue masks are probability maps from SPM.
            # We use linear interpolation during resampling.
            prob = resample_mask_to_img(load_img(p), img, order=1)
            regions[region] = prob > tissue_threshold

    return regions


# =============================================================================
# B1 homogeneity analysis
# =============================================================================

def compute_b1_homogeneity(
    overview: pd.DataFrame,
    manifest: pd.DataFrame,
    tissue_threshold: float,
) -> pd.DataFrame:
    """Compute B1-specific homogeneity metrics for STX/PTX maps.

    Metrics are computed in brain/WM/GM/CSF regions. The most useful ones are:
      - b1_cv: SD / mean
      - b1_robust_cv: 1.4826*MAD / median
      - b1_p95_p05: robust spatial range
      - b1_homogeneity_index: 1 - SD/mean
      - b1_robust_homogeneity_index: 1 - (p95-p05)/(2*median)
      - b1_left_right_asymmetry: abs(left_mean-right_mean)/whole_region_mean

    These are intended mainly for STX/PTX comparisons, not pass/fail thresholds yet.
    """
    rows = []

    for _, row in tqdm(overview.iterrows(), total=len(overview), desc="Computing B1 homogeneity"):
        subject = get_subject(row)
        session = get_session(row)
        scanT = str_or_empty(row.get("scanT"))
        scanner = str_or_empty(row.get("scanner"))

        row_manifest = manifest[(manifest["subject"] == subject) & (manifest["session"] == session)]
        b1_entries = row_manifest[
            (row_manifest["category"] == "map")
            & (row_manifest["name"].isin(["B1_STX", "B1_PTX", "TB1_relax"]))
            & (row_manifest["path"].astype(str) != "")
        ]

        for _, ent in b1_entries.iterrows():
            b1_name = ent["name"]
            p = Path(ent["path"])
            if not p.exists():
                continue

            try:
                img = load_img(p)
                data = finite_data(img)
                if data.ndim != 3:
                    continue
                regions = load_regions_for_image(img, row_manifest, tissue_threshold)
            except Exception as exc:
                sys.stderr.write(f"WARNING: B1 homogeneity failed loading {p}: {exc}\\n")
                continue

            for region_name, region_mask in regions.items():
                if region_mask.shape != data.shape:
                    continue

                valid_mask = region_mask & np.isfinite(data) & (data > 0)
                values = data[valid_mask]
                if values.size == 0:
                    continue

                # Classical statistics.
                mean = float(np.nanmean(values))
                median = float(np.nanmedian(values))
                sd = float(np.nanstd(values))
                mad = robust_mad(values)
                robust_sigma = 1.4826 * mad
                p05, p95 = np.nanpercentile(values, [5, 95])
                p95_p05 = float(p95 - p05)

                b1_cv = float(sd / mean) if mean > 0 else np.nan
                b1_robust_cv = float(robust_sigma / median) if median > 0 else np.nan
                b1_hi = float(1.0 - b1_cv) if np.isfinite(b1_cv) else np.nan
                b1_hi_robust = float(1.0 - (p95_p05 / (2.0 * median))) if median > 0 else np.nan

                # Approximate left/right split in voxel x-dimension.
                #
                # NOTE:
                # This assumes images are approximately neurologically aligned.
                # This is intended as a QC heuristic, not a formal anatomical metric.
                # This assumes images are not wildly reoriented; it is QC only.
                x_mid = data.shape[0] // 2
                left_mask = valid_mask.copy()
                right_mask = valid_mask.copy()
                left_mask[x_mid:, :, :] = False
                right_mask[:x_mid, :, :] = False

                if np.sum(left_mask) > 0 and np.sum(right_mask) > 0 and mean > 0:
                    left_mean = float(np.nanmean(data[left_mask]))
                    right_mean = float(np.nanmean(data[right_mask]))
                    lr_asym = float(abs(left_mean - right_mean) / mean)
                else:
                    left_mean = right_mean = lr_asym = np.nan

                # AP and IS gradients: slope normalized by mean across voxel index.
                # These are crude but useful for systematic B1 shading.
                # Estimate systematic spatial gradients across the brain.
                #
                # Example:
                #   inferior-superior shading
                #   anterior-posterior RF inhomogeneity
                #
                # We fit a simple linear slope across slice means.
                def normalized_axis_slope(axis: int) -> float:
                    prof = []
                    idxs = []
                    for i in range(data.shape[axis]):
                        slicer = [slice(None)] * 3
                        slicer[axis] = i
                        sl_mask = valid_mask[tuple(slicer)]
                        if np.sum(sl_mask) < 10:
                            continue
                        sl_data = data[tuple(slicer)]
                        prof.append(float(np.nanmean(sl_data[sl_mask])))
                        idxs.append(i)
                    if len(prof) < 5 or mean <= 0:
                        return np.nan
                    x = np.asarray(idxs, dtype=float)
                    y = np.asarray(prof, dtype=float)
                    x = (x - np.nanmean(x)) / max(np.nanstd(x), 1e-8)
                    slope = np.polyfit(x, y, 1)[0]
                    return float(slope / mean)

                rows.append({
                    "subject": subject,
                    "session": session,
                    "scanT": scanT,
                    "scanner": scanner,
                    "b1_map": b1_name,
                    "region": region_name,
                    "path": str(p),
                    "n_vox": int(values.size),
                    "mean": mean,
                    "median": median,
                    "sd": sd,
                    "mad": mad,
                    "p05": float(p05),
                    "p95": float(p95),
                    "p95_p05": p95_p05,
                    "b1_cv": b1_cv,
                    "b1_robust_cv": b1_robust_cv,
                    "b1_homogeneity_index": b1_hi,
                    "b1_robust_homogeneity_index": b1_hi_robust,
                    "left_mean": left_mean,
                    "right_mean": right_mean,
                    "b1_left_right_asymmetry": lr_asym,
                    "b1_ap_gradient_norm": normalized_axis_slope(axis=1),
                    "b1_is_gradient_norm": normalized_axis_slope(axis=2),
                })

    return pd.DataFrame(rows)


def plot_b1_homogeneity(b1: pd.DataFrame, outdir: Path) -> None:
    """Create B1 homogeneity plots separately for each field strength.

    B1 physics differs strongly between 3T and 7T, and STX/PTX comparisons are
    mainly meaningful within the same scanT group. This avoids mixed 3T/7T plots.
    """
    if b1.empty:
        return
    outdir.mkdir(parents=True, exist_ok=True)

    metrics = [
        "b1_cv",
        "b1_robust_cv",
        "b1_homogeneity_index",
        "b1_robust_homogeneity_index",
        "b1_left_right_asymmetry",
        "b1_ap_gradient_norm",
        "b1_is_gradient_norm",
    ]

    for scanT_lab, scan_df in iter_scanT_groups(b1):
        scan_dir = outdir / scanT_lab
        scan_dir.mkdir(parents=True, exist_ok=True)

        for region in ["brain", "WM", "GM", "CSF"]:
            sub = scan_df[scan_df["region"] == region]
            if sub.empty:
                continue

            for metric in metrics:
                if metric not in sub.columns:
                    continue
                vals = sub[np.isfinite(sub[metric])]
                if vals.empty:
                    continue

                labels = sorted(vals["b1_map"].dropna().unique())
                if not labels:
                    continue

                fig, ax = plt.subplots(figsize=(max(6, 1.5 * len(labels)), 4))
                data = [vals.loc[vals["b1_map"] == lab, metric].dropna().values for lab in labels]
                if all(len(d) == 0 for d in data):
                    plt.close(fig)
                    continue

                ax.boxplot(data, tick_labels=labels, showfliers=True)
                ax.set_title(f"B1 {metric}, region={region}, {scanT_lab}")
                ax.set_ylabel(metric)
                ax.tick_params(axis="x", rotation=30)
                fig.tight_layout()
                fig.savefig(scan_dir / f"b1_{metric}_{region}_{scanT_lab}.png", dpi=140)
                plt.close(fig)

        # STX vs PTX paired scatter within scanT only.
        for region in ["brain", "WM"]:
            sub = scan_df[scan_df["region"] == region]
            for metric in ["b1_cv", "b1_robust_cv", "b1_left_right_asymmetry", "b1_homogeneity_index"]:
                if metric not in sub.columns:
                    continue

                pivot = sub.pivot_table(
                    index=["subject", "session"],
                    columns="b1_map",
                    values=metric,
                    aggfunc="median",
                )
                if "B1_STX" not in pivot.columns or "B1_PTX" not in pivot.columns:
                    continue
                paired = pivot[["B1_STX", "B1_PTX"]].dropna()
                if paired.empty:
                    continue

                fig, ax = plt.subplots(figsize=(4.5, 4.5))
                ax.scatter(paired["B1_STX"], paired["B1_PTX"])
                mn = float(np.nanmin(paired.values))
                mx = float(np.nanmax(paired.values))
                if np.isfinite(mn) and np.isfinite(mx) and mx > mn:
                    ax.plot([mn, mx], [mn, mx], linestyle="--", linewidth=1)
                ax.set_xlabel(f"STX {metric}")
                ax.set_ylabel(f"PTX {metric}")
                ax.set_title(f"STX vs PTX {metric}, {region}, {scanT_lab}")
                fig.tight_layout()
                fig.savefig(scan_dir / f"b1_STX_vs_PTX_{metric}_{region}_{scanT_lab}.png", dpi=140)
                plt.close(fig)




# =============================================================================
# Quantitative QC metric computation
# =============================================================================

def process_metrics(
    overview: pd.DataFrame,
    manifest: pd.DataFrame,
    tissue_threshold: float,
    histogram_bins: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    hist_rows = []
    error_rows = []

    for _, row in tqdm(overview.iterrows(), total=len(overview), desc="Computing metrics"):
        subject = get_subject(row)
        session = get_session(row)
        scanT = str_or_empty(row.get("scanT"))
        scanner = str_or_empty(row.get("scanner"))

        row_manifest = manifest[(manifest["subject"] == subject) & (manifest["session"] == session)]
        image_entries = row_manifest[
            (row_manifest["category"].isin(["map", "error"]))
            & (row_manifest["path"].astype(str) != "")
        ]

        # Iterate through all discovered maps and error maps.
        for _, ent in image_entries.iterrows():
            path = Path(ent["path"])
            if not path.exists():
                continue

            try:
                img = load_img(path)
                data = finite_data(img)
                if data.ndim != 3:
                    continue
                regions = load_regions_for_image(img, row_manifest, tissue_threshold)
            except Exception as exc:
                sys.stderr.write(f"WARNING: failed loading {path}: {exc}\n")
                continue

            category = ent["category"]
            map_name = ent["name"]

            for region_name, region_mask in regions.items():
                if region_mask.shape != data.shape:
                    continue
                # Extract voxel values inside the current region.
                values = data[region_mask & np.isfinite(data)]
                value_range = MAP_VALUE_RANGES.get(map_name)
                stats = safe_stats(
                    values,
                    value_range=value_range if category == "map" else None,
                    signed_or_zero_centered=map_name in SIGNED_OR_ZERO_CENTERED_MAPS,
                )
                common = {
                    "subject": subject,
                    "session": session,
                    "scanT": scanT,
                    "scanner": scanner,
                    "category": category,
                    "map": map_name,
                    "region": region_name,
                    "path": str(path),
                }
                out = {**common, **stats}

                if category == "error":
                    error_rows.append(out)
                else:
                    metric_rows.append(out)

                centers, density = make_histogram(
                    values,
                    bins=histogram_bins,
                    value_range=value_range if category == "map" else None,
                )
                hist_rows.append({
                    **common,
                    "bin_centers": json.dumps(np.round(centers.astype(float), 8).tolist()),
                    "density": json.dumps(np.round(density.astype(float), 8).tolist()),
                })

    return pd.DataFrame(metric_rows), pd.DataFrame(error_rows), pd.DataFrame(hist_rows)


def add_outlier_flags(metrics: pd.DataFrame, z_thresh: float = 3.0) -> pd.DataFrame:
    if metrics.empty:
        return metrics
    out = metrics.copy()
    candidate_cols = [
        "median",
        "cv",
        "robust_cv",
        "p95",
        "sd_norm_range",
        "out_of_range_fraction",
        "relative_spread",
    ]
    cols = [c for c in candidate_cols if c in out.columns]

    for col in cols:
        out[f"{col}_z"] = np.nan
        out[f"{col}_outlier"] = False

    # Crucial: outliers are computed within scanT, not across mixed 3T/7T.
    group_cols = ["category", "map", "region", "scanT"]
    for _, idx in out.groupby(group_cols, dropna=False).groups.items():
        idx = list(idx)
        for col in cols:
            vals = out.loc[idx, col].astype(float)
            sd = vals.std(skipna=True)
            if sd and np.isfinite(sd):
                z = (vals - vals.mean(skipna=True)) / sd
                out.loc[idx, f"{col}_z"] = z
                out.loc[idx, f"{col}_outlier"] = np.abs(z) >= z_thresh
    return out


def make_segmentation_overlays(
    overview: pd.DataFrame,
    manifest: pd.DataFrame,
    outdir: Path,
    tissue_threshold: float,
    n_slices: int,
) -> None:
    """
    Create segmentation overlay screenshots.

    Overlays:
        - brain mask
        - GM
        - WM
        - CSF

    Background:
        R1 map

    Display:
        red overlay at 50% opacity

    Purpose:
        Quick visual QC of segmentation quality and alignment.
    """
    overlay_dir = outdir / "figures" / "segmentation_overlays"
    for _, row in tqdm(overview.iterrows(), total=len(overview), desc="Making segmentation overlays"):
        subject = get_subject(row)
        session = get_session(row)
        row_manifest = manifest[(manifest["subject"] == subject) & (manifest["session"] == session)]

        r1 = row_manifest[
            (row_manifest["category"] == "map")
            & (row_manifest["name"] == "R1")
            & (row_manifest["path"].astype(str) != "")
        ]
        if r1.empty:
            continue

        r1_path = Path(r1.iloc[0]["path"])
        if not r1_path.exists():
            continue

        try:
            r1_img = load_img(r1_path)
            regions = load_regions_for_image(r1_img, row_manifest, tissue_threshold)
        except Exception as exc:
            sys.stderr.write(f"WARNING: could not load R1/regions for overlays {subject} {session}: {exc}\\n")
            continue

        for region_name in ["brain", "GM", "WM", "CSF"]:
            region = regions.get(region_name)
            if region is None:
                continue
            out_png = overlay_dir / subject / session / f"{subject}_{session}_R1_overlay_{region_name}_red50.png"
            try:
                save_overlay_screenshot(
                    r1_path,
                    region.astype(bool),
                    out_png,
                    title=f"{subject} {session} R1 + {region_name} mask, red 50%",
                    n_slices=n_slices,
                    bg_fixed_range=MAP_VALUE_RANGES.get("R1"),
                    alpha=0.5,
                )
            except Exception as exc:
                sys.stderr.write(f"WARNING: overlay failed for {subject} {session} {region_name}: {exc}\\n")


def make_screenshots(
    overview: pd.DataFrame,
    manifest: pd.DataFrame,
    outdir: Path,
    maps_to_plot: Sequence[str],
    tissue_threshold: float,
    n_slices: int,
) -> None:
    shot_dir = outdir / "figures" / "screenshots"
    for _, row in tqdm(overview.iterrows(), total=len(overview), desc="Making screenshots"):
        subject = get_subject(row)
        session = get_session(row)
        row_manifest = manifest[(manifest["subject"] == subject) & (manifest["session"] == session)]

        for map_name in maps_to_plot:
            sub = row_manifest[
                (row_manifest["category"] == "map")
                & (row_manifest["name"] == map_name)
                & (row_manifest["path"].astype(str) != "")
            ]
            if sub.empty:
                continue
            p = Path(sub.iloc[0]["path"])
            if not p.exists():
                continue
            try:
                img = load_img(p)
                regions = load_regions_for_image(img, row_manifest, tissue_threshold)
                brain = regions.get("brain")
                out_png = shot_dir / subject / session / f"{subject}_{session}_{map_name}.png"
                save_coronal_screenshot(
                    p,
                    out_png,
                    title=f"{subject} {session} {map_name}",
                    mask_data=brain.astype(float) if brain is not None else None,
                    n_slices=n_slices,
                    fixed_range=MAP_VALUE_RANGES.get(map_name),
                    colorbar_label=map_name,
                )
            except Exception as exc:
                sys.stderr.write(f"WARNING: screenshot failed for {p}: {exc}\n")

        # Error map screenshots too, because they are low-hanging fruit.
        errors = row_manifest[
            (row_manifest["category"] == "error")
            & (row_manifest["path"].astype(str) != "")
        ]
        for _, ent in errors.iterrows():
            p = Path(ent["path"])
            if not p.exists():
                continue
            try:
                img = load_img(p)
                regions = load_regions_for_image(img, row_manifest, tissue_threshold)
                brain = regions.get("brain")
                safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", ent["name"])
                out_png = shot_dir / subject / session / f"{subject}_{session}_{safe_name}.png"
                save_coronal_screenshot(
                    p,
                    out_png,
                    title=f"{subject} {session} {ent['name']}",
                    mask_data=brain.astype(float) if brain is not None else None,
                    n_slices=n_slices,
                    fixed_range=None,
                    colorbar_label=str(ent["name"]),
                )
            except Exception as exc:
                sys.stderr.write(f"WARNING: error screenshot failed for {p}: {exc}\n")


def make_histogram_figures(hist: pd.DataFrame, outdir: Path) -> None:
    fig_dir = outdir / "figures" / "histograms"
    if hist.empty:
        return

    # Main cohort overlays for brain region, separately by map.
    for map_name in sorted(hist["map"].dropna().unique()):
        sub = hist[(hist["map"] == map_name) & (hist["region"] == "brain") & (hist["category"] == "map")]
        if sub.empty:
            continue
        plot_histogram_png(sub, fig_dir / f"hist_overlay_brain_{map_name}.png", f"{map_name} brain histogram overlay")

        # 3T/7T separated.
        for scanT in sorted(sub["scanT"].dropna().astype(str).unique()):
            st = sub[sub["scanT"].astype(str) == scanT]
            if len(st) > 0:
                plot_histogram_png(st, fig_dir / f"hist_overlay_brain_{map_name}_{scanT}T.png",
                                   f"{map_name} brain histogram overlay, {scanT}T")

    # Error histograms, separated by scanT.
    for err_name in sorted(hist.loc[hist["category"] == "error", "map"].dropna().unique()):
        sub = hist[(hist["map"] == err_name) & (hist["region"] == "brain") & (hist["category"] == "error")]
        if sub.empty:
            continue
        for scanT_lab, st in iter_scanT_groups(sub):
            if st.empty:
                continue
            plot_histogram_png(
                st,
                fig_dir / f"hist_overlay_error_brain_{err_name}_{scanT_lab}.png",
                f"{err_name} brain error histogram overlay, {scanT_lab}",
            )



# =============================================================================
# Integrated per-file FSLeyes scene browser
# =============================================================================

def valid_path_string(x) -> str:
    """Return usable path string or empty string."""
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    s = str(x)
    if s.lower() in {"nan", "none", "nat", ""}:
        return ""
    return s


def make_safe_filename(text: str, max_len: int = 180) -> str:
    """Make a compact filename-safe string."""
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")
    safe = re.sub(r"_+", "_", safe)
    return safe[:max_len] if len(safe) > max_len else safe


def quote_cmd_arg(x) -> str:
    """Shell-quote one command argument for display/copying."""
    return shlex.quote(str(x))


def fsleyes_args_for_overlay(name: str, path: Path, *, as_context: bool = False) -> List[str]:
    """Return FSLeyes arguments for exactly one grayscale image.

    Manual QC browser rule:
        - one scene = one file
        - no masks
        - no tissue overlays
        - no context overlays
        - grayscale only

    R2star is deliberately left without a fixed display range here because
    3T/7T and slice-dependent contrast made the fixed 0-100 range unhelpful.
    """
    args = [
        str(path),
        "--name", name,
        "--overlayType", "volume",
        "--cmap", "greyscale",
    ]

    # Apply only robustly useful ranges. Leave R2star autoscaled in FSLeyes.
    browser_ranges = {
        "R1": (0.0, 2.0),
        "PD": (0.0, 120.0),
        "MTsat": (0.0, 2.0),
        "Chimap": (-0.5, 0.5),
        "T2": (0.0, 1.0),
        "R2": (0.0, 100.0),
        "R2prime": (-50.0, 100.0),
        "B1_STX": (0.0, 200.0),
        "B1_PTX": (0.0, 200.0),
        "TB1_relax": (0.0, 200.0),
    }

    if name in browser_ranges:
        lo, hi = browser_ranges[name]
        args += ["--displayRange", str(lo), str(hi)]

    return args



def browser_overlay_name(category: str, name: str) -> str:
    """Convert manifest category/name to a concise FSLeyes overlay label."""
    if category == "weighted":
        return name.replace("_undistorted", "")
    if category == "mask":
        return f"mask_{name}"
    return name


def manifest_valid_file_rows(manifest: pd.DataFrame) -> pd.DataFrame:
    """Return manifest rows that correspond to actual files on disk."""
    rows = manifest.copy()
    if "path" not in rows.columns:
        return rows.iloc[0:0].copy()

    rows["path"] = rows["path"].apply(valid_path_string)
    rows = rows[rows["path"] != ""].copy()
    rows = rows[rows["path"].apply(lambda p: Path(p).exists())].copy()
    # Manual browser should contain only actual image targets.
    # No masks/tissue-class segmentations are included here because each generated
    # FSLeyes script should open exactly one baseline image in grayscale.
    rows = rows[rows["category"].isin(["map", "error", "weighted"])].copy()

    return rows


def context_masks_for_row(row: pd.Series, manifest: pd.DataFrame) -> List[Tuple[str, Path]]:
    """Find same-subject/session masks to add as context overlays."""
    subject = str(row["subject"])
    session = str(row["session"])

    mask_rows = manifest[
        (manifest["subject"].astype(str) == subject)
        & (manifest["session"].astype(str) == session)
        & (manifest["category"].astype(str) == "mask")
    ].copy()

    out = []
    for _, m in mask_rows.iterrows():
        p_str = valid_path_string(m.get("path", ""))
        if not p_str:
            continue
        p = Path(p_str)
        if not p.exists():
            continue
        label = browser_overlay_name("mask", str(m.get("name", "")))
        out.append((label, p))
    return out


def create_single_file_fsleyes_script(
    row: pd.Series,
    manifest: pd.DataFrame,
    out_script: Path,
    *,
    include_context_masks: bool = False,
) -> None:
    """Create one FSLeyes launcher script for one manifest file."""
    subject = str(row["subject"])
    session = str(row["session"])
    scanT = str(row.get("scanT", ""))
    scanner = str(row.get("scanner", ""))
    category = str(row["category"])
    name = str(row["name"])
    path = Path(valid_path_string(row["path"]))
    overlay_name = browser_overlay_name(category, name)

    args = ["fsleyes", "--scene", "ortho"]
    args.extend(fsleyes_args_for_overlay(overlay_name, path))

    # No extra overlays. This script opens exactly one image.

    cmd = " ".join(quote_cmd_arg(a) for a in args)

    out_script.parent.mkdir(parents=True, exist_ok=True)
    script_text = f'''#!/bin/bash
# Auto-generated single-file HistoPark FSLeyes scene
# Subject: {subject}
# Session: {session}
# Field strength: {scanT}
# Scanner: {scanner}
# Category: {category}
# Name: {name}
# File: {path}
#
# Run:
#   bash {out_script.name}

set -euo pipefail

{cmd}
'''
    out_script.write_text(script_text)
    out_script.chmod(out_script.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)


def write_tk_file_browser_script(tk_script_path: Path, scene_manifest_path: Path) -> None:
    """Write standalone Tkinter launcher for generated FSLeyes scripts."""
    browser_code = '''#!/usr/bin/env python3
"""
Standalone HistoPark file-scene browser.

Shows one long searchable list of FSLeyes scene scripts.
Double-click or press "Run selected" to launch the selected script.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox
import pandas as pd


MANIFEST = Path(__file__).resolve().parent / "__SCENE_MANIFEST_NAME__"


class FileSceneBrowser:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("HistoPark FSLeyes File Scene Browser")
        self.root.geometry("1350x780")

        self.df = pd.read_csv(MANIFEST)
        self.df["search"] = self.df["search"].fillna("").astype(str)
        self.df["display"] = self.df["display"].fillna("").astype(str)
        self.filtered_idx = list(range(len(self.df)))
        self.current_idx = 0

        self._build_ui()
        self._refresh()

    def _build_ui(self):
        top = ttk.Frame(self.root, padding=8)
        top.pack(fill="x")

        ttk.Label(top, text="Search:").pack(side="left")
        self.search_var = tk.StringVar()
        self.search_var.trace_add("write", lambda *_: self._apply_filter())
        search = ttk.Entry(top, textvariable=self.search_var, width=85)
        search.pack(side="left", padx=6)

        ttk.Button(top, text="Previous", command=self.previous_item).pack(side="left", padx=4)
        ttk.Button(top, text="Next", command=self.next_item).pack(side="left", padx=4)
        ttk.Button(top, text="Run selected", command=self.run_selected).pack(side="left", padx=12)
        ttk.Button(top, text="Open script folder", command=self.open_script_folder).pack(side="left", padx=4)

        body = ttk.Frame(self.root, padding=(8, 0, 8, 8))
        body.pack(fill="both", expand=True)

        self.listbox = tk.Listbox(
            body,
            selectmode=tk.EXTENDED,
            exportselection=False,
            font=("DejaVu Sans Mono", 10),
        )
        self.listbox.pack(side="left", fill="both", expand=True)
        self.listbox.bind("<Double-Button-1>", lambda _e: self.run_selected())
        self.listbox.bind("<<ListboxSelect>>", lambda _e: self._on_select())

        scrollbar = ttk.Scrollbar(body, orient="vertical", command=self.listbox.yview)
        scrollbar.pack(side="right", fill="y")
        self.listbox.configure(yscrollcommand=scrollbar.set)

        self.detail_var = tk.StringVar(value="")
        ttk.Label(self.root, textvariable=self.detail_var, anchor="w", relief="sunken").pack(fill="x")

        self.root.bind("<Return>", lambda _e: self.run_selected())
        self.root.bind("<Up>", lambda _e: self.previous_item())
        self.root.bind("<Down>", lambda _e: self.next_item())

    def _apply_filter(self):
        query = self.search_var.get().lower().strip()
        if not query:
            self.filtered_idx = list(range(len(self.df)))
        else:
            tokens = query.split()
            self.filtered_idx = [
                i for i, s in enumerate(self.df["search"].tolist())
                if all(tok in s for tok in tokens)
            ]
        self.current_idx = 0
        self._refresh()

    def _refresh(self):
        self.listbox.delete(0, tk.END)
        for i in self.filtered_idx:
            self.listbox.insert(tk.END, self.df.iloc[i]["display"])
        if self.filtered_idx:
            self.listbox.selection_set(0)
            self.listbox.see(0)
            self._on_select()
        else:
            self.detail_var.set("No matching scene scripts")

    def _selected_manifest_rows(self):
        selection = self.listbox.curselection()
        rows = []
        for pos in selection:
            if 0 <= pos < len(self.filtered_idx):
                rows.append(self.df.iloc[self.filtered_idx[pos]])
        return rows

    def _on_select(self):
        rows = self._selected_manifest_rows()
        if not rows:
            self.detail_var.set("")
            return
        if len(rows) == 1:
            row = rows[0]
            self.current_idx = int(self.listbox.curselection()[0])
            self.detail_var.set(f"{row['display']}  |  {row['file']}")
        else:
            self.detail_var.set(f"{len(rows)} scene scripts selected")

    def previous_item(self):
        if not self.filtered_idx:
            return
        self.current_idx = max(0, self.current_idx - 1)
        self.listbox.selection_clear(0, tk.END)
        self.listbox.selection_set(self.current_idx)
        self.listbox.see(self.current_idx)
        self._on_select()

    def next_item(self):
        if not self.filtered_idx:
            return
        self.current_idx = min(len(self.filtered_idx) - 1, self.current_idx + 1)
        self.listbox.selection_clear(0, tk.END)
        self.listbox.selection_set(self.current_idx)
        self.listbox.see(self.current_idx)
        self._on_select()

    def run_selected(self):
        rows = self._selected_manifest_rows()
        if not rows:
            messagebox.showinfo("Nothing selected", "Select one or more scene scripts first.")
            return

        for row in rows:
            script = Path(row["scene_script"])
            if not script.exists():
                messagebox.showerror("Missing script", str(script))
                continue
            subprocess.Popen(["bash", str(script)])

    def open_script_folder(self):
        rows = self._selected_manifest_rows()
        folder = Path(rows[0]["scene_script"]).parent if rows else MANIFEST.parent
        try:
            subprocess.Popen(["xdg-open", str(folder)])
        except FileNotFoundError:
            subprocess.Popen(["open", str(folder)])


def main():
    root = tk.Tk()
    FileSceneBrowser(root)
    root.mainloop()


if __name__ == "__main__":
    main()
'''
    browser_code = browser_code.replace("__SCENE_MANIFEST_NAME__", scene_manifest_path.name)
    tk_script_path.write_text(browser_code)
    tk_script_path.chmod(tk_script_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)


def create_file_scene_browser(
    manifest: pd.DataFrame,
    outdir: Path,
    *,
    skip_existing: bool = True,
    include_context_masks: bool = False,
) -> Tuple[Path, pd.DataFrame]:
    """Generate one FSLeyes script per important file plus a Tk launcher script."""
    browser_dir = outdir / "file_scene_browser"
    scenes_dir = browser_dir / "scenes"
    browser_dir.mkdir(parents=True, exist_ok=True)
    scenes_dir.mkdir(parents=True, exist_ok=True)

    valid_rows = manifest_valid_file_rows(manifest)

    scene_rows = []
    for _, row in valid_rows.iterrows():
        subject = str(row["subject"])
        session = str(row["session"])
        scanT = str(row.get("scanT", ""))
        category = str(row["category"])
        name = str(row["name"])
        path = Path(valid_path_string(row["path"]))

        t_label = f"{scanT}T" if scanT and scanT.lower() not in {"nan", "none"} else "unknownT"
        fname_base = make_safe_filename(
            f"{subject}_{session}_{t_label}_{category}_{name}_{path.stem}"
        )
        script_path = scenes_dir / f"{fname_base}.sh"

        if (not skip_existing) or (not script_path.exists()):
            create_single_file_fsleyes_script(
                row,
                valid_rows,
                script_path,
                include_context_masks=include_context_masks,
            )

        search_text = " ".join([
            subject,
            session,
            t_label,
            str(row.get("scanner", "")),
            category,
            name,
            path.name,
            str(path),
        ])

        scene_rows.append({
            "subject": subject,
            "session": session,
            "scanT": scanT,
            "scanner": str(row.get("scanner", "")),
            "category": category,
            "name": name,
            "file": str(path),
            "scene_script": str(script_path),
            "display": f"{subject} | {session} | {t_label} | {category} | {name} | {path.name}",
            "search": search_text.lower(),
        })

    scene_manifest = pd.DataFrame(scene_rows)
    scene_manifest_path = browser_dir / "file_scene_manifest.csv"
    if (not skip_existing) or (not scene_manifest_path.exists()):
        scene_manifest.to_csv(scene_manifest_path, index=False)

    tk_script_path = browser_dir / "tk_file_browser.py"
    if (not skip_existing) or (not tk_script_path.exists()):
        write_tk_file_browser_script(tk_script_path, scene_manifest_path)

    return tk_script_path, scene_manifest


def launch_file_scene_browser(scene_browser_script: Path) -> None:
    """Launch the generated Tk file-scene browser."""
    subprocess.Popen([sys.executable, str(scene_browser_script)])


# =============================================================================
# Main execution entry point
# =============================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="HistoPark MRI QC v1")
    parser.add_argument("--overview", required=True, help="CSV table with subject/session and derivative paths.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--path-replace", nargs=2, action="append", default=[],
                        metavar=("OLD", "NEW"),
                        help="Replace OLD with NEW in all paths from the overview table. Can be used multiple times.")
    parser.add_argument("--tissue-threshold", type=float, default=0.7,
                        help="Threshold for SPM c1/c2/c3 tissue probability maps. Default: 0.7")
    parser.add_argument("--histogram-bins", type=int, default=80)
    parser.add_argument("--z-threshold", type=float, default=3.0)
    parser.add_argument("--skip-screenshots", action="store_true")
    parser.add_argument("--make-file-browser", action="store_true",
                        help="Create one-file grayscale FSLeyes scene scripts and a one-column Tk browser.")
    parser.add_argument("--launch-tk-browser", action="store_true",
                        help="Create and launch the one-column grayscale file browser after QC.")
    # Browser scenes intentionally do not include context masks/tissue overlays.
    parser.add_argument("--n-slices", type=int, default=5)
    parser.add_argument("--screenshot-maps", nargs="*", default=DEFAULT_SCREENSHOT_MAPS,
                        help="Map names to screenshot. Default: selected core maps.")
    parser.add_argument("--max-rows", type=int, default=None,
                        help="Optional debugging limit.")
    args = parser.parse_args(argv)

    overview_path = Path(args.overview)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load cohort overview table.
    #
    # Expected columns typically include:
    #   bids_subjID
    #   bids_sesID
    #   qMRI
    #   QSMxT
    #   relax_R2
    #   scanT
    #
    # Additional columns are tolerated.
    overview = pd.read_csv(overview_path)
    if args.max_rows:
        overview = overview.iloc[:args.max_rows].copy()

    # Normalize core columns to strings where useful.
    for col in ["bids_subjID", "bids_sesID", "scanner", "scanT"]:
        if col in overview.columns:
            overview[col] = overview[col].fillna("").astype(str)

    replacements = [(old, new) for old, new in args.path_replace]

    # -------------------------------------------------------------------------
    # Step 1: File discovery
    # -------------------------------------------------------------------------
    #
    # Build a manifest containing:
    #   - all discovered quantitative maps
    #   - error maps
    #   - masks
    #   - weighted images
    #
    # This makes the pipeline auditable and helps debugging missing files.
    manifest_rows: List[FoundFile] = []
    for _, row in tqdm(overview.iterrows(), total=len(overview), desc="Finding files"):
        found, _ = find_maps_for_row(row, replacements)
        manifest_rows.extend(found)

    manifest = pd.DataFrame([r.__dict__ for r in manifest_rows])
    manifest.to_csv(outdir / "qc_file_manifest.csv", index=False)

    # -------------------------------------------------------------------------
    # Step 2: Quantitative metrics
    # -------------------------------------------------------------------------
    metrics, error_metrics, hist = process_metrics(
        overview,
        manifest,
        tissue_threshold=args.tissue_threshold,
        histogram_bins=args.histogram_bins,
    )

    metrics_flagged = add_outlier_flags(metrics, z_thresh=args.z_threshold)
    error_flagged = add_outlier_flags(error_metrics, z_thresh=args.z_threshold)

    # -------------------------------------------------------------------------
    # Step 3: B1 homogeneity analysis
    # -------------------------------------------------------------------------
    b1_homogeneity = compute_b1_homogeneity(
        overview,
        manifest,
        tissue_threshold=args.tissue_threshold,
    )

    metrics_flagged.to_csv(outdir / "qc_metrics.csv", index=False)
    error_flagged.to_csv(outdir / "qc_error_metrics.csv", index=False)
    hist.to_csv(outdir / "qc_histograms.csv", index=False)
    b1_homogeneity.to_csv(outdir / "qc_b1_homogeneity.csv", index=False)

    # Combined outlier table.
    combined = pd.concat([metrics_flagged, error_flagged], ignore_index=True)
    if not combined.empty:
        outlier_cols = [c for c in combined.columns if c.endswith("_outlier")]
        if outlier_cols:
            outliers = combined[np.any(combined[outlier_cols].fillna(False).values, axis=1)]
            outliers.to_csv(outdir / "qc_outlier_flags.csv", index=False)
        else:
            pd.DataFrame().to_csv(outdir / "qc_outlier_flags.csv", index=False)

    # -------------------------------------------------------------------------
    # Step 4: Cohort plots and screenshots
    # -------------------------------------------------------------------------
    plot_metric_boxplots(metrics_flagged, outdir / "figures" / "cohort_plots")
    plot_metric_boxplots(error_flagged, outdir / "figures" / "error_plots")
    plot_b1_homogeneity(b1_homogeneity, outdir / "figures" / "b1_homogeneity")
    plot_heatmap(metrics_flagged, outdir / "figures" / "cohort_plots" / "heatmap_brain_cv_z.png",
                 region="brain", metric="cv")
    plot_heatmap(metrics_flagged, outdir / "figures" / "cohort_plots" / "heatmap_brain_sd_norm_range_z.png",
                 region="brain", metric="sd_norm_range")
    plot_heatmap(metrics_flagged, outdir / "figures" / "cohort_plots" / "heatmap_brain_out_of_range_fraction_z.png",
                 region="brain", metric="out_of_range_fraction")
    make_histogram_figures(hist, outdir)

    if not args.skip_screenshots:
        make_screenshots(
            overview,
            manifest,
            outdir,
            maps_to_plot=args.screenshot_maps,
            tissue_threshold=args.tissue_threshold,
            n_slices=args.n_slices,
        )
        make_segmentation_overlays(
            overview,
            manifest,
            outdir,
            tissue_threshold=args.tissue_threshold,
            n_slices=args.n_slices,
        )

    if args.make_file_browser or args.launch_tk_browser:
        browser_script, browser_manifest = create_file_scene_browser(
            manifest,
            outdir,
            skip_existing=True,
            include_context_masks=False,
        )
        print(f"File-scene browser written to: {browser_script.parent}")
        print(f"Run browser with: python {browser_script}")
        if args.launch_tk_browser:
            launch_file_scene_browser(browser_script)

    summary = {
        "overview": str(overview_path),
        "n_rows_in_overview": int(len(overview)),
        "n_manifest_rows": int(len(manifest)),
        "n_metric_rows": int(len(metrics_flagged)),
        "n_error_metric_rows": int(len(error_flagged)),
        "n_histogram_rows": int(len(hist)),
        "n_b1_homogeneity_rows": int(len(b1_homogeneity)),
        "outputs": [
            "qc_file_manifest.csv",
            "qc_metrics.csv",
            "qc_error_metrics.csv",
            "qc_histograms.csv",
            "qc_b1_homogeneity.csv",
            "qc_outlier_flags.csv",
            "figures/",
        ],
        "notes": [
            "No voxelwise 3T-vs-7T comparison is performed in v1.",
            "Masks are resampled into each image grid.",
            "QSM uses Chimap only and excludes singlepass/minIP/SWI/R2star.",
            "SPM tissue masks are thresholded according to --tissue-threshold.",
            "Screenshots use fixed display ranges for core maps.",
            "Chimap CV is set to NaN; use relative_spread/robust_cv instead.",
            "Histograms for core maps are clipped to fixed plausible ranges.",
            "R1 overlays are generated for brain/GM/WM/CSF masks in red at 50% opacity.",
            "B1 homogeneity metrics are written to qc_b1_homogeneity.csv.",
            "One-file grayscale FSLeyes scene browser can be generated with --make-file-browser.",
            "Error plots, cohort plots, and B1 homogeneity plots are separated by scanT.",
        ],
    }
    with open(outdir / "qc_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\nQC complete.")
    print(f"Output directory: {outdir}")
    print(f"Manifest rows: {len(manifest)}")
    print(f"Metric rows: {len(metrics_flagged)}")
    print(f"Error metric rows: {len(error_flagged)}")
    print(f"Histogram rows: {len(hist)}")
    print(f"B1 homogeneity rows: {len(b1_homogeneity)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main([

        "--overview", "/data/u_krohn_software/my_projects/IronSleep_Histopark/pipeline_ironsleep/overview_table.csv",

        "--outdir", "/data/u_krohn_software/my_projects/IronSleep_Histopark/pipeline_ironsleep/histopark_QC_v1",

        "--make-file-browser"

    ]))
