"""Verification and validation of the lid-driven cavity runs.

Reads the centerline profiles written by OpenFOAM (runs/Re*_N*/postProcessing/centerlines/<time>/),
compares them with Ghia, Ghia & Shin (1982), estimates discretization error with the Grid
Convergence Index (Celik et al., 2008) and writes figures/ and results/.

Usage: python3 scripts/postprocess.py   (run from the repository root)
"""
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RUNS, FIG, RES = ROOT / "runs", ROOT / "figures", ROOT / "results"
RE_COLUMNS = {100: 1, 1000: 3}  # column of each Re in the Ghia tables (after the coordinate)


def load_ghia(name, re_number):
    rows = [line.split() for line in (ROOT / "reference" / name).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    data = np.array([[float(r[0]), float(r[RE_COLUMNS[re_number]])] for r in rows])
    return data[np.argsort(data[:, 0])]


def load_extrema(re_number):
    for line in (ROOT / "reference" / "high_accuracy_extrema.txt").read_text().splitlines():
        parts = line.split()
        if parts and not line.startswith("#") and int(parts[0]) == re_number:
            return dict(zip(("u_min", "v_max", "v_min"), map(float, parts[1:])))
    raise KeyError(re_number)


def load_profile(case, line, time_dir):
    data = np.loadtxt(case / "postProcessing" / "centerlines" / time_dir / f"{line}_U.xy")
    return data[:, 0], data[:, 1], data[:, 2]


def write_times(case):
    times = (case / "postProcessing" / "centerlines").iterdir()
    return sorted((t.name for t in times), key=float)


def extremum(x, f, kind):
    # Parabola through the discrete extremum and its neighbours (sub-sample location).
    i = int(np.argmin(f) if kind == "min" else np.argmax(f))
    i = min(max(i, 1), len(f) - 2)
    a, b, c = np.polyfit(x[i - 1:i + 2], f[i - 1:i + 2], 2)
    xe = -b / (2 * a)
    return float(a * xe ** 2 + b * xe + c), float(xe)


def gci(phi1, phi2, phi3, r=2.0):
    """Celik et al. (2008), constant refinement ratio; phi1 is the finest mesh."""
    e21, e32 = phi2 - phi1, phi3 - phi2
    if e21 == 0 or e32 == 0:
        return {"note": "no change between meshes"}
    ratio = e32 / e21
    p = abs(math.log(abs(ratio))) / math.log(r)
    extrapolated = (r ** p * phi1 - phi2) / (r ** p - 1)
    ea21 = abs(e21 / phi1)
    ea32 = abs(e32 / phi2)
    gci21 = 1.25 * ea21 / (r ** p - 1)
    gci32 = 1.25 * ea32 / (r ** p - 1)
    return {
        "observed_order": p,
        "monotonic": ratio > 0,
        "extrapolated": extrapolated,
        "gci_fine_percent": 100 * gci21,
        "asymptotic_ratio": gci32 / (r ** p * gci21),
    }


def analyse():
    FIG.mkdir(exist_ok=True)
    (RES / "profiles").mkdir(parents=True, exist_ok=True)
    summary = {}
    for re_number in (100, 1000):
        ghia_u = load_ghia("ghia1982_table1_u.txt", re_number)
        ghia_v = load_ghia("ghia1982_table2_v.txt", re_number)
        cases = sorted(RUNS.glob(f"Re{re_number}_N*"), key=lambda p: int(re.search(r"N(\d+)", p.name)[1]))
        rows, profiles = [], {}
        for case in cases:
            n = int(re.search(r"N(\d+)", case.name)[1])
            times = write_times(case)
            last, prev = times[-1], times[-2]
            y, u, _ = load_profile(case, "vertical", last)
            x, _, v = load_profile(case, "horizontal", last)
            _, u_prev, _ = load_profile(case, "vertical", prev)
            _, _, v_prev = load_profile(case, "horizontal", prev)
            inner_u = (ghia_u[:, 0] > 0) & (ghia_u[:, 0] < 1)   # wall values are boundary conditions
            inner_v = (ghia_v[:, 0] > 0) & (ghia_v[:, 0] < 1)
            du = np.interp(ghia_u[inner_u, 0], y, u) - ghia_u[inner_u, 1]
            dv = np.interp(ghia_v[inner_v, 0], x, v) - ghia_v[inner_v, 1]
            u_min, y_umin = extremum(y, u, "min")
            v_max, x_vmax = extremum(x, v, "max")
            v_min, x_vmin = extremum(x, v, "min")
            rows.append({
                "N": n, "cells": n * n, "final_time": float(last),
                "steady_change": float(max(np.abs(u - u_prev).max(), np.abs(v - v_prev).max())),
                "u_rms_vs_ghia": float(np.sqrt(np.mean(du ** 2))), "u_max_vs_ghia": float(np.abs(du).max()),
                "v_rms_vs_ghia": float(np.sqrt(np.mean(dv ** 2))), "v_max_vs_ghia": float(np.abs(dv).max()),
                "u_min": u_min, "y_at_u_min": y_umin, "v_max": v_max, "x_at_v_max": x_vmax,
                "v_min": v_min, "x_at_v_min": x_vmin,
            })
            profiles[n] = (y, u, x, v)
            # Keep the final centerline data in the repository (runs/ is not versioned).
            np.savetxt(RES / "profiles" / f"Re{re_number}_N{n}_vertical.csv", np.column_stack([y, u]),
                       delimiter=",", header="y,u", comments="", fmt="%.8g")
            np.savetxt(RES / "profiles" / f"Re{re_number}_N{n}_horizontal.csv", np.column_stack([x, v]),
                       delimiter=",", header="x,v", comments="", fmt="%.8g")
        fine = rows[-3:]
        convergence = {q: gci(fine[2][q], fine[1][q], fine[0][q]) for q in ("u_min", "v_max", "v_min")}
        convergence["meshes"] = [r["N"] for r in fine]
        reference = load_extrema(re_number)
        for q in ("u_min", "v_max", "v_min"):
            convergence[q]["reference"] = reference[q]
            convergence[q]["finest"] = fine[2][q]
        summary[re_number] = {"cases": rows, "gci": convergence}
        plot_profiles(re_number, profiles, ghia_u, ghia_v)
    plot_errors(summary)
    (RES / "summary.json").write_text(json.dumps(summary, indent=2))
    (RES / "summary.md").write_text(markdown(summary))
    return summary


def plot_profiles(re_number, profiles, ghia_u, ghia_v):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.4))
    for n, (y, u, x, v) in profiles.items():
        ax1.plot(u, y, lw=1.2, label=f"OpenFOAM {n}×{n}")
        ax2.plot(x, v, lw=1.2, label=f"OpenFOAM {n}×{n}")
    ax1.plot(ghia_u[:, 1], ghia_u[:, 0], "ko", mfc="none", label="Ghia et al. (1982)")
    ax2.plot(ghia_v[:, 0], ghia_v[:, 1], "ko", mfc="none", label="Ghia et al. (1982)")
    ax1.set(xlabel="u / U", ylabel="y / L", title=f"u on vertical centerline, Re = {re_number}")
    ax2.set(xlabel="x / L", ylabel="v / U", title=f"v on horizontal centerline, Re = {re_number}")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / f"centerlines_Re{re_number}.png", dpi=150)
    plt.close(fig)


def plot_errors(summary):
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for re_number, data in summary.items():
        n = [r["N"] for r in data["cases"]]
        ax.loglog(n, [r["u_rms_vs_ghia"] for r in data["cases"]], "o-", label=f"u, Re = {re_number}")
        ax.loglog(n, [r["v_rms_vs_ghia"] for r in data["cases"]], "s--", label=f"v, Re = {re_number}")
    ax.set(xlabel="cells per side N", ylabel="RMS difference vs Ghia et al.",
           title="Mesh refinement vs benchmark")
    meshes = sorted({r["N"] for data in summary.values() for r in data["cases"]})
    ax.set_xticks(meshes, [str(n) for n in meshes])
    ax.minorticks_off()
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "error_vs_mesh.png", dpi=150)
    plt.close(fig)


def markdown(summary):
    out = ["# Results\n", "Generated by `scripts/postprocess.py`. RMS/max differences use the interior Ghia points "
           "(walls excluded); `steady change` is the largest centerline difference between the last two writes "
           "(10 s apart).\n"]
    for re_number, data in summary.items():
        out.append(f"\n## Re = {re_number}\n")
        out.append("| Mesh | Steady change | u RMS | u max | v RMS | v max | u_min | v_max | v_min |")
        out.append("|---|---|---|---|---|---|---|---|---|")
        for r in data["cases"]:
            out.append(f"| {r['N']}×{r['N']} | {r['steady_change']:.1e} | {r['u_rms_vs_ghia']:.4f} | "
                       f"{r['u_max_vs_ghia']:.4f} | {r['v_rms_vs_ghia']:.4f} | {r['v_max_vs_ghia']:.4f} | "
                       f"{r['u_min']:.5f} | {r['v_max']:.5f} | {r['v_min']:.5f} |")
        g = data["gci"]
        out.append(f"\nGrid Convergence Index on meshes {g['meshes']} (r = 2):\n")
        out.append("| Quantity | Observed order p | Monotonic | Finest mesh | GCI fine (%) | Asymptotic ratio "
                   "| Extrapolated | High-accuracy reference | Extrapolated vs reference (%) |")
        out.append("|---|---|---|---|---|---|---|---|---|")
        for q in ("u_min", "v_max", "v_min"):
            c = g[q]
            if "observed_order" not in c:
                out.append(f"| {q} | – | – | {c['finest']:.5f} | – | – | – | {c['reference']:.7f} | {c.get('note', '')} |")
            else:
                diff = 100 * abs(c["extrapolated"] - c["reference"]) / abs(c["reference"])
                out.append(f"| {q} | {c['observed_order']:.2f} | {'yes' if c['monotonic'] else 'no'} | "
                           f"{c['finest']:.5f} | {c['gci_fine_percent']:.3f} | {c['asymptotic_ratio']:.3f} | "
                           f"{c['extrapolated']:.5f} | {c['reference']:.7f} | {diff:.3f} |")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    analyse()
    print((RES / "summary.md").read_text())
