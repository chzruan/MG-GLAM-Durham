#!/usr/bin/env python3
"""Measure plot-ready summaries of a pair of GLAM runs' on-the-fly products with glamio.

usage: measure_products.py {fast,cat,planes} [--lcdm RUN1] [--fr RUN1] [--out DIR]

Stages:
  fast    schedule, P(k) tables of both runs, mode cubes (delta and momentum spectra,
          band-limited real-space projection), lens planes (one plane, variances,
          PDFs), subsamples (velocity rms, PDFs, a thin slab of particles)
  cat     halo catalogues of both runs (cumulative mass functions) and f_R at the
          halo centres of the f(R) run
  planes  per-plane rms of Sigma/Sigma_bar on 8 Mpc/h pixels at every output
Each stage writes products_<stage>.npz into --out (default: the current directory) with
a JSON provenance string under the key "provenance".  plot_products.py turns them
into figures.  The defaults point at the node-2 DEGRACE companion runs (LCDM
confirmation run and f(R) pilot); --lcdm/--fr accept any two Run1 directories with
the same output schedule (the f(R) one needs fR_halo files for the cat stage).
Light enough for a login node: about one minute and < 4 GB per stage.
"""
import argparse
import datetime
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))
import glamio as G  # noqa: E402

RUNS = {"lcdm": ROOT + "/DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048_confirm20261018/Run1",
        "fr": ROOT + "/DEGRACE-companion/runs/Node_002_FR_L1024Np1024Ng2048/Run1"}
JOBS = {"lcdm": 12185817, "fr": 12180110}
DEFAULT_RUNS = dict(RUNS)
OUT = "."
BOX = 1024.0
Z_CUBES = (0.0, 1.0, 2.0)          # epochs for the mode-cube spectra, subsample PDFs, plane PDFs
Z_CATS = (0.0, 0.5, 1.0, 2.0)      # epochs for the halo catalogues and f_R at halo centres
SLAB_PART = 16.0                   # Mpc/h, thickness of the particle slab shown in the field figure
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f} s] {msg}", flush=True)


def provenance(stage, extra):
    return json.dumps({"stage": stage, "created": datetime.datetime.now().isoformat(timespec="seconds"),
                       "script": os.path.abspath(__file__), "glamio": G.__version__, "runs": RUNS,
                       "jobs": JOBS if RUNS == DEFAULT_RUNS else {}, "box": BOX, **extra})


def step_at(S, z):
    return int(S.row_at_z(z)["step"])


def shell_power(slab_fn, n, kf):
    """Shell-averaged Box^3 <|F_k|^2> over the integer shells 1..n of a half cube.

    slab_fn(iz, IX, IY) returns |F_k|^2 on the (2n+1)^2 slab iz; modes with iz > 0
    stand for a Hermitian pair and get weight 2.  Returns (k_mean, P) per shell.
    """
    idx = np.arange(-n, n + 1)
    IX, IY = np.meshgrid(idx, idx, indexing="ij")
    sum_w = np.zeros(n + 2)
    sum_p = np.zeros(n + 2)
    sum_k = np.zeros(n + 2)
    for iz in range(0, n + 1):
        r = np.sqrt(IX ** 2 + IY ** 2 + iz ** 2)
        s = np.rint(r).astype(int)
        m = (s >= 1) & (s <= n)
        w = 2.0 if iz > 0 else 1.0
        p2 = slab_fn(iz, IX, IY)
        sum_p += np.bincount(s[m], weights=w * p2[m], minlength=n + 2)
        sum_w += np.bincount(s[m], weights=np.full(m.sum(), w), minlength=n + 2)
        sum_k += np.bincount(s[m], weights=w * r[m], minlength=n + 2)
    sh = np.arange(1, n + 1)
    return kf * sum_k[sh] / sum_w[sh], BOX ** 3 * sum_p[sh] / sum_w[sh]


def stage_fast():
    out = {}
    # ---- schedule and P(k) tables of both runs
    for tag, run in RUNS.items():
        S = G.read_schedule(run + "/OutputSchedule.dat")
        out[f"sched_{tag}_step"] = np.asarray(S.steps, int)
        out[f"sched_{tag}_a"] = np.asarray(S.a, float)
        out[f"sched_{tag}_z"] = np.asarray(S.z, float)
        out[f"sched_{tag}_chi"] = np.asarray(S.chi, float)
        out[f"sched_{tag}_class"] = np.asarray(S.classes, int)
        ks, ps, ns, aa = [], [], [], []
        for st in out[f"sched_{tag}_step"]:
            P = G.read_powerdm(f"{run}/PowerDM.R.{st:04d}.0001.dat", box=BOX)
            ks.append(P.k); ps.append(P.pk); ns.append(P.nmodes); aa.append(P.a)
        nb = min(len(k) for k in ks)
        assert all(np.array_equal(k[:nb], ks[0][:nb]) for k in ks), "P(k) binning differs between steps"
        out[f"pk_{tag}_k"] = ks[0][:nb]
        out[f"pk_{tag}_nmodes"] = ns[0][:nb]
        out[f"pk_{tag}_p"] = np.array([p[:nb] for p in ps])
        out[f"pk_{tag}_a"] = np.array(aa)
        out[f"pk_{tag}_knyq"] = P.k_nyq
        log(f"{tag}: schedule {len(S.steps)} rows, P(k) tables {len(ps)} x {nb} bins")
    # linear P(k) of the input table (5 header lines, then k [h/Mpc] P [(Mpc/h)^3])
    tab = np.loadtxt(os.path.dirname(RUNS["lcdm"]) + "/PkTable.dat", skiprows=5)
    out["pk_lin_k"], out["pk_lin_p"] = tab[:, 0], tab[:, 1]
    out["box"] = BOX
    # ---- mode cubes of the LCDM run at Z_CUBES
    S = G.read_schedule(RUNS["lcdm"] + "/OutputSchedule.dat")
    meta = G.read_meta(RUNS["lcdm"] + "/deltam.0181.0001.meta")
    out["omega_m"] = float(meta["Omega_m"].split()[0]); out["omega_l"] = float(meta["Omega_L"].split()[0])
    out["hubble"] = float(meta["hubble"].split()[0])
    cube_rows = []
    for z in Z_CUBES:
        st = step_at(S, z)
        D = G.read_modes(f"{RUNS['lcdm']}/deltam.{st:04d}.0001.meta", mmap=True)
        k, pk, _ = D.power()
        out[f"cube_dd_k_z{z:g}"], out[f"cube_dd_p_z{z:g}"] = k, pk
        mom = [G.read_modes(f"{RUNS['lcdm']}/mom_{ax}.{st:04d}.0001.meta") for ax in "xyz"]
        n, n0, kf = mom[0].n, D.n, D.k_fund
        assert all(m.n == n for m in mom) and n <= n0

        def theta2(iz, IX, IY, mom=mom, kf=kf):
            px, py, pz = (m.data[:, :, iz].astype(np.complex128) for m in mom)
            return np.abs(kf * (IX * px + IY * py + iz * pz)) ** 2

        def delta2(iz, IX, IY, D=D, n=n, n0=n0):
            return np.abs(D.data[n0 - n:n0 + n + 1, n0 - n:n0 + n + 1, iz].astype(np.complex128)) ** 2

        kt, ptt = shell_power(theta2, n, kf)
        kd, pdd = shell_power(delta2, n, kf)
        assert np.allclose(kt, kd)
        out[f"cube_tt_k_z{z:g}"], out[f"cube_tt_p_z{z:g}"], out[f"cube_tt_pdd_z{z:g}"] = kt, ptt, pdd
        cube_rows.append((z, st, float(D.a), float(mom[0].a_v), n0, n))
        log(f"cubes z={z:g} step {st}: deltam n={n0} ({len(k)} shells), momentum n={n} -> P_tt computed")
        if z == 0.0:
            # band-limited real-space field (k <= k_c) on a (2n0+1)^3 grid and its projection over the first z-slab
            N = 2 * n0 + 1
            order = np.fft.fftfreq(N, 1.0 / N).astype(int)
            A = D.data[np.ix_(order + n0, order + n0, np.arange(n0 + 1))]
            parseval = float(np.sum(np.abs(A[:, :, 0]) ** 2) + 2.0 * np.sum(np.abs(A[:, :, 1:]) ** 2))
            field = np.fft.irfftn(A, s=(N, N, N)) * N ** 3
            del A
            nz = int(round(256.0 / BOX * N))
            out["field_z0_proj"] = field[:, :, :nz].mean(axis=2)
            out["field_z0_depth"] = nz * BOX / N
            out["field_z0_var"] = float(field.var())
            out["field_z0_parseval"] = parseval
            out["field_z0_ngrid"] = N
            log(f"real-space field: N={N}, var {field.var():.5f} vs Parseval {parseval:.5f}, projection depth {nz * BOX / N:.1f}")
            del field
    out["cube_rows"] = np.array(cube_rows)
    # ---- lens planes: z=0 plane, per-plane std at all 42 outputs, pixel PDFs at Z_CUBES
    L = G.read_lensplanes(f"{RUNS['lcdm']}/LensPlanes.0181.0001.meta", mmap=True)
    P = L.plane("z", 1)
    out["plane_z0"] = P.reshape(512, 4, 512, 4).mean(axis=(1, 3))
    out["plane_depth"] = float(L.depth)
    stds = []
    for st in out["sched_lcdm_step"]:
        Lp = G.read_lensplanes(f"{RUNS['lcdm']}/LensPlanes.{st:04d}.0001.meta", mmap=True)
        stds.append([float(Lp.planes[i].astype(np.float64).std()) for i in range(Lp.n_planes)])
    out["plane_std"] = np.array(stds)
    log(f"planes: std of Sigma/Sigma_bar at {len(stds)} outputs, z=0 plane stored at 512^2")
    edges = np.linspace(-1.5, 2.0, 141)
    out["plane_pdf_edges"] = edges
    for z in Z_CUBES:
        st = step_at(S, z)
        Lp = G.read_lensplanes(f"{RUNS['lcdm']}/LensPlanes.{st:04d}.0001.meta", mmap=True)
        v = np.log10(np.clip(Lp.planes[8:12].astype(np.float64), 1e-6, None)).ravel()
        out[f"plane_pdf_z{z:g}"] = np.histogram(v, bins=edges)[0] / v.size / np.diff(edges)
    # ---- subsamples at the 17 primaries
    prim = [int(s) for s, c in zip(out["sched_lcdm_step"], out["sched_lcdm_class"]) if c == 1]
    vrms, va, vav, vz = [], [], [], []
    vedges = np.arange(-2500.0, 2500.1, 25.0)
    out["vel_pdf_edges"] = vedges
    for st in prim:
        Sb = G.read_subsample(f"{RUNS['lcdm']}/Subsample.{st:04d}.0001.meta")
        v = Sb.velocities()
        vrms.append(v.std(axis=0)); va.append(Sb.a); vav.append(Sb.a_v); vz.append(Sb.z)
        if any(abs(Sb.z - z) < 1e-3 for z in Z_CUBES):
            out[f"vel_pdf_z{Sb.z:g}"] = np.histogram(v[:, 0], bins=vedges)[0] / len(v) / 25.0
        if abs(Sb.z) < 1e-6:
            x = Sb.positions()
            m = x[:, 2] < SLAB_PART
            out["part_z0_xy"] = x[m][:, :2].astype(np.float32)
            out["part_z0_vz"] = v[m][:, 2].astype(np.float32)
            out["part_slab"] = SLAB_PART
            out["sub_m"], out["sub_nsel"], out["sub_np"] = Sb.m, Sb.n_sel, Sb.n_p
        del v
    out["sub_step"], out["sub_a"], out["sub_av"], out["sub_z"] = np.array(prim), np.array(va), np.array(vav), np.array(vz)
    out["sub_vrms"] = np.array(vrms)
    log(f"subsamples: {len(prim)} primaries, 1-D rms at z=0 {vrms[-1].mean():.1f} km/s, slab particles {m.sum()}")
    out["provenance"] = provenance("fast", {"z_cubes": Z_CUBES, "slab_part": SLAB_PART})
    np.savez_compressed(os.path.join(OUT, "products_fast.npz"), **out)
    log("wrote products_fast.npz")


def stage_cat():
    out = {}
    medges = np.logspace(12.0, 15.5, 36)
    out["hmf_edges"] = medges
    fedges = np.arange(12.0, 15.21, 0.2)
    out["fr_mass_edges"] = fedges
    for tag, run in RUNS.items():
        S = G.read_schedule(run + "/OutputSchedule.dat")
        for z in Z_CATS:
            st = step_at(S, z)
            stem = f"{run}/CATALOGS/CatshortV.{st:04d}.0001"
            if tag == "fr":
                F = G.read_frhalo(f"{run}/fR_halo.{st:04d}.0001.meta")
                cat = F.catalogue(prefer="text")
            else:
                cat = G.read_catalogue(stem)
            mass = cat["Mtot"].astype(np.float64)
            out[f"hmf_{tag}_z{z:g}"] = np.array([(mass > e).sum() for e in medges]) / BOX ** 3
            out[f"hmf_{tag}_n_z{z:g}"] = len(cat)
            out[f"hmf_{tag}_a_z{z:g}"] = float(cat.a)
            log(f"{tag} z={z:g} step {st}: {len(cat)} haloes, M_min {mass.min():.3e}, M_max {mass.max():.3e}")
            if tag == "fr":
                for col in ("fR_centre", "fR_mean"):
                    r = np.log10(np.clip(np.abs(F[col].astype(np.float64) / F.fR_bg_old), 1e-8, None))
                    lm = np.log10(mass)
                    rows = []
                    for lo, hi in zip(fedges[:-1], fedges[1:]):
                        sel = (lm >= lo) & (lm < hi)
                        rows.append([sel.sum()] + (list(np.percentile(r[sel], [16, 50, 84])) if sel.sum() >= 20 else [np.nan] * 3))
                    out[f"fr_{col}_z{z:g}"] = np.array(rows)
                out[f"fr_bg_z{z:g}"] = float(F.fR_bg_old)
                out[f"fr_aold_z{z:g}"] = float(F.a_old)
                out["fr_fR0"] = float(F.fR0)
                out["fr_n"] = int(F.n_fr)
                log(f"   f_R at {len(F)} halo centres: fR_bg(a_old) = {F.fR_bg_old:.3e}")
    out["provenance"] = provenance("cat", {"z_cats": Z_CATS, "mass_column": "Mtot"})
    np.savez_compressed(os.path.join(OUT, "products_cat.npz"), **out)
    log("wrote products_cat.npz")


def stage_planes():
    """Per-plane std of Sigma/Sigma_bar after block-averaging to 8 Mpc/h pixels (discreteness noise < 1 %)."""
    out = {}
    S = G.read_schedule(RUNS["lcdm"] + "/OutputSchedule.dat")
    b = 16
    stds = []
    for st in S.steps:
        Lp = G.read_lensplanes(f"{RUNS['lcdm']}/LensPlanes.{st:04d}.0001.meta", mmap=True)
        n = Lp.n_lens // b
        stds.append([float(Lp.planes[i].astype(np.float64).reshape(n, b, n, b).mean(axis=(1, 3)).std()) for i in range(Lp.n_planes)])
    out["plane_std8"] = np.array(stds)
    out["pixel8"] = b * BOX / Lp.n_lens
    out["provenance"] = provenance("planes", {"block": b})
    np.savez_compressed(os.path.join(OUT, "products_planes.npz"), **out)
    log(f"wrote products_planes.npz ({len(stds)} outputs, {out['pixel8']:g} Mpc/h pixels)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=("fast", "cat", "planes"))
    ap.add_argument("--lcdm", default=RUNS["lcdm"], help="Run1 directory of the LCDM run")
    ap.add_argument("--fr", default=RUNS["fr"], help="Run1 directory of the f(R) run")
    ap.add_argument("--out", default=".", help="directory for products_<stage>.npz")
    args = ap.parse_args()
    RUNS = {"lcdm": os.path.abspath(args.lcdm), "fr": os.path.abspath(args.fr)}
    OUT = os.path.abspath(args.out)
    os.makedirs(OUT, exist_ok=True)
    {"fast": stage_fast, "cat": stage_cat, "planes": stage_planes}[args.stage]()
