#!/usr/bin/env python3
"""Smoothed overdensity and tidal eigenvalues at every halo of one GLAM output.

Reads the ``deltam`` mode cube and the BDM catalogue of one output through
``glamio``, builds the band-limited matter field on a periodic grid, and
interpolates at each halo position

* ``delta_R``          the overdensity smoothed on the scale R,
* ``lambda_R``         the three eigenvalues of the tidal tensor
                       T_ij = d_i d_j phi with nabla^2 phi = delta_R
                       (descending, so that lambda_1 >= lambda_2 >= lambda_3
                       and lambda_1 + lambda_2 + lambda_3 = delta_R),
* ``alpha_R``          the tidal anisotropy sqrt(q^2) / (1 + delta_R) with
                       q^2 = [(l1-l2)^2 + (l2-l3)^2 + (l3-l1)^2] / 2
                       (Paranjape, Hahn & Sheth 2018),
* ``nweb_R``           the T-web class: number of eigenvalues above
                       ``--lambda-th`` (0 void, 1 sheet, 2 filament, 3 knot).

Recipe.  delta(k) W(kR) / W_CIC(k) on an N^3 grid (zero beyond the cube's
k_c), inverse FFT; T_ij(k) = (k_i k_j / k^2) delta(k) W(kR), six inverse FFTs;
cubic-spline interpolation with periodic wrap at the halo positions (grid node
j at x = j L/N, the cube's convention; ``--interp linear`` is faster but
underestimates peaks by 10-40 % at R = 2 on a 2 Mpc/h grid).  W is a Gaussian exp(-k^2 R^2 / 2) or a
top-hat 3 (sin x - x cos x) / x^3; the Gaussian is band-limited to the cube
for R >= 2 Mpc/h (leakage exp(-k_c^2 R^2 / 2) = 1 % at R = 2 for k_c = 1.5).

Cost (node-2 run, 777,663 haloes, one login-node core).  Seven N^3 inverse
FFTs per R: 40 s and 4 GB per R at N = 512 (2 Mpc/h cells, 1-3 % low at the
peaks for R = 2, converged for R >= 3), 100 s and about 11 GB at N = 768.
Use ``nice -n 19`` and keep to one output on a login node; the outputs of a
suite go to cosma8-serial.

Usage
    micromamba run -n cosemu python3 python/examples/halo_environment.py RUN1 \\
        --z 0 --R 2 5 --filter gaussian --out /cosma8/data/dp203/dc-ruan1/scratch/env-20261020/env_z0.npz

(the output is about 150 bytes per halo, 117 MB for the node-2 run at z = 0:
keep it on the project filesystem, not in /tmp.)

Output: an ``.npz`` with the catalogue columns copied for every halo (``row``,
positions, masses, ``Cvir`` ...), one set of ``delta_R<R>``, ``lambda_R<R>``
(n_halo, 3), ``alpha_R<R>``, ``nweb_R<R>`` per R, grid-level checks, and JSON
provenance under ``provenance``.
"""

import argparse
import datetime
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import map_coordinates

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import glamio as G  # noqa: E402

COPY_COLUMNS = ("row", "x", "y", "z", "Mtot", "Mbound", "Rvir", "Vcirc", "Vrms", "Cvir", "Npart",
                "Xoff", "virial_ratio", "lambda", "RadRMS", "b_a", "c_a")

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f} s] {msg}", flush=True)


def window(k, R, kind):
    """Smoothing kernel W(kR); W(0) = 1."""
    x = k * R
    if kind == "gaussian":
        return np.exp(-0.5 * x * x)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3
    return np.where(x < 1e-6, 1.0, w)


def cic_window(kx, ky, kz, k_nyq):
    """CIC assignment window as three 1-D factors (sinc^2 per axis); k_nyq = pi N_dep / L."""
    return tuple((np.sinc(ki / (2.0 * k_nyq)) ** 2).astype(np.float32) for ki in (kx, ky, kz))


def place_cube(D, N):
    """The cube's modes on the rfftn layout of an N^3 grid (complex64), zero outside k_c."""
    n = D.n
    if N < 2 * n + 2:
        raise SystemExit(f"--ngrid {N} too small for the cube: need at least {2 * n + 2}")
    A = np.zeros((N, N, N // 2 + 1), dtype=np.complex64)
    idx = np.arange(-n, n + 1)
    A[np.ix_(idx % N, idx % N, np.arange(n + 1))] = D.data[:, :, :]
    A[0, 0, 0] = 0.0
    return A


def rfft_wavenumbers(N, box):
    kf = 2.0 * np.pi / box
    kx = (np.fft.fftfreq(N, 1.0 / N) * kf)[:, None, None]
    ky = kx.reshape(1, N, 1)
    kz = (np.arange(N // 2 + 1) * kf)[None, None, :]
    return kx, ky, kz


def parseval(A):
    """Grid variance implied by a half-spectrum on the rfftn layout (nz = 0 plane once, nz > 0 twice)."""
    a2 = np.abs(A) ** 2
    return float(a2[:, :, 0].sum() + 2.0 * a2[:, :, 1:-1].sum() + a2[:, :, -1].sum())


def to_real(A, N):
    """f(x) = sum_k delta_k exp(i k x); numpy's irfftn divides by N^3."""
    return (np.fft.irfftn(A, s=(N, N, N), axes=(0, 1, 2)) * N ** 3).astype(np.float32)


ORDER = {"linear": 1, "cubic": 3}


def bspline_prefilter(N):
    """Fourier-space cubic B-spline prefilter of a periodic N^3 grid on the rfftn layout.

    The cubic B-spline sampled at the nodes is (1/6, 2/3, 1/6), so the spline
    coefficients of periodic data are its DFT divided by 2/3 + cos(2 pi m/N)/3
    per axis.  Applying this before the inverse FFT replaces scipy's recursive
    prefilter (about 60 s per 512^3 field with mode='grid-wrap') by a multiply.
    """
    m = np.fft.fftfreq(N, 1.0 / N)
    b = (2.0 / 3.0 + np.cos(2.0 * np.pi * m / N) / 3.0).astype(np.float32)
    return b[:, None, None], b[None, :, None], b[None, None, : N // 2 + 1]


def sample(field, coords, interp):
    """Periodic interpolation at fractional grid coordinates (3, n); fields are prefiltered already."""
    return map_coordinates(field, coords, order=ORDER[interp], mode="grid-wrap", prefilter=False).astype(np.float64)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", help="a Run1 directory holding OutputSchedule.dat, deltam.*.meta and CATALOGS/")
    when = ap.add_mutually_exclusive_group()
    when.add_argument("--z", type=float, default=0.0, help="output redshift (resolved through the schedule; default 0)")
    when.add_argument("--step", type=int, help="output step instead of --z")
    ap.add_argument("--box", type=int, default=1, help="realization number of the product files (default 1)")
    ap.add_argument("--R", type=float, nargs="+", default=[2.0], help="smoothing scales in Mpc/h (default 2)")
    ap.add_argument("--filter", choices=("gaussian", "tophat"), default="gaussian")
    ap.add_argument("--ngrid", type=int, default=512, help="grid side for the real-space fields (default 512)")
    ap.add_argument("--lambda-th", type=float, default=0.0, help="T-web eigenvalue threshold (default 0)")
    ap.add_argument("--interp", choices=tuple(ORDER), default="cubic",
                    help="interpolation at the halo positions (default cubic; linear underestimates peaks on coarse grids)")
    ap.add_argument("--no-deconvolve", action="store_true", help="keep the CIC window of the cube")
    ap.add_argument("--out", help="output .npz (default halo_environment_<step>.npz in the working directory)")
    args = ap.parse_args()

    run = os.path.abspath(args.run)
    if args.step is None:
        S = G.read_schedule(G.product_path(run, "schedule", 0))
        step = int(S.row_at_z(args.z)["step"])
    else:
        step = args.step
    cube_path = G.product_path(run, "deltam", step, args.box) + ".meta"
    cat_stem = G.product_path(run, "catalogue", step, args.box)
    out_path = args.out or f"halo_environment_{step:04d}.npz"

    D = G.read_modes(cube_path, mmap=True, expect_product="deltam")
    cat = G.read_catalogue(cat_stem)
    if cat.step != D.step:
        raise SystemExit(f"catalogue step {cat.step} != cube step {D.step}")
    N, L = args.ngrid, float(D.box)
    nh = len(cat)
    log(f"step {step} (z = {D.z:.4f}): cube n = {D.n} (k_c = {D.k_c} h/Mpc, N_dep = {D.n_dep}), "
        f"{nh} haloes, grid {N}^3 ({L / N:.2f} Mpc/h cells)")

    kx, ky, kz = rfft_wavenumbers(N, L)
    kmag = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2).astype(np.float32)  # dense float32; with inv_k2 the only full-size k arrays
    with np.errstate(divide="ignore"):
        inv_k2 = np.where(kmag > 0, 1.0 / kmag ** 2, 0.0).astype(np.float32)
    pos = cat.positions().astype(np.float64)
    coords = (pos.T / L * N) % N
    base = place_cube(D, N)
    if not args.no_deconvolve:
        for w in cic_window(kx, ky, kz, np.pi * D.n_dep / L):
            base /= w
    log("cube placed on the grid" + ("" if args.no_deconvolve else " and CIC-deconvolved"))
    prefilter = bspline_prefilter(N) if args.interp == "cubic" else ()

    out = {c: cat[c].copy() for c in COPY_COLUMNS if c in cat.columns}
    checks = {}
    pairs = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    kvec = tuple(k.astype(np.float32) for k in (kx, ky, kz))

    for R in args.R:
        tag = f"R{R:g}"
        if L / N > 0.67 * R:
            log(f"{tag}: cells of {L / N:.2f} Mpc/h are coarse for R = {R:g}: cubic interpolation is 1-3 % low at peaks "
                f"(0.3 % with --ngrid {int(np.ceil(1.5 * L / R / 64)) * 64})")
        if args.filter == "tophat":
            log(f"{tag}: the top-hat is not band-limited: |W(k_c R)| = {abs(window(D.k_c, R, 'tophat')):.3f} at the cube "
                f"edge, so the field is the cube's truncation of the top-hat")
        A = base * window(kmag, R, args.filter).astype(np.float32)
        var_k = parseval(A)
        for b in prefilter:  # spline coefficients instead of node values, see bspline_prefilter
            A /= b
        dens = to_real(A, N)
        dR = sample(dens, coords, args.interp)
        del dens
        log(f"{tag}: delta_R on the grid, sigma = {np.sqrt(var_k):.4f} (Parseval); "
            f"at haloes mean {dR.mean():+.4f}, min {dR.min():+.4f}, max {dR.max():+.3f}")

        T = np.empty((nh, 3, 3))
        for i, j in pairs:
            comp = to_real(A * (kvec[i] * kvec[j] * inv_k2), N)
            T[:, i, j] = T[:, j, i] = sample(comp, coords, args.interp)
            del comp
        del A
        lam = np.linalg.eigvalsh(T)[:, ::-1]  # descending
        trace_err = float(np.abs(lam.sum(axis=1) - dR).max())
        q2 = 0.5 * ((lam[:, 0] - lam[:, 1]) ** 2 + (lam[:, 1] - lam[:, 2]) ** 2 + (lam[:, 2] - lam[:, 0]) ** 2)
        unphysical = 1.0 + dR <= 0.0  # possible only through the k_c truncation or the CIC deconvolution
        with np.errstate(divide="ignore", invalid="ignore"):
            alpha = np.where(unphysical, np.nan, np.sqrt(q2) / (1.0 + dR))
        nweb = (lam > args.lambda_th).sum(axis=1).astype(np.int8)
        log(f"{tag}: eigenvalues done, max |trace - delta_R| = {trace_err:.2e}; "
            f"alpha_R median {np.nanmedian(alpha):.3f} ({int(unphysical.sum())} haloes with 1 + delta_R <= 0 -> NaN); "
            f"T-web void/sheet/filament/knot = " + "/".join(str(int((nweb == c).sum())) for c in range(4)))
        out[f"delta_{tag}"] = dR
        out[f"lambda_{tag}"] = lam
        out[f"alpha_{tag}"] = alpha
        out[f"nweb_{tag}"] = nweb
        checks[tag] = {"sigma_parseval": float(np.sqrt(var_k)),
                       "max_trace_error": trace_err, "mean_delta_at_haloes": float(dR.mean()),
                       "n_alpha_nan": int(unphysical.sum())}

    out["provenance"] = json.dumps({
        "script": os.path.abspath(__file__), "glamio": G.__version__, "date": datetime.datetime.now().isoformat(),
        "run": run, "step": step, "realization": args.box, "z": float(D.z), "a": float(D.a), "box": L, "cube": cube_path,
        "catalogue": str(cat.path), "n_halo": nh,
        "cube_n": int(D.n), "k_c": float(D.k_c), "N_dep": int(D.n_dep), "ngrid": N, "R": list(map(float, args.R)),
        "filter": args.filter, "interp": args.interp, "cic_deconvolved": not args.no_deconvolve, "lambda_th": args.lambda_th,
        "columns_copied": [c for c in COPY_COLUMNS if c in cat.columns], "checks": checks,
        "conventions": "delta_R = sum of eigenvalues; lambda_R descending; alpha_R = sqrt(q^2)/(1+delta_R), NaN where "
                       "1+delta_R <= 0; "
                       "nweb_R = number of eigenvalues above lambda_th; positions as in the catalogue (Mpc/h)"})
    np.savez(out_path, **out)
    log(f"wrote {out_path} ({os.path.getsize(out_path) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
