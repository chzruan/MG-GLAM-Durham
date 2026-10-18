#!/usr/bin/env python3
"""Generate MG-GLAM inputs (PkTable.dat, Init.dat, provenance) for a DEGRACE
companion run: either a companion of one of Baojiu Li's ECOSMOG nodes (same
cosmology, same 42 output epochs, his CAMB linear spectrum) or a new cosmology
with the 17 primary redshifts plus 25 stations chi = 256 k Mpc/h computed
here (matter + Lambda only, the convention of the ECOSMOG suite).

Item I3 of docs/degrace_companion_implementation_plan.md.

Examples
  # companion of ECOSMOG node 2, LCDM, Ng = 2048 on Baojiu's 1024^3 ICs
  python3 gen_inputs.py --node 2 --gravity lcdm --out runs/Node_002_LCDM
  # f(R) twin (fR0 from the design, n = 1)
  python3 gen_inputs.py --node 2 --gravity fr   --out runs/Node_002_FR
  # new cosmology with its own CAMB z=0 spectrum
  python3 gen_inputs.py --cosmology 0.3089 0.6774 0.8159 --ombh2 0.02237 \
      --pk pofk_z0.dat --out runs/fid_box1
  # verification S2: regenerate every ECOSMOG list and the fiducial Slack list
  python3 gen_inputs.py --check

The ICs themselves are ingested separately with ic2pm (default `half`).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
ECOSMOG = REPO / "docs" / "degrace_ecosmog"
BL267_ICS = Path("/cosma8/data/dp203/bl267/Projects/Ongoing/DEGRACE/ICs/IC_data/L1024")
C_H0 = 2997.92458  # Mpc/h
PRIMARY_Z = np.arange(0.0, 2.0001, 0.125)
STATION_CHI = 256.0 * np.arange(1, 26)
# Baojiu's 42 scale factors for the Planck fiducial (Slack, slack-chat-degrace-prep.txt)
FIDUCIAL_A = sorted([
    0.09885636873988636, 0.11436339923360468, 0.13101192462219735, 0.14880957763634076,
    0.16776699787816512, 0.18789864408246054, 0.2092238939643731, 0.2317683231209623,
    0.25556525683790143, 0.2806576241072646, 0.30710027529154543, 0.3333333333333333,
    0.33496251467163335, 0.34782608695652173, 0.36363636363636365, 0.3643315849360659,
    0.38095238095238093, 0.3953165138959508, 0.4, 0.42105263157894735, 0.42805303395318656,
    0.4444444444444444, 0.4627096181121927, 0.47058823529411764, 0.49949502290382997, 0.5,
    0.5333333333333333, 0.5386678234812446, 0.5714285714285714, 0.5805485755215236,
    0.6153846153846154, 0.6255358719332695, 0.6666666666666666, 0.6741267910946328,
    0.7269451849577723, 0.7272727272727273, 0.7847794661808203, 0.8, 0.8486352893278191,
    0.8888888888888888, 0.9198100312830662, 1.0])


# ----------------------------------------------------------------- epochs
def chi(a: float, om: float) -> float:
    f = lambda x: 1.0 / (x * x * np.sqrt(om / x**3 + 1.0 - om))  # noqa: E731
    return C_H0 * quad(f, a, 1.0, epsabs=0, epsrel=1e-12)[0]


def output_epochs(om: float) -> np.ndarray:
    """The 42 expansion factors: 17 primary + 25 stations, ascending in a."""
    stations = [brentq(lambda x: chi(x, om) - c, 0.005, 0.9999999, xtol=1e-14) for c in STATION_CHI]
    return np.sort(np.concatenate([1.0 / (1.0 + PRIMARY_Z), stations]))


def check_lists() -> int:
    """S2: reproduce the fiducial Slack list and all 96 ECOSMOG lists."""
    rc = 0
    pred = output_epochs(0.3089)
    d = np.max(np.abs(pred / np.array(FIDUCIAL_A) - 1.0))
    print(f"fiducial (Om=0.3089) vs Slack list: max |da/a| = {d:.2e}  {'OK' if d < 1e-7 else 'FAIL'}")
    rc |= d >= 1e-7
    aout = json.load(open(ECOSMOG / "aout.json"))
    worst = 0.0
    for r in csv.DictReader(open(ECOSMOG / "nodes.csv")):
        pred = output_epochs(float(r["Omega_m"]))
        d = np.max(np.abs(pred / np.sort(aout[r["node"]]) - 1.0))
        worst = max(worst, d)
    print(f"96 ECOSMOG nodes vs aout.json: max |da/a| = {worst:.2e}  {'OK' if worst < 1e-6 else 'FAIL'}")
    rc |= worst >= 1e-6
    return int(rc)


# ----------------------------------------------------------------- files
def pktable(ombh2: float, omch2: float, om: float, sigma8: float, pk: np.ndarray) -> str:
    """5-line header (PMP2init derives hubble = sqrt((ombh2+omch2)/Om)) + body."""
    head = (f"#Om_b h^2             =  {ombh2}\n"
            f"#Om_c h^2             =  {omch2}\n"
            f"#Om_Lambda            =  {1.0 - om}\n"
            f"#Om_m (1-Om_K-Om_L)   =  {om}\n"
            f"#at z                 =  0.00E+00  sigma8 (all matter)=  {sigma8:.6f}\n")
    body = "".join(f"    {k:.5E}    {p:.5E}\n".replace("E-0", "E-").replace("E+0", "E+") for k, p in pk)
    return head + body


def read_pk(path: Path) -> np.ndarray:
    pk = np.loadtxt(path, comments="#")
    if pk.ndim != 2 or pk.shape[1] < 2:
        sys.exit(f"{path}: expected two columns k/h, P(k)")
    pk = pk[:, :2]
    r = pk[1:, 0] / pk[:-1, 0]
    if np.any(np.abs(r / r.mean() - 1.0) > 2e-2):
        sys.exit(f"{path}: k must be uniformly log-spaced (PMP2init requirement)")
    return pk


def is_primary(a: float) -> bool:
    return bool(np.any(np.isclose(a, 1.0 / (1.0 + PRIMARY_Z), rtol=0, atol=1e-9)))


def init_dat(a: dict, epochs: np.ndarray) -> str:
    z = 1.0 / epochs[::-1] - 1.0           # descending in z, i.e. ascending in a
    z[z < 1e-12] = 0.0
    lines = [" ".join(f"{v:10.7f}" for v in z[i:i + 6]) for i in range(0, len(z), 6)]
    if a["classes"] == "primary":          # same order as the redshift list
        cl = [1 if is_primary(e) else 2 for e in epochs[::-1]]
        class_block = (f"Output classes            = {len(cl):9d}  0-all class 1; 1 full, 2 station per redshift\n"
                       + "\n".join(" ".join(f"{c:2d}" for c in cl[i:i + 20]) for i in range(0, len(cl), 20)) + "\n")
    else:
        class_block = "Output classes            =         0  0-all class 1\n"
    extra_block = (class_block
                   + f"Lens slabs per axis       = {a['lens_slabs']:9d}  0-off\n"
                   + f"Lens plane grid           = {a['lens_grid']:9d}  0-Ngrid\n"
                   + f"Extra outputs             = {a['extra']:9d}  0-off 1-low-k modes, subsample, fR at halos\n"
                   + f"kc_delta                  = {a['kc_delta']:9.3f}  h/Mpc\n"
                   + f"kc_mom                    = {a['kc_mom']:9.3f}  h/Mpc\n"
                   + f"N_dep_mom                 = {a['ndep_mom']:9d}  0-Ngrid/2\n"
                   + f"Subsample M               = {a['subsample_m']:9d}  0-off; keep 1 in M particles\n"
                   + f"fR at halos               = {a['fr_halos']:9d}  0-off 1-on (f(R) runs)\n")
    params = "".join("Parametr =     0.000\n" for _ in range(9))
    fr0 = "0.0" if a["mg_flag"] == 0 else f"{a['fR0']:.4E}"
    return (
        f"Box      =   {a['box']:.2f}\n"
        f"Nrow     =      {a['nrow']}\n"
        f"Ngrid    =     {a['ngrid']}\n"
        f"sig8     =    {a['sigma8']:.5f}\n"
        f"z_init   =   {a['zinit']:.3f}\n"
        f"step da  =  {a['da']:.4E}\n"
        f"z_final  =     0.000\n"
        f"#outputs =       {-len(z):3d}  (<0: hit these redshifts exactly)\n"
        + "\n".join(lines) + "\n"
        f"dens_thr =    30.000\n"
        f"Vrms     =     0.000\n"
        f"#Params  =        10\n"
        f"{params}"
        f"Pk tune  =     1.000\n"
        f"Steps between checkpoints =       {a['ncheckpoint']:3d}\n"
        f"Save snapshots            =         0  0-no, 1-yes, 2-only last\n"
        f"DM power spectrum         =         1  0-no\n"
        f"Redshift distortions      =         0  0-no\n"
        f"DM PDF                    =         0  0-no\n"
        f"Biasing model             =         0  0-no 1-particles 2-density\n"
        f"Write Mocks               =         0  0-no 1-yes\n"
        f"Dump random fraction      =         0  0-no 1-yes\n"
        f"Find BDM halos            =         {a['bdm']}  0-no 1-yes\n"
        f"MG_flag                   =         {a['mg_flag']}  \n"
        f"MG_test                   =         0  \n"
        f"MG_model                  =         1  \n"
        f"Vcycle                    =         1  \n"
        f"Fcycle                    =         0  \n"
        f"Wcycle                    =         0  \n"
        f"iter_count_max            =         2  \n"
        f"res_conv                  = 1.0000E-08 \n"
        f"fR_n                      =         1  \n"
        f"fR0                       = {fr0}\n"
        f"N_branch                  =         0  \n"
        f"S_branch                  =         0  \n"
        f"H0rc                      = 1.000  \n"
        f"sym_astar                 = 0.50000  \n"
        f"sym_xi                    = 0.00100  \n"
        f"sym_beta                  = 0.10000  \n"
        f"power_law                 =         0  \n"
        f"born_infeld               =         0  \n"
        f"kmf_n                     =         3  \n"
        f"kmf_K0                    = 1.00000 \n"
        f"kmf_beta                  = 0.10000  \n"
        f"csf_coupling              =         1\n"
        f"csf_potential             =         1\n"
        f"csf_alpha                 = 0.10000\n"
        f"csf_beta                  = -0.20000\n"
        + extra_block
    )


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_hash() -> str:
    try:
        return subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


# ----------------------------------------------------------------- main
def bdm_config(catalogue_format=None):
    """BDM.config of the companion runs: the finder settings used by the pilots
    (test/BDM.config) plus, when requested, the catalogue file format (text,
    binary or both; the key is omitted otherwise so that older builds accept the file)."""
    return (
        "! ------------- BDM configuration file. Spaces between entries are significant\n"
        "!                   These are main parameters:\n"
        "          iVirial   =      1           ! switch: 1 =Virial,  0=200 overdensity\n"
        "          NradP     =     100           ! Number of shells for halo profiles\n"
        "          Rext      =  0.6             ! extr radius shift at m=1e15\n"
        "          SlopeR    =  1.667           ! slope for extr radius shift\n"
        "          MassMin   =  1.000E+11       ! Minimum halo mass\n"
        "          dLogR     =  2.000E-02       ! size of log binning for potential\n"
        "          dLogP     =  1.000E-02       ! size of log binning for profiles\n"
        + (f"          CatalogueFormat = {catalogue_format}   ! binary (.bin + .meta, default), text (CatshortV.*.DAT) or both\n"
           if catalogue_format else "")
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--node", type=int, help="ECOSMOG node 1..96 (companion run)")
    ap.add_argument("--cosmology", nargs=3, type=float, metavar=("OMEGA_M", "H", "SIGMA8"), help="new cosmology")
    ap.add_argument("--ombh2", type=float, default=0.02237, help="baryon density for new cosmologies")
    ap.add_argument("--pk", type=Path, help="CAMB z=0 matter power (k/h, P) for new cosmologies")
    ap.add_argument("--fR0", type=float, default=None, help="f_R0 (negative) for a new f(R) cosmology")
    ap.add_argument("--gravity", choices=["lcdm", "fr"], default="lcdm")
    ap.add_argument("--out", type=Path, help="run directory to create (holds Init.dat, PkTable.dat, Run1/)")
    ap.add_argument("--box", type=float, default=1024.0)
    ap.add_argument("--nrow", type=int, default=1024, help="particles per side (1024 for Baojiu's ICs)")
    ap.add_argument("--ngrid", type=int, default=2048)
    ap.add_argument("--zinit", type=float, default=49.0)
    ap.add_argument("--da", type=float, default=6e-4, help="initial step da0 (6e-4: legacy late-time ladder from z=49)")
    ap.add_argument("--ncheckpoint", type=int, default=50)
    ap.add_argument("--no-bdm", action="store_true")
    ap.add_argument("--catalogue-format", choices=["text", "binary", "both"], default=None,
                    help="BDM catalogue files (BDM.config CatalogueFormat, needs a PMP2main built from cz after 2026-10-18): "
                         "binary (.bin + .meta, the compiled default since 2026-10-18), text (legacy CatshortV.*.DAT) or both. "
                         "Not given: the key is not written and any build accepts the file; a build from 2026-10-18 on "
                         "then writes the binary only, an older build the text. Pass both (or text) only while a "
                         "consumer still needs the .DAT, e.g. halocat before it imports glamio.")
    ap.add_argument("--classes", choices=["all", "primary"], default="all",
                    help="per-output classes: 'all' = every output class 1 (BDM everywhere); "
                         "'primary' = the 17 primary redshifts class 1, the 25 stations class 2 (no BDM/extras)")
    ap.add_argument("--lens-slabs", type=int, default=0, help="lens planes: slabs per axis (0 = off; 4 = 256 Mpc/h at L=1024)")
    ap.add_argument("--lens-grid", type=int, default=0, help="lens planes: plane grid (0 = Ngrid)")
    ap.add_argument("--extra", action="store_true", help="extra outputs at class-1 outputs: low-k density and momentum modes (+ subsample, fR at halos if set)")
    ap.add_argument("--kc-delta", type=float, default=1.5, help="k_c [h/Mpc] of the density mode dump")
    ap.add_argument("--kc-mom", type=float, default=0.75, help="k_c [h/Mpc] of the momentum mode dumps; must stay below k_Nyquist(N_dep_mom)/4 (0.785 at N_dep = 1024, L = 1024), PMP2init refuses more")
    ap.add_argument("--ndep-mom", type=int, default=0, help="momentum deposit grid (0 = Ngrid/2)")
    ap.add_argument("--subsample-m", type=int, default=0, help="keep 1 in M particles (0 = off)")
    ap.add_argument("--fr-halos", action="store_true", help="f_R at halo centres (f(R) runs)")
    ap.add_argument("--check", action="store_true", help="S2: verify the epoch generator against the ECOSMOG lists")
    args = ap.parse_args()

    if args.check:
        return check_lists()
    if not args.out:
        sys.exit("--out is required")

    prov = {"generator": "DEGRACE-companion/gen_inputs.py", "git_hash": git_hash()}
    if args.node:
        rows = {r["node"]: r for r in csv.DictReader(open(ECOSMOG / "nodes.csv"))}
        node = f"Node_{args.node:03d}"
        if node not in rows:
            sys.exit(f"{node} not in nodes.csv")
        r = rows[node]
        om, h = float(r["Omega_m"]), float(r["h"])
        sigma8 = float(r["sigma8_ic"])           # what 2LPTic actually normalised to
        fr0 = float(r["fR0"])
        epochs = np.array(sorted(json.load(open(ECOSMOG / "aout.json"))[node]))
        pkfile = BL267_ICS / node / f"PK_CAMB_z00_{args.node:03d}.dat"
        prov.update(node=node, design={k: r[k] for k in ("Omega_m", "h", "S_8", "n_s", "log10_fR0", "sigma8_ic", "Omega_m_camb")},
                    ics=str(BL267_ICS / node / "ics."), epochs_source=str(ECOSMOG / "aout.json"),
                    pk_caveat="CAMB was run with omch2 = Omega_m h^2 (total Omega_m too high by Omega_b+Omega_nu); see docs/degrace_ecosmog/README.md section 3")
    else:
        if not (args.cosmology and args.pk):
            sys.exit("give --node, or --cosmology OMEGA_M H SIGMA8 with --pk")
        om, h, sigma8 = args.cosmology
        fr0 = args.fR0 if args.fR0 is not None else 0.0
        if args.gravity == "fr" and fr0 >= 0:
            sys.exit("--gravity fr needs --fR0 < 0")
        epochs = output_epochs(om)
        pkfile = args.pk
        prov.update(cosmology=dict(Omega_m=om, h=h, sigma8=sigma8, ombh2=args.ombh2, fR0=fr0),
                    epochs_source="computed: 17 primary z + 25 stations chi=256k Mpc/h, matter+Lambda")
    ombh2 = args.ombh2
    omch2 = om * h * h - ombh2                   # so that PMP2init's hubble = sqrt((ombh2+omch2)/Om) = h
    pk = read_pk(pkfile)
    a = dict(box=args.box, nrow=args.nrow, ngrid=args.ngrid, sigma8=sigma8, zinit=args.zinit, da=args.da,
             ncheckpoint=args.ncheckpoint, bdm=0 if args.no_bdm else 1,
             mg_flag=1 if args.gravity == "fr" else 0, fR0=fr0,
             classes=args.classes, lens_slabs=args.lens_slabs, lens_grid=args.lens_grid,
             extra=1 if args.extra else 0, kc_delta=args.kc_delta, kc_mom=args.kc_mom, ndep_mom=args.ndep_mom,
             subsample_m=args.subsample_m, fr_halos=1 if (args.fr_halos and args.gravity == "fr") else 0,
             catalogue_format=args.catalogue_format or "default (key not written: binary with a build from 2026-10-18 on, text before)")
    if args.lens_slabs and args.ngrid % args.lens_slabs:
        sys.exit("--lens-slabs must divide --ngrid")
    if args.lens_grid and args.ngrid % args.lens_grid:
        sys.exit("--lens-grid must divide --ngrid")
    if any(1.0 / e - 1.0 >= args.zinit for e in epochs):
        sys.exit("an output epoch is at or before z_init")
    out = args.out
    (out / "Run1" / "CATALOGS").mkdir(parents=True, exist_ok=True)   # BDM writes its catalogues and logs here
    (out / "PkTable.dat").write_text(pktable(ombh2, omch2, om, sigma8, pk))
    (out / "Init.dat").write_text(init_dat(a, epochs))
    (out / "Run1" / "BDM.config").write_text(bdm_config(args.catalogue_format))
    prov.update(pk_file=str(pkfile), pk_md5=md5(pkfile), gravity=args.gravity, init=a,
                epochs=[float(e) for e in epochs], hubble_expected=h, omch2=omch2, ombh2=ombh2)
    (out / "inputs_provenance.json").write_text(json.dumps(prov, indent=1))
    print(f"wrote {out}/PkTable.dat ({len(pk)} rows), Init.dat ({len(epochs)} exact outputs), "
          f"Run1/BDM.config (CatalogueFormat {args.catalogue_format or 'not written: compiled default binary'}), inputs_provenance.json")
    print(f"  Om={om} h={h} sigma8={sigma8} gravity={args.gravity} fR0={fr0 if a['mg_flag'] else 'n/a'}; "
          f"Nrow={args.nrow} Ngrid={args.ngrid} z_init={args.zinit} da0={args.da}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
