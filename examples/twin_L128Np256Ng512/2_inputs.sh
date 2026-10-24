#!/bin/bash
# Step 2: the input files of both twins with DEGRACE-companion/gen_inputs.py, plus a frozen copy of the executables.
# Writes $WORK/<twin>/: Init.dat, PkTable.dat, inputs_provenance.json, Run<N>/BDM.config and bin/ (the executables,
# md5sum.txt, git_hash.txt, main_exe.txt). The twins differ only in gravity: MG_flag, fR0 and f_R at haloes.
set -euo pipefail
source "$(dirname "$0")/config.sh"
for exe in PMP2init.exe ic2pm.exe PMP2main.exe PMP2MG.exe; do
    [ -x "$GLAM_BIN/$exe" ] || { echo "$GLAM_BIN/$exe not found: build GLAM first (see README)"; exit 1; }
done
common=(--cosmology "$OMEGA_M" "$HUBBLE" "$SIGMA8" --ombh2 "$OMBH2" --pk "$PK_Z0"
        --box "$BOX" --nrow "$NP" --ngrid "$NG" --zinit "$ZINIT" --run "$REALISATION"
        --redshifts $FULL_Z --lens-slabs 1 --extra --subsample-m 100)
[ -n "${STATION_Z// /}" ] && common+=(--stations $STATION_Z)
for twin in $TWINS; do
    run=$WORK/$twin
    if [ -f "$run/Init.dat" ]; then echo "$run/Init.dat exists: skipped (to change the settings, use a new WORK)"; continue; fi
    case $twin in
        LCDM) gravity=(--gravity lcdm);                       exe=PMP2main.exe ;;
        fR5)  gravity=(--gravity fr --fR0 "$FR0" --fr-halos); exe=PMP2MG.exe ;;
        *)    echo "unknown twin $twin (expected LCDM or fR5)"; exit 1 ;;
    esac
    $PYTHON "$REPO/DEGRACE-companion/gen_inputs.py" "${common[@]}" "${gravity[@]}" --out "$run"
    mkdir -p "$run/bin"
    cp "$GLAM_BIN/PMP2init.exe" "$GLAM_BIN/ic2pm.exe" "$GLAM_BIN/$exe" "$run/bin/"
    (cd "$run/bin" && md5sum *.exe > md5sum.txt)
    echo "$(git -C "$REPO" rev-parse --short HEAD) (repository HEAD when the executables were copied from $GLAM_BIN; md5sum.txt identifies the build)" \
        > "$run/bin/git_hash.txt"
    echo "$exe" > "$run/bin/main_exe.txt"
done
echo "inputs written: $(for t in $TWINS; do printf '%s ' "$WORK/$t"; done)"
