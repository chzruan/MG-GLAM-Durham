"""Generate the convergence-only native IC variant without editing production source."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise ValueError(f"Expected exactly one source anchor: {old!r}")
    return text.replace(old, new)


def generate(root: Path) -> str:
    text = (root / "PMP2start.f90").read_text()
    text = text[:text.index("!-------------------------------------------------------------------------\nPROGRAM PMstartMp")]
    text = once(text, "    use Tools\n", "    use Tools\n    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite\n")
    text = once(text, "    Real*8 :: SKINE\n", """    Real*8 :: SKINE
    ! Explicit, mandatory campaign controls; production PMP2start is unchanged.
    integer :: ic_master_nrow = 1024, ic_origin_ngrid = 2048
    real*4 :: ic_alpha = -1.
    real*4 :: AEXP0 = 0., AU0 = 0.
    logical :: ic_normalize_only = .false.
    real*8 :: ic_spectrum_sum = 0.d0
""")
    text = once(text, "    Integer*4, PARAMETER :: nbyteword = 4", "    Integer*4, PARAMETER :: nbyteword = 1")
    text = once(text, "        if (.not. exst) Call SetSeeds", "        if (.not. exst) error stop 'Provide explicit ../TableSeeds.dat'")
    text = once(text, "Contains\n", "Contains\n" + (HERE / "controls.f90").read_text())
    text = once(text, "        gSet = 0.\n", "        gSet = 0.\n        iFlag = 0 ! Defined first Box-Muller state for every k plane.\n")
    text = once(text, "        myMemory = Memory(9_8*Nparticles)\n", "        if (.not. ic_normalize_only) then\n        myMemory = Memory(9_8*Nparticles)\n")
    text = once(text, "        ALLOCATE (GRZ(NROW, NROW, NROW))     !\n", "        ALLOCATE (GRZ(NROW, NROW, NROW))     !\n        end if\n")
    text = once(text, "        INTEGER*4, allocatable, DIMENSION(:) :: mapz, map3\n", """        INTEGER*4, allocatable, DIMENSION(:) :: mapz, map3
        real*8, allocatable :: plane_sum(:)
        real*8 :: local_sum
""")
    text = once(text, "        allocate (mapz(NROW), map3(NROW))\n", "        allocate (mapz(NROW), map3(NROW), plane_sum(NROW))\n        plane_sum = 0.d0\n")
    text = once(text, "!$OMP PRIVATE(Mk3,Mk2,Mk1)\n        DO Mk3", "!$OMP PRIVATE(Mk3,Mk2,Mk1)\n        DO Mk3")
    text = once(text, "!$OMP PARALLEL DO DEFAULT(SHARED)  &\n!$OMP PRIVATE(Mk3,Mk2,Mk1)", "        if (.not. ic_normalize_only) then\n!$OMP PARALLEL DO DEFAULT(SHARED)  &\n!$OMP PRIVATE(Mk3,Mk2,Mk1)")
    text = once(text, "!                                      Set Spectrum\n", "        end if\n!                                      Set Spectrum\n")
    text = once(text, "!$OMP PRIVATE(WD,Wk,TS,TRX,m,NRAND,gSet,iFlag)           &\n!$OMP REDUCTION(+:SUMM)", "!$OMP PRIVATE(WD,Wk,TS,TRX,m,NRAND,gSet,iFlag,local_sum)")
    text = once(text, "        DO k = 1, NROW\n            i24", "        DO k = 1, NROW\n            local_sum = 0.d0\n            i24")
    text = once(text, "                DO i = 1, NROW\n                    Wi3 = map3(i)**2\n                    TS = GAUSS3(gSet, iFlag)\n", """                DO i = 1, ic_master_nrow
                    ! Preserve the master luxury draw at each packed (i,j,k),
                    ! including skipped coefficients and the terminal column.
                    TS = GAUSS3(gSet, iFlag)
                    if (i >= NROW .or. j == NROW .or. k == NROW) cycle
                    Wi3 = map3(i)**2
""")
    text = once(text, "                        GRX(1, 1, 1) = 0.\n                        GRY(1, 1, 1) = 0.\n                        GRZ(1, 1, 1) = 0.\n", "                        cycle ! Uniform mode is zero; arrays were cleared.\n")
    text = once(text, "                        GRX(mapz(i), j, k) = TRX*map3(i)\n", "                        if (.not. ic_normalize_only) then\n                        GRX(mapz(i), j, k) = TRX*map3(i)\n")
    text = once(text, "                        SUMM = SUMM + TS**2\n", "                        end if\n                        local_sum = local_sum + TS**2\n")
    text = once(text, "        END DO               ! k\n", "            plane_sum(k) = local_sum\n        END DO               ! k\n        do k = 1, NROW\n            SUMM = SUMM + plane_sum(k)\n        end do\n        ic_spectrum_sum = SUMM\n")
    text = once(text, "        ALPHA = AMPLT/SQRT(SUMM)*sqrt(8.)\n", """        if (SUMM <= 0.d0 .or. .not. ieee_is_finite(SUMM)) error stop 'Invalid spectrum sum'
        if (ic_alpha < 0.) then
            ALPHA = AMPLT/SQRT(SUMM)*sqrt(8.)
        else
            ALPHA = ic_alpha
        end if
        if (.not. ieee_is_finite(ALPHA) .or. ALPHA <= 0.) error stop 'Invalid alpha'
        if (ic_normalize_only) then
            deallocate(mapz, map3, plane_sum)
            CALL Timing(3, 1)
            return
        end if
""")
    text = once(text, "        deallocate (mapz, map3)\n", "        deallocate (mapz, map3, plane_sum)\n")
    text = once(text, "        xShift = 0.5\n", "        xShift = 0.5*real(NGRID)/real(ic_origin_ngrid)\n")
    text = once(text, "            sqrt(sDispl/max(Icurrent, 1))\n", "            sqrt(sDispl/max(Nparticles, 1_8))\n")
    text = once(text, "        lenr = Ngrid\n", "        lenr = NROW ! Actual transform length, independent of PM mesh.\n")
    return text + (HERE / "entry.f90").read_text()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    text = generate(args.repo)
    args.output.write_text(text)
    receipt = {
        "generator": "native-luxury-master-stride-v1",
        "production_source_sha256": hashlib.sha256((args.repo / "PMP2start.f90").read_bytes()).hexdigest(),
        "generated_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "template_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (HERE / "generate.py", HERE / "controls.f90", HERE / "entry.f90")},
    }
    args.output.with_suffix(".generation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
