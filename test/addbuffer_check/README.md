# AddBuffer bit-identity check

`run.sh` builds `addbuffer_check.f90` against the objects of a full build in the
repository root and runs it at 1 and 4 threads. It fills 300,000 pseudo-random
particles (plus edge values that round to the box size), builds the halo
finder's periodically buffered analysis set along the historical path
(`RescaleCoords(1)` then `AddBuffer`) and along the inline finder's path since
2026-10-21 (`RetainParticles` then `AddBuffer(rescaleRetained=.true.)`, no
rescaled copy; memory audit cut 4), and requires: identical bits in every
buffered array and in `OriginalParticleId`, exact restoration of the
simulation arrays by `RemoveBuffer`, and a memory counter back at its start.
Expected output: `ADDBUFFER_CHECK PASS` twice.
