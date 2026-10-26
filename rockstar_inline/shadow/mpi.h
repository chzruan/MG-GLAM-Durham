/* shadow/mpi.h: the in-line Rockstar adapter needs no MPI.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
 *
 * Rockstar's error.h includes <mpi.h> for nothing it uses on the kernel path.
 * This header (first on the include path of the Rockstar sources) declares the
 * two types only and defines no MPI symbol: the GLAM *MPI executables link the
 * same objects next to a real MPI library, which must stay the only definer. */
#ifndef RS_SHADOW_MPI_H
#define RS_SHADOW_MPI_H
typedef int MPI_Comm;
typedef int MPI_Datatype;
#endif
