/*
 * gorilla_cuda.h - frozen C ABI between the GORILLA Fortran batch host and the
 * CUDA device backend.
 *
 * Layout rules
 *   - All arrays are plain contiguous buffers, Fortran column major. A Fortran
 *     array mesh(ntetr,NGORILLA_MESH_NFIELD) is addressed in C as
 *     mesh[field*ntetr + cell] with 0-based cell and field indices.
 *   - No Fortran derived type and no native Fortran LOGICAL crosses this ABI.
 *   - Particle and topology indices keep the Fortran 1-based values used inside
 *     GORILLA. A negative neighbour tetrahedron denotes a physical boundary.
 *   - Every function returns an int error code, 0 (NGORILLA_ERR_OK) on success.
 */
#ifndef GORILLA_CUDA_H
#define GORILLA_CUDA_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define NGORILLA_CUDA_ABI_VERSION 1

/* Number of double fields per mesh cell in the packed mesh buffer. */
#define NGORILLA_MESH_NFIELD 59
/* Number of integer fields per mesh cell in the packed topology buffer. */
#define NGORILLA_TOPO_NFIELD 16
/* Number of double fields per marker in the packed run state buffer. */
#define NGORILLA_RSTATE_NFIELD 9
/* Number of integer fields per marker in the packed run state buffer. */
#define NGORILLA_ISTATE_NFIELD 4
/* Number of double entries of the timing profile, in milliseconds. */
#define NGORILLA_PROFILE_NFIELD 5

/* Field offsets into the packed mesh buffer (double). */
#define NGORILLA_MESH_X1          0  /* x1[3]      coordinates of first vertex */
#define NGORILLA_MESH_ANORM       3  /* anorm[3*4] face-major inward normals   */
#define NGORILLA_MESH_DIST_REF    15 /* dist_ref   |n_1|*distance(v1,face1)    */
#define NGORILLA_MESH_DT_DTAU     16 /* dt_dtau_const (0th order time tracing) */
#define NGORILLA_MESH_B1          17 /* bmod1      |B| at first vertex          */
#define NGORILLA_MESH_GB          18 /* gB[3]      gradient of |B|              */
#define NGORILLA_MESH_PHI1        21 /* Phi1       electrostatic potential      */
#define NGORILLA_MESH_GPHI        22 /* gPhi[3]    gradient of Phi             */
#define NGORILLA_MESH_CURLH       25 /* curlh[3]                                */
#define NGORILLA_MESH_GBXH1       28 /* gBxh1[3]                                */
#define NGORILLA_MESH_GPHIXH1     31 /* gPhixh1[3]                              */
#define NGORILLA_MESH_GBXCURLA    34 /* gBxcurlA   scalar                       */
#define NGORILLA_MESH_PHIXCURLA   35 /* gPhixcurlA scalar                       */
#define NGORILLA_MESH_ALPMAT      36 /* alpmat[3*3] column major              */
#define NGORILLA_MESH_BETMAT      45 /* betmat[3*3] column major             */
#define NGORILLA_MESH_SPALPMAT    54 /* spalpmat                              */
#define NGORILLA_MESH_SPBETMAT    55 /* spbetmat                              */
#define NGORILLA_MESH_CURLA       56 /* curlA[3]                              */

/* Field offsets into the packed topology buffer (int64_t, 1-based indices). */
#define NGORILLA_TOPO_NEIGHBOUR_TETR 0  /* neighbour_tetr[4]        */
#define NGORILLA_TOPO_NEIGHBOUR_FACE 4  /* neighbour_face[4]        */
#define NGORILLA_TOPO_PERBOU_PHI     8  /* neighbour_perbou_phi[4]  */
#define NGORILLA_TOPO_PERBOU_THETA   12 /* neighbour_perbou_theta[4]*/

/* Field offsets into the packed marker state buffer (double). */
#define NGORILLA_RSTATE_X        0 /* x[3]        position in coord_system      */
#define NGORILLA_RSTATE_VPAR     3 /* vpar        parallel velocity             */
#define NGORILLA_RSTATE_PERPINV  4 /* perpinv = -0.5*vperp^2/B, set ONCE per
                                       orbital time step, never recomputed          */
#define NGORILLA_RSTATE_TREMAIN  5 /* signed remaining time of the time step    */
#define NGORILLA_RSTATE_LAST_Z   6 /* last_z[3] relative position in last cell  */

/* Field offsets into the packed marker state buffer (int64_t). */
#define NGORILLA_ISTATE_TETR      0 /* current (or last) tetrahedron, 1-based   */
#define NGORILLA_ISTATE_IFACE     1 /* entry face, 0 if the step finished inside */
#define NGORILLA_ISTATE_LAST_TETR 2 /* cell of last_z, kept also for losses      */
#define NGORILLA_ISTATE_STATUS    3 /* NGORILLA_STATUS_*                          */

/* Marker status, identical to SRC/gorilla_status_mod.f90. */
enum {
    NGORILLA_STATUS_ACTIVE            = 0,
    NGORILLA_STATUS_COMPLETE          = 1,
    NGORILLA_STATUS_PHYSICAL_LOSS     = 2,
    NGORILLA_STATUS_REPLAY_CELL       = 3,
    NGORILLA_STATUS_NUMERICAL_FAILURE = 4
};

/* Error codes returned by all three ABI functions. */
enum {
    NGORILLA_ERR_OK          = 0,
    NGORILLA_ERR_BAD_ARG     = 1, /* null pointer, wrong length, n < 0        */
    NGORILLA_ERR_UNSUPPORTED = 2, /* settings outside the supported device mode */
    NGORILLA_ERR_NO_MESH     = 3, /* create called without a loaded mesh       */
    NGORILLA_ERR_NO_CUDA     = 4, /* backend requested but not compiled in     */
    NGORILLA_ERR_LAUNCH      = 5, /* device allocation or kernel launch failed  */
    NGORILLA_ERR_STALE_MESH  = 6  /* handle does not match the live mesh       */
};

/*
 * Creates a device backend handle and uploads the packed mesh.
 *
 *   ntetr        number of tetrahedra, length of mesh and topology buffers
 *   mesh         double[NGORILLA_MESH_NFIELD*ntetr]
 *   topology     int64_t[NGORILLA_TOPO_NFIELD*ntetr]
 *   coord_system 1 (R,phi,z) or 2 (s,theta,phi)
 *   nperiods     n_field_periods
 *   sign_sqg     sign of the metric determinant sqrt(g)
 *   cm_over_e    mass-to-charge ratio of the species divided by e
 *   clight       speed of light in Gaussian units
 *   eps          root-solver tolerance from the constants module
 *   handle       receives the opaque handle, set to NULL on failure
 */
int gorilla_cuda_create(int64_t ntetr, const double *mesh, const int64_t *topology,
                        int coord_system, int nperiods, int sign_sqg,
                        double cm_over_e, double clight, double eps, void **handle);

/*
 * Advances n markers for as many cell crossings as allowed.
 *
 *   n                number of markers
 *   rstate           double[NGORILLA_RSTATE_NFIELD*n], in and out
 *   istate           int64_t[NGORILLA_ISTATE_NFIELD*n], in and out
 *   crossing_counts  int64_t[n], accepted cell crossings per marker, in and out
 *   max_crossings    upper bound of crossings per call; markers that reach it stay
 *                    NGORILLA_STATUS_ACTIVE and are resumed by the next call, they
 *                    are never reported as complete or lost
 *   profile_ms       double[NGORILLA_PROFILE_NFIELD]: allocate, upload, kernel,
 *                    download, total for this device call; host replay is timed separately
 *
 * A marker whose cell cannot be decided on the device is returned with
 * NGORILLA_STATUS_REPLAY_CELL and its entry state untouched; the host replays
 * that cell on the CPU and resumes the marker.
 */
int gorilla_cuda_advance(void *handle, int64_t n, double *rstate, int64_t *istate,
                         int64_t *crossing_counts, int64_t max_crossings,
                         double *profile_ms);

/* Destroys the handle and frees all device memory. NULL is accepted. */
int gorilla_cuda_destroy(void **handle);

#ifdef __cplusplus
}
#endif

#endif /* GORILLA_CUDA_H */
