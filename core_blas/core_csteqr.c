/**
 * @file core_csteqr.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gr�goire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated c Fri Apr  1 11:02:33 2016
 *
 **/
#include "common.h"
#include <lapacke.h>

/**
 *******************************************************************************
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_csteqr - solves the symmetric tridiagonal eigensystem using QR
 *
 *******************************************************************************
 *
 * @param[in] compz
 *          = PlasmaNoVec: computes eigenvalues only.
 *          = PlasmaVec: computes eigenpairs of the original symmetric
 *            matrix. On entry, Z must contain the orthogonal matrix used to
 *            reduce the original matrix to tridiagonal form.
 *          = PlasmaIVec: computes eigenpairs of the tridiagonal matrix. Z is
 *            initialized to the Identity Matrix.
 *
 * @param[in] n
 *          n specifies the order of the matrix. N >= 0
 *
 * @param[in,out] D
 *          On entry, D contains the diagonal elements of the tridiagonal matrix.
 *          On exit, D contains the eigenvalues sorted into increasing order.
 *
 * @param[in] E
 *          On entry, E contains the extra-diagonal elements of the tridiagonal
 *          matrix.
 *          On exit, E is destroyed.
 *
 * @param[in,out] Z
 *          On entry, Z has to be set to the Identity matrix.
 *          On exit, Z contains the eigenvectors.
 *
 * @param[in] LDZ
 *          LDZ specifies the leading direction of Z
 *
 * @param[in,out] WORK
 *          Workspace.
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_csteqr = PCORE_csteqr
#define CORE_csteqr PCORE_csteqr
#endif
int CORE_csteqr( PLASMA_enum compz, int n,
                 float *D, float *E,
                 PLASMA_Complex32_t *Z, int LDZ,
                 float *WORK )
{
    int info;

    if (WORK == NULL){
        info = LAPACKE_csteqr( LAPACK_COL_MAJOR, lapack_const(compz),
                               n, D, E, Z, LDZ );
    } else {
        info = LAPACKE_csteqr_work( LAPACK_COL_MAJOR, lapack_const(compz),
                                    n, D, E, Z, LDZ, WORK );
    }

    assert(!info);
    return info;
}
