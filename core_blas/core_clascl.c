/**
 * @file core_clascl.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated c Fri Apr  1 11:02:32 2016
 *
 **/
#include <lapacke.h>
#include "common.h"

#define IN_LAPACKE

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_clascl scales all or part of a two-dimensional matrix A.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Specifies the type of the matrix A.
 *            = PlasmaGeneral : A is a general matrix
 *	      = PlasmaLowerTriangular : A is a lower triangular matrix
 *	      = PlasmaUpperTriangular : A is an upper triangular
 *	      matrix
 *	      = PlasmaUpperHessenberg : A is an upper Hessenberg
 *	      matrix
 *	      = PlasmaSymetricBandLowerStored : A is a symmetric band
 *	      matrix with lower bandwidth KL and upper bandwidth KU
 *	      and with the only the lower half stored
 *	      = PlasmaSymetricBandUpperStored : A is a symmetric band
 *	      matrix with lower bandwidth KL and upper bandwidth KU
 *	      and with the only the upper half stored
 *	      = PlasmaBand : A is a band matrix with lower bandwidth
 *	      KL and upper bandwidth KU. See ZGBTRF for storage
 *	      details.
 *
 * @param[in] kl is the lower bandwidth of A. Referenced only if type =
 *            PlasmaSymetricBandLowerStored, PlasmaSymetricBandUpperStored or
 *            PlasmaBand.
 *
 * @param[in] ku is the upper bandwidth of A. Referenced only if type =
 *            PlasmaSymetricBandLowerStored, PlasmaSymetricBandUpperStored or
 *            PlasmaBand.
 *
 * @param[in] cfrom is real
 *
 * @param[in] cto is real
 *            The matrix A is multiplied bt cto/cfrom. cfrom must be
 *            nonzero. The final result ctot*A(i,j)/cfrom is computed
 *            without over/underflow
 *
 * @param[in] m is the number of rows of the matrix A. m >= 0
 *
 * @param[in] n is the number of columns of the matrix A. n >= 0
 *
 * @param[in,out] A is the matrix to be multiplied by cto/cfrom
 *
 * @param[in] lda is the leading dimension of the array A. lda >= max(1,m).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int
CORE_clascl( PLASMA_enum type, int kl, int ku, float cfrom, float cto,
             int m, int n, PLASMA_Complex32_t *A, int lda )
{
    int info;
    if ( (cfrom == 0.0) || (cfrom != cfrom) ) {
        printf("error scale with %f\n", (float)cfrom );
        coreblas_error(-1, "error lascl\n");
    }

#if defined(IN_LAPACKE)
    /* To be used when integrated into lapacke */
    info = LAPACKE_clascl_work( LAPACK_COL_MAJOR, lapack_const(type), kl, ku, cfrom, cto,
                                m, n, A, lda );
#else
    PLASMA_FCALL(clascl, CLASCL)( &(lapack_const(type)), &kl, &ku, &cfrom, &cto,
                                  &m, &n, A, &lda, &info );
#endif

    if (info != 0){
        coreblas_error(info, "numerical error in clascl\n");
    }

    return info;
}
