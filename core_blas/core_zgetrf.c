/**
 *
 * @file core_zgetrf.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zgetrf - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm with partial tile pivoting with row interchanges.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrix A. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A. n >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] lda
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations.
 *
 * @param[out] info
 *          - 0 on successful exit
 *          - <0 if -i, the i-th argument had an illegal value
 *          - >0 if i, U(i,i) is exactly zero. The factorization has been
 *            completed, but the factor U is exactly singular, and division by
 *            zero will occur if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgetrf = PCORE_zgetrf
#define CORE_zgetrf PCORE_zgetrf
#endif
int CORE_zgetrf(int m, int n,
                 PLASMA_Complex64_t *A, int lda,
                 int *IPIV, int *info)
{
    *info = LAPACKE_zgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV );
    return PLASMA_SUCCESS;
}
