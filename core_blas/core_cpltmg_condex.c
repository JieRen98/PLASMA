/**
 *
 * @file core_cpltmg_condex.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:30 2016
 *
 **/
#include <math.h>
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cpltmg_condexq generates the Q used in condex matrix generation
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999898
 *  gallery('condex',n,4,100)
 *
 *  Returns a "counter-example" matrix to a condition estimator. It has order n
 *  and scalar parameter theta (default 100).
 *
 *  LAPACK (RCOND): It is the inverse of this matrix that is a counter-example.
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the matrix Q used in condex generation. M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A to be generated. N >= 0.
 *
 * @param[out] Q
 *         On entry, the M-by-3 matrix to be initialized.
 *         On exit, the housholder reflectors required for condex generation.
 *
 * @param[in] LDQ
 *         The leading dimension of the matrix Q. LDQ >= max(1,M).
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cpltmg_condexq = PCORE_cpltmg_condexq
#define CORE_cpltmg_condexq PCORE_cpltmg_condexq
#endif
void CORE_cpltmg_condexq( int M, int N, PLASMA_Complex32_t *Q, int LDQ )
{
    PLASMA_Complex32_t tau[3];
    PLASMA_Complex32_t *tQ = Q;
    int i;

    /* First column is [ 1 ... 1 ] */
    for( i=0; i < M; i++, tQ++ )
        *tQ = (PLASMA_Complex32_t)1.0;

    /* Second column is [1 0 0 ... 0] */
    tQ = Q + LDQ;
    *tQ = (PLASMA_Complex32_t)1.;
    tQ++;
    memset( tQ, 0, (M-1) * sizeof(PLASMA_Complex32_t) );

    /* Third column is  (-1)^i * (1. + i / (N-1)) */
    tQ = Q + 2 * LDQ;
    for( i=0; i<M; i++, tQ++ )
        *tQ = (PLASMA_Complex32_t)( cpowf( -1.0, (float)i ) * (1.0 + (float)i/(N-1) ) );

    /* Generate orthogonal projector */
    LAPACKE_cgeqrf( LAPACK_COL_MAJOR, M, 3,    Q, LDQ, tau );
    LAPACKE_cungqr( LAPACK_COL_MAJOR, M, 3, 3, Q, LDQ, tau );

    return;
}
