/**
 * @file core_slaed3_computeW.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:33 2016
 *
 **/
#include "common.h"
#include <math.h>
#include <lapacke.h>

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_slaed2_computeK - Computes the number of deflated eigenvalues
 *  of a symmetric tridiagonal matrix.
 *  A eigenvalue is deflated if it's close enough to another eigenvalue
 *  or if the vector Z contains a tiny entry.
 *  This operation is useful when two subproblems, which were previously solved,
 *  are merge into a larger problem.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          n specifies the dimension of the symmetric tridiagonal matrix
 *
 * @param[in] K
 *          K specifies the number of non-deflated eigenvalues
 *
 * @param[in] Q
 *          On entry, Q contains the eigenvectors after slaed4
 *
 * @param[in] LDQ
 *          LDQ specifies the leading direction of Q
 *
 * @param[in] DLAMBDA
 *          DLAMBDA contains a copy of the first K eigenvalues
 *
 * @param[out] W
 *          W
 *
 * @param[out] INDX
 *          The permutation used to sort the contents of DLAMBDA into ascending order
 *
 * @param[in] start
 *          start specifies the first column index to be considered by this kernel
 *
 * @param[in] end
 *          end specifies the last column index to be considered by this kernel
 *          note that this index can be > n, it is supported
 *
 *******************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed3_computeW = PCORE_slaed3_computeW
#define CORE_slaed3_computeW PCORE_slaed3_computeW
#endif
void CORE_slaed3_computeW(int n, int K,
                          const float *Q, int LDQ,
                          const float *DLAMBDA, float *W,
                          const int *INDX,
                          int start, int end)
{
    int i, j, is;

    LAPACKE_slaset_work(LAPACK_COL_MAJOR, 'A', K, 1, 1., 1., W, 1 );

    if ( K <= 2 ){
        return;
    }

    /* Compute the subset of W */
    end = min(end, n);
    for (i=start; i<end; i++){
        is = INDX[i];
        if ( is < K ){
            for (j=0; j<K; j++){
                if ( is != j ){
                    W[j] = W[j] * Q[LDQ*is+j] / (DLAMBDA[j] - DLAMBDA[is]);
                }
            }
        }
    }
}
