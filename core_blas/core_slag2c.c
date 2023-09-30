/**
 * @file core_slag2c.c
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
#include "common.h"
#include <math.h>

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_slag2c - Copy a real matrix R into the real part of a complex matrix Z
 *
 *******************************************************************************
 *
 * @param[in] m
 *          m specifies the number of rows of the matrices R and Z.
 *
 * @param[in] n
 *          n specifies the number of columns of the matrices R and Z.
 *
 * @param[in] R
 *          R contains the real coefficient to be copied into the complex
 *          matrix. R is of size ldr -by- n.
 *
 * @param[in] ldr
 *          ldr specifies the leading dimension of R. ldr >= max(1,m).
 *
 * @param[out] Z
 *          On exit, the real part of the values of Z will be set to the ones
 *          stored in R, and the imagfinary part set to 0.
 *          Z is of size ldz -by- n.
 *
 * @param[in] ldz
 *          ldz specifies the leading dimension of Z.  ldz >= max(1,m).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slag2c = PCORE_slag2c
#define CORE_slag2c PCORE_slag2c
#endif
int CORE_slag2c( int m, int n,
                 const float *R, int ldr,
                 PLASMA_Complex32_t *Z, int ldz )
{
    int i, j;

    if (m < 0) {
        coreblas_error(1, "Illegal value of m");
        return -1;
    }
    if (n < 0) {
        coreblas_error(2, "Illegal value of n");
        return -2;
    }
    if ( (ldr < max(1,m)) && (m > 0) ) {
        coreblas_error(4, "Illegal value of ldr");
        return -4;
    }
    if ( (ldz < max(1,m)) && (m > 0) ) {
        coreblas_error(7, "Illegal value of ldz");
        return -7;
    }

    for (j=0; j<n; j++){
        for (i=0; i<m; i++, Z++, R++){
            *Z = *R + I * 0.;
        }
        Z += ldz - i;
        R += ldr - i;
    }

    return 0;
}
