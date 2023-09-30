/**
 *
 * @file core_dsetvar.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:31 2016
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 * CORE_dsetvar sets a single variable, x := alpha.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *         Scalar to set x to, passed by pointer so it can depend on runtime value.
 *
 * @param[out] x
 *         On exit, x = alpha.
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsetvar = PCORE_dsetvar
#define CORE_dsetvar PCORE_dsetvar
#endif
void CORE_dsetvar(const double *alpha, double *x)
{
    *x = *alpha;
}
