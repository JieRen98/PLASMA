/**
 *
 * @file coreblas.c
 *
 *  PLASMA core_blas tracing kernels
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * This file provides the init and conclude functions for
 * EZTrace module.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/

#define _GNU_SOURCE 1 /* or _BSD_SOURCE or _SVID_SOURCE */
#define _REENTRANT

#include <stdlib.h>
#include <unistd.h>
#include <sys/syscall.h>

#include <eztrace.h>
#include <ev_codes.h>

#include "common.h"
#include "coreblas_ev_codes.h"
#include "coreblas_macros.h"

FUNCTION_VOID( CORE_pivot_update, PIVOT_UPDATE, void ,
               (int m, int n, int *ipiv, int *indices, int offset, int init),
               (m, n, ipiv, indices, offset, init) );

static void __coreblas_init (void) __attribute__ ((constructor));
static void
__coreblas_init (void)
{
#ifdef EZTRACE_AUTOSTART
  eztrace_start ();
#endif
}

static void __coreblas_conclude (void) __attribute__ ((destructor));
static void
__coreblas_conclude (void)
{
  eztrace_stop ();
}
