/**
 *
 * @file qwrapper_slaed3_computevectors.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include <stdlib.h>
#include <lapacke.h>
#include "common.h"

void
CORE_slaed3_compvec_quark(Quark *quark);

void
CORE_slaed3_compvec_ws3_quark(Quark *quark);

void
CORE_slaed3_wscopy_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed3_computevectors(Quark *quark, Quark_Task_Flags *task_flags,
                                      int wsmode, int n, const int *K,
                      const int *il_nondef, const int *iu_nondef,
                                      float *Q, int LDQ, float *W, const int *INDXC,
                                      float **WSglobal, float **WSlocal,
                                      int start, int end )
{
    plasma_profile_by_kernel( task_flags, LAED3_COMPVEC );

    if ( wsmode == 3 ) {
        /* Use ws as the Q pointer */
        QUARK_Insert_Task(quark, CORE_slaed3_compvec_ws3_quark, task_flags,
                          sizeof(int),       K,              NODEP,
                          sizeof(int),       il_nondef,      NODEP,
                          sizeof(int),       iu_nondef,      NODEP,
                          sizeof(float),    WSglobal,       NODEP,
                          sizeof(float),    W,                  INPUT,
                          sizeof(float)*n,  NULL,               SCRATCH,
                          sizeof(int),       INDXC,          NODEP,
                          sizeof(int),      &start,          VALUE,
                          sizeof(int),      &end,            VALUE,
                          sizeof(float),    Q+start*LDQ,        INOUT,
                          sizeof(float),    WSlocal,            OUTPUT,
                          0);
    }
    else {
        QUARK_Insert_Task(quark, CORE_slaed3_compvec_quark, task_flags,
                          sizeof(int),       K,                  INPUT,
                          sizeof(int),       il_nondef,          INPUT,
                          sizeof(int),       iu_nondef,          INPUT,
                          sizeof(float),    Q,                  INPUT,
                          sizeof(int),      &LDQ,            VALUE,
                          sizeof(float),    W,                  INPUT,
                          sizeof(float)*n,  NULL,               SCRATCH,
                          sizeof(int),       INDXC,              NODEP,
                          sizeof(int),      &start,          VALUE,
                          sizeof(int),      &end,            VALUE,
                          sizeof(float),    Q+start*LDQ,        INOUT,
                          sizeof(float),    WSlocal,            ((wsmode==0) ? NODEP : OUTPUT),
                          0);
    }
}
/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed3_compvec_ws3_quark = PCORE_slaed3_compvec_ws3_quark
#define CORE_slaed3_compvec_ws3_quark PCORE_slaed3_compvec_ws3_quark
#endif
void CORE_slaed3_compvec_ws3_quark(Quark *quark)
{
    const int *K;
    const int *il_nondef;
    const int *iu_nondef;
    float **WSg;
    float *W;
    float *S;
    const int *INDXC;
    int start;
    int end;
    void *fake, *fake2;
    int il, iu;

    quark_unpack_args_11(quark, K, il_nondef, iu_nondef,
                         WSg, W, S, INDXC,
                         start, end, fake, fake2 );

    il = ( il_nondef == NULL ) ? 0  : *il_nondef;
    iu = ( iu_nondef == NULL ) ? *K : *iu_nondef;
    CORE_slaed3_computevectors(*K, il, iu,
                               (*WSg), *K, W, S, INDXC,
                               start, end);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed3_compvec_quark = PCORE_slaed3_compvec_quark
#define CORE_slaed3_compvec_quark PCORE_slaed3_compvec_quark
#endif
void CORE_slaed3_compvec_quark(Quark *quark)
{
    const int *K;
    const int *il_nondef;
    const int *iu_nondef;
    float *Q;
    int LDQ;
    float *W;
    float *S;
    const int *INDXC;
    int start;
    int end;
    void *fake, *fake2;
    int il, iu;

    quark_unpack_args_12(quark, K, il_nondef, iu_nondef,
                         Q, LDQ, W, S, INDXC,
                         start, end, fake, fake2 );

    il = ( il_nondef == NULL ) ? 0  : *il_nondef;
    iu = ( iu_nondef == NULL ) ? *K : *iu_nondef;
    CORE_slaed3_computevectors(*K, il, iu,
                               Q, LDQ, W, S, INDXC,
                               start, end);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed3_wscopy( Quark *quark, Quark_Task_Flags *task_flags,
                               const int *K, const int *il_nondef, const int *iu_nondef,
                               const float *Q, int LDQ, float **WORK,
                               int start, int end )
{
    QUARK_Insert_Task(quark, CORE_slaed3_wscopy_quark, task_flags,
        sizeof(int),     K,               INPUT,
        sizeof(int),     il_nondef,       NODEP,
        sizeof(int),     iu_nondef,       NODEP,
        sizeof(float),  Q,               NODEP,
        sizeof(int),    &LDQ,           VALUE,
        sizeof(float),  WORK,            INOUT,
        sizeof(int),    &start,         VALUE,
        sizeof(int),    &end,           VALUE,
        sizeof(float),  Q+start*LDQ,     INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed3_wscopy_quark = PCORE_slaed3_wscopy_quark
#define CORE_slaed3_wscopy_quark PCORE_slaed3_wscopy_quark
#endif
void CORE_slaed3_wscopy_quark(Quark *quark)
{
    const int *K;
    const int *il_nondef;
    const int *iu_nondef;
    const float *Q;
    int LDQ;
    float **WORK;
    int start;
    int end;
    int size;
    void *fake;
    int il, iu;

    quark_unpack_args_9(quark, K, il_nondef, iu_nondef,
                        Q, LDQ, WORK,
                        start, end, fake);

    il = ( il_nondef == NULL ) ? 0  : *il_nondef;
    iu = ( iu_nondef == NULL ) ? *K : *iu_nondef;

    /* Compute the size */
    start = max( 0,   max(  start, il ) );
    end   = min( end, min( *K,     iu ) );
    size  = max( 0, end-start );

    if ((size > 0) && (*K > 0)) {
        *WORK = malloc( size * (*K) * sizeof(float) );
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,
                            lapack_const(PlasmaUpperLower), *K, size,
                            Q+start*LDQ, LDQ, *WORK, *K);
    } else {
        *WORK = NULL;
    }
}
