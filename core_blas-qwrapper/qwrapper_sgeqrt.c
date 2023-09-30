/**
 *
 * @file qwrapper_sgeqrt.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:39 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt)
{
    plasma_profile_by_kernel( task_flags, GEQRT );

    QUARK_Insert_Task(quark, CORE_sgeqrt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(float)*lda*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(float)*ldt*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(float)*nb,       NULL,          SCRATCH,
        sizeof(float)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgeqrt_quark = PCORE_sgeqrt_quark
#define CORE_sgeqrt_quark PCORE_sgeqrt_quark
#endif
void CORE_sgeqrt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    float *A;
    int lda;
    float *T;
    int ldt;
    float *TAU;
    float *WORK;

    quark_unpack_args_9(quark, m, n, ib, A, lda, T, ldt, TAU, WORK);
    CORE_sgeqrt(m, n, ib, A, lda, T, ldt, TAU, WORK);
}
