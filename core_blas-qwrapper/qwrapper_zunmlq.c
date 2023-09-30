/**
 *
 * @file qwrapper_zunmlq.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zunmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const PLASMA_Complex64_t *A, int lda,
                       const PLASMA_Complex64_t *T, int ldt,
                       PLASMA_Complex64_t *C, int ldc)
{
    plasma_profile_by_kernel( task_flags, UNMLQ );

    QUARK_Insert_Task(quark, CORE_zunmlq_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,             INPUT | QUARK_REGION_U,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    C,             INOUT,
        sizeof(int),                        &ldc,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &nb,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zunmlq_quark = PCORE_zunmlq_quark
#define CORE_zunmlq_quark PCORE_zunmlq_quark
#endif
void CORE_zunmlq_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum trans;
    int m;
    int n;
    int k;
    int ib;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *T;
    int ldt;
    PLASMA_Complex64_t *C;
    int ldc;
    PLASMA_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_14(quark, side, trans, m, n, k, ib,
                         A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_zunmlq(side, trans, m, n, k, ib,
                A, lda, T, ldt, C, ldc, WORK, ldwork);
}
