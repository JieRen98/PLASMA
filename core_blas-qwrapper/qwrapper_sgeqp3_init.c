/**
 *
 * @file qwrapper_sgeqp3_init.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgeqp3_init( Quark *quark, Quark_Task_Flags *task_flags,
                             int n, int *jpvt )
{
    Quark_Task *task;

    DAG_SET_PROPERTIES("init", "brown");
    task = QUARK_Task_Init( quark, CORE_sgeqp3_init_quark, task_flags );

    QUARK_Task_Pack_Arg( quark, task, sizeof(int),    &n,   VALUE  );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int)*n,  jpvt, OUTPUT );

    QUARK_Insert_Task_Packed( quark, task );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgeqp3_init_quark = PCORE_sgeqp3_init_quark
#define CORE_sgeqp3_init_quark PCORE_sgeqp3_init_quark
#endif
void CORE_sgeqp3_init_quark( Quark *quark )
{
    int n;
    int *jpvt;

    quark_unpack_args_2( quark, n, jpvt );
    CORE_sgeqp3_init(           n, jpvt );
}
