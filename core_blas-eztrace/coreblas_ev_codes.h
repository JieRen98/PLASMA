/**
 *
 * @file coreblas_ev_codes.h
 *
 *  PLASMA core_blas tracing kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#ifndef EVCODES_COREBLAS_H
#define EVCODES_COREBLAS_H

#define COREBLAS_EVENTS_ID    USER_MODULE_ID(0x01)
#define COREBLAS_PREFIX       GENERATE_USER_MODULE_PREFIX(COREBLAS_EVENTS_ID)

#define PREFIX_MASK  (((1 << NB_BITS_PREFIX) -1) << NB_BITS_EVENTS)
#define EVENTS_MASK  (~PREFIX_MASK)

#define IS_A_COREBLAS_EV(ev)  ((LITL_READ_GET_CODE(ev) & PREFIX_MASK) == COREBLAS_PREFIX)
#define COREBLAS_GET_CODE(ev) (LITL_READ_GET_CODE(ev) & EVENTS_MASK)
#define FUT_COREBLAS(event)   (COREBLAS_PREFIX | COREBLAS_##event )

enum coreblas_ev_code_e {
    COREBLAS_STOP,
    COREBLAS_FREE,
    COREBLAS_FOO,
    COREBLAS_GEMM,
    COREBLAS_HERK,
    COREBLAS_SYRK,
    COREBLAS_HEMM,
    COREBLAS_SYMM,
    COREBLAS_TRMM,
    COREBLAS_TRSM,
    COREBLAS_HER2K,
    COREBLAS_SYR2K,
    COREBLAS_GEMV,
    COREBLAS_GBMV,
    COREBLAS_HEMV,
    COREBLAS_HBMV,
    COREBLAS_HPMV,
    COREBLAS_SYMV,
    COREBLAS_SBMV,
    COREBLAS_SPMV,
    COREBLAS_TRMV,
    COREBLAS_TBMV,
    COREBLAS_TPMV,
    COREBLAS_TRSV,
    COREBLAS_TBSV,
    COREBLAS_TPSV,
    COREBLAS_GER,
    COREBLAS_GERU,
    COREBLAS_GERC,
    COREBLAS_HER,
    COREBLAS_HPR,
    COREBLAS_HER2,
    COREBLAS_HPR2,
    COREBLAS_SYR,
    COREBLAS_SPR,
    COREBLAS_SYR2,
    COREBLAS_SPR2,
    COREBLAS_ROTG,
    COREBLAS_ROTMG,
    COREBLAS_ROT,
    COREBLAS_ROTM,
    COREBLAS_SWAP,
    COREBLAS_SCAL,
    COREBLAS_COPY,
    COREBLAS_AXPY,
    COREBLAS_GEADD,
    COREBLAS_DOT,
    COREBLAS_DOTU,
    COREBLAS_DOTC,
    COREBLAS_XDOT,
    COREBLAS_NRM2,
    COREBLAS_ASUM,
    COREBLAS_AMAX,
    COREBLAS_LACPY,
    COREBLAS_LANGE,
    COREBLAS_LANHE,
    COREBLAS_LANSY,
    COREBLAS_LARFB,
    COREBLAS_LARFT,
    COREBLAS_LASWP,
    COREBLAS_LAUUM,
    COREBLAS_POTRF,
    COREBLAS_TRTRI,
    COREBLAS_LASET,
    COREBLAS_LASSQ,
    COREBLAS_GELQT,
    COREBLAS_GEQRT,
    COREBLAS_GESSM,
    COREBLAS_GETRF,
    COREBLAS_LATRO,
    COREBLAS_SSSSM,
    COREBLAS_TITRO,
    COREBLAS_TRBMM,
    COREBLAS_TRGMM,
    COREBLAS_TSLQT,
    COREBLAS_TSMLQ,
    COREBLAS_TSMQR,
    COREBLAS_TSQRT,
    COREBLAS_TSRFB,
    COREBLAS_TSTRF,
    COREBLAS_TTLQT,
    COREBLAS_TTMLQ,
    COREBLAS_TTMQR,
    COREBLAS_TTQRT,
    COREBLAS_TTRFB,
    COREBLAS_UNMLQ,
    COREBLAS_UNMQR,
    COREBLAS_GETRIP,
    COREBLAS_PLGHE,
    COREBLAS_PLGSY,
    COREBLAS_SHIFT,
    COREBLAS_SHIFTW,
    COREBLAS_SWPAB,
    COREBLAS_PLRNT,
    COREBLAS_PEMV,
    COREBLAS_BRDALG,
    COREBLAS_TRDALG,
    COREBLAS_HEGST,
    COREBLAS_SYGST,
    COREBLAS_HERFB,
    COREBLAS_SYRFB,
    COREBLAS_LARFG,
    COREBLAS_GEQP3_INIT,
    COREBLAS_GEQP3_NORMS,
    COREBLAS_GEQP3_PIVOT,
    COREBLAS_GEQP3_UPDATE,
    COREBLAS_PIVOT_UPDATE,
    COREBLAS_SETVAR,
    COREBLAS_STEDC,
    COREBLAS_STEQR,
    COREBLAS_LASCL,
    COREBLAS_LAED2_COMPRESSQ,
    COREBLAS_LAED2_COMPUTEK,
    COREBLAS_LAED2_COPYDEF,
    COREBLAS_LAED3_COMPVEC,
    COREBLAS_LAED3_WSCOPY,
    COREBLAS_LAED3_COMPW,
    COREBLAS_LAED3_REDUCEW,
    COREBLAS_LAED3_UPDATEVECTORS,
    COREBLAS_LAED4,
    COREBLAS_LAG2Z,
    COREBLAS_TYPE1,
    COREBLAS_TYPE2,
    COREBLAS_TYPE3,
    COREBLAS_NBMAX_EVENTS,
};

#endif /* COREBLAS_CODES_H */
