/**
 *
 * @file qwrapper_z.c
 *
 *  PLASMA core_blas tracing kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This file provides the wrapper for each function of the
 *  core_blas library which will generate an event before and
 *  after the execution of the kernel.
 *  This file is automatically generated with convert2eztrace.pl
 *  script. DO NOT MANUALLY EDIT THIS FILE.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <eztrace.h>
#include <ev_codes.h>
#include "common.h"
#include "coreblas_ev_codes.h"
#include "coreblas_macros.h"
#undef REAL
#define COMPLEX
#ifdef REAL
FUNCTION_QUARK( CORE_dlaed0_lascl_quark, LASCL )
FUNCTION_QUARK( CORE_dlaed2_compressq_quark, LAED2_COMPRESSQ )
FUNCTION_QUARK( CORE_dlaed2_computeK_quark, LAED2_COMPUTEK )
FUNCTION_QUARK( CORE_dlaed2_copydef_quark, LAED2_COPYDEF )
FUNCTION_QUARK( CORE_dlaed3_compvec_ws3_quark, LAED3_COMPVEC )
FUNCTION_QUARK( CORE_dlaed3_compvec_quark, LAED3_COMPVEC )
FUNCTION_QUARK( CORE_dlaed3_wscopy_quark, LAED3_WSCOPY )
FUNCTION_QUARK( CORE_dlaed3_compW_p2f1_quark, LAED3_COMPW )
FUNCTION_QUARK( CORE_dlaed3_reduceW_quark, LAED3_REDUCEW )
FUNCTION_QUARK( CORE_dlaed3_reduceW_p2_quark, LAED3_REDUCEW )
FUNCTION_QUARK( CORE_dlaed3_updatevectors_quark, LAED3_UPDATEVECTORS )
FUNCTION_QUARK( CORE_dlaed4_p2f1_quark, LAED4 )
FUNCTION_QUARK( CORE_dlag2z_quark, LAG2Z )
#endif
FUNCTION_QUARK( CORE_dzasum_quark, ASUM )
FUNCTION_QUARK( CORE_dzasum_f1_quark, ASUM )
FUNCTION_QUARK( CORE_zbrdalg1_quark, BRDALG )
FUNCTION_QUARK( CORE_zgeadd_quark, GEADD )
FUNCTION_QUARK( CORE_zgelqt_quark, GELQT )
FUNCTION_QUARK( CORE_zgemm_quark, GEMM )
FUNCTION_QUARK( CORE_zgemm_f2_quark, GEMM )
FUNCTION_QUARK( CORE_zgemm_p2_quark, GEMM )
FUNCTION_QUARK( CORE_zgemm_p3_quark, GEMM )
FUNCTION_QUARK( CORE_zgemm_p2f1_quark, GEMM )
FUNCTION_QUARK( CORE_zgemm_tile_quark, GEMM )
FUNCTION_QUARK( CORE_zgemv_quark, GEMV )
FUNCTION_QUARK( CORE_zgemv_tile_quark, GEMV )
FUNCTION_QUARK( CORE_zgeqp3_init_quark, GEQP3_INIT )
FUNCTION_QUARK( CORE_zgeqp3_larfg_quark, LARFG )
FUNCTION_QUARK( CORE_zgeqp3_norms_quark, GEQP3_NORMS )
FUNCTION_QUARK( CORE_zgeqp3_pivot_quark, GEQP3_PIVOT )
FUNCTION_QUARK( CORE_zgeqp3_tntpiv_quark, GEQRT )
FUNCTION_QUARK( CORE_zgeqp3_update_quark, GEQP3_UPDATE )
FUNCTION_QUARK( CORE_zgeqrt_quark, GEQRT )
FUNCTION_QUARK( CORE_zgessm_quark, GESSM )
FUNCTION_QUARK( CORE_zgessq_quark, LASSQ )
FUNCTION_QUARK( CORE_zgessq_f1_quark, LASSQ )
FUNCTION_QUARK( CORE_zgetrf_quark, GETRF )
FUNCTION_QUARK( CORE_zgetrf_incpiv_quark, GETRF )
FUNCTION_QUARK( CORE_zgetrf_nopiv_quark, GETRF )
FUNCTION_QUARK( CORE_zgetrf_reclap_quark, GETRF )
FUNCTION_QUARK( CORE_zgetrf_rectil_quark, GETRF )
FUNCTION_QUARK( CORE_zgetrip_quark, GETRIP )
FUNCTION_QUARK( CORE_zgetrip_f1_quark, GETRIP )
FUNCTION_QUARK( CORE_zgetrip_f2_quark, GETRIP )
FUNCTION_QUARK( CORE_zhegst_quark, HEGST )
#ifdef COMPLEX
FUNCTION_QUARK( CORE_zhemm_quark, HEMM )
#endif
#ifdef COMPLEX
FUNCTION_QUARK( CORE_zher2k_quark, HER2K )
#endif
FUNCTION_QUARK( CORE_zherfb_quark, HERFB )
#ifdef COMPLEX
FUNCTION_QUARK( CORE_zherk_quark, HERK )
#endif
#ifdef COMPLEX
FUNCTION_QUARK( CORE_zhessq_quark, LASSQ )
FUNCTION_QUARK( CORE_zhessq_f1_quark, LASSQ )
#endif
FUNCTION_QUARK( CORE_zlacpy_quark, LACPY )
FUNCTION_QUARK( CORE_zlacpy_f1_quark, LACPY )
FUNCTION_QUARK( CORE_zlacpy_pivot_quark, LACPY )
FUNCTION_QUARK( CORE_zlange_quark, LANGE )
FUNCTION_QUARK( CORE_zlange_f1_quark, LANGE )
#ifdef COMPLEX
FUNCTION_QUARK( CORE_zlanhe_quark, LANHE )
FUNCTION_QUARK( CORE_zlanhe_f1_quark, LANHE )
#endif
FUNCTION_QUARK( CORE_zlansy_quark, LANSY )
FUNCTION_QUARK( CORE_zlansy_f1_quark, LANSY )
FUNCTION_QUARK( CORE_zlantr_quark, LANGE )
FUNCTION_QUARK( CORE_zlantr_f1_quark, LANGE )
FUNCTION_QUARK( CORE_zlascl_quark, LASCL )
FUNCTION_QUARK( CORE_zlascl_p2f1_quark, LASCL )
FUNCTION_QUARK( CORE_zlaset2_quark, LASET )
FUNCTION_QUARK( CORE_zlaset_quark, LASET )
FUNCTION_QUARK( CORE_zlaset_identity_quark, LASET )
FUNCTION_QUARK( CORE_zlaswp_quark, LASWP )
FUNCTION_QUARK( CORE_zlaswp_f2_quark, LASWP )
FUNCTION_QUARK( CORE_zlaswp_ontile_quark, LASWP )
FUNCTION_QUARK( CORE_zlaswp_ontile_f2_quark, LASWP )
FUNCTION_QUARK( CORE_zswptr_ontile_quark, TRSM )
FUNCTION_QUARK( CORE_zlaswpc_ontile_quark, LASWP )
FUNCTION_QUARK( CORE_zlatro_quark, LATRO )
FUNCTION_QUARK( CORE_zlatro_f1_quark, LATRO )
FUNCTION_QUARK( CORE_zlauum_quark, LAUUM )
#ifdef COMPLEX
FUNCTION_QUARK( CORE_zplghe_quark, PLGHE )
#endif
FUNCTION_QUARK( CORE_zplgsy_quark, PLGSY )
FUNCTION_QUARK( CORE_zplrnt_quark, PLRNT )
FUNCTION_QUARK( CORE_zplssq_quark, LASSQ )
FUNCTION_QUARK( CORE_zpltmg_quark, PLRNT )
FUNCTION_QUARK( CORE_zpltmg_chebvand_quark, PLRNT )
FUNCTION_QUARK( CORE_zpltmg_circul_quark, PLRNT )
FUNCTION_QUARK( CORE_zpltmg_fiedler_quark, PLRNT )
FUNCTION_QUARK( CORE_zpltmg_hankel_quark, PLRNT )
FUNCTION_QUARK( CORE_zpltmg_toeppd1_quark, PLRNT )
FUNCTION_QUARK( CORE_zpltmg_toeppd2_quark, PLRNT )
FUNCTION_QUARK( CORE_zpotrf_quark, POTRF )
FUNCTION_QUARK( CORE_zsetvar_quark, SETVAR )
FUNCTION_QUARK( CORE_zshiftw_quark, SHIFTW )
FUNCTION_QUARK( CORE_zshift_quark, SHIFT )
FUNCTION_QUARK( CORE_zssssm_quark, SSSSM )
FUNCTION_QUARK( CORE_zstedc_quark, STEDC )
FUNCTION_QUARK( CORE_zstedc_f2_quark, STEDC )
FUNCTION_QUARK( CORE_zsteqr_quark, STEQR )
FUNCTION_QUARK( CORE_zswap_quark, SWAP )
FUNCTION_QUARK( CORE_zswpab_quark, SWPAB )
FUNCTION_QUARK( CORE_zsymm_quark, SYMM )
FUNCTION_QUARK( CORE_zsyr2k_quark, SYR2K )
FUNCTION_QUARK( CORE_zsyrk_quark, SYRK )
FUNCTION_QUARK( CORE_zsyssq_quark, LASSQ )
FUNCTION_QUARK( CORE_zsyssq_f1_quark, LASSQ )
FUNCTION_QUARK( CORE_ztradd_quark, GEADD )
FUNCTION_QUARK( CORE_ztrasm_quark, ASUM )
FUNCTION_QUARK( CORE_ztrasm_f1_quark, ASUM )
FUNCTION_QUARK( CORE_ztrdalg1_quark, TRDALG )
FUNCTION_QUARK( CORE_ztrmm_quark, TRMM )
FUNCTION_QUARK( CORE_ztrmm_p2_quark, TRMM )
FUNCTION_QUARK( CORE_ztrsm_quark, TRSM )
FUNCTION_QUARK( CORE_ztrssq_quark, LASSQ )
FUNCTION_QUARK( CORE_ztrssq_f1_quark, LASSQ )
FUNCTION_QUARK( CORE_ztrtri_quark, TRTRI )
FUNCTION_QUARK( CORE_ztslqt_quark, TSLQT )
FUNCTION_QUARK( CORE_ztsmlq_quark, TSMLQ )
FUNCTION_QUARK( CORE_ztsmlq_corner_quark, TSMLQ )
FUNCTION_QUARK( CORE_ztsmlq_hetra1_quark, TSMLQ )
FUNCTION_QUARK( CORE_ztsmqr_quark, TSMQR )
FUNCTION_QUARK( CORE_ztsmqr_corner_quark, TSMQR )
FUNCTION_QUARK( CORE_ztsmqr_hetra1_quark, TSMQR )
FUNCTION_QUARK( CORE_ztsqrt_quark, TSQRT )
FUNCTION_QUARK( CORE_ztstrf_quark, TSTRF )
FUNCTION_QUARK( CORE_zttlqt_quark, TTLQT )
FUNCTION_QUARK( CORE_zttmlq_quark, TTMLQ )
FUNCTION_QUARK( CORE_zttmqr_quark, TTMQR )
FUNCTION_QUARK( CORE_zttqrt_quark, TTQRT )
FUNCTION_QUARK( CORE_zunmlq_quark, UNMLQ )
FUNCTION_QUARK( CORE_zunmqr_quark, UNMQR )

