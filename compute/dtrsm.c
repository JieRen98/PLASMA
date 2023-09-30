/**
 *
 * @file dtrsm.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:54 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dtrsm - Computes triangular solve A*X = B or X*A = B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  A*X = B
 *          = PlasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or ugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_dtrsm_Tile
 * @sa PLASMA_dtrsm_Tile_Async
 * @sa PLASMA_ctrsm
 * @sa PLASMA_dtrsm
 * @sa PLASMA_strsm
 *
 ******************************************************************************/
int PLASMA_dtrsm(PLASMA_enum side, PLASMA_enum uplo,
                 PLASMA_enum transA, PLASMA_enum diag,
                 int N, int NRHS, double alpha,
                 double *A, int LDA,
                 double *B, int LDB)
{
    int NB, NA;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtrsm", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    if (side == PlasmaLeft) {
      NA = N;
    } else {
      NA = NRHS;
    }

    /* Check input arguments */
    if (side != PlasmaLeft && side != PlasmaRight) {
        plasma_error("PLASMA_dtrsm", "illegal value of side");
        return -1;
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_dtrsm", "illegal value of uplo");
        return -2;
    }
    if (transA != PlasmaTrans && transA != PlasmaNoTrans && transA != PlasmaTrans ) {
        plasma_error("PLASMA_dtrsm", "illegal value of transA");
        return -3;
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        plasma_error("PLASMA_dtrsm", "illegal value of diag");
        return -4;
    }
    if (N < 0) {
        plasma_error("PLASMA_dtrsm", "illegal value of N");
        return -5;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_dtrsm", "illegal value of NRHS");
        return -6;
    }
    if (LDA < max(1, NA)) {
        plasma_error("PLASMA_dtrsm", "illegal value of LDA");
        return -8;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_dtrsm", "illegal value of LDB");
        return -10;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_DPOSV, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dtrsm", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, NA,   0, 0, NA, NA,   sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_dooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N,  NRHS, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, NA,   0, 0, NA, NA,  
                            sequence, &request);
        plasma_diplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N,  NRHS,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_dtrsm_Tile_Async(
        side, uplo, transA, diag, alpha, &descA, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, NA,    sequence, &request);
        plasma_diptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}
/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dtrsm_Tile - Computes triangular solve.
 *  Tile equivalent of PLASMA_dtrsm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  A*X = B
 *          = PlasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or ugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dtrsm
 * @sa PLASMA_dtrsm_Tile_Async
 * @sa PLASMA_ctrsm_Tile
 * @sa PLASMA_dtrsm_Tile
 * @sa PLASMA_strsm_Tile
 *
 ******************************************************************************/
int PLASMA_dtrsm_Tile(PLASMA_enum side, PLASMA_enum uplo,
                      PLASMA_enum transA, PLASMA_enum diag,
                      double alpha, PLASMA_desc *A, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtrsm_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dtrsm_Tile_Async(side, uplo, transA, diag, alpha, A, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dtrsm_Tile_Async - Computes triangular solve.
 *  Non-blocking equivalent of PLASMA_dtrsm_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_dtrsm
 * @sa PLASMA_dtrsm_Tile
 * @sa PLASMA_ctrsm_Tile_Async
 * @sa PLASMA_dtrsm_Tile_Async
 * @sa PLASMA_strsm_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dtrsm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo,
                            PLASMA_enum transA, PLASMA_enum diag,
                            double alpha, PLASMA_desc *A, PLASMA_desc *B,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtrsm_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dtrsm_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dtrsm_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dtrsm_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dtrsm_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_dtrsm_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (side != PlasmaLeft && side != PlasmaRight) {
        plasma_error("PLASMA_dtrsm_Tile", "illegal value of side");
        return plasma_request_fail(sequence, request, -1);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_dtrsm_Tile", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -2);
    }
    if (transA != PlasmaTrans && transA != PlasmaNoTrans && transA != PlasmaTrans) {
        plasma_error("PLASMA_dtrsm_Tile", "illegal value of transA");
        return plasma_request_fail(sequence, request, -3);
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        plasma_error("PLASMA_dtrsm_Tile", "illegal value of diag");
        return plasma_request_fail(sequence, request, -4);
    }

    /* Quick return */
    plasma_parallel_call_9(plasma_pdtrsm,
        PLASMA_enum, side,
        PLASMA_enum, uplo,
        PLASMA_enum, transA,
        PLASMA_enum, diag,
        double, alpha,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
