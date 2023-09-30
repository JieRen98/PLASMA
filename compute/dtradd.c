/**
 *
 * @file dtradd.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:54 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dtradd - Performs a matrix addition similarly to the pdtradd()
 *  function from the PBLAS library:
 *
 *    \f[ C = \alpha op( A ) + \beta B \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = g( X' )
 *
 *  alpha and beta are scalars, and A, and B are two trapezoidal matrices, with
 *  op( A ) and B two m by n matrices.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A and B matrices:
 *          = PlasmaUpperLower: A and B are general matrices.
 *          = PlasmaUpper: op(A) and B are upper trapezoidal matrices.
 *          = PlasmaLower: op(A) and B are lower trapezoidal matrices.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          ugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix B. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( A ) and of the matrix B. N >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is N when trans = PlasmaNoTrans,
 *          and is M otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,K), where K is M
 *          when trans = PlasmaNoTrans, and is N when otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] B
 *          B is a LDB-by-N matrix.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dtradd_Tile
 * @sa PLASMA_ctradd
 * @sa PLASMA_dtradd
 * @sa PLASMA_stradd
 *
 ******************************************************************************/
int PLASMA_dtradd(PLASMA_enum uplo, PLASMA_enum trans, int M, int N,
                  double alpha, double *A, int LDA,
                  double beta,  double *B, int LDB)
{
    int NB;
    int Am, An;
    int status;
    PLASMA_desc descA, descB;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtradd", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != PlasmaUpperLower) && (uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        plasma_error("PLASMA_dtradd", "illegal value of uplo");
        return -1;
    }
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans) && (trans != PlasmaTrans)) {
        plasma_error("PLASMA_dtradd", "illegal value of trans");
        return -2;
    }
    if ( trans == PlasmaNoTrans ) {
        Am = M; An = N;
    } else {
        Am = N; An = M;
    }
    if (M < 0) {
        plasma_error("PLASMA_dtradd", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_dtradd", "illegal value of N");
        return -4;
    }
    if (LDA < max(1, Am)) {
        plasma_error("PLASMA_dtradd", "illegal value of LDA");
        return -7;
    }
    if (LDB < max(1, M)) {
        plasma_error("PLASMA_dtradd", "illegal value of LDB");
        return -10;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (double)0.0) && beta == (double)1.0))
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_DGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dtradd", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_dooplap2tile( descB, B, NB, NB, LDB, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An,
                            sequence, &request);
        plasma_diplap2tile( descB, B, NB, NB, LDB, N, 0, 0, M, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_dtradd_Tile_Async(
        uplo, trans, alpha, &descA, beta, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descB, B, NB, NB, LDB, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, An, sequence, &request);
        plasma_diptile2lap( descB, B, NB, NB, LDB, N,  sequence, &request);
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
 *  PLASMA_dtradd_Tile - Performs a matrix addition similarly to the pdtradd()
 *  function from the PBLAS library.
 *  Tile equivalent of PLASMA_dtradd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A and B matrices:
 *          = PlasmaUpperLower: A and B are general matrices.
 *          = PlasmaUpper: op(A) and B are upper trapezoidal matrices.
 *          = PlasmaLower: op(A) and B are lower trapezoidal matrices.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          ugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is N when trans = PlasmaNoTrans,
 *          and is M otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] B
 *          B is a LDB-by-N matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dtradd
 * @sa PLASMA_dtradd_Tile_Async
 * @sa PLASMA_ctradd_Tile
 * @sa PLASMA_dtradd_Tile
 * @sa PLASMA_stradd_Tile
 *
 ******************************************************************************/
int PLASMA_dtradd_Tile(PLASMA_enum uplo, PLASMA_enum trans,
                       double alpha, PLASMA_desc *A,
                       double beta,  PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtradd_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dtradd_Tile_Async(uplo, trans, alpha, A, beta, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dtradd_Tile_Async - Performs a matrix addition similarly to the
 *  pdtradd() function from the PBLAS library.
 *  Non-blocking equivalent of PLASMA_dtradd_Tile().
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
 * @sa PLASMA_dtradd
 * @sa PLASMA_dtradd_Tile
 * @sa PLASMA_ctradd_Tile_Async
 * @sa PLASMA_dtradd_Tile_Async
 * @sa PLASMA_stradd_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dtradd_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans,
                             double alpha, PLASMA_desc *A,
                             double beta,  PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descB;
    int M, N;
    int Am, An, Ai, Aj, Amb, Anb;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtradd_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dtradd_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dtradd_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dtradd_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dtradd_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans) && (trans != PlasmaTrans)) {
        plasma_error("PLASMA_dtradd_Tile_Async", "illegal value of trans");
        return plasma_request_fail(sequence, request, -1);
    }

    if ( trans == PlasmaNoTrans ) {
        Am  = descA.m;
        An  = descA.n;
        Amb = descA.mb;
        Anb = descA.nb;
        Ai  = descA.i;
        Aj  = descA.j;
    } else {
        Am  = descA.n;
        An  = descA.m;
        Amb = descA.nb;
        Anb = descA.mb;
        Ai  = descA.j;
        Aj  = descA.i;
    }

    if ( (Amb != descB.mb) || (Anb != descB.nb) ) {
        plasma_error("PLASMA_dtradd_Tile_Async", "tile sizes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (Am != descB.m) || (An != descB.n) ) {
        plasma_error("PLASMA_dtradd_Tile_Async", "sizes of matrices have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (Ai%Amb != descB.i%descB.mb) ||
         (Aj%Anb != descB.j%descB.nb) )
    {
        plasma_error("PLASMA_dtradd_Tile_Async", "start indexes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    M = descB.m;
    N = descB.n;

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (double)0.0) && beta == (double)1.0))
        return PLASMA_SUCCESS;

    plasma_parallel_call_8(plasma_pdtradd,
        PLASMA_enum, uplo,
        PLASMA_enum, trans,
        double, alpha,
        PLASMA_desc, descA,
        double, beta,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
