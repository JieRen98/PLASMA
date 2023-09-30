/**
 *
 * @file pstrmm.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:58 2016
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, float, m, n)
#define B(m,n) BLKADDR(B, float, m, n)

#if 0
/***************************************************************************//**
 *  Parallel tile triangular matrix-matrix multiplication - static scheduling
 **/
void plasma_pstrmm(plasma_context_t *plasma)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    PLASMA_enum trans;
    PLASMA_enum diag;
    float alpha;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_m;
    int next_n;
    int ldam, ldan, ldak, ldbm, ldbk;
    int tempkm, tempkn, tempmm, tempnn;

    float *lB;
    float zone  = (float)1.0;

    plasma_unpack_args_9(side, uplo, trans, diag, alpha, A, B, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= B.mt && n < B.nt) {
        n++;
        m = m-B.mt;
    }

    while (n < B.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= B.mt && next_n < B.nt) {
            next_n++;
            next_m = next_m - B.mt;
        }

        tempmm = m == B.mt-1 ? B.m - m * B.mb : B.mb;
        tempnn = n == B.nt-1 ? B.n - n * B.nb : B.nb;
        ldbm = BLKLDD(B, m);
        lB   = B(m, n);

        if ( side == PlasmaLeft ) {
            if ( uplo == PlasmaUpper ) {
                if ( trans == PlasmaNoTrans ) {
                    ldam = BLKLDD(A, m);

                    CORE_strmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(m, m), ldam,  /* lda * tempkm */
                               lB,      ldbk); /* ldb * tempnn */

                    for (k = m+1; k < A.mt; k++) {
                        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;
                        ldbk = BLKLDD(B, k);

                        CORE_sgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkn,
                            alpha, A(m, k), ldam,
                                   B(k, n), ldbk,
                            zone,  lB,      ldbm);
                    }
                }
                /*
                 *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
                 */
                else {
                    ldam = BLKLDD(A, m);

                    CORE_strmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(m, m), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = 0; k < m; k++) {
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        CORE_sgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            alpha, A(k, m), ldak,
                                   B(k, n), ldbk,
                            zone,  lB,      ldbm);
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaLower / PlasmaNoTrans
             */
            else {
                if ( trans == PlasmaNoTrans ) {
                    ldam = BLKLDD(A, m);

                    CORE_strmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(m, m), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = 0; k < m; k++) {
                        ldbk = BLKLDD(B, k);
                        CORE_sgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            alpha, A(m, k), ldam,
                                   B(k, n), ldbk,
                            zone,  lB,      ldbm);
                    }
                }
                /*
                 *  PlasmaLeft / PlasmaLower / Plasma[Conj]Trans
                 */
                else {
                    lda = BLKLDD(A, m);
                    CORE_strmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(m, m), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = m+1; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m - k * A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);

                        CORE_sgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            alpha, A(k, m), ldak,
                                   B(k, n), ldbk,
                            zone,  lB,      ldbm);
                    }
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
         */
        else {
            if (uplo == PlasmaUpper) {
                if ( trans == PlasmaNoTrans ) {
                    ldan = BLKLDD(A, n);

                    CORE_strmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(n, n), ldan,  /* lda * tempkm */
                               lB,      ldbm); /* ldb * tempnn */

                    for (k = 0; k < n; k++) {
                        ldak = BLKLDD(A, k);
                        CORE_sgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, B.mb,
                            alpha, B(m, k), ldbm,
                                   A(k, n), ldak,
                            zone,  lB,      ldbm);
                    }
                }
                /*
                 *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
                 */
                else {
                    ldan = BLKLDD(A, n);

                    CORE_strmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(n, n), ldan,  /* lda * tempkm */
                               lB,      ldbm); /* ldb * tempnn */

                    for (k = n+1; k < A.mt; k++) {
                        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;

                        CORE_sgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, tempkn,
                            alpha, B(m, k), ldbm,
                                   A(n, k), ldan,
                            zone,  lB,      ldbm);
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower / PlasmaNoTrans
             */
            else {
                if (trans == PlasmaNoTrans) {
                    ldan = BLKLDD(A, n);

                    CORE_strmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(n, n), ldan,  /* lda * tempkm */
                               lB,      ldbm); /* ldb * tempnn */

                    for (k = n+1; k < A.mt; k++) {
                        ldak = BLKLDD(A, k);
                        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;

                        CORE_sgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, tempkn,
                            alpha, B(m, k), ldbm,
                                   A(k, n), ldak,
                            zone,  lB,      ldbm);
                    }
                }
                /*
                 *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
                 */
                else {
                    ldan = BLKLDD(A, n);

                    CORE_strmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(n, n), ldan,  /* lda * tempkm */
                               lB,      ldbm); /* ldb * tempnn */

                    for (k = 0; k < n; k++) {
                        CORE_sgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, B.mb,
                            alpha, B(m, k), ldbm,
                                   A(n, k), ldan,
                            zone,  lB,      ldbm);
                    }
                }
            }
        }

        m = next_m;
        n = next_n;
    }
}
#endif

/***************************************************************************//**
 *  Parallel tile triangular matrix-matrix multiplication - dynamic scheduling
 **/
void plasma_pstrmm_quark(PLASMA_enum side, PLASMA_enum uplo,
                         PLASMA_enum trans, PLASMA_enum diag,
                         float alpha, PLASMA_desc A, PLASMA_desc B,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldam, ldan, ldbk, ldbm;
    int tempkm, tempkn, tempmm, tempnn;

    float zone = (float)1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, TRMM );

    /*
     *  PlasmaLeft / PlasmaUpper / PlasmaNoTrans
     */
    if (side == PlasmaLeft) {
        if (uplo == PlasmaUpper) {
            if (trans == PlasmaNoTrans) {
                for (m = 0; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    ldam = BLKLDD(A, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_strmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(m, m), ldam,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = m+1; k < A.mt; k++) {
                            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                            ldbk = BLKLDD(B, k);
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, A(m, k), ldam,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (m = B.mt-1; m > -1; m--) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    ldam = BLKLDD(A, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_strmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(m, m), ldam,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = 0; k < m; k++) {
                            ldbk = BLKLDD(B, k);
                            ldak = BLKLDD(A, k);
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaLeft / PlasmaLower / PlasmaNoTrans
         */
        else {
            if (trans == PlasmaNoTrans) {
                for (m = B.mt-1; m > -1; m--) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    ldam = BLKLDD(A, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_strmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(m, m), ldam,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = 0; k < m; k++) {
                            ldbk = BLKLDD(B, k);
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                alpha, A(m, k), ldam,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (m = 0; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    ldam = BLKLDD(A, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_strmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(m, m), ldam,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = m+1; k < A.mt; k++) {
                            tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                            ldak = BLKLDD(A, k);
                            ldbk = BLKLDD(B, k);
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm, A.mb,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
        }
    }
    /*
     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
     */
    else {
        if (uplo == PlasmaUpper) {
            if (trans == PlasmaNoTrans) {
                for (n = B.nt-1; n > -1; n--) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_strmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(n, n), ldan,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = 0; k < n; k++) {
                            ldak = BLKLDD(A, k);
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, B.mb, A.mb,
                                alpha, B(m, k), ldbm,
                                       A(k, n), ldak,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_strmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(n, n), ldan,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = n+1; k < A.mt; k++) {
                            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, B(m, k), ldbm,
                                       A(n, k), ldan,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaLower / PlasmaNoTrans
         */
        else {
            if (trans == PlasmaNoTrans) {
                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_strmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(n, n), ldan,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = n+1; k < A.mt; k++) {
                            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                            ldak = BLKLDD(A, k);
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, B(m, k), ldbm,
                                       A(k, n), ldak,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (n = B.nt-1; n > -1; n--) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_strmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(n, n), ldan,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = 0; k < n; k++) {
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, B.mb, A.mb,
                                alpha, B(m, k), ldbm,
                                       A(n, k), ldan,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
        }
    }
}
