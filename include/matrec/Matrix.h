#ifndef MATREC_MATRIX_H
#define MATREC_MATRIX_H

#include "Shared.h"
#include <limits.h>

#ifdef __cplusplus
extern "C"{
#endif

///A standard Compressed Storage Row (CSR) matrix
typedef struct {
    MATREC_matrix_size numRows;
    MATREC_matrix_size numColumns;
    MATREC_matrix_size numNonzeros;
    MATREC_matrix_size *firstRowIndex;
    MATREC_matrix_size *entryColumns;
    int *entryValues;
} MATRECCSMatrixInt;

MATREC_ERROR MATRECcreateIntMatrix(MATREC *env,
                                   MATRECCSMatrixInt **mat,
                                   MATREC_matrix_size numRows,
                                   MATREC_matrix_size numColumns,
                                   MATREC_matrix_size numNonzeros
                           );

void MATRECfreeIntMatrix(MATREC *env,
                         MATRECCSMatrixInt **mat
                   );

typedef struct {
    MATREC_row row;
    MATREC_col column;
    int value;
} MATRECIntMatrixTriplet;

MATREC_ERROR MATRECcreateIntMatrixWithNonzeros(MATREC *env,
                                               MATRECCSMatrixInt **rowMat,
                                               MATREC_matrix_size numRows,
                                               MATREC_matrix_size numColumns,
                                               MATREC_matrix_size numNonzeros,
                                               const MATRECIntMatrixTriplet *nonzeros
);

void MATRECwriteIntMatrixToStream(const MATRECCSMatrixInt *matrix,
                                  FILE * file
                            );

MATREC_ERROR MATRECreadIntMatrixFromStream(MATREC* env,
                                           MATRECCSMatrixInt ** presult,
                                           FILE* stream
                                   );

MATREC_matrix_size MATRECintMatrixNumRows(MATRECCSMatrixInt * matrix);
MATREC_matrix_size MATRECintMatrixNumColumns(MATRECCSMatrixInt * matrix);
MATREC_matrix_size MATRECintMatrixNumNonzeros(MATRECCSMatrixInt * matrix);
MATREC_matrix_size MATRECintMatrixRowNumNonzeros(MATRECCSMatrixInt * matrix,
                                                 MATREC_row row
                                    );
MATREC_matrix_size * MATRECintMatrixRowColumnIndices(MATRECCSMatrixInt * matrix,
                                                     MATREC_row row
                                        );

int * MATRECintMatrixRowColumnValues(MATRECCSMatrixInt * matrix,
                                     MATREC_row row
                               );

bool MATRECintMatrixIsBinary(const MATRECCSMatrixInt * matrix);

bool MATRECintMatrixIsTernary(const MATRECCSMatrixInt * matrix);


MATREC_ERROR MATRECtransposeIntMatrix(MATREC *env,
                                      MATRECCSMatrixInt *in,
                                      MATRECCSMatrixInt **pOut
);


typedef struct {
    MATREC_matrix_size numRows;
    MATREC_matrix_size numColumns;
    MATREC_matrix_size numNonzeros;
    MATREC_matrix_size *firstRowIndex;
    MATREC_matrix_size *entryColumns;
    double *entryValues;
} MATRECCSMatrixDouble;

MATREC_ERROR MATRECcreateDoubleMatrix(MATREC *env,
                                      MATRECCSMatrixDouble **mat,
                                      MATREC_matrix_size numRows,
                                      MATREC_matrix_size numColumns,
                                      MATREC_matrix_size numNonzeros
                              );

void MATRECfreeDoubleMatrix(MATREC *env,
                            MATRECCSMatrixDouble **mat
                      );

typedef struct {
    MATREC_row row;
    MATREC_col column;
    double value;
} MATRECMatrixTripletDouble;

MATREC_ERROR MATRECcreateDoubleMatrixWithNonzeros(MATREC *env,
                                                  MATRECCSMatrixDouble **rowMat,
                                                  MATREC_matrix_size numRows,
                                                  MATREC_matrix_size numColumns,
                                                  MATREC_matrix_size numNonzeros,
                                                  const MATRECMatrixTripletDouble *nonzeros
                                      );

void MATRECwriteDoubleMatrixToStream(const MATRECCSMatrixDouble * matrix,
                                     FILE * file
                               );

MATREC_ERROR MATRECreadDoubleMatrixFromStream(MATREC* env,
                                              MATRECCSMatrixDouble ** presult,
                                              FILE* stream
                                      );

MATREC_matrix_size MATRECdoubleMatrixNumRows(MATRECCSMatrixDouble * matrix);
MATREC_matrix_size MATRECdoubleMatrixNumColumns(MATRECCSMatrixDouble * matrix);
MATREC_matrix_size MATRECdoubleMatrixNumNonzeros(MATRECCSMatrixDouble * matrix);
MATREC_matrix_size MATRECdoubleMatrixRowNumNonzeros(MATRECCSMatrixDouble * matrix,
                                                    MATREC_row row
);
MATREC_matrix_size * MATRECdoubleMatrixRowColumnIndices(MATRECCSMatrixDouble * matrix,
                                                        MATREC_row row
);

double * MATRECdoubleMatrixRowColumnValues(MATRECCSMatrixDouble * matrix,
                                           MATREC_row row
);

MATREC_ERROR MATRECtransposeDoubleMatrix(MATREC *env,
                                         MATRECCSMatrixDouble *in,
                                         MATRECCSMatrixDouble**pOut
                                 );

typedef struct {
    MATRECCSMatrixInt *rowMat;
    MATRECCSMatrixInt *colMat;
} MATRECCompressedSparseMatrixPairInt;

/**
 * Note this passes ownership of the given input matrix to the matrix pair!
 */
MATREC_ERROR MATRECcreateIntMatrixPair(MATREC *env,
                                       MATRECCSMatrixInt *rowMatrix,
                                       MATRECCompressedSparseMatrixPairInt **matrixPair
                               );

void MATRECfreeIntMatrixPair(MATREC *env, MATRECCompressedSparseMatrixPairInt **matrixPair);


typedef struct {
    MATRECCSMatrixDouble *rowMat;
    MATRECCSMatrixDouble *colMat;
} MATRECCompressedSparseMatrixPairDouble;

/**
 * Note this passes ownership of the given input matrix to the matrix pair!
 */
MATREC_ERROR MATRECcreateMatrixPairFromDoubleRowMatrix(MATREC *env,
                                                       MATRECCSMatrixDouble *rowMatrix,
                                                       MATRECCompressedSparseMatrixPairDouble **matrixPair
                                               );

void MATRECfreeMatrixPairDouble(MATREC *env,
                                MATRECCompressedSparseMatrixPairDouble **matrixPair
                          );

///Represents a submatrix by keeping two arrays with the rows and columns in the submatrix
typedef struct{
    MATREC_matrix_size numRows;     /**< \brief Number of rows. */
    MATREC_row *rows;       /**< \brief Array with row indices. */
    MATREC_matrix_size numColumns;  /**< \brief Number of columns. */
    MATREC_col *columns;    /**< \brief Array with column indices. */
} MATRECSubMatrix;

/**
 * \brief Creates a submatrix of given size.
 *
 * Only allocates the memory. Use rows and columns attributes of *\p psubmatrix to actually set the row and column
 * indices, respectively.
 */
MATREC_ERROR MATRECcreateSubMatrix(
        MATREC* env,                   /**< MATREC environment. */
        MATREC_matrix_size numRows,         /**< Number of rows */
        MATREC_matrix_size numColumns,      /**< Number of columns */
        MATRECSubMatrix ** psubmatrix      /**< Pointer to where the submatrix is to be stored. */
);

/**
 * \brief Frees a submatrix.
 */
void MATRECfreeSubMatrix(
        MATREC *env,               /**< \ref environment. */
        MATRECSubMatrix **psubmatrix /**< Pointer to submatrix. */
);

/**
 * \brief Writes the submatrix \p submatrix to the file \p stream by means of lists of row and column indices.
 */

void MATRECwriteSubMatrixToStream(
        MATRECSubMatrix *submatrix,  /**< Reference submatrix. */
        MATREC_matrix_size numRows,         /**< Number of rows of original matrix. */
        MATREC_matrix_size numColumns,      /**< Number of columns of original matrix. */
        FILE *stream            /**< File stream to save submatrix to. */
);

/**
 * \brief Reads the submatrix \p *psubmatrix from the file \p stream.
 */
MATREC_ERROR MATRECreadSubMatrixFromStream(
        MATREC *env,                   /**< \ref environment. */
        MATRECSubMatrix **psubmatrix,    /**< Pointer for storing the submatrix. */
        MATREC_matrix_size *pnumMatrixRows,     /**< Pointer for storing the number of rows of the original matrix; may be \c NULL. */
        MATREC_matrix_size *pnumMatrixColumns,  /**< Pointer for storing the number of rows of the original matrix; may be \c NULL. */
        FILE *stream                /**< File stream to save submatrix to. */
);

MATREC_matrix_size MATRECcountIntSubMatrixNonzeros(const MATRECSubMatrix* subMatrix,
                                                   const MATRECCSMatrixInt *rowMatrix
);

MATREC_matrix_size MATRECcountDoubleSubMatrixNonzeros(const MATRECSubMatrix* subMatrix,
                                                      const MATRECCSMatrixDouble * rowMatrix
                                         );

void MATRECtransposeSubmatrix(MATRECSubMatrix * submatrix);

#ifdef __cplusplus
}
#endif
#endif //MATREC_MATRIX_H
