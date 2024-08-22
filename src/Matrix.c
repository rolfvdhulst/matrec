#include "matrec/Matrix.h"

bool MATRECrowIsInvalid(MATREC_row row){
    return row == MATREC_INVALID_ROW;
}
bool MATRECrowIsValid(MATREC_row row){
    return !MATRECrowIsInvalid(row);
}
bool MATRECcolIsInvalid(MATREC_col col){
    return col == MATREC_INVALID_COL;
}
bool MATRECcolIsValid(MATREC_col col){
    return !MATRECcolIsInvalid(col);
}

MATREC_ERROR MATRECcreateSubMatrix(
        MATREC* env,
        MATREC_matrix_size numRows,
        MATREC_matrix_size numColumns,
        MATRECSubMatrix ** psubmatrix
){
    assert(psubmatrix);

    MATREC_CALL( MATRECallocBlock(env, psubmatrix) );
    MATRECSubMatrix * submatrix = *psubmatrix;

    submatrix->numRows = numRows;
    submatrix->numColumns = numColumns;
    submatrix->rows = NULL;
    submatrix->columns = NULL;
    MATREC_CALL( MATRECallocBlockArray(env, &submatrix->rows, numRows) );
    MATREC_CALL( MATRECallocBlockArray(env, &submatrix->columns, numColumns) );


    return MATREC_OKAY;
}

void MATRECfreeSubMatrix(
        MATREC *env,
        MATRECSubMatrix **psubmatrix
){
    assert(psubmatrix);
    MATRECSubMatrix * submatrix = *psubmatrix;
    if (!submatrix){
        return;
    }

    if (submatrix->rows)
        MATRECfreeBlockArray(env, &submatrix->rows);
    if (submatrix->columns)
        MATRECfreeBlockArray(env, &submatrix->columns);
    MATRECfreeBlock(env, psubmatrix);
}

void MATRECwriteSubMatrixToStream(
        MATRECSubMatrix *submatrix,
        MATREC_matrix_size numRows,
        MATREC_matrix_size numColumns,
        FILE *stream
){
    assert(submatrix);
    assert(stream);

    fprintf(stream, "%lu %lu %lu %lu\n", numRows, numColumns, submatrix->numRows, submatrix->numColumns);
    for (MATREC_matrix_size row = 0; row < submatrix->numRows; ++row)
        fprintf(stream, "%lu ", submatrix->rows[row] + 1);
    fputc('\n', stream);
    for (MATREC_matrix_size column = 0; column < submatrix->numColumns; ++column)
        fprintf(stream, "%lu ", submatrix->columns[column] + 1);
    fputc('\n', stream);
}

MATREC_ERROR MATRECreadSubMatrixFromStream(
        MATREC *env,
        MATRECSubMatrix **psubmatrix,
        MATREC_matrix_size *pnumMatrixRows,
        MATREC_matrix_size *pnumMatrixColumns,
        FILE *stream
){
    assert(env);
    assert(psubmatrix);
    assert(stream);

    MATREC_matrix_size numOriginalRows;
    MATREC_matrix_size numOriginalColumns;
    MATREC_matrix_size numRows;
    MATREC_matrix_size numColumns;
    if (fscanf(stream, "%lu %lu %lu %lu", &numOriginalRows, &numOriginalColumns, &numRows, &numColumns) != 4)
        return MATREC_ERROR_INPUT;

    if (numRows > numOriginalRows || numColumns > numOriginalColumns)
        return MATREC_ERROR_INPUT;

    MATREC_CALL(MATRECcreateSubMatrix(env, numRows, numColumns, psubmatrix) );
    MATRECSubMatrix * submatrix = *psubmatrix;
    for (MATREC_matrix_size r = 0; r < numRows; ++r)
    {
        MATREC_matrix_size row;
        int numRead = fscanf(stream, "%lu", &row);
        if (numRead != 1)
            return MATREC_ERROR_INPUT;
        submatrix->rows[r] = row - 1;
    }
    for (MATREC_matrix_size c = 0; c < numColumns; ++c)
    {
        MATREC_matrix_size col;
        int numRead = fscanf(stream, "%lu", &col);
        if (numRead != 1)
            return MATREC_ERROR_INPUT;
        submatrix->columns[c] = col - 1;
    }

    if (pnumMatrixRows)
        *pnumMatrixRows = numOriginalRows;
    if (pnumMatrixColumns)
        *pnumMatrixColumns = numOriginalColumns;

    return MATREC_OKAY;
}

MATREC_ERROR MATRECcreateIntMatrix(MATREC * env, MATRECCSMatrixInt ** pmat,
                                   MATREC_matrix_size numRows, MATREC_matrix_size numColumns, MATREC_matrix_size numNonzeros){
    assert(env);
    assert(pmat);

    MATREC_CALL(MATRECallocBlock(env, pmat));
    MATRECCSMatrixInt * mat =  *pmat;
    mat->numRows = numRows;
    mat->numColumns = numColumns;
    mat->numNonzeros = numNonzeros;
    MATREC_CALL(MATRECallocBlockArray(env, &mat->firstRowIndex, numRows + 1));
    if(numNonzeros > 0){
        MATREC_CALL(MATRECallocBlockArray(env, &mat->entryColumns, numNonzeros));
        MATREC_CALL(MATRECallocBlockArray(env,&mat->entryValues,numNonzeros));
    }else{
        mat->entryColumns = NULL;
        mat->entryValues = NULL;
    }
    return MATREC_OKAY;
}

MATREC_ERROR MATRECcreateIntMatrixWithNonzeros(MATREC * env, MATRECCSMatrixInt **rowMat,
                                               MATREC_matrix_size numRows, MATREC_matrix_size numColumns, MATREC_matrix_size numNonzeros, const MATRECIntMatrixTriplet * nonzeros){
    MATREC_CALL(MATRECcreateIntMatrix(env, rowMat, numRows, numColumns, numNonzeros));
#ifndef NDEBUG
    //assert that the nonzero array is sorted
    if(numNonzeros > 0){
        for (MATREC_matrix_size i = 0; i < numNonzeros - 1; ++i) {
            assert( (nonzeros[i+1].row > nonzeros[i].row) ||
                    (nonzeros[i+1].row == nonzeros[i].row && nonzeros[i+1].column > nonzeros[i].column));
        }
    }
#endif
    MATREC_row currentRowIdx = 0;
    MATRECCSMatrixInt * matrix = *rowMat;

    MATREC_matrix_size index = 0;
    while(currentRowIdx < numRows){
        matrix->firstRowIndex[currentRowIdx] = index;
        while( index < numNonzeros && nonzeros[index].row == currentRowIdx){
            assert(matrix->firstRowIndex[currentRowIdx] == index ||
                   (index == 0 || nonzeros[index].column > nonzeros[index-1].column));//Each row must have increasing column indices
            matrix->entryColumns[index] = nonzeros[index].column;
            matrix->entryValues[index] = nonzeros[index].value;
            index++;
        }
        currentRowIdx++;
    }
    assert(index == numNonzeros);
    matrix->firstRowIndex[numRows] = index;

    return MATREC_OKAY;
}

void MATRECfreeIntMatrix(MATREC * env, MATRECCSMatrixInt ** mat){
    MATRECCSMatrixInt * matrix = *mat;
    MATRECfreeBlockArray(env,&matrix->entryValues);
    MATRECfreeBlockArray(env,&matrix->entryColumns);
    MATRECfreeBlockArray(env,&matrix->firstRowIndex);
    MATRECfreeBlock(env,mat);
}

MATREC_ERROR MATRECtransposeIntMatrix(MATREC * env, MATRECCSMatrixInt * in, MATRECCSMatrixInt ** pOut){
    assert(env);
    assert(in);
    assert(pOut);
    MATREC_CALL(MATRECcreateIntMatrix(env, pOut, in->numColumns, in->numRows, in->numNonzeros));
    MATRECCSMatrixInt * outMat = *pOut;

    //Read documentation for this function as if we are going from row -> column

    //count number of nonzeros in each column, storing it in the next entry
    for (MATREC_matrix_size i = 0; i <= in->numColumns; ++i) {
        outMat->firstRowIndex[i] = 0;
    }
    for (MATREC_matrix_size i = 0; i < in->numNonzeros; ++i) {
        ++(outMat->firstRowIndex[in->entryColumns[i] + 1]);
    }
    for (MATREC_matrix_size i = 1; i < in->numColumns; ++i) {
        outMat->firstRowIndex[i] += outMat->firstRowIndex[i - 1];
    }
    //Start indices for each column
    for (MATREC_matrix_size row = 0; row < in->numRows; ++row) {
        MATREC_matrix_size first = in->firstRowIndex[row];
        MATREC_matrix_size beyond = in->firstRowIndex[row + 1];
        for (MATREC_matrix_size entry = first; entry < beyond; ++entry) {
            MATREC_matrix_size column = in->entryColumns[entry];
            MATREC_matrix_size transIndex = outMat->firstRowIndex[column];
            outMat->entryColumns[transIndex] = row;
            outMat->entryValues[transIndex] = in->entryValues[entry];
            ++(outMat->firstRowIndex[column]);
        }
    }
    /* We shifted rowSlice of result, so we shift it back. */
    for (MATREC_matrix_size c = outMat->numRows; c > 0; --c){
        outMat->firstRowIndex[c] = outMat->firstRowIndex[c - 1];
    }
    outMat->firstRowIndex[0] = 0;

    return MATREC_OKAY;
}
MATREC_ERROR MATRECcreateIntMatrixPair(MATREC * env, MATRECCSMatrixInt * rowMatrix, MATRECCompressedSparseMatrixPairInt ** matrixPair){
    assert(env);
    assert(rowMatrix);
    assert(matrixPair);
    MATREC_CALL(MATRECallocBlock(env,matrixPair));
    MATRECCompressedSparseMatrixPairInt * pair = *matrixPair;
    pair->rowMat = rowMatrix;
    pair->colMat = NULL;
    MATREC_CALL(MATRECtransposeIntMatrix(env, rowMatrix, &pair->colMat));
    return MATREC_OKAY;
}

void MATRECfreeIntMatrixPair(MATREC * env, MATRECCompressedSparseMatrixPairInt ** matrixPair){
    MATRECCompressedSparseMatrixPairInt  * pair = *matrixPair;
    MATRECfreeIntMatrix(env, &pair->rowMat);
    MATRECfreeIntMatrix(env, &pair->colMat);
    MATRECfreeBlock(env,matrixPair);
}

bool MATRECintMatrixIsTernary(const MATRECCSMatrixInt * matrix){
    for(MATREC_matrix_size i = 0; i < matrix->numNonzeros; i++){
        int value = matrix->entryValues[i];
        if(!(value == 1 || value == -1)){
            return false;
        }
    }
    return true;
}

bool MATRECintMatrixIsBinary(const MATRECCSMatrixInt * matrix){
    for(MATREC_matrix_size i = 0; i < matrix->numNonzeros; i++){
        int value = matrix->entryValues[i];
        if(value != 1){
            return false;
        }
    }
    return true;
}

void MATRECwriteIntMatrixToStream(const MATRECCSMatrixInt *matrix, FILE * file){
    assert(matrix);
    assert(file);

    fprintf(file, "%lu %lu %lu\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
    for (MATREC_matrix_size row = 0; row < matrix->numRows; ++row)
    {
        MATREC_matrix_size first = matrix->firstRowIndex[row];
        MATREC_matrix_size beyond = matrix->firstRowIndex[row + 1];
        for (MATREC_matrix_size entry = first; entry < beyond; ++entry){
            fprintf(file, "%lu %lu %d\n", row+1, matrix->entryColumns[entry] + 1, matrix->entryValues[entry]);
        }
    }
}
static
int compareIntNonzeros(const void* pa, const void* pb)
{
    MATREC_matrix_size  aRow = ((const MATRECIntMatrixTriplet *)pa)->row;
    MATREC_matrix_size  bRow = ((const MATRECIntMatrixTriplet *)pb)->row;
    if (aRow < bRow)
        return -1;
    else if (aRow > bRow)
        return +1;

    MATREC_matrix_size  aColumn = ((const MATRECIntMatrixTriplet*)pa)->column;
    MATREC_matrix_size  bColumn = ((const MATRECIntMatrixTriplet*)pb)->column;
    if(aColumn < bColumn){
        return -1;
    }else if (aColumn == bColumn){
        return 0;
    }
    return 1;
}

MATREC_ERROR MATRECreadIntMatrixFromStream(MATREC* env, MATRECCSMatrixInt ** presult, FILE* stream)
{
    assert(env);
    assert(presult);
    assert(stream);

    MATREC_matrix_size numRows, numColumns, numNonzeros;
    int numRead = fscanf(stream, "%lu %lu %lu", &numRows, &numColumns, &numNonzeros);
    if (numRead < 3)
    {
        return MATREC_ERROR_INPUT;
    }

    /* Read all nonzeros. */

    MATRECIntMatrixTriplet * nonzeros = NULL;


    MATREC_CALL( MATRECallocBlockArray(env, &nonzeros, numNonzeros) );
    MATREC_matrix_size entry = 0;
    for (MATREC_matrix_size i = 0; i < numNonzeros; ++i)
    {
        MATREC_matrix_size row;
        MATREC_matrix_size column;
        int value;
        numRead = fscanf(stream, "%lu %lu %d", &row, &column, &value);
        if (numRead < 3 || row == 0 || column == 0 || row > numRows || column > numColumns)
        {
            MATRECfreeBlockArray(env,&numNonzeros);
            return MATREC_ERROR_INPUT;
        }
        if (value != 0)
        {
            nonzeros[entry].row = row - 1;
            nonzeros[entry].column = column - 1;
            nonzeros[entry].value = value;
            ++entry;
        }
    }
    numNonzeros = entry;

    /* We sort all nonzeros by row and then by column. */
    qsort(nonzeros, numNonzeros, sizeof(MATRECIntMatrixTriplet), compareIntNonzeros);

    MATREC_CALL(MATRECcreateIntMatrixWithNonzeros(env, presult, numRows, numColumns, numNonzeros, nonzeros));

    MATRECfreeBlockArray(env,&nonzeros);
    return MATREC_OKAY;
}

MATREC_matrix_size MATRECintMatrixNumRows(MATRECCSMatrixInt * matrix){
    assert(matrix);
    return matrix->numRows;
}
MATREC_matrix_size MATRECintMatrixNumColumns(MATRECCSMatrixInt * matrix){
    return matrix->numColumns;
}
MATREC_matrix_size MATRECintMatrixNumNonzeros(MATRECCSMatrixInt * matrix){
    return matrix->numNonzeros;
}
MATREC_matrix_size MATRECintMatrixRowNumNonzeros(MATRECCSMatrixInt * matrix,
                                                 MATREC_row row
){
    assert(row < matrix->numRows);
    return matrix->firstRowIndex[row+1] - matrix->firstRowIndex[row];
}

MATREC_matrix_size * MATRECintMatrixRowColumnIndices(MATRECCSMatrixInt * matrix,
                                                     MATREC_row row
){

    return &matrix->entryColumns[matrix->firstRowIndex[row]];
}

int * MATRECintMatrixRowColumnValues(MATRECCSMatrixInt * matrix,
                                     MATREC_row row
){
    return &matrix->entryValues[matrix->firstRowIndex[row]];
}

MATREC_matrix_size MATRECcountIntSubMatrixNonzeros(const MATRECSubMatrix* subMatrix, const MATRECCSMatrixInt *rowMatrix){
    MATREC_matrix_size nonZeros = 0;
    for (MATREC_matrix_size i = 0; i < subMatrix->numRows; ++i) {
        MATREC_row row = (subMatrix->rows[i] - 1); //submatrix rows are stored as 1 higher!!
        MATREC_matrix_size matrixIndex = rowMatrix->firstRowIndex[row];
        MATREC_matrix_size matrixRowEnd = rowMatrix->firstRowIndex[row + 1];
        assert(matrixRowEnd <= rowMatrix->numNonzeros);
        for(MATREC_matrix_size subMatColIndex = 0; subMatColIndex < subMatrix->numColumns; subMatColIndex++){
            MATREC_col subMatCol = subMatrix->columns[subMatColIndex] - 1; //submatrix columns are stored as 1 higher!!
            while(matrixIndex != matrixRowEnd &&
                  rowMatrix->entryColumns[matrixIndex] < subMatCol){
                matrixIndex++;
            }
            if(matrixIndex == matrixRowEnd){
                break;
            }
            if(rowMatrix->entryColumns[matrixIndex] == subMatCol){
                nonZeros++;
            }
        }
    }
    return nonZeros;
}

MATREC_matrix_size MATRECcountDoubleSubMatrixNonzeros(const MATRECSubMatrix* subMatrix, const MATRECCSMatrixDouble * rowMatrix){
    MATREC_matrix_size nonZeros = 0;
    for (MATREC_matrix_size i = 0; i < subMatrix->numRows; ++i) {
        MATREC_row row = (subMatrix->rows[i] - 1); //submatrix rows are stored as 1 higher!!
        MATREC_matrix_size matrixIndex = rowMatrix->firstRowIndex[row];
        MATREC_matrix_size matrixRowEnd = rowMatrix->firstRowIndex[row + 1];
        assert(matrixRowEnd <= rowMatrix->numNonzeros);
        for(MATREC_matrix_size subMatColIndex = 0; subMatColIndex < subMatrix->numColumns; subMatColIndex++){
            MATREC_col subMatCol = subMatrix->columns[subMatColIndex] - 1; //submatrix columns are stored as 1 higher!!
            while(matrixIndex != matrixRowEnd &&
                  rowMatrix->entryColumns[matrixIndex] < subMatCol){
                matrixIndex++;
            }
            if(matrixIndex == matrixRowEnd){
                break;
            }
            if(rowMatrix->entryColumns[matrixIndex] == subMatCol){
                nonZeros++;
            }
        }
    }
    return nonZeros;
}

MATREC_ERROR MATRECcreateDoubleMatrix(MATREC *env, MATRECCSMatrixDouble **pmat,
                                      MATREC_matrix_size numRows, MATREC_matrix_size numColumns, MATREC_matrix_size numNonzeros){
    assert(env);
    assert(pmat);
    MATREC_CALL(MATRECallocBlock(env, pmat));
    MATRECCSMatrixDouble * mat = *pmat;
    mat->numRows = numRows;
    mat->numColumns = numColumns;
    mat->numNonzeros = numNonzeros;
    MATREC_CALL(MATRECallocBlockArray(env, &mat->firstRowIndex, numRows + 1));
    if(numNonzeros > 0){
        MATREC_CALL(MATRECallocBlockArray(env,& mat->entryColumns, numNonzeros));
        MATREC_CALL(MATRECallocBlockArray(env,& mat->entryValues,numNonzeros));
    }else{
        mat->entryColumns = NULL;
        mat->entryValues = NULL;
    }
    return MATREC_OKAY;
}


MATREC_ERROR MATRECcreateDoubleMatrixWithNonzeros(MATREC *env, MATRECCSMatrixDouble **rowMat,
                                                  MATREC_matrix_size numRows, MATREC_matrix_size numColumns, MATREC_matrix_size numNonzeros,
                                                  const MATRECMatrixTripletDouble *nonzeros){
    MATREC_CALL(MATRECcreateDoubleMatrix(env, rowMat, numRows, numColumns, numNonzeros));
#ifndef NDEBUG
    //assert that the nonzero array is sorted
    if(numNonzeros > 0) {
        for (MATREC_matrix_size i = 0; i < numNonzeros - 1; ++i) {
            assert((nonzeros[i + 1].row > nonzeros[i].row) ||
                   (nonzeros[i + 1].row == nonzeros[i].row && nonzeros[i + 1].column > nonzeros[i].column));
        }
    }
#endif
    MATREC_row currentRowIdx = 0;
    MATRECCSMatrixDouble * matrix = *rowMat;

    MATREC_matrix_size index = 0;
    while(currentRowIdx < numRows){
        matrix->firstRowIndex[currentRowIdx] = index;
        while( index < numNonzeros && nonzeros[index].row == currentRowIdx){
            assert(matrix->firstRowIndex[currentRowIdx] == index ||
                   (index == 0 || nonzeros[index].column > nonzeros[index-1].column));//Each row must have increasing column indices
            matrix->entryColumns[index] = nonzeros[index].column;
            matrix->entryValues[index] = nonzeros[index].value;
            index++;
        }
        currentRowIdx++;
    }
    assert(index == numNonzeros);
    matrix->firstRowIndex[numRows] = index;

    return MATREC_OKAY;
}
void MATRECfreeDoubleMatrix(MATREC *env, MATRECCSMatrixDouble **mat){
    MATRECCSMatrixDouble * matrix = *mat;
    MATRECfreeBlockArray(env,&matrix->entryValues);
    MATRECfreeBlockArray(env,&matrix->entryColumns);
    MATRECfreeBlockArray(env,&matrix->firstRowIndex);
    MATRECfreeBlock(env,mat);
}

MATREC_ERROR MATRECtransposeDoubleMatrix(MATREC *env, MATRECCSMatrixDouble *in, MATRECCSMatrixDouble**pOut){
    assert(env);
    assert(in);
    assert(pOut);
    assert(*pOut == NULL);
    MATREC_CALL(MATRECcreateDoubleMatrix(env, pOut, in->numColumns, in->numRows, in->numNonzeros));
    MATRECCSMatrixDouble * outMat = *pOut;

    //Read documentation for this function as if we are going from row -> column

    //count number of nonzeros in each column, storing it in the next entry
    for (MATREC_matrix_size i = 0; i <= in->numColumns; ++i) {
        outMat->firstRowIndex[i] = 0;
    }
    for (MATREC_matrix_size i = 0; i < in->numNonzeros; ++i) {
        ++(outMat->firstRowIndex[in->entryColumns[i] + 1]);
    }
    for (MATREC_matrix_size i = 1; i < in->numColumns; ++i) {
        outMat->firstRowIndex[i] += outMat->firstRowIndex[i - 1];
    }
    //Start indices for each column
    for (MATREC_matrix_size row = 0; row < in->numRows; ++row) {
        MATREC_matrix_size first = in->firstRowIndex[row];
        MATREC_matrix_size beyond = in->firstRowIndex[row + 1];
        for (MATREC_matrix_size entry = first; entry < beyond; ++entry) {
            MATREC_matrix_size column = in->entryColumns[entry];
            MATREC_matrix_size transIndex = outMat->firstRowIndex[column];
            outMat->entryColumns[transIndex] = row;
            outMat->entryValues[transIndex] = in->entryValues[entry];
            ++(outMat->firstRowIndex[column]);
        }
    }
    /* We shifted rowSlice of result, so we shift it back. */
    for (MATREC_matrix_size c = outMat->numRows; c > 0; --c){
        outMat->firstRowIndex[c] = outMat->firstRowIndex[c - 1];
    }
    outMat->firstRowIndex[0] = 0;
    return MATREC_OKAY;
}

/**
 * Note this passes ownership of the given input matrix to the matrix pair!
 */
MATREC_ERROR MATRECcreateMatrixPairFromDoubleRowMatrix(MATREC *env, MATRECCSMatrixDouble *rowMatrix,
                                                       MATRECCompressedSparseMatrixPairDouble **matrixPair){
    assert(env);
    assert(rowMatrix);
    assert(matrixPair);
    assert(*matrixPair == NULL);
    MATREC_CALL(MATRECallocBlock(env,matrixPair));
    MATRECCompressedSparseMatrixPairDouble * pair = *matrixPair;
    pair->rowMat = rowMatrix;
    pair->colMat = NULL;
    MATREC_CALL(MATRECtransposeDoubleMatrix(env, rowMatrix, &pair->colMat));
    return MATREC_OKAY;
}

void MATRECfreeMatrixPairDouble(MATREC *env, MATRECCompressedSparseMatrixPairDouble **matrixPair) {
    MATRECCompressedSparseMatrixPairDouble *pair = *matrixPair;
    MATRECfreeDoubleMatrix(env, &pair->rowMat);
    MATRECfreeDoubleMatrix(env, &pair->colMat);
    MATRECfreeBlock(env, matrixPair);
}

void MATRECwriteDoubleMatrixToStream(const MATRECCSMatrixDouble * matrix, FILE * file){
    assert(matrix);
    assert(file);

    fprintf(file, "%lu %lu %lu\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
    for (MATREC_matrix_size row = 0; row < matrix->numRows; ++row)
    {
        MATREC_matrix_size first = matrix->firstRowIndex[row];
        MATREC_matrix_size beyond = matrix->firstRowIndex[row + 1];
        for (MATREC_matrix_size entry = first; entry < beyond; ++entry){
            fprintf(file, "%lu %lu %lf\n", row+1, matrix->entryColumns[entry] + 1, matrix->entryValues[entry]);
        }
    }
}

MATREC_matrix_size MATRECdoubleMatrixNumRows(MATRECCSMatrixDouble * matrix){
    return matrix->numRows;
}
MATREC_matrix_size MATRECdoubleMatrixNumColumns(MATRECCSMatrixDouble * matrix){
    return matrix->numColumns;
}
MATREC_matrix_size MATRECdoubleMatrixNumNonzeros(MATRECCSMatrixDouble * matrix){
    return matrix->numNonzeros;
}
MATREC_matrix_size MATRECdoubleMatrixRowNumNonzeros(MATRECCSMatrixDouble * matrix,
                                                    MATREC_row row
){
    return matrix->firstRowIndex[row+1] - matrix->firstRowIndex[row];
}

MATREC_matrix_size * MATRECdoubleMatrixRowColumnIndices(MATRECCSMatrixDouble * matrix,
                                                        MATREC_row row
){
    return &matrix->entryColumns[matrix->firstRowIndex[row]];
}

double * MATRECdoubleMatrixRowColumnValues(MATRECCSMatrixDouble * matrix,
                                           MATREC_row row
){
    return &matrix->entryValues[matrix->firstRowIndex[row]];
}


MATREC_ERROR MATRECreadDoubleMatrixFromStream(MATREC* env, MATRECCSMatrixDouble ** presult, FILE* stream){
    assert(env);
    assert(presult);
    assert(stream);

    MATREC_matrix_size numRows, numColumns, numNonzeros;
    int numRead = fscanf(stream, "%lu %lu %lu", &numRows, &numColumns, &numNonzeros);
    if (numRead < 3)
    {
        return MATREC_ERROR_INPUT;
    }

    /* Read all nonzeros. */

    MATRECMatrixTripletDouble * nonzeros = NULL;


    MATREC_CALL( MATRECallocBlockArray(env, &nonzeros, numNonzeros) );
    char rowString[128];
    char colString[128];
    char valueString[1024];
    MATREC_matrix_size entry = 0;
    for (MATREC_matrix_size i = 0; i < numNonzeros; ++i)
    {
        MATREC_matrix_size row;
        MATREC_matrix_size column;
        double value;
        numRead = fscanf(stream, "%s %s %s", rowString, colString, valueString);
        numRead += sscanf(rowString,"%lu",&row);
        numRead += sscanf(colString,"%lu",&column);
        value =  strtod(valueString,NULL);

        if (numRead < 5  || row == 0 || column == 0 || row > numRows || column > numColumns)
        {
            MATRECfreeBlockArray(env,&numNonzeros);
            return MATREC_ERROR_INPUT;
        }
        if (value != 0)
        {
            nonzeros[entry].row = row - 1;
            nonzeros[entry].column = column - 1;
            nonzeros[entry].value = value;
            ++entry;
        }
    }
    numNonzeros = entry;

    /* We sort all nonzeros by row and then by column. */
    qsort(nonzeros, numNonzeros, sizeof(MATRECIntMatrixTriplet), compareIntNonzeros);

    MATREC_CALL(MATRECcreateDoubleMatrixWithNonzeros(env, presult, numRows, numColumns, numNonzeros, nonzeros));

    MATRECfreeBlockArray(env,&nonzeros);
    return MATREC_OKAY;
}

void MATRECtransposeSubmatrix(MATRECSubMatrix *submatrix) {
    MATREC_matrix_size temp = submatrix->numRows;
    submatrix->numRows = submatrix->numColumns;
    submatrix->numColumns = temp;

    MATREC_matrix_size * tempPointer = submatrix->rows;
    submatrix->rows = submatrix->columns;
    submatrix->columns = tempPointer;

}
