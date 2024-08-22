#ifndef MATREC_GRAPHIC_H
#define MATREC_GRAPHIC_H

#include "Shared.h"

#ifdef __cplusplus
extern "C"{
#endif

/**
 * This class stores the SPQR decomposition
 */
typedef struct MATRECGraphicDecompositionImpl MATRECGraphicDecomposition;


MATREC_ERROR MATRECGraphicDecompositionCreate(MATREC * env, MATRECGraphicDecomposition **pDecomposition, int numRows, int numColumns);

void MATRECGraphicDecompositionFree(MATRECGraphicDecomposition **pDecomposition);

/**
 * Returns if the MATREC decomposition contains the given row
 */
bool MATRECGraphicDecompositionContainsRow(const MATRECGraphicDecomposition * decomposition, MATREC_row row);
/**
 * Returns if the MATREC decomposition contains the given column
 */
bool MATRECGraphicDecompositionContainsColumn(const MATRECGraphicDecomposition *decomposition, MATREC_col column);

/**
 * Checks if the MATREC decomposition of the graph is minimal
 */
bool MATRECGraphicDecompositionIsMinimal(const MATRECGraphicDecomposition * decomposition);

//TODO: method to convert decomposition into a graphic realization
//TODO: method to remove complete components of the SPQR tree

/**
 * A method to check if the cycle stored in the SPQR cycle matches the given array. Mostly useful in testing.
 */
bool MATRECGraphicDecompositionVerifyCycle(const MATRECGraphicDecomposition * dec, MATREC_col column, MATREC_row * column_rows,
                                           int num_rows, MATREC_row * computed_column_storage);

/**
 * This class stores all data for performing sequential column additions to a matrix and checking if it is graphic or not.
 */
typedef struct MATRECGraphicColumnAdditionImpl MATRECGraphicColumnAddition;

/**
 * @brief Creates the data structure for managing column-addition for an SPQR decomposition
 */
MATREC_ERROR MATRECcreateGraphicColumnAddition(MATREC* env, MATRECGraphicColumnAddition** pNewCol );
/**
 * @brief Destroys the data structure for managing column-addition for SPQR decomposition
 */
void MATRECfreeGraphicColumnAddition(MATREC* env, MATRECGraphicColumnAddition ** pNewCol);

/**
 * Checks if adding a column of the given matrix creates a graphic SPQR decomposition.
 * Adding a column which is already in the decomposition is undefined behavior and not checked for.
 * @param dec Current SPQR-decomposition
 * @param newRow Data structure to store information on how to add the new column (if applicable).
 * @param column The index of the column to be added
 * @param rows An array with the row indices of the nonzero entries of the column.
 * @param numRows The number of nonzero entries of the column
 */
MATREC_ERROR MATRECGraphicColumnAdditionCheck(MATRECGraphicDecomposition * dec, MATRECGraphicColumnAddition * newCol, MATREC_col column, const MATREC_row * rows, size_t numRows);
/**
 * @brief Adds the most recently checked column from checkNewRow() to the Decomposition.
 * //TODO: specify (and implement) behavior in special cases (e.g. zero columns, columns with a single entry)
 * In Debug mode, adding a column for which MATRECGraphicColumnAdditionRemainsGraphic() returns false will exit the program.
 * In Release mode, adding a column for which MATRECGraphicColumnAdditionRemainsGraphic() return false is undefined behavior
 * @param dec Current SPQR-decomposition
 * @param newRow Data structure containing information on how to add the new column.
 */
MATREC_ERROR MATRECGraphicColumnAdditionAdd(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol);

/**
 * TODO: specify (and implement) behavior in special cases (e.g. zero columns, columns with a single entry)
 * @param newColumn
 * @return True if the most recently checked column is addable to the SPQR decomposition passed to it, i.e. the submatrix
 * given by both remains graphic.
 */
bool MATRECGraphicColumnAdditionRemainsGraphic(MATRECGraphicColumnAddition *newCol);


/**
 * This class stores all data for performing sequential row-additions to a matrix and checking if it is graphic or not.
 */
typedef struct MATRECGraphicRowAdditionImpl MATRECGraphicRowAddition;

/**
 * @brief Creates the data structure for managing row-addition for an SPQR decomposition
 */
MATREC_ERROR MATRECcreateGraphicRowAddition(MATREC* env, MATRECGraphicRowAddition** pNewRow );
/**
 * @brief Destroys the data structure for managing row-addition for SPQR decomposition
 */
void MATRECfreeGraphicRowAddition(MATREC* env, MATRECGraphicRowAddition ** pNewRow);

/**
 * Checks if adding a row of the given matrix creates a graphic SPQR decomposition.
 * Adding a row which is already in the decomposition is undefined behavior and not checked for.
 * @param dec Current SPQR-decomposition
 * @param newRow Data structure to store information on how to add the new row (if applicable).
 * @param row The index of the row to be added
 * @param columns An array with the column indices of the nonzero entries of the row.
 * @param numColumns The number of nonzero entries of the row
 */
MATREC_ERROR MATRECGraphicRowAdditionCheck(MATRECGraphicDecomposition * dec, MATRECGraphicRowAddition * newRow, MATREC_row row, const MATREC_col * columns, size_t numColumns);
/**
 * @brief Adds the most recently checked column from checkNewRow() to the Decomposition.
 * //TODO: specify (and implement) behavior in special cases (e.g. zero rows, rows with a single entry?)
 * In Debug mode, adding a column for which rowAdditionRemainsGraphic() returns false will exit the program.
 * In Release mode, adding a column for which rowAdditionRemainsGraphic() return false is undefined behavior
 * @param dec Current SPQR-decomposition
 * @param newRow Data structure containing information on how to add the new row.
 */
MATREC_ERROR MATRECGraphicRowAdditionAdd(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow);

/**
 * TODO: specify (and implement) behavior in special cases (e.g. zero rows, rows with a single entry?)
 * @param newRow
 * @return True if the most recently checked row is addable to the SPQR decomposition passed to it, i.e. the submatrix
 * given by both remains graphic.
 */
bool MATRECGraphicRowAdditionRemainsGraphic(const MATRECGraphicRowAddition *newRow);


#ifdef __cplusplus
}
#endif

#endif //MATREC_GRAPHIC_H
