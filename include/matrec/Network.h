#ifndef MATREC_NETWORK_H
#define MATREC_NETWORK_H

#ifdef __cplusplus
extern "C"{
#endif

#include "Shared.h"

/**
 * This class stores the SPQR decomposition
 */
typedef struct MATRECNetworkDecompositionImpl MATRECNetworkDecomposition;


MATREC_ERROR MATRECNetworkDecompositionCreate(MATREC * env, MATRECNetworkDecomposition **pDecomposition, int numRows, int numColumns);

void MATRECNetworkDecompositionFree(MATRECNetworkDecomposition **pDecomposition);

/**
 * Returns if the MATREC decomposition contains the given row
 */
bool MATRECNetworkDecompositionContainsRow(const MATRECNetworkDecomposition * decomposition, MATREC_row row);
/**
 * Returns if the MATREC decomposition contains the given column
 */
bool MATRECNetworkDecompositionContainsColumn(const MATRECNetworkDecomposition *decomposition, MATREC_col column);

/**
 * Checks if the MATREC decomposition of the graph is minimal
 */
bool MATRECNetworkDecompositionIsMinimal(const MATRECNetworkDecomposition * decomposition);

//TODO: method to convert decomposition into a realization
//TODO: method to remove complete components of the MATREC tree

/**
 * A method to check if the cycle stored in the MATREC cycle matches the given array. Mostly useful in testing.
 */
bool MATRECNetworkDecompositionVerifyCycle(const MATRECNetworkDecomposition * dec, MATREC_col column, MATREC_row * column_rows,
                                           double * column_values, int num_rows,
                                           MATREC_row * computed_column_storage,
                                           bool * computedSignStorage);

/**
 * This class stores all data for performing sequential column additions to a matrix and checking if it is network or not.
 */
typedef struct MATRECNetworkColumnAdditionImpl MATRECNetworkColumnAddition;

/**
 * @brief Creates the data structure for managing column-addition for an MATREC decomposition
 */
MATREC_ERROR MATRECcreateNetworkColumnAddition(MATREC* env, MATRECNetworkColumnAddition** pNewCol );
/**
 * @brief Destroys the data structure for managing column-addition for MATREC decomposition
 */
void MATRECfreeNetworkColumnAddition(MATREC* env, MATRECNetworkColumnAddition ** pNewCol);

/**
 * Checks if adding a column of the given matrix creates a network MATREC decomposition.
 * Adding a column which is already in the decomposition is undefined behavior and not checked for.
 * @param dec Current MATREC-decomposition
 * @param newRow Data structure to store information on how to add the new column (if applicable).
 * @param column The index of the column to be added
 * @param rows An array with the row indices of the nonzero entries of the column.
 * @param numRows The number of nonzero entries of the column
 */
MATREC_ERROR MATRECNetworkColumnAdditionCheck(MATRECNetworkDecomposition * dec, MATRECNetworkColumnAddition * newCol, MATREC_col column,
                                              const MATREC_row * nonzeroRows, const double * nonzeroValues, size_t numNonzeros);
/**
 * @brief Adds the most recently checked column from checkNewRow() to the Decomposition.
 * //TODO: specify (and implement) behavior in special cases (e.g. zero columns, columns with a single entry)
 * In Debug mode, adding a column for which MATRECNetworkColumnAdditionRemainsNetwork() returns false will exit the program.
 * In Release mode, adding a column for which MATRECNetworkColumnAdditionRemainsNetwork() return false is undefined behavior
 * @param dec Current MATREC-decomposition
 * @param newRow Data structure containing information on how to add the new column.
 */
MATREC_ERROR MATRECNetworkColumnAdditionAdd(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol);

/**
 * TODO: specify (and implement) behavior in special cases (e.g. zero columns, columns with a single entry)
 * @param newColumn
 * @return True if the most recently checked column is addable to the MATREC decomposition passed to it, i.e. the submatrix
 * given by both remains network.
 */
bool MATRECNetworkColumnAdditionRemainsNetwork(MATRECNetworkColumnAddition *newCol);


/**
 * This class stores all data for performing sequential row-additions to a matrix and checking if it is network or not.
 */
typedef struct MATRECNetworkRowAdditionImpl MATRECNetworkRowAddition;

/**
 * @brief Creates the data structure for managing row-addition for an MATREC decomposition
 */
MATREC_ERROR MATRECcreateNetworkRowAddition(MATREC* env, MATRECNetworkRowAddition** pNewRow );
/**
 * @brief Destroys the data structure for managing row-addition for MATREC decomposition
 */
void MATRECfreeNetworkRowAddition(MATREC* env, MATRECNetworkRowAddition ** pNewRow);

/**
 * Checks if adding a row of the given matrix creates a network MATREC decomposition.
 * Adding a row which is already in the decomposition is undefined behavior and not checked for.
 * @param dec Current MATREC-decomposition
 * @param newRow Data structure to store information on how to add the new row (if applicable).
 * @param row The index of the row to be added
 * @param columns An array with the column indices of the nonzero entries of the row.
 * @param numColumns The number of nonzero entries of the row
 */
MATREC_ERROR MATRECNetworkRowAdditionCheck(MATRECNetworkDecomposition * dec, MATRECNetworkRowAddition * newRow, MATREC_row row,
                                           const MATREC_col * nonzeroCols, const double * nonzeroValues, size_t numNonzeros);
/**
 * @brief Adds the most recently checked column from checkNewRow() to the Decomposition.
 * //TODO: specify (and implement) behavior in special cases (e.g. zero rows, rows with a single entry?)
 * In Debug mode, adding a column for which rowAdditionRemainsNetwork() returns false will exit the program.
 * In Release mode, adding a column for which rowAdditionRemainsNetwork() return false is undefined behavior
 * @param dec Current MATREC-decomposition
 * @param newRow Data structure containing information on how to add the new row.
 */
MATREC_ERROR MATRECNetworkRowAdditionAdd(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow);

/**
 * TODO: specify (and implement) behavior in special cases (e.g. zero rows, rows with a single entry?)
 * @param newRow
 * @return True if the most recently checked row is addable to the MATREC decomposition passed to it, i.e. the submatrix
 * given by both remains network.
 */
bool MATRECNetworkRowAdditionRemainsNetwork(const MATRECNetworkRowAddition *newRow);


#ifdef __cplusplus
}
#endif

#endif //MATREC_NETWORK_H
