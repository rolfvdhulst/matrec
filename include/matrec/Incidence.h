#ifndef MATREC_INCIDENCEADDITION_H
#define MATREC_INCIDENCEADDITION_H

#include "Shared.h"
#include "Matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MATRECIncidenceAdditionImpl MATRECIncidenceAddition;
typedef enum{
    MATREC_INIT_NONE,
    MATREC_INIT_ALL_ROWS, //Start with all rows already in the incidence (sub)-matrix
    MATREC_INIT_ALL_COLUMNS, //Start with all columns already in the incidence (sub)-matrix
} MATRECIncidenceAdditionInit;

MATREC_ERROR MATRECcreateIncidenceAddition(MATREC* env,
                                           MATRECIncidenceAddition ** pIncidenceAddition,
                                           MATREC_matrix_size numRows,
                                           MATREC_matrix_size numColumns,
                                           MATRECIncidenceAdditionInit init
                                    );

void MATRECfreeIncidenceAddition(MATREC* env,
                                 MATRECIncidenceAddition ** pIncidenceAddition
                           );

bool MATRECincidenceAdditionAddRow(MATRECIncidenceAddition *addition,
                                   MATREC_row row,
                                   MATREC_matrix_size nRowNonzeros,
                                   const MATREC_col *entryColumns,
                                   const int *entryValues
                             );

bool MATRECincidenceAdditionAddColumn(MATRECIncidenceAddition * addition,
                                      MATREC_col column,
                                      MATREC_matrix_size nColumnNonzeros,
                                      const MATREC_row * columnRows,
                                      const int * columnValues
                                );

bool MATRECincidenceContainsNonemptyColumn(MATRECIncidenceAddition * addition, MATREC_col column);
bool MATRECincidenceContainsColumn(MATRECIncidenceAddition * addition, MATREC_col column);
bool MATRECincidenceContainsRow(MATRECIncidenceAddition * addition, MATREC_row row);

///Returns 1 if the row is not in the incidence submatrix
int MATRECincidenceRowSign(MATRECIncidenceAddition * addition, MATREC_row row);
#ifdef __cplusplus
}
#endif

#endif //MATREC_INCIDENCEADDITION_H
