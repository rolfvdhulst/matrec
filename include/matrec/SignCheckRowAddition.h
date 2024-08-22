#ifndef MATREC_SIGNEDCHECK_H
#define MATREC_SIGNEDCHECK_H
#ifdef __cplusplus
extern "C"{
#endif

#include "Shared.h"
#include "Matrix.h"

typedef struct MATRECSignCheckRowAdditionImpl MATRECSignCheckRowAddition;

/**
 * @brief Don't modify matrix at point during the call of this algorithm.
 * Modification of any rows which are already added by SignCheckAddition may result in the wrong answers or undefined behavior.
 */
MATREC_ERROR MATRECcreateSignCheckRowAddition(MATREC *env, MATRECSignCheckRowAddition **sca, MATRECCompressedSparseMatrixPairInt *matrix);

void MATRECfreeSignCheckRowAddition(MATREC *env, MATRECSignCheckRowAddition ** sca);

bool MATRECcheckSigningNewRow(MATRECSignCheckRowAddition *sca, MATREC_row row);

MATREC_ERROR MATRECaddSigningNewRow(MATRECSignCheckRowAddition *sca, MATREC_row row);

#ifdef __cplusplus
}
#endif

#endif //MATREC_SIGNEDCHECK_H
