#include "matrec/SignCheckRowAddition.h"
typedef struct {
    size_t predecessorNode;
    int predecessorValue;
    int targetValue;
    char status;
}ScaGraphNode;

struct MATRECSignCheckRowAdditionImpl{
    MATRECCompressedSparseMatrixPairInt * matrix;
    ScaGraphNode * graphNodes;
    size_t * bfsCallStack;
    bool * isRowAdded;
};

MATREC_ERROR MATRECcreateSignCheckRowAddition(MATREC * env, MATRECSignCheckRowAddition ** psca, MATRECCompressedSparseMatrixPairInt * matrix){
    assert(env);
    assert(psca);
    assert(*psca == NULL);
    assert(matrix);
    MATREC_CALL(MATRECallocBlock(env, psca));
    MATRECSignCheckRowAddition * sca = *psca;
    sca->matrix = matrix;
    assert(MATRECintMatrixIsTernary(matrix->rowMat) && MATRECintMatrixIsTernary(matrix->colMat));
    MATREC_CALL(MATRECallocBlockArray(env,&sca->graphNodes, matrix->rowMat->numRows + matrix->rowMat->numColumns));
    MATREC_CALL(MATRECallocBlockArray(env,&sca->bfsCallStack, matrix->rowMat->numRows + matrix->rowMat->numColumns));
    MATREC_CALL(MATRECallocBlockArray(env,&sca->isRowAdded,matrix->rowMat->numRows));
    for (size_t i = 0; i < matrix->rowMat->numRows; ++i) {
        sca->isRowAdded[i] = false;
    }
    return MATREC_OKAY;
}
void MATRECfreeSignCheckRowAddition(MATREC * env, MATRECSignCheckRowAddition ** psca){
    MATRECSignCheckRowAddition * sca = *psca;
    MATRECfreeBlockArray(env,&sca->isRowAdded);
    MATRECfreeBlockArray(env,&sca->bfsCallStack);
    MATRECfreeBlockArray(env,&sca->graphNodes);
    MATRECfreeBlock(env,&sca);
}
bool MATRECcheckSigningNewRow(MATRECSignCheckRowAddition * sca, MATREC_row row){
    assert(!sca->isRowAdded[row]); //If we already added the row, it is nonsensical to check its signing

    MATRECCSMatrixInt * rowMat = sca->matrix->rowMat;
    MATRECCSMatrixInt * colMat = sca->matrix->colMat;

    size_t first = rowMat->firstRowIndex[row];
    size_t beyond = rowMat->firstRowIndex[row + 1];
    if(first == beyond){
        //row is empty
        return true;
    }
    for (size_t i = 0; i < rowMat->numRows + rowMat->numColumns; ++i) {
        sca->graphNodes[i].targetValue = 0;
        sca->graphNodes[i].predecessorNode = MATREC_INVALID_ROW;
        sca->graphNodes[i].status = 0;
    }
    const size_t firstRowGraphNode = rowMat->numColumns;

    for (size_t i = first; i < beyond; ++i) {
        assert(abs(rowMat->entryValues[i]) == 1);
        sca->graphNodes[rowMat->entryColumns[i]].targetValue = rowMat->entryValues[i];
    }

    for (size_t rootIndex = first; rootIndex < beyond; ++rootIndex) {
        size_t rootSearchNode = rowMat->entryColumns[rootIndex];
        //Bfs marks the predecessors; if it is not marked, The entry is in a different component of the matrix
        //So we need to do BFS with that entry as a root again.
        if(sca->graphNodes[rootSearchNode].status != 0) continue;

        sca->bfsCallStack[0] = rootSearchNode;
        sca->graphNodes[rootSearchNode].status = 1;
        int bfsBegin = 0;
        int bfsEnd = 1;
        while(bfsBegin < bfsEnd){
            size_t currentNode = sca->bfsCallStack[bfsBegin];
            assert(sca->graphNodes[currentNode].status == 1);
            sca->graphNodes[currentNode].status = 2;
            if(currentNode >= firstRowGraphNode){
                size_t nodeRow = currentNode-firstRowGraphNode;

                size_t firstRow = rowMat->firstRowIndex[nodeRow];
                size_t beyondRow = rowMat->firstRowIndex[nodeRow + 1];
                for (size_t i = firstRow; i < beyondRow ; ++i) {
                    size_t entryCol = rowMat->entryColumns[i];
                    if(sca->graphNodes[entryCol].status == 0){
                        //If column is new, push it onto the bfs stack
                        sca->graphNodes[entryCol].status = 1;
                        sca->graphNodes[entryCol].predecessorNode = currentNode;
                        sca->graphNodes[entryCol].predecessorValue = rowMat->entryValues[i];

                        sca->bfsCallStack[bfsEnd] = entryCol;
                        bfsEnd++;

                        //If we reach a target column for the first time, trace back to the previous target column
                        if(sca->graphNodes[entryCol].targetValue != 0){
                            int sum = sca->graphNodes[entryCol].targetValue;
                            size_t pathNode = entryCol;
                            do{
                                sum += sca->graphNodes[pathNode].predecessorValue;
                                pathNode = sca->graphNodes[pathNode].predecessorNode;
                            }while(sca->graphNodes[pathNode].targetValue == 0);
                            sum += sca->graphNodes[pathNode].targetValue;

                            //By adding the first and final value we simplified this check a bit
                            if(sum % 4 != 0){
                                assert(sum % 4 == -2 || sum % 4 == 2);
                                return false;
                            }
                        }
                    }
                }
            }
            else{
                size_t nodeColumn = currentNode;

                //Iterate over outgoing edges (rows of the column)
                size_t firstCol = colMat->firstRowIndex[nodeColumn];
                size_t beyondCol = colMat->firstRowIndex[nodeColumn + 1];
                for (size_t i = firstCol; i < beyondCol ; ++i) {
                    size_t entryRow = colMat->entryColumns[i];
                    if(!sca->isRowAdded[entryRow]) continue;

                    size_t nodeIndex = firstRowGraphNode + entryRow;
                    //If row is new, push it onto the bfs stack
                    if(sca->graphNodes[nodeIndex].status == 0){
                        sca->graphNodes[nodeIndex].status = 1;
                        sca->graphNodes[nodeIndex].predecessorNode = currentNode;
                        sca->graphNodes[nodeIndex].predecessorValue = colMat->entryValues[i];
                        sca->bfsCallStack[bfsEnd] = nodeIndex;
                        bfsEnd++;
                    }
                }
            }
            ++bfsBegin;
        }
    }

    return true;
}
MATREC_ERROR MATRECaddSigningNewRow(MATRECSignCheckRowAddition * sca, MATREC_row row){
    assert(MATRECcheckSigningNewRow(sca, row)); // we double-check that the user has actually checked that the row can be added
    sca->isRowAdded[row] = true;

    return MATREC_OKAY;
}

