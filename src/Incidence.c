#include "matrec/Incidence.h"

typedef struct {
    bool picked;
    int numEntries;
    int lastEntryRow;
    int lastEntrySign;
} ColumnInfo;

typedef struct {
    int edgeSign;
    int representative;
} SignedUnionFindInfo;

typedef struct {
    int numNormal;
    int numReflected;
} RowComponentInfo;
struct MATRECIncidenceAdditionImpl{
    ColumnInfo * columnInfo;
    bool * rowPicked;
    SignedUnionFindInfo * unionFind; //first the rows, then the columns

    //Temporary storage. We allocate this once to prevent reallocations for every call
    RowComponentInfo * rowComponentInfo;
    int * componentRepresentatives;
    int numComponentRepresentatives;
};

MATREC_ERROR MATRECcreateIncidenceAddition(MATREC* env,
                                           MATRECIncidenceAddition ** pIncidenceAddition,
                                           MATREC_matrix_size numRows,
                                           MATREC_matrix_size numColumns,
                                           MATRECIncidenceAdditionInit init
){
    MATREC_CALL(MATRECallocBlock(env,pIncidenceAddition));
    MATRECIncidenceAddition* incidenceAddition = *pIncidenceAddition;
    MATREC_CALL(MATRECallocBlockArray(env,&incidenceAddition->columnInfo,numColumns));
    bool initialCols = init == MATREC_INIT_ALL_COLUMNS ? true : false;

    for (MATREC_matrix_size i = 0; i < numColumns; ++i) {
        incidenceAddition->columnInfo[i].picked = initialCols;
        incidenceAddition->columnInfo[i].numEntries = 0;
        incidenceAddition->columnInfo[i].lastEntryRow = -1;
        incidenceAddition->columnInfo[i].lastEntrySign = 0;
    }
    MATREC_CALL(MATRECallocBlockArray(env,&incidenceAddition->rowPicked,numRows));

    bool initialRows = init == MATREC_INIT_ALL_ROWS ? true : false;
    for(MATREC_matrix_size i = 0; i < numRows; ++i){
        incidenceAddition->rowPicked[i] = initialRows;
    }


    MATREC_CALL(MATRECallocBlockArray(env,&incidenceAddition->unionFind,numRows));
    for (MATREC_matrix_size i = 0; i < numRows; ++i) {
        incidenceAddition->unionFind[i].representative = -1;
        incidenceAddition->unionFind[i].edgeSign = 1;
    }

    MATREC_CALL(MATRECallocBlockArray(env,&incidenceAddition->rowComponentInfo,numRows));
    for (MATREC_matrix_size i = 0; i < numRows; ++i) {
        incidenceAddition->rowComponentInfo[i].numNormal = 0;
        incidenceAddition->rowComponentInfo[i].numReflected = 0;
    }
    MATREC_CALL(MATRECallocBlockArray(env,&incidenceAddition->componentRepresentatives,numRows));
    //Memory initialized because all entries are set before accessing during algorithm

    incidenceAddition->numComponentRepresentatives = 0;

    return MATREC_OKAY;
}

void MATRECfreeIncidenceAddition(MATREC* env,
                                 MATRECIncidenceAddition ** pIncidenceAddition
){
    MATRECIncidenceAddition * incidenceAddition = *pIncidenceAddition;
    MATRECfreeBlockArray(env,&incidenceAddition->componentRepresentatives);
    MATRECfreeBlockArray(env,&incidenceAddition->rowComponentInfo);
    MATRECfreeBlockArray(env,&incidenceAddition->rowPicked);
    MATRECfreeBlockArray(env,&incidenceAddition->columnInfo);
    MATRECfreeBlockArray(env,&incidenceAddition->unionFind);
    MATRECfreeBlock(env,pIncidenceAddition);
}
bool isNegative(int row){
    return row < 0;
}

void findRowRepresentative(MATRECIncidenceAddition * addition, int row,
                           int * representative, int * sign){
    assert(addition);

    int current = row;
    int next;

    int totalSign = 1;
    //traverse down tree to find the root
    while (!isNegative(next = addition->unionFind[current].representative)) {
        totalSign *= addition->unionFind[current].edgeSign;
        current = next;
    }

    int root = current;
    current = row;

    int currentSign = totalSign;
    //update all pointers along path to point to root, flattening the tree
    while (!isNegative(next = addition->unionFind[current].representative)) {
        addition->unionFind[current].representative = root;
        int previousSign = addition->unionFind[current].edgeSign;
        addition->unionFind[current].edgeSign = currentSign;
        currentSign *= previousSign;
        current = next;
    }
    *representative = root;
    *sign = totalSign;
}

int mergeRowRepresentatives(MATRECIncidenceAddition * addition,
                            int first,
                            int second,
                            int sign) {
    assert(first != second); //We cannot merge a member into itself

    //The rank is stored as a negative number: we decrement it making the negative number larger.
    // We want the new root to be the one with 'largest' rank, so smallest number. If they are equal, we decrement.
    //This 'rank storing' scheme is necessary to ensure the inverse ackermann time complexity
    int firstRank = addition->unionFind[first].representative;
    int secondRank = addition->unionFind[second].representative;
    if (firstRank > secondRank) {
        int temp = first;
        first = second;
        second = temp;
    }
    assert(addition->unionFind[first].edgeSign == 1 && addition->unionFind[second].edgeSign == 1);
    addition->unionFind[second].representative = first;
    addition->unionFind[second].edgeSign = sign;
    if (firstRank == secondRank) {
        --addition->unionFind[first].representative;
    }
    return first;
}

void cleanUpRowTemporaryStorage(MATRECIncidenceAddition *addition){
    for (int i = 0; i < addition->numComponentRepresentatives; ++i) {
        int row = addition->componentRepresentatives[i];
        addition->rowComponentInfo[row].numNormal = 0;
        addition->rowComponentInfo[row].numReflected = 0;
    }
    //TODO: assert that all component info is set to zero
    addition->numComponentRepresentatives = 0;
}

bool MATRECincidenceAdditionAddRow(MATRECIncidenceAddition *addition,
                                   MATREC_row row,
                                   MATREC_matrix_size nRowNonzeros,
                                   const MATREC_col *entryColumns,
                                   const int *entryValues
){

    assert(addition->numComponentRepresentatives == 0);

    for(MATREC_matrix_size i = 0; i < nRowNonzeros; ++i){
        MATREC_col column = entryColumns[i];
        if(!addition->columnInfo[column].picked) continue;
        if(addition->columnInfo[column].numEntries >= 2){
            cleanUpRowTemporaryStorage(addition);
            return false;
        }
        //Find to which component the columns belong (along with their signs)
        assert((addition->columnInfo[column].numEntries == 0 ) ==
        (addition->columnInfo[column].lastEntryRow == -1));
        int componentRow = addition->columnInfo[column].lastEntryRow;
        if(componentRow != -1){
            int rowSign;
            int representative;
            findRowRepresentative(addition,componentRow,&representative,&rowSign);
            if(addition->rowComponentInfo[representative].numNormal == 0 &&
            addition->rowComponentInfo[representative].numReflected == 0){
                addition->componentRepresentatives[addition->numComponentRepresentatives] = representative;
                ++addition->numComponentRepresentatives;
            }

            int otherValue = addition->columnInfo[column].lastEntrySign * rowSign;
            assert(otherValue == 1 || otherValue == -1);
            assert(entryValues[i] == 1 || entryValues[i] == -1);
            if(entryValues[i] == otherValue){
                ++addition->rowComponentInfo[representative].numReflected;
            }else{
                ++addition->rowComponentInfo[representative].numNormal;
            }
        }
    }
    for (int i = 0; i < addition->numComponentRepresentatives; ++i) {
        int representative = addition->componentRepresentatives[i];
        assert(addition->rowComponentInfo[representative].numNormal > 0 ||
        addition->rowComponentInfo[representative].numReflected > 0);
        if(addition->rowComponentInfo[representative].numNormal != 0 &&
        addition->rowComponentInfo[representative].numReflected != 0){
            cleanUpRowTemporaryStorage(addition);
            return false;
        }
    }
    //Otherwise, we are good
    addition->rowPicked[row] = true;

    //If the component with the root in it is multiplied by -1, we need to change the multiplier
    //of the other components during merging. We do this using the following multiplier
    int rootMovedSign = 1;

    int thisRepresentative = (int) row;
    for (int i = 0; i < addition->numComponentRepresentatives; ++i) {
        int representative = addition->componentRepresentatives[i];
        assert(addition->rowComponentInfo[representative].numReflected == 0 || addition->rowComponentInfo[representative].numNormal == 0);
        int sign = addition->rowComponentInfo[representative].numNormal == 0 ? -1 : 1;
        sign *= rootMovedSign;
        int newRepresentative = mergeRowRepresentatives(addition,representative,thisRepresentative,sign);
        if(newRepresentative != thisRepresentative){
            rootMovedSign *= addition->unionFind[thisRepresentative].edgeSign;
        }
        thisRepresentative = newRepresentative;
    }
    int thisRow = (int) row;
    //If the row could be added;
    //fixup column information
    for(MATREC_matrix_size i = 0; i < nRowNonzeros; ++i){
        MATREC_col column = entryColumns[i];
        if(!addition->columnInfo[column].picked) continue;
        if(addition->columnInfo[column].numEntries == 0){
            addition->columnInfo[column].lastEntryRow = thisRow;
            addition->columnInfo[column].lastEntrySign = entryValues[i];
        }
        ++addition->columnInfo[column].numEntries;
    }
    cleanUpRowTemporaryStorage(addition);
    return true;
}
bool MATRECincidenceAdditionAddColumn(MATRECIncidenceAddition * addition,
                                      MATREC_col column,
                                      MATREC_matrix_size nColumnNonzeros,
                                      const MATREC_row * columnRows,
                                      const int * columnValues) {
    int entryRowRepresentative[2];
    int entryRow = -1;
    int entrySign = 0;

    int signSum = 0;
    int numEntries = 0;
    for (MATREC_matrix_size i = 0; i < nColumnNonzeros; ++i) {
        int row = (int) columnRows[i];
        if(!addition->rowPicked[row]) continue;

        if(numEntries >= 2){
            return false; //We found a third entry in the column; No way to add it.
        }
        int rowSign;
        findRowRepresentative(addition,row,&entryRowRepresentative[numEntries],&rowSign);
        assert(entryRowRepresentative[numEntries] >= 0);
        assert(rowSign == 1 || rowSign == -1);
        assert(columnValues[i] == 1 || columnValues[i] == -1);
        signSum += rowSign * columnValues[i];
        entrySign = columnValues[i];
        entryRow = row;

        ++numEntries;
    }
    if(numEntries <= 1){
        addition->columnInfo[column].picked = true;
        addition->columnInfo[column].numEntries = numEntries;
        if(numEntries != 0){
            addition->columnInfo[column].lastEntryRow = entryRow;
            addition->columnInfo[column].lastEntrySign = entrySign;
        }
        return true;
    }
    assert(numEntries == 2);
    if(entryRowRepresentative[0] == entryRowRepresentative[1]){
        //Columns are in the same component. Then, we can only add it if the signs differ correctly
        if(signSum == 0){
            addition->columnInfo[column].picked = true;
            addition->columnInfo[column].lastEntryRow = entryRow;
            addition->columnInfo[column].lastEntrySign = entrySign;
            addition->columnInfo[column].numEntries = 2;
            return true;
        }
        return false;
    }
    assert(entryRowRepresentative[0] != entryRowRepresentative[1]);
    //merge the two components into one, optionally reflecting one if necessary
    //if signSum == 0, we do not need to change the reflection, otherwise we do need to do so
    int sign = signSum == 0 ? 1 : -1;
    mergeRowRepresentatives(addition, entryRowRepresentative[0], entryRowRepresentative[1], sign);
    addition->columnInfo[column].picked = true;
    addition->columnInfo[column].lastEntryRow = entryRow;
    addition->columnInfo[column].lastEntrySign = entrySign;
    addition->columnInfo[column].numEntries = 2;
    return true;
}

bool MATRECincidenceContainsNonemptyColumn(MATRECIncidenceAddition * addition, MATREC_col column){
    assert(addition);
    return addition->columnInfo[column].picked && addition->columnInfo[column].numEntries > 0;
}
bool MATRECincidenceContainsColumn(MATRECIncidenceAddition * addition, MATREC_col column){
    assert(addition);
    return addition->columnInfo[column].picked;
}
bool MATRECincidenceContainsRow(MATRECIncidenceAddition * addition, MATREC_row row){
    assert(addition);
    return addition->rowPicked[row];
}

///Returns 1 if the row is not in the incidence submatrix
int MATRECincidenceRowSign(MATRECIncidenceAddition * addition, MATREC_row row){
    int sign, representative;
    findRowRepresentative(addition,(int) row,&representative,&sign);
    return sign;
}
