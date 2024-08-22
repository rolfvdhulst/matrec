#include <gtest/gtest.h>
#include "TestHelpers.h"
#include <matrec/Network.h>
#include <matrec/Graphic.h>
#include <matrec/SignCheckRowAddition.h>

MATREC_ERROR runGraphicCheck(MATREC * env,
        const DirectedColTestCase& testCase,
        bool& isGraphic){

    MATRECGraphicDecomposition * dec = NULL;
    MATREC_CALL(MATRECGraphicDecompositionCreate(env,&dec,testCase.rows,testCase.cols));
    MATRECGraphicColumnAddition * addition = NULL;
    MATREC_CALL(MATRECcreateGraphicColumnAddition(env,&addition));
    std::vector<MATREC_row> rows;

    bool runSuccessful = true;
    for (std::size_t col = 0; col < testCase.cols; ++col) {
        rows.clear();
        for(const auto& nonz : testCase.matrix[col]){
            rows.push_back(nonz.index);
        }
        MATREC_CALL(MATRECGraphicColumnAdditionCheck(dec, addition, col, rows.data(),rows.size()));
        if (MATRECGraphicColumnAdditionRemainsGraphic(addition)) {
            MATREC_CALL(MATRECGraphicColumnAdditionAdd(dec, addition));
        } else {
            runSuccessful = false;
            break;
        }
    }
    isGraphic = runSuccessful;
    MATRECfreeGraphicColumnAddition(env,&addition);
    MATRECGraphicDecompositionFree(&dec);
    return MATREC_OKAY;
}
MATREC_ERROR runSignedCheck(MATREC * env,
                           const DirectedColTestCase& testCase,
                           bool& isSigned){

    std::vector<MATRECIntMatrixTriplet> triplets;
    for(std::size_t i = 0 ; i < testCase.cols; ++i){
        for(const auto& entry : testCase.matrix[i]){
            triplets.push_back(MATRECIntMatrixTriplet{.row = entry.index,.column = i,.value = (entry.value > 0.0 ) ? 1 : -1});
        }
    }
    std::sort(triplets.begin(),triplets.end(),
    [](const MATRECIntMatrixTriplet &a , const MATRECIntMatrixTriplet & b){
        if(a.row == b.row){
            return a.column < b.column;
        }
        return a.row < b.row;
    });
    MATRECCSMatrixInt * rowMatrix = NULL;
    MATREC_CALL(MATRECcreateIntMatrixWithNonzeros(env, &rowMatrix, testCase.rows, testCase.cols, triplets.size(),
                                                  triplets.data()));

    MATRECCompressedSparseMatrixPairInt * matrixPair = NULL;
    MATREC_CALL(MATRECcreateIntMatrixPair(env, rowMatrix, &matrixPair));

    MATRECSignCheckRowAddition * signcheck = NULL;
    MATREC_CALL(MATRECcreateSignCheckRowAddition(env, &signcheck, matrixPair));
    bool runSuccessful = true;
    for (std::size_t i = 0; i < testCase.rows; ++i) {
        if(!MATRECcheckSigningNewRow(signcheck, i)){
            runSuccessful = false;
            break;
        }
        MATREC_CALL(MATRECaddSigningNewRow(signcheck, i));
    }
    isSigned = runSuccessful;
    MATRECfreeSignCheckRowAddition(env, &signcheck);
    MATRECfreeIntMatrixPair(env, &matrixPair);
    return MATREC_OKAY;
}

MATREC_ERROR runRowGraphicCheck(MATREC * env,
                           const DirectedTestCase& testCase,
                           bool& isGraphic){

    MATRECGraphicDecomposition * dec = NULL;
    MATREC_CALL(MATRECGraphicDecompositionCreate(env,&dec,testCase.rows,testCase.cols));
    MATRECGraphicRowAddition * addition = NULL;
    MATREC_CALL(MATRECcreateGraphicRowAddition(env,&addition));
    std::vector<MATREC_col> cols;

    bool runSuccessful = true;
    for (std::size_t row = 0; row < testCase.rows; ++row) {
        cols.clear();
        for(const auto& nonz : testCase.matrix[row]){
            cols.push_back(nonz.index);
        }
        MATREC_CALL(MATRECGraphicRowAdditionCheck(dec, addition, row, cols.data(),cols.size()));
        if (MATRECGraphicRowAdditionRemainsGraphic(addition)) {
            MATREC_CALL(MATRECGraphicRowAdditionAdd(dec, addition));
        } else {
            runSuccessful = false;
            break;
        }
    }
    isGraphic = runSuccessful;
    MATRECfreeGraphicRowAddition(env,&addition);
    MATRECGraphicDecompositionFree(&dec);
    return MATREC_OKAY;
}
MATREC_ERROR runRowSignedCheck(MATREC * env,
                          const DirectedTestCase& testCase,
                          bool& isSigned){

    std::vector<MATRECIntMatrixTriplet> triplets;
    for(std::size_t i = 0 ; i < testCase.rows; ++i){
        for(const auto& entry : testCase.matrix[i]){
            triplets.push_back(MATRECIntMatrixTriplet{.row = i,.column = entry.index,
                                                .value = ( entry.value > 0.0 ) ? 1 : -1});
        }
    }
    std::sort(triplets.begin(),triplets.end(),
              [](const MATRECIntMatrixTriplet &a , const MATRECIntMatrixTriplet & b){
                  if(a.row == b.row){
                      return a.column < b.column;
                  }
                  return a.row < b.row;
              });
    MATRECCSMatrixInt * rowMatrix = NULL;
    MATREC_CALL(MATRECcreateIntMatrixWithNonzeros(env, &rowMatrix, testCase.rows, testCase.cols, triplets.size(),
                                                  triplets.data()));

    MATRECCompressedSparseMatrixPairInt * matrixPair = NULL;
    MATREC_CALL(MATRECcreateIntMatrixPair(env, rowMatrix, &matrixPair));

    MATRECSignCheckRowAddition * signcheck = NULL;
    MATREC_CALL(MATRECcreateSignCheckRowAddition(env, &signcheck, matrixPair));
    bool runSuccessful = true;
    for (std::size_t i = 0; i < testCase.rows; ++i) {
        if(!MATRECcheckSigningNewRow(signcheck, i)){
            runSuccessful = false;
            break;
        }
        MATREC_CALL(MATRECaddSigningNewRow(signcheck, i));
    }
    isSigned = runSuccessful;
    MATRECfreeSignCheckRowAddition(env, &signcheck);
    MATRECfreeIntMatrixPair(env, &matrixPair);
    return MATREC_OKAY;
}
MATREC_ERROR runNetworkDecomposition(DirectedColTestCase& testCase,
                                   bool& isGood,
                                   std::optional<bool> expectNetwork){

    MATREC *env = NULL;
    MATREC_CALL(MATRECcreateEnvironment(&env));
    bool * signStorage = NULL;
    MATREC_CALL(MATRECallocBlockArray(env,&signStorage,testCase.rows));
    std::vector<MATREC_row> rowStorage(testCase.rows,MATREC_INVALID_ROW);
    MATRECNetworkDecomposition *dec = NULL;
    MATREC_CALL(MATRECNetworkDecompositionCreate(env, &dec, testCase.rows, testCase.cols));
    MATRECNetworkColumnAddition *newCol = NULL;
    MATREC_CALL(MATRECcreateNetworkColumnAddition(env, &newCol));
    {
        std::vector<MATREC_row> rows;
        std::vector<double> values;
        bool isNetwork = true;
        for (std::size_t col = 0; col < testCase.cols; ++col) {
            rows.clear();
            values.clear();
            for(const auto& nonz : testCase.matrix[col]){
                rows.push_back(nonz.index);
                values.push_back(nonz.value);
            }
            MATREC_CALL(MATRECNetworkColumnAdditionCheck(dec, newCol, col, rows.data(),values.data(),rows.size()));
            if (MATRECNetworkColumnAdditionRemainsNetwork(newCol)) {
                MATREC_CALL(MATRECNetworkColumnAdditionAdd(dec, newCol));
            } else {
                isNetwork = false;
                break;
            }

            bool isMinimal = MATRECNetworkDecompositionIsMinimal(dec);
            EXPECT_TRUE(isMinimal); //Check that there are no series-series and bond-bond connections
            //TODO: fix
            bool fundamental_cycles_good = true;
            std::vector<MATREC_row> colRows;
            std::vector<double> colValues;
            for (std::size_t check_col = 0; check_col <= col; ++check_col) {
                colRows.clear();
                colValues.clear();
                for(const auto& elem : testCase.matrix[check_col]){
                    colRows.push_back(elem.index);
                    colValues.push_back(elem.value);
                }
                bool correct = MATRECNetworkDecompositionVerifyCycle(dec, check_col, colRows.data(),colValues.data(),colRows.size(),
                                                      rowStorage.data(),signStorage);
                EXPECT_TRUE(correct);
                if (!correct) {
                    fundamental_cycles_good = false;
                }
            }
            if (!fundamental_cycles_good) {
                isGood = false;
                break;
            }
        }
        if(expectNetwork.has_value()){
            EXPECT_EQ(expectNetwork.value(),isNetwork);
        }else{
            bool isGraphic = false;
            bool isSigned = false;
            if(isNetwork){
               MATREC_CALL(runGraphicCheck(env,testCase,isGraphic));
               EXPECT_TRUE(isGraphic);
               MATREC_CALL(runSignedCheck(env,testCase,isSigned));
               EXPECT_TRUE(isSigned);
            }else{
                MATREC_CALL(runGraphicCheck(env,testCase,isGraphic));
                if(isGraphic){
                    MATREC_CALL(runSignedCheck(env,testCase,isSigned));
                    EXPECT_FALSE(isSigned);
                }
            }
        }

    }
    MATRECfreeNetworkColumnAddition(env, &newCol);
    MATRECNetworkDecompositionFree(&dec);
    MATREC_CALL(MATRECfreeEnvironment(&env));
    return MATREC_OKAY;
}

void runTestCase(const DirectedTestCase &testCase, std::optional<bool> expectNetwork = std::nullopt) {
    ASSERT_EQ(testCase.rows, testCase.matrix.size());
    DirectedColTestCase colTestCase(testCase);

    bool isGood = true;

    MATREC_ERROR error = runNetworkDecomposition(colTestCase, isGood, expectNetwork);
    EXPECT_EQ(error, MATREC_OKAY);
    EXPECT_TRUE(isGood);
}

MATREC_ERROR runRowNetworkDecomposition(const DirectedTestCase& testCase,
                                      bool& isGood,
                                      std::optional<bool> expectNetwork = std::nullopt){
    MATREC *env = NULL;
    MATREC_CALL(MATRECcreateEnvironment(&env));
    bool * signStorage = NULL;
    MATREC_CALL(MATRECallocBlockArray(env,&signStorage,testCase.rows));
    std::vector<MATREC_row> rowStorage(testCase.rows,MATREC_INVALID_ROW);
    MATRECNetworkDecomposition *dec = NULL;
    MATREC_CALL(MATRECNetworkDecompositionCreate(env, &dec, testCase.rows, testCase.cols));
    MATRECNetworkRowAddition *newRow = NULL;
    MATREC_CALL(MATRECcreateNetworkRowAddition(env, &newRow));
    std::vector<std::vector<MATREC_row>> columnCycles(testCase.cols,std::vector<MATREC_col>());
    std::vector<std::vector<double>> columnCycleValues(testCase.cols,std::vector<double>());
    {
        std::vector<MATREC_col> cols;
        std::vector<double> values;
        bool isNetwork = true;
        for (std::size_t row = 0; row < testCase.rows; ++row) {
            cols.clear();
            values.clear();
            for(const auto& nonz : testCase.matrix[row]){
                cols.push_back(nonz.index);
                values.push_back(nonz.value);
            }
            MATREC_CALL(MATRECNetworkRowAdditionCheck(dec, newRow, row, cols.data(),values.data(),cols.size()));
            if (MATRECNetworkRowAdditionRemainsNetwork(newRow)) {
                MATREC_CALL(MATRECNetworkRowAdditionAdd(dec, newRow));
            } else {
                isNetwork = false;
                break;
            }

            for(const auto& nonz : testCase.matrix[row]){
                columnCycles[nonz.index].push_back(row);
                columnCycleValues[nonz.index].push_back(nonz.value);
            }
            bool isMinimal = MATRECNetworkDecompositionIsMinimal(dec);
            EXPECT_TRUE(isMinimal); //Check that there are no series-series and bond-bond connections
            //TODO: fix
            bool fundamental_cycles_good = true;
            for (std::size_t check_col = 0; check_col < columnCycles.size(); ++check_col) {
                bool correct = MATRECNetworkDecompositionVerifyCycle(dec, check_col,
                                                                   columnCycles[check_col].data(),
                                                                   columnCycleValues[check_col].data(),
                                                                   columnCycles[check_col].size(),
                                                                   rowStorage.data(),signStorage);
                EXPECT_TRUE(correct);
                if (!correct) {
                    fundamental_cycles_good = false;
                }
            }
            if (!fundamental_cycles_good) {
                isGood = false;
                break;
            }
        }
        if(expectNetwork.has_value()){
            EXPECT_EQ(expectNetwork.value(),isNetwork);
        }else{
            bool isGraphic = false;
            bool isSigned = false;
            if(isNetwork){
                MATREC_CALL(runRowGraphicCheck(env,testCase,isGraphic));
                EXPECT_TRUE(isGraphic);
                MATREC_CALL(runRowSignedCheck(env,testCase,isSigned));
                EXPECT_TRUE(isSigned);
            }else{
                MATREC_CALL(runRowGraphicCheck(env,testCase,isGraphic));
                if(isGraphic){
                    MATREC_CALL(runRowSignedCheck(env,testCase,isSigned));
                    EXPECT_FALSE(isSigned);
                }
            }
        }

    }
    MATRECfreeNetworkRowAddition(env, &newRow);
    MATRECNetworkDecompositionFree(&dec);
    MATREC_CALL(MATRECfreeEnvironment(&env));
    return MATREC_OKAY;
}
void runRowTestCase(const DirectedTestCase& testCase, std::optional<bool> expectNetwork = std::nullopt){
    bool isGood = true;

    MATREC_ERROR error = runRowNetworkDecomposition(testCase, isGood, expectNetwork);
    EXPECT_EQ(error, MATREC_OKAY);
    EXPECT_TRUE(isGood);
}

MATREC_ERROR runInterleavedNetworkDecomposition(const DirectedTestCase& testCase,
                                      bool& isGood, std::size_t seed,
                                      std::optional<bool> expectNetwork = std::nullopt){
    MATREC *env = NULL;
    MATREC_CALL(MATRECcreateEnvironment(&env));
    bool * signStorage = NULL;
    MATREC_CALL(MATRECallocBlockArray(env,&signStorage,testCase.rows));
    std::vector<MATREC_row> rowStorage(testCase.rows,MATREC_INVALID_ROW);
    MATRECNetworkDecomposition *dec = NULL;
    MATREC_CALL(MATRECNetworkDecompositionCreate(env, &dec, testCase.rows, testCase.cols));
    MATRECNetworkRowAddition *newRow = NULL;
    MATREC_CALL(MATRECcreateNetworkRowAddition(env, &newRow));
    MATRECNetworkColumnAddition *newCol = NULL;
    MATREC_CALL(MATRECcreateNetworkColumnAddition(env,&newCol));
    std::vector<std::vector<MATREC_row>> columnCycles(testCase.cols,std::vector<MATREC_col>());
    std::vector<std::vector<double>> columnCycleValues(testCase.cols,std::vector<double>());

    const DirectedColTestCase colTestCase(testCase);
    {
        std::vector<int> ordering;
        for (int i = 0; i < testCase.rows; ++i) {
            ordering.push_back(i);
        }
        for (int i = 0; i < testCase.cols; ++i) {
            ordering.push_back(-i - 1);
        }

        auto rng = std::mt19937(seed);
        std::shuffle(ordering.begin(),ordering.end(),rng);

        std::vector<MATREC_col> effectiveCols;
        std::vector<MATREC_row> effectiveRows;
        std::vector<double> effectiveValues;
        bool isNetwork = true;
        for (int orderIndex: ordering) {
            if (orderIndex < 0) {
                int column = -orderIndex - 1;
                effectiveRows.clear();
                effectiveValues.clear();
                for (auto pair: colTestCase.matrix[column]) {
                    if (MATRECNetworkDecompositionContainsRow(dec, pair.index)) {
                        effectiveRows.push_back(pair.index);
                        effectiveValues.push_back(pair.value);
                    }
                }
                MATREC_CALL(MATRECNetworkColumnAdditionCheck(dec, newCol, column, effectiveRows.data(),
                                                         effectiveValues.data(),effectiveRows.size()));
                if (MATRECNetworkColumnAdditionRemainsNetwork(newCol)) {
                    MATREC_CALL(MATRECNetworkColumnAdditionAdd(dec, newCol));
                } else {
                    isNetwork = false;
                    break;
                }
                assert(columnCycles[column].empty());
                for (MATREC_row row: effectiveRows) {
                    columnCycles[column].push_back(row);
                }
                for(double value : effectiveValues){
                    columnCycleValues[column].push_back(value);
                }
            } else {
                int row = orderIndex;
                effectiveCols.clear();
                effectiveValues.clear();
                for (auto pair: testCase.matrix[row]) {
                    if (MATRECNetworkDecompositionContainsColumn(dec, pair.index)) {
                        effectiveCols.push_back(pair.index);
                        effectiveValues.push_back(pair.value);
                    }
                }

                MATREC_CALL(MATRECNetworkRowAdditionCheck(dec, newRow, row, effectiveCols.data(),
                                                      effectiveValues.data(),effectiveCols.size()));
                if (MATRECNetworkRowAdditionRemainsNetwork(newRow)) {
                    MATREC_CALL(MATRECNetworkRowAdditionAdd(dec, newRow));
                } else {
                    isNetwork = false;
                    break;
                }
                for(std::size_t i = 0; i < effectiveCols.size(); ++i){
                    MATREC_col col = effectiveCols[i];
                    columnCycles[col].push_back(row);
                    columnCycleValues[col].push_back(effectiveValues[i]);
                }
            }
            bool isMinimal = MATRECNetworkDecompositionIsMinimal(dec);
            EXPECT_TRUE(isMinimal); //Check that there are no series-series and bond-bond connections

            bool fundamental_cycles_good = true;
            for (std::size_t check_col = 0; check_col < columnCycles.size(); ++check_col) {
                bool correct = MATRECNetworkDecompositionVerifyCycle(dec, check_col,
                                                                   columnCycles[check_col].data(),
                                                                   columnCycleValues[check_col].data(),
                                                                   columnCycles[check_col].size(),
                                                                   rowStorage.data(),signStorage);
                EXPECT_TRUE(correct);
                if (!correct) {
                    fundamental_cycles_good = false;
                }
            }
            if (!fundamental_cycles_good) {
                isGood = false;
                break;
            }
        }
        if(expectNetwork.has_value()){
            EXPECT_EQ(expectNetwork.value(),isNetwork);
        }else{
            bool isGraphic = false;
            bool isSigned = false;
            if(isNetwork){
                MATREC_CALL(runRowGraphicCheck(env,testCase,isGraphic));
                EXPECT_TRUE(isGraphic);
                MATREC_CALL(runRowSignedCheck(env,testCase,isSigned));
                EXPECT_TRUE(isSigned);
            }else{
                MATREC_CALL(runRowGraphicCheck(env,testCase,isGraphic));
                if(isGraphic){
                    MATREC_CALL(runRowSignedCheck(env,testCase,isSigned));
                    EXPECT_FALSE(isSigned);
                }
            }
        }

    }
    MATRECfreeNetworkColumnAddition(env,&newCol);
    MATRECfreeNetworkRowAddition(env, &newRow);
    MATRECNetworkDecompositionFree(&dec);
    MATREC_CALL(MATRECfreeEnvironment(&env));
    return MATREC_OKAY;
}


void runInterleavedTestCase(const DirectedTestCase& testCase, std::size_t seed,
                            std::optional<bool> expectNetwork =  std::nullopt){
    bool isGood = true;

    MATREC_ERROR error = runInterleavedNetworkDecomposition(testCase, isGood,seed, expectNetwork);
    EXPECT_EQ(error, MATREC_OKAY);
    EXPECT_TRUE(isGood);
}
namespace NetworkAdditionTest{

    TEST(NetworkColAddition, SingleColumn){
        auto testCase = stringToDirectedTestCase(
                "+1 "
                "+1 "
                "-1 ",
                3,1);

        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, DoubleColumnInvalidSign){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                "+1  0 "
                "-1 +1 ",
                3,2);

        runTestCase(testCase,false);
    }
    TEST(NetworkColAddition, DoubleColumnInvalidSign2){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                "+1  0 "
                "-1 -1 ",
                3,2);

        runTestCase(testCase,false);

    }

    //Create and convert both false;
    TEST(NetworkColAddition, SplitSeries1){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                "+1  0 "
                " 0  0 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries2){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                "+1  0 "
                " 0  0 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries1r){
        auto testCase = stringToDirectedTestCase(
                "-1 -1 "
                "+1  0 "
                " 0  0 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries2r){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 "
                "+1  0 "
                " 0  0 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries3){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                "+1  0 "
                " 0 +1 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries4){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                "+1  0 "
                " 0 +1 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries3r){
        auto testCase = stringToDirectedTestCase(
                "-1 -1 "
                "+1  0 "
                " 0 +1 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries4r){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 "
                "+1  0 "
                " 0 +1 ",
                3,2);
        runTestCase(testCase,true);
    }

//    //Create = false, convert = true
    TEST(NetworkColAddition, SplitSeries5){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                " 0  0 "
                " 0  0 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries6){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                " 0  0 "
                " 0  0 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries7){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                " 0  0 "
                " 0  +1 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries8){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                " 0  0 "
                " 0  +1 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries5r){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 "
                " 0  0 "
                " 0  0 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries6r){
        auto testCase = stringToDirectedTestCase(
                "-1 -1 "
                " 0  0 "
                " 0  0 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries7r){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 "
                " 0  0 "
                " 0  +1 ",
                3,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries8r){
        auto testCase = stringToDirectedTestCase(
                "-1 -1 "
                " 0  0 "
                " 0  +1 ",
                3,2);
        runTestCase(testCase,true);
    }

    //Create = true, convert = false
    TEST(NetworkColAddition, SplitSeries9){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                "+1  0 "
                "-1  +1 "
                " 0  0 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries10){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                "+1  0 "
                "-1  -1 "
                " 0  0 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries11){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                "+1  0 "
                "-1  -1 "
                " 0  +1 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries12){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                "+1  0 "
                "-1  +1 "
                " 0  +1 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries9r){
        auto testCase = stringToDirectedTestCase(
                "-1 -1 "
                "-1  0 "
                "+1  +1 "
                " 0  0 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries10r){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 "
                "-1  0 "
                "+1  -1 "
                " 0  0 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries11r){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 "
                "-1  0 "
                "+1  -1 "
                " 0  +1 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries12r){
        auto testCase = stringToDirectedTestCase(
                "-1 -1 "
                "-1  0 "
                "+1  +1 "
                " 0  +1 ",
                4,2);
        runTestCase(testCase,true);
    }

    //Create == convert == true
    TEST(NetworkColAddition, SplitSeries13){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                "+1 +1 "
                "-1 -1 "
                " 0  0 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries14){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                "+1 -1 "
                "-1 +1 "
                " 0  0 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries15){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 "
                "+1 +1 "
                "-1 -1 "
                " 0  +1 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries16){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                "+1 -1 "
                "-1 +1 "
                " 0  +1 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries13r){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 "
                "-1 +1 "
                "+1 -1 "
                " 0  0 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries14r){
        auto testCase = stringToDirectedTestCase(
                "-1 -1 "
                "-1 -1 "
                "+1 +1 "
                " 0  0 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries15r){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 "
                "-1 +1 "
                "+1 -1 "
                " 0  +1 ",
                4,2);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition, SplitSeries16r){
        auto testCase = stringToDirectedTestCase(
                "-1 -1 "
                "-1 -1 "
                "+1 +1 "
                " 0  +1 ",
                4,2);
        runTestCase(testCase,true);
    }

    TEST(NetworkColAddition,ParallelSimple){
        auto testCase = stringToDirectedTestCase(
                "1 1 1 -1 ",
                1,4);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition,ParallelSimple2){
        auto testCase = stringToDirectedTestCase(
                "1 1 1 -1 "
                "0 0 -1 0 ",
                2,4);
        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition,ParallelSimple3){
        auto testCase = stringToDirectedTestCase(
                "1 1 1 1 "
                "0 0 1 0 ",
                2,4);
        runTestCase(testCase,true);
    }

    TEST(NetworkColAddition,Components1){
        auto testCase = stringToDirectedTestCase(
                "0 -1 -1 "
                "1 0 1 ",
                2,3);

        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition,Components2){
        auto testCase = stringToDirectedTestCase(
                "0 1 -1 "
                "1  0 1 ",
                2,3);

        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition,Components3){
        auto testCase = stringToDirectedTestCase(
                "0 1  1 "
                "-1  0 1 ",
                2,3);

        runTestCase(testCase,true);
    }
    TEST(NetworkColAddition,ThreeByThree1){
        auto testCase = seedToDirectedTestCase(1,3,3);

        runTestCase(testCase,false);
    }
    TEST(NetworkColAddition,ThreeByThree2){
        auto testCase = seedToDirectedTestCase(2,3,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree3){
        auto testCase = seedToDirectedTestCase(3,3,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree4){
        auto testCase = seedToDirectedTestCase(8,3,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree5){
        auto testCase = seedToDirectedTestCase(14,3,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree6){
        auto testCase = seedToDirectedTestCase(41,3,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree7){
        auto testCase = seedToDirectedTestCase(59,3,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree8){
        auto testCase = seedToDirectedTestCase(29,3,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree9){
        auto testCase = seedToDirectedTestCase(43,3,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree10){
        auto testCase = seedToDirectedTestCase(90,3,3);
        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByThree11){
        auto testCase = seedToDirectedTestCase(110,3,3);
        runTestCase(testCase);
    }
    TEST(NetworkColAddition,SixByThree1){
        auto testCase = seedToDirectedTestCase(24,6,3);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,SixByThree2){
        auto testCase = seedToDirectedTestCase(8,6,3);

        runTestCase(testCase);
    }

    TEST(NetworkColAddition,ThreeByFour1){
        auto testCase = seedToDirectedTestCase(3,3,4);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,ThreeByFive1){
        auto testCase = seedToDirectedTestCase(14,3,5);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FourByEight1){
        auto testCase = seedToDirectedTestCase(0,4,8);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FourByEight2){
        auto testCase = seedToDirectedTestCase(30,4,8);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FourByEight3){
        auto testCase = seedToDirectedTestCase(65,4,8);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FourByEight4){
        auto testCase = seedToDirectedTestCase(243,4,8);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FourByEight5){
        auto testCase = seedToDirectedTestCase(381,4,8);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FourByEight6){
        auto testCase = seedToDirectedTestCase(3804,4,8);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FourByEight7){
        auto testCase = seedToDirectedTestCase(27,4,8);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FourByFour1){
        auto testCase = stringToDirectedTestCase(
                "-1 +1 -1 0 "
                "-1 0 -1 0 "
                "0 0 -1 +1 "
                "-1 +1 0 -1",
                4,4);
        runTestCase(testCase);

    }
    TEST(NetworkColAddition,FiveByFive1){
        auto testCase = seedToDirectedTestCase(326,5,5);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FiveByFive2){
        auto testCase = seedToDirectedTestCase(2099,5,5);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FiveByFive3){
        auto testCase = seedToDirectedTestCase(6083,5,5);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FiveByFive4){
        auto testCase = seedToDirectedTestCase(6649,5,5);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FiveByFive5){
        auto testCase = seedToDirectedTestCase(84293,5,5);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FiveByFive6){
        auto testCase = seedToDirectedTestCase(1332,5,5);

        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FiveByFive7){
        auto testCase = seedToDirectedTestCase(3474,5,5);
        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FiveByFive8){
        auto testCase = seedToDirectedTestCase(201847,5,5);
        runTestCase(testCase);
    }
    TEST(NetworkColAddition,FiveByFive9){
        auto testCase = seedToDirectedTestCase(481421,5,5);
        runTestCase(testCase);
    }
    TEST(NetworkColAddition,Random6By2){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,6,2);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random3By3){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,3,3);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random6By3){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,6,3);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random2By6){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,2,6);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random3By5){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,3,5);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random3By8){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,3,8);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random4By8){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,4,8);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random5By5){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,5,5);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random8By8){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,8,8);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,Random4By10){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,4,10);

            runTestCase(testCase);
        }
    }
    TEST(NetworkColAddition,ER10D50){
        for(std::size_t seed = 0; seed < 100; ++seed){
            auto testCase = erdosRenyiDirectedTestCase(10,0.25,seed);
            runTestCase(testCase,true);
        }
    }
    TEST(NetworkColAddition,ER20D50){
        for(std::size_t seed = 0; seed < 100; ++seed){
            auto testCase = erdosRenyiDirectedTestCase(20,0.5,seed);
            runTestCase(testCase,true);
        }
    }
    TEST(NetworkColAddition,ER100D3){
        for(std::size_t seed = 0; seed < 100; ++seed){
            auto testCase = erdosRenyiDirectedTestCase(100,0.03,seed);
            runTestCase(testCase,true);
        }
    }

    TEST(NetworkRowAddition,OneByTwo1){
        auto testCase = stringToDirectedTestCase(
                "+1 0 ",
                1,2);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,OneByTwo2){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 ",
                1,2);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,OneByTwo3){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 ",
                1,2);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,TwoByThree1){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 +1 "
                "-1 +1 -1 ",
                2,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,TwoByThree2){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 +1 "
                "-1 +1 +1 ",
                2,3);
        runRowTestCase(testCase,false);
    }
    TEST(NetworkRowAddition,TwoByThree3){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 +1 "
                "+1 0 +1 ",
                2,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,TwoByThree4){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 0 "
                "+1 0  0 ",
                2,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,TwoByThree5){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 0 "
                "+0 +1  0 ",
                2,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,TwoByThree6){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 0 "
                "+0 -1  0 ",
                2,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,TwoByThree7){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 0 "
                "+0 +1 +1 ",
                2,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,TwoByThree8){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 0 "
                "+0 -1 +1 ",
                2,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,ThreeBySix1){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 0 0 0 0 "
                "0 0 +1 -1 0 0 "
                "-1 +1 -1 0 0 0 ",
                3,6);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,ThreeBySix2){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 0 0 0 0 "
                "0 0 +1 -1 0 0 "
                "-1 +1 -1 0 0 +1 ",
                3,6);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,ThreeByTwo1){
        auto testCase = stringToDirectedTestCase(
                "+1 -1 "
                "-1 +1 "
                "+1 -1 ",
                3,2);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,ThreeByOne1){
        auto testCase = stringToDirectedTestCase(
                "+1 "
                "-1 "
                "+1 ",
                3,1);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,ThreeByThree1){
        auto testCase = seedToDirectedTestCase(8,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeByThree2){
        auto testCase = seedToDirectedTestCase(14,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeByThree3){
        auto testCase = seedToDirectedTestCase(16,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeByThree4){
        auto testCase = seedToDirectedTestCase(19,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeByThree5){
        auto testCase = seedToDirectedTestCase(52,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeByThree6){
        auto testCase = seedToDirectedTestCase(67,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeByThree7){
        auto testCase = seedToDirectedTestCase(77,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeByThree8){
        auto testCase = seedToDirectedTestCase(87,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeByThree9){
        auto testCase = seedToDirectedTestCase(135,3,3);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeBySix3){
        auto testCase = seedToDirectedTestCase(73,3,6);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,ThreeBySix4){
        auto testCase = seedToDirectedTestCase(195,3,6);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour1){
        auto testCase = seedToDirectedTestCase(25,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour2){
        auto testCase = seedToDirectedTestCase(40,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour3){
        auto testCase = seedToDirectedTestCase(146,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour4){
        auto testCase = seedToDirectedTestCase(175,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour5){
        auto testCase = seedToDirectedTestCase(210,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour6){
        auto testCase = seedToDirectedTestCase(160,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour7){
        auto testCase = seedToDirectedTestCase(160,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour8){
        auto testCase = seedToDirectedTestCase(74,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour9){
        auto testCase = seedToDirectedTestCase(37,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour10){
        auto testCase = seedToDirectedTestCase(217,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour11){
        auto testCase = seedToDirectedTestCase(460,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour12){
        auto testCase = seedToDirectedTestCase(473,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour13){
        auto testCase = seedToDirectedTestCase(490,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour14){
        auto testCase = seedToDirectedTestCase(735,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour15){
        auto testCase = seedToDirectedTestCase(929,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour16){
        auto testCase = seedToDirectedTestCase(1300,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour17){
        auto testCase = seedToDirectedTestCase(1400,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FourByFour18){
        auto testCase = seedToDirectedTestCase(2884,4,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FiveByFive1){
        auto testCase = seedToDirectedTestCase(8,5,5);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FiveByFive2){
        auto testCase = seedToDirectedTestCase(7070,5,5);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FiveByFive3){
        auto testCase = seedToDirectedTestCase(14940,5,5);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FiveByFive4){
        auto testCase = seedToDirectedTestCase(39220,5,5);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,FiveByFive5){
        auto testCase = seedToDirectedTestCase(47629,5,5);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,EightByFour1){
        auto testCase = seedToDirectedTestCase(1722,8,4);

        runRowTestCase(testCase);
    }
    TEST(NetworkRowAddition,SingleRigid1){
        auto testCase = stringToDirectedTestCase(
                "+1 0 +1 "
                "+1 +1 0 "
                "0 -1 +1 "
                "+1 +1 0 ",
                4,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,SingleRigid2){
        auto testCase = stringToDirectedTestCase(
                "+1 0 +1 "
                "+1 +1 0 "
                "0 -1 +1 "
                "-1 -1 0 ",
                4,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,SingleRigid3){
        auto testCase = stringToDirectedTestCase(
                "+1 0 +1 "
                "+1 +1 0 "
                "0 -1 +1 "
                "-1 +1 0 ",
                4,3);
        runRowTestCase(testCase,false);
    }
    TEST(NetworkRowAddition,SingleRigid4){
        auto testCase = stringToDirectedTestCase(
                "+1 0 +1 "
                "-1 -1 -1 "
                "0 +1 +1 "
                "+1 +1 +1 ",
                4,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,SingleRigid5){
        auto testCase = stringToDirectedTestCase(
                "+1 0 +1 "
                "-1 -1 -1 "
                "0 +1 +1 "
                "-1 -1 -1 ",
                4,3);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,SingleRigid6){
        auto testCase = stringToDirectedTestCase(
                "+1 0 +1 "
                "-1 -1 -1 "
                "0 +1 +1 "
                "-1 +1 -1 ",
                4,3);
        runRowTestCase(testCase,false);
    }
    TEST(NetworkRowAddition,SingleRigid7){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 0 0 +1 "
                "+1 0 +1 0 0 "
                "0 -1 +1 +1 -1 "
                "0 0 0 -1 +1 "
                "+1 +1 0 0 0 ",
                5,5);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,SingleRigid8){
        auto testCase = stringToDirectedTestCase(
                "+1 +1 0 0 +1 "
                "+1 0 +1 0 0 "
                "0 -1 +1 +1 -1 "
                "0 0 0 -1 +1 "
                "+1 +1 0 0 0 "
                "+1 0 +1 +1 0 ",
                6,5);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,RandomTwoBySix){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,2,6);

            runRowTestCase(testCase);
        }
    }
    TEST(NetworkRowAddition,RandomThreeByThree){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,3,3);

            runRowTestCase(testCase);
        }
    }
    TEST(NetworkRowAddition,RandomThreeBySix){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,3,6);

            runRowTestCase(testCase);
        }
    }
    TEST(NetworkRowAddition,RandomFourByFour){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,4,4);

            runRowTestCase(testCase);
        }
    }

    TEST(NetworkRowAddition,RandomFiveByFive){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,5,5);

            runRowTestCase(testCase);
        }
    }
    TEST(NetworkRowAddition,RandomEightByEight){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,8,8);

            runRowTestCase(testCase);
        }
    }
    TEST(NetworkRowAddition,ER10D50){
        for(std::size_t seed = 0; seed < 100; ++seed){
            auto testCase = erdosRenyiDirectedTestCase(10,0.25,seed);
            runRowTestCase(testCase,true);
        }
    }
    TEST(NetworkRowAddition,ER10D50Case1){
        auto testCase = erdosRenyiDirectedTestCase(10,0.25,2733);
        runRowTestCase(testCase,true);
    }
    TEST(NetworkRowAddition,ER20D50){
        for(std::size_t seed = 0; seed < 100; ++seed){
            auto testCase = erdosRenyiDirectedTestCase(20,0.5,seed);
            runRowTestCase(testCase,true);
        }
    }
    TEST(NetworkRowAddition,ER100D3){
        for(std::size_t seed = 0; seed < 100; ++seed){
            auto testCase = erdosRenyiDirectedTestCase(100,0.03,seed);
            runRowTestCase(testCase,true);
        }
    }
    TEST(NetworkInterleavedAddition,FourByFour){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,4,4);

            runInterleavedTestCase(testCase,i);
        }
    }
    TEST(NetworkInterleavedAddition,FourByFour1){
        std::size_t i = 2;
        auto testCase = seedToDirectedTestCase(i,4,4);
        runInterleavedTestCase(testCase,i);
    }
    TEST(NetworkInterleavedAddition,FourByFour2){
        std::size_t i = 13;
        auto testCase = seedToDirectedTestCase(i,4,4);
        runInterleavedTestCase(testCase,i);
    }
    TEST(NetworkInterleavedAddition,FourByFour3){
        std::size_t i = 10;
        auto testCase = seedToDirectedTestCase(i,4,4);
        runInterleavedTestCase(testCase,i);
    }
    TEST(NetworkInterleavedAddition,FourByFour4){
        std::size_t i = 5;
        auto testCase = seedToDirectedTestCase(i,4,4);
        runInterleavedTestCase(testCase,i);
    }
    TEST(NetworkInterleavedAddition,FourByFour5){
        std::size_t i = 28;
        auto testCase = seedToDirectedTestCase(i,4,4);
        runInterleavedTestCase(testCase,i);
    }
    TEST(NetworkInterleavedAddition,FourByFour6){
        std::size_t i = 3026;
        auto testCase = seedToDirectedTestCase(i,4,4);
        runInterleavedTestCase(testCase,i);
    }
    TEST(NetworkInterleavedAddition,FiveByFive){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,5,5);

            runInterleavedTestCase(testCase,i);
        }
    }
    TEST(NetworkInterleavedAddition,EightByEight){
        for(std::size_t i = 0; i < 10'000; ++i){
            auto testCase = seedToDirectedTestCase(i,8,8);

            runInterleavedTestCase(testCase,i);
        }
    }
    TEST(NetworkInterleavedAddition,ER10D50){
        for(std::size_t seed = 0; seed < 100; ++seed){
            auto testCase = erdosRenyiDirectedTestCase(10,0.25,seed);
            runInterleavedTestCase(testCase,seed,true);
        }
    }
}