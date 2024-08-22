#include "matrec/Incidence.h"
#include "matrec/Network.h"
#include <gtest/gtest.h>


MATREC_ERROR stringToMatrix(MATREC * env,
                            MATRECCSMatrixInt ** rowMatrix,
                            std::string string){
    FILE * file = fmemopen(string.data(),string.size(),"r");
    MATREC_CALL(MATRECreadIntMatrixFromStream(env, rowMatrix, file));
    fclose(file);
    return MATREC_OKAY;
}

MATREC_ERROR runIncidenceAlgorithms(MATREC * env, MATRECCSMatrixInt * rowMatrix, MATRECCSMatrixInt * transpose,
                                    std::optional<bool> mustBeResult){

    bool rowWiseIncidence = true;
    {
        MATRECIncidenceAddition *addition;
        MATREC_CALL(MATRECcreateIncidenceAddition(env, &addition, rowMatrix->numRows, rowMatrix->numColumns, MATREC_INIT_NONE));
        for (MATREC_matrix_size i = 0; i < rowMatrix->numColumns; ++i) {
            bool added = MATRECincidenceAdditionAddColumn(addition, i, MATRECintMatrixRowNumNonzeros(transpose, i),
                                                          MATRECintMatrixRowColumnIndices(transpose, i),
                                                          MATRECintMatrixRowColumnValues(transpose, i));
            EXPECT_TRUE(added);
        }
        for (MATREC_matrix_size i = 0; i < rowMatrix->numRows; ++i) {
            bool added = MATRECincidenceAdditionAddRow(addition, i, MATRECintMatrixRowNumNonzeros(rowMatrix, i),
                                                       MATRECintMatrixRowColumnIndices(rowMatrix, i),
                                                       MATRECintMatrixRowColumnValues(rowMatrix, i));
            if(!added){
                rowWiseIncidence = false;
                break;
            }
        }
        for (int i = 0; i < rowMatrix->numColumns; ++i) {
            int sum = 0;
            MATREC_matrix_size nnonz = MATRECintMatrixRowNumNonzeros(transpose, i);
            MATREC_row * rows = MATRECintMatrixRowColumnIndices(transpose, i);
            int * values = MATRECintMatrixRowColumnValues(transpose, i);
            for (MATREC_matrix_size j = 0; j < nnonz; ++j) {
                if(!MATRECincidenceContainsRow(addition, rows[j])) continue;
                sum += MATRECincidenceRowSign(addition, rows[j]) * values[j];
            }
            EXPECT_LE(abs(sum),1);
        }

        MATRECfreeIncidenceAddition(env, &addition);
    }
    bool colWiseIncidence = true;
    {
        MATRECIncidenceAddition *addition;
        MATREC_CALL(MATRECcreateIncidenceAddition(env, &addition, rowMatrix->numRows, rowMatrix->numColumns, MATREC_INIT_NONE));
        for (MATREC_matrix_size i = 0; i < rowMatrix->numRows; ++i) {
            bool added = MATRECincidenceAdditionAddRow(addition, i, MATRECintMatrixRowNumNonzeros(rowMatrix, i),
                                                       MATRECintMatrixRowColumnIndices(rowMatrix, i),
                                                       MATRECintMatrixRowColumnValues(rowMatrix, i));

            EXPECT_TRUE(added);
        }

        for (MATREC_matrix_size i = 0; i < rowMatrix->numColumns; ++i) {
            bool added = MATRECincidenceAdditionAddColumn(addition, i, MATRECintMatrixRowNumNonzeros(transpose, i),
                                                          MATRECintMatrixRowColumnIndices(transpose, i),
                                                          MATRECintMatrixRowColumnValues(transpose, i));
            if(!added){
                colWiseIncidence = false;
                break;
            }
        }
        for (int i = 0; i < rowMatrix->numColumns; ++i) {
            if(!MATRECincidenceContainsColumn(addition, i)) continue;
            int sum = 0;
            MATREC_matrix_size nnonz = MATRECintMatrixRowNumNonzeros(transpose, i);
            MATREC_row * rows = MATRECintMatrixRowColumnIndices(transpose, i);
            int * values = MATRECintMatrixRowColumnValues(transpose, i);
            for (MATREC_matrix_size j = 0; j < nnonz; ++j) {
                sum += MATRECincidenceRowSign(addition, rows[j]) * values[j];
            }
            EXPECT_LE(abs(sum),1);
        }
        MATRECfreeIncidenceAddition(env, &addition);
    }
    EXPECT_EQ(rowWiseIncidence,colWiseIncidence);
    if(rowWiseIncidence || colWiseIncidence){
        //TODO: fix
//        bool isNetwork = true;
//        MATRECNetworkDecomposition * dec = NULL;
//        MATREC_CALL(MATRECNetworkDecompositionCreate(env,&dec,rowMatrix->numRows,rowMatrix->numColumns));
//
//        MATRECNetworkRowAddition * nra = NULL;
//        MATREC_CALL(MATRECcreateNetworkRowAddition(env,&nra));
//
//        for (int i = 0; i < rowMatrix->numRows; ++i) {
//            MATREC_CALL(MATRECNetworkRowAdditionCheck(dec,nra,i,
//                                                      intMatr(rowMatrix,i),));
//            if(!MATRECNetworkRowAdditionRemainsNetwork(nra)){
//                isNetwork = false;
//                break;
//            }
//            MATREC_CALL(MATRECNetworkRowAdditionAdd(dec,nra));
//        }
//        MATRECfreeNetworkRowAddition(env,&nra);
//        MATRECNetworkDecompositionFree(&dec);
//        EXPECT_TRUE(isNetwork);
//        for (int i = 0; i < rowMatrix->numColumns; ++i) {
//            EXPECT_LE(intMatrixRowNumNonzeros(transpose,i),2);
//        }
    }
    if(mustBeResult.has_value()){
        EXPECT_EQ(mustBeResult.value(),rowWiseIncidence);
    }

    return MATREC_OKAY;
}
MATREC_ERROR runTestCase(std::string string,std::optional<bool> guaranteedResult = std::nullopt){
    MATREC * env;
    MATREC_CALL(MATRECcreateEnvironment(&env));
    MATRECCSMatrixInt * matrix;
    MATRECCSMatrixInt * transpose;
    MATREC_CALL(stringToMatrix(env,&matrix,string));
    MATREC_CALL(MATRECtransposeIntMatrix(env, matrix, &transpose));
    MATREC_CALL(runIncidenceAlgorithms(env,matrix,transpose,guaranteedResult));
    MATRECfreeIntMatrix(env, &matrix);
    MATRECfreeIntMatrix(env, &transpose);
    MATREC_CALL(MATRECfreeEnvironment(&env));
    return MATREC_OKAY;
}
TEST(IncidenceAddition,simple){
    std::string testCase = "3 3 5\n"
                           "1 1 1\n"
                           "1 2 1\n"
                           "2 2 1\n"
                           "2 3 -1\n"
                           "3 3 1";
    MATREC_ERROR error = runTestCase(testCase,true);
    EXPECT_EQ(error, MATREC_OKAY);
}

TEST(IncidenceAddition,simple2){
    std::string testCase = "3 3 6\n"
                           "1 1 1\n"
                           "1 2 1\n"
                           "2 2 1\n"
                           "2 3 -1\n"
                           "3 1 1\n"
                           "3 3 1";
    MATREC_ERROR error = runTestCase(testCase,true);
    EXPECT_EQ(error, MATREC_OKAY);
}
TEST(IncidenceAddition,simple3){
    std::string testCase = "3 3 6\n"
                           "1 1 1\n"
                           "1 2 1\n"
                           "2 2 1\n"
                           "2 3 -1\n"
                           "3 1 1\n"
                           "3 3 -1";
    MATREC_ERROR error = runTestCase(testCase,false);
    EXPECT_EQ(error, MATREC_OKAY);
}

TEST(IncidenceAddition,badColumn){
    std::string testCase = "3 3 6\n"
                           "1 1 1\n"
                           "1 2 1\n"
                           "1 3 -1\n"
                           "2 2 1\n"
                           "2 3 -1\n"
                           "3 3 1";
    MATREC_ERROR error = runTestCase(testCase,false);
    EXPECT_EQ(error, MATREC_OKAY);
}
TEST(IncidenceAddition,multipleComponents){
    std::string testCase ="7 9 18\n"
                          "1 1 1\n"
                          "1 2 1\n"
                          "2 2 -1\n"
                          "2 3 -1\n"
                          "3 4 1\n"
                          "3 5 1\n"
                          "4 5 1\n"
                          "4 6 -1\n"
                          "5 7 1\n"
                          "5 9 1\n"
                          "6 7 1\n"
                          "6 8 -1\n"
                          "6 9 1\n"
                          "7 1 -1\n"
                          "7 3 1\n"
                          "7 4 1\n"
                          "7 6 1\n"
                          "7 8 -1";
    MATREC_ERROR error = runTestCase(testCase,true);
    EXPECT_EQ(error, MATREC_OKAY);

}

TEST(IncidenceAddition,multipleComponents2){
    std::string testCase ="7 9 18\n"
                          "1 1 1\n"
                          "1 2 1\n"
                          "2 2 -1\n"
                          "2 3 -1\n"
                          "3 4 1\n"
                          "3 5 1\n"
                          "4 5 1\n"
                          "4 6 -1\n"
                          "5 7 1\n"
                          "5 9 1\n"
                          "6 7 1\n"
                          "6 8 -1\n"
                          "6 9 1\n"
                          "7 1 -1\n"
                          "7 3 1\n"
                          "7 4 -1\n"
                          "7 6 -1\n"
                          "7 8 -1";
    MATREC_ERROR error = runTestCase(testCase,true);
    EXPECT_EQ(error, MATREC_OKAY);

}
TEST(IncidenceAddition,multipleComponents3){
    std::string testCase ="4 5 10\n"
                          "1 1 1\n"
                          "1 2 -1\n"
                          "1 5 1\n"
                          "2 1 -1\n"
                          "2 2 1\n"
                          "3 3 -1\n"
                          "3 4 1\n"
                          "4 3 -1\n"
                          "4 4 1\n"
                          "4 5 1\n";
    MATREC_ERROR error = runTestCase(testCase,true);
    EXPECT_EQ(error, MATREC_OKAY);

}