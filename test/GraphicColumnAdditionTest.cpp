#include "TestHelpers.h"
#include <matrec/Graphic.h>
#include <cmr/graphic.h>
#include <gtest/gtest.h>
namespace ColAdditionTest {
    CMR_ERROR cmrTestGraphicness(bool &result, const ColTestCase &testCase) {
        CMR *cmr = NULL;

        CMR_CALL(CMRcreateEnvironment(&cmr));
        CMR_CHRMAT *matrix = NULL;
        std::string filename = "temptestfile.sparse";
        {
            std::ofstream writeStream(filename, std::ofstream::out);
            std::size_t nonzeros = std::accumulate(testCase.matrix.begin(), testCase.matrix.end(),
                                                   0ul, [](std::size_t b,
                                                           const std::vector<MATREC_row> &colData) -> std::size_t {
                        return b + colData.size();
                    });
            writeStream << testCase.rows << " " << testCase.cols << " " << nonzeros << "\n";

            //TODO: what type of format does CMR expect, row or column wise?
            for (std::size_t col_idx = 0; col_idx < testCase.cols; ++col_idx) {
                for (const auto &row_idx: testCase.matrix[col_idx]) {
                    writeStream << row_idx + 1 << " " << col_idx + 1 << " 1\n";
                }
            }
        }

        FILE *readFile = fopen(filename.c_str(), "r");
        CMR_CALL(CMRchrmatCreateFromSparseStream(cmr, readFile, &matrix));
        fclose(readFile);
        unlink(filename.c_str());//Destroy the file again.

        CMR_GRAPHIC_STATISTICS stats;
        CMR_CALL(CMRgraphicStatsInit(&stats));
        bool isGraphic;
        CMR_CALL(CMRgraphicTestMatrix(cmr, matrix, &isGraphic, NULL, NULL, NULL, NULL, &stats, std::numeric_limits<double>::infinity()));

        result = isGraphic;

        CMR_CALL(CMRchrmatFree(cmr, &matrix));
        CMR_CALL(CMRfreeEnvironment(&cmr));
        return CMR_OKAY;
    }

    bool CMRisGraphic(const ColTestCase &testCase) {
        bool result;
        CMR_ERROR code = cmrTestGraphicness(result, testCase);
        EXPECT_EQ(code, CMR_OKAY);
        return result;
    }

//TODO: should probably write this using a googletest fixture to prevent reallocating all the time
    MATREC_ERROR runDecomposition(ColTestCase &testCase, std::vector<std::vector<MATREC_col>> &column_cycles, bool &isGood,
                                bool guaranteedGraphic) {
        MATREC *env = NULL;
        MATREC_CALL(MATRECcreateEnvironment(&env));
        MATRECGraphicDecomposition *dec = NULL;
        MATREC_CALL(MATRECGraphicDecompositionCreate(env, &dec, testCase.rows, testCase.cols));
        MATRECGraphicColumnAddition *newCol = NULL;
        MATREC_CALL(MATRECcreateGraphicColumnAddition(env, &newCol));
        {
            std::vector<MATREC_row> column_storage(testCase.rows, -1);// A buffer to hold the computed cycles in

            bool isGraphic = true;
            for (std::size_t col = 0; col < testCase.cols; ++col) {
                MATREC_CALL(MATRECGraphicColumnAdditionCheck(dec, newCol, col, testCase.matrix[col].data(),
                                                         testCase.matrix[col].size()));
                if (MATRECGraphicColumnAdditionRemainsGraphic(newCol)) {
                    MATREC_CALL(MATRECGraphicColumnAdditionAdd(dec, newCol));
                } else {
                    isGraphic = false;
                    break;
                }

                bool isMinimal = MATRECGraphicDecompositionIsMinimal(dec);
                EXPECT_TRUE(isMinimal); //Check that there are no series-series and bond-bond connections
                if (!isMinimal) {
                    std::cout << "Seed: " << testCase.seed << " is not minimal!\n";
                }
                bool fundamental_cycles_good = true;
                for (std::size_t check_col = 0; check_col <= col; ++check_col) {
                    bool correct = MATRECGraphicDecompositionVerifyCycle(dec, check_col, testCase.matrix[check_col].data(),
                                                          testCase.matrix[check_col].size(), column_storage.data());
                    EXPECT_TRUE(correct);
                    if (!correct) {
                        fundamental_cycles_good = false;
                        std::cout << "seed: " << testCase.seed << " has mismatching fundamental cycle for " << check_col
                                  << "\n";
                    }
                }
                if (!fundamental_cycles_good) {
                    isGood = false;
                    break;
                }
            }

            if (!isGraphic) {
                if (guaranteedGraphic) {
                    EXPECT_FALSE(true);
                }
                bool cmrSaysGraphic = CMRisGraphic(testCase);
                EXPECT_FALSE(cmrSaysGraphic);
                if (cmrSaysGraphic) {
                    std::cout << "seed: " << testCase.seed << " is graphic according to CMR!\n";
                }
            }
        }
        MATRECfreeGraphicColumnAddition(env, &newCol);
        MATRECGraphicDecompositionFree(&dec);
        MATREC_CALL(MATRECfreeEnvironment(&env));
        return MATREC_OKAY;
    }

    void runTestCase(const TestCase &testCase, bool guaranteedGraphic = false) {
        ASSERT_EQ(testCase.rows, testCase.matrix.size());
        ColTestCase colTestCase(testCase);

        std::vector<std::vector<MATREC_col>> column_cycles(testCase.cols);
        bool isGood = true;

        MATREC_ERROR error = runDecomposition(colTestCase, column_cycles, isGood, guaranteedGraphic);
        EXPECT_EQ(error, MATREC_OKAY);
        EXPECT_TRUE(isGood);
    }

    void runAllMatrices(uint64_t rows, uint64_t cols) {
        uint64_t numElems = rows * cols;
        if (numElems >= 64ul) {
            std::cout << "cannot run on matrices this large...\n";
            EXPECT_TRUE(false);
            return;
        }

        uint64_t max = 1ul << numElems;
        float minPercent = 0.05;
        for (uint64_t i = 0; i < max; ++i) { //TODO: change back to zero
            TestCase testCase = seedToTestCase(i, rows, cols);
            runTestCase(testCase);

            float percentage = float(i) / float(max);
            if (percentage > minPercent) {
                std::cout << std::round(minPercent * 100.0) << "% complete" << std::endl;
                minPercent += 0.05;
            }
        }
    }

    void randomlySample(uint64_t rows, uint64_t cols, int numRepetitions, uint64_t firstSeed) {
        uint64_t numElems = rows * cols;
        if (numElems >= 64ul) {
            std::cout << "cannot run on matrices this large...\n";
            EXPECT_TRUE(false);
            return;
        }
        std::cout << "running with first seed: " << firstSeed << "\n";
        uint64_t max = 1ul << numElems;

        auto dist = std::uniform_int_distribution<uint64_t>(0, max);
        std::mt19937_64 gen(firstSeed);

        float minPercent = 0.05;
        for (int i = 0; i < numRepetitions; ++i) {
            uint64_t num = dist(gen);
            TestCase testCase = seedToTestCase(num, rows, cols);
            runTestCase(testCase);

            float percentage = float(i) / float(numRepetitions);
            if (percentage > minPercent) {
                std::cout << std::round(minPercent * 100.0) << "% complete" << std::endl;
                minPercent += 0.05;
            }
        }
    }

    TEST(ColAddition, SingleColumn){
        TestCase testCase({{0},{0},{},{0}},4,1);
        runTestCase(testCase);
    }
    TEST(ColAddition, SingleCycle1){
        TestCase testCase({{0,1},{0},{0}},3,2);
        runTestCase(testCase);
    }
    TEST(ColAddition, SingleCycle2){
        TestCase testCase({{0,1},{0,1},{0}},3,2);
        runTestCase(testCase);
    }
    TEST(ColAddition, SingleCycle3){
        TestCase testCase({{0,1},{0,1},{0,1}},3,2);
        runTestCase(testCase);
    }
    TEST(ColAddition,SimpleColEmptyCol){
        TestCase testCase({{0,1},{},{},{},{},{}},6,2);
        runTestCase(testCase);
    }
    TEST(ColAddition,test){
        TestCase testCase = seedToTestCase(3, 2, 6);
        runTestCase(testCase);
    }
    TEST(ColAddition,ReducedRootMovingProblem){
        TestCase testCase = seedToTestCase(391, 2, 6);
        runTestCase(testCase);
    }
    TEST(ColAddition,NonMinimalCase1){
        TestCase testCase = seedToTestCase(455,2,6);
        runTestCase(testCase);
    }
    TEST(ColAddition,NonMinimalCase2){
        TestCase testCase = seedToTestCase(911,2,6);
        runTestCase(testCase);
    }
    TEST(ColAddition,InfiniteLoop){
        TestCase testCase = seedToTestCase(219,2,6);
        runTestCase(testCase);
    }
    TEST(ColAddition,test1){
        TestCase testCase = seedToTestCase(119,3,3);
        runTestCase(testCase);
    }
    TEST(ColAddition,misplacedEdge){
        TestCase testCase = seedToTestCase(238,3,3);
        runTestCase(testCase);
    }
    TEST(ColAddition, MultipleComponents1){
        TestCase testCase({{1,2},
                           {0},
                           {0,1}}, 3, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition, Triconnected1){
        TestCase testCase({{0,1},
                           {1, 2},
                           {0,2}}, 3, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition, Triconnected2){
        TestCase testCase({{0,1},
                           {1, 2},
                           {0,1,2}}, 3, 3);
        runTestCase(testCase);
    }

    TEST(ColAddition, Triconnected3){
        TestCase testCase({{0,2},
                           {0,1,2},
                           {0,1}}, 3, 3);
        runTestCase(testCase);
    }

    TEST(ColAddition, PropagatingToOne1) {
        TestCase testCase({{0,1},
                           {1},
                           {0, 1}}, 3, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition, PropagatingToOne2) {
        TestCase testCase({{1},
                           {0, 1},
                           {0, 1}}, 3, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition, PropagatingToOne3) {
        TestCase testCase({{0,1},
                           {0},
                           {0, 1}}, 3, 3);
        runTestCase(testCase);
    }

    TEST(ColAddition,TriconnectedSplit3){
        TestCase testCase({{1,2},
                           {0, 2},
                           {0, 1},
                           {0,1}}, 4, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedSplit4){
        TestCase testCase({{0,2},
                           {0,1,2},
                           {0, 2},
                           {0,1}}, 4, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedSplit5){
        TestCase testCase({{0,2},
                           {1, 2},
                           {0, 1},
                           {0,1}}, 4, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedArticulationNode1){
        TestCase testCase({{1,2},
                           {0, 2},
                           {0, 1},
                           {0,1,2}}, 4, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedArticulationNode2){
        TestCase testCase({{0,1,2},
                           {0, 2},
                           {0, 1},
                           {0,1,2}}, 4, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedArticulationNode3){
        TestCase testCase({{0,1,2},
                           {1, 2},
                           {0, 1},
                           {0,1,2}}, 4, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedMerging2){
        TestCase testCase({{0,2},
                           {0,1},
                           {0,1,2},
                           {0,1,2}}, 4, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedSplit1){
        TestCase testCase({{1,2},
                           {0, 2},
                           {0, 1},
                           {0}}, 4, 3);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedSplit2){
        TestCase testCase({{0,2},
                           {0, 1},
                           {1, 2},
                           {0}}, 4, 3);
        runTestCase(testCase);
    }


    TEST(ColAddition,ThreeByFour1){
        TestCase testCase = seedToTestCase(862, 3, 4);
        runTestCase(testCase);
    }
    TEST(ColAddition,ThreeByFour2){
        TestCase testCase = seedToTestCase(1964, 3, 4);
        runTestCase(testCase);
    }


    TEST(ColAddition,FourByFourInfiniteLoop){
        TestCase testCase = seedToTestCase(13722, 4, 4);
        runTestCase(testCase);
    }
    TEST(ColAddition,FourByFour1){
        TestCase testCase = seedToTestCase(13740, 4, 4);
        runTestCase(testCase);
    }
    TEST(ColAddition,FourByFour2){
        TestCase testCase = seedToTestCase(13743, 4, 4);
        runTestCase(testCase);
    }
    TEST(ColAddition,FourByFour3){
        TestCase testCase = seedToTestCase(13996, 4, 4);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedArticulationNode4){
        TestCase testCase({{0,1,3},
                           {1, 2},
                           {0, 2},
                           {0,1,}}, 4, 4);
        runTestCase(testCase);
    }
    TEST(ColAddition,TriconnectedMerging1){
        TestCase testCase({{1,2,3},
                           {0,3},
                           {0,2},
                           {0,1}}, 4, 4);
        runTestCase(testCase);
    }

    TEST(ColAddition,Triconnected4){
        TestCase testCase({{1,2,3},
                           {0, 2},
                           {0,1}}, 3, 5);
        runTestCase(testCase);
    }
    TEST(ColAddition,Triconnected5){
        TestCase testCase({{1,2,3,4},
                           {0, 3},
                           {0,1,2}}, 3, 5);
        runTestCase(testCase);
    }
    TEST(ColAddition,Triconnected6){
        TestCase testCase({{0,1,3},
                           {0, 1,2},
                           {0,1,2,3}}, 3, 5);
        runTestCase(testCase);
    }


    TEST(ColAddition,AllSixByTwo){
        runAllMatrices(6,2);
    }
    TEST(ColAddition,AllTwoBySix){
        //This should cover most, if not all, cases with two rows.
        //Note that these matrices are always graphic by definition
        runAllMatrices(2,6);
    }
    TEST(ColAddition,AllThreeByThree){
        runAllMatrices(3,3);
    }
    TEST(ColAddition,AllSixByThree){
        runAllMatrices(6,3);
    }
    TEST(ColAddition,AllThreeByFour){
        runAllMatrices(3,4);
    }
    TEST(ColAddition,AllFourByFour){
        runAllMatrices(4,4);
    }
//    TEST(ColAddition,AllSixByFour){
//        runAllMatrices(6,4);
//    }
//    TEST(ColAddition,AllThreeBySeven){
//        runAllMatrices(3,7);
//    }
//    TEST(CollAddition,AllFourBySix){
//        runAllMatrices(4,6);
//    }
//    TEST(ColAddition,AllFiveByFive){
//        runAllMatrices(5,5);
//    }

    TEST(ColAddition,NotMinimal){
        runTestCase(seedToTestCase(14041,5,3));
    }
    TEST(ColAddition,CutEdgeNoMemory){
        runTestCase(seedToTestCase(29653,5,3));
    }
    TEST(ColAddition,reorderingBug){
        runTestCase(seedToTestCase(1652,4,3));
    }
    TEST(ColAddition,RigidNoCuts){
        runTestCase(seedToTestCase(111290,4,5));
    }
    TEST(ColAddition,WrongSplit){
        runTestCase(seedToTestCase(239554,4,5));
    }
    TEST(ColAddition,CMRBug){
        runTestCase(seedToTestCase(241337,4,5));
    }
    TEST(ColAddition,ArticulationEdgeRigidIntoParent){
        runTestCase(seedToTestCase(509643,4,5));
    }
    TEST(ColAddition,RigidMergingArticulation){
        runTestCase(seedToTestCase(1019318,4,5));
    }
    TEST(ColAddition,RigidMergingArticulation2){
        runTestCase(seedToTestCase(1019319,4,5));
    }
    TEST(ColAddition,FiveByFourNotGraphic1){
        runTestCase(seedToTestCase(206238,5,4));
    }
    TEST(ColAddition,FiveByFourCase2){
        runTestCase(seedToTestCase(214474,5,4));
    }
    TEST(ColAddition,FiveByFourIncorrectMerge){
        runTestCase(seedToTestCase(472493,5,4));
    }
    TEST(ColAddition,FiveByFourBadColoring){
        //case where coloring for articulation node is not bipartite
        runTestCase(seedToTestCase(472508,5,4));
    }

    TEST(ColAddition,FiveByFourUnknown){
        runTestCase(seedToTestCase(472523,5,4));
    }
    TEST(ColAddition,ArticulationPointSearchBug1){
        runTestCase(seedToTestCase(996783,5,4));
    }
    TEST(ColAddition,ArticulationPointSearchBug2){
        runTestCase(seedToTestCase(996797,5,4));
    }
    TEST(ColAddition,SixBySix1){
        runTestCase(seedToTestCase(47091484421,6,6));
    }
    TEST(ColAddition,SixBySixThreeArticulationPoints){
        runTestCase(seedToTestCase(16273736387,6,6));
    }
    TEST(ColAddition,SixByFourThreeArticulationPoints){
        runTestCase(seedToTestCase(7559087,6,4));
    }
    TEST(ColAddition,FiveByFiveOne){
        runTestCase(seedToTestCase(31804393,5,5));
    }
    TEST(ColAddition,FiveByFiveParallelSplitting){
        runTestCase(seedToTestCase(33223995,5,5));
    }
    TEST(ColAddition,FiveByFiveWrong){
        runTestCase(seedToTestCase(13616375,5,5));
    }
    TEST(ColAddition,FiveByFiveTwo){
        runTestCase(seedToTestCase(3510,5,5));
    }
    TEST(ColAddition,FiveByFiveThree){
        runTestCase(seedToTestCase(112282,5,5));
    }
    TEST(ColAddition,FiveByFiveFour){
        runTestCase(seedToTestCase(3324504,5,5));
    }
    TEST(ColAddition,FiveByFiveFive){
        runTestCase(seedToTestCase(3324702,5,5));
    }
    TEST(ColAddition,FiveByFiveSix){
        runTestCase(seedToTestCase(7645847,5,5));
    }
    TEST(ColAddition,FourBySixOne){
        runTestCase(seedToTestCase(824988,4,6));
    }
    TEST(ColAddition,cmrOne){
        auto testCase = stringToTestCase(    "1 0 1 1 "
                                             "1 1 0 0 "
                                             "0 1 1 1 ",
                                             3,4);
        runTestCase(testCase);
    }
    TEST(ColAddition,cmrTwo){
        auto testCase = stringToTestCase(    "1 0 0 1 1 "
                                             "1 1 0 0 1 "
                                             "0 1 1 0 1 "
                                             "0 1 0 1 0 "
                                             "0 0 0 1 1 ",
                                             5,5);
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr3){
        auto testCase = stringToTestCase(    "1 1 0 1 0 "
                                             "1 0 0 0 1 "
                                             "1 0 1 1 0 "
                                             "0 0 0 1 1 "
                                             "0 1 1 0 0 ",
                                             5,5);
        runTestCase(testCase);
    }

    TEST(ColAddition,cmr4){
        auto testCase = stringToTestCase(            "1 0 1 1 1 0 1 0 0 0 "
                                                     "0 0 1 0 0 1 1 0 0 0 "
                                                     "0 1 0 1 0 0 1 1 0 0 "
                                                     "0 0 0 0 1 0 1 1 0 0 "
                                                     "0 0 1 0 0 0 1 0 1 1 ",
                                                     5, 10
                                             );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr5){
        auto testCase = stringToTestCase(    "0 0 0 1 1 1 1 "
                                             "0 1 1 1 1 0 0 "
                                             "0 0 0 1 1 0 1 "
                                             "0 1 1 1 0 0 1 "
                                             "0 1 0 0 0 1 1 "
                                             ,5,7
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr6){
        auto testCase = stringToTestCase(                "0 1 1 1 0 1 "
                                                         "1 1 1 1 1 1 "
                                                         "0 1 0 0 1 0 "
                                                         "0 0 1 0 1 1 "
                                                         "1 1 1 0 1 0 ",
                                                     5, 6
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr7){
        auto testCase = stringToTestCase(    "1 1 0 0 0 1 0 "
                                             "1 1 0 1 1 0 0 "
                                             "1 0 0 1 0 0 0 "
                                             "0 0 0 0 1 1 0 "
                                             "1 1 0 1 1 1 0 ",
                                                     5, 7
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr8){
        auto testCase = stringToTestCase("0 1 1 0 1 0 "
                                         "0 1 1 1 1 1 "
                                         "1 0 0 1 0 0 "
                                         "0 0 1 1 0 0 "
                                         "1 0 1 1 0 0 "
                                         "1 1 0 1 0 0 ",
                                                     6,6
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr9){
        auto testCase = stringToTestCase(    "1 1 1 0 0 0 1 1 "
                                             "0 0 1 1 1 1 0 0 "
                                             "0 0 0 0 0 1 1 0 "
                                             "0 1 0 1 0 0 1 0 "
                                             "0 0 0 0 1 0 1 0 ",
                                                     5, 8
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr10){
        auto testCase = stringToTestCase(        "0 0 0 0 0 1 "
                                                 "1 1 0 1 1 0 "
                                                 "1 1 1 0 0 0 "
                                                 "0 1 0 1 1 0 "
                                                 "0 1 1 0 1 0 ",
                                             5, 6
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr11){
        auto testCase = stringToTestCase(    "0 0 1 0 1 1 "
                                             "0 1 1 0 0 1 "
                                             "0 0 0 1 0 0 "
                                             "0 1 1 1 0 0 "
                                             "1 0 0 1 1 1 ",
                                             5, 6
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr12){
        auto testCase = stringToTestCase(    "0 0 0 1 0 1 "
                                             "1 1 1 0 1 1 "
                                             "1 0 0 1 1 1 "
                                             "1 1 0 0 0 1 "
                                             "0 0 0 1 1 0 ",
                                             5, 6
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr13){
        auto testCase = stringToTestCase(    "0 1 1 1 1 1 "
                                             "0 0 0 0 1 1 "
                                             "0 1 1 0 1 1 "
                                             "0 1 0 0 0 1 "
                                             "1 0 1 1 0 1 ",
                                             5, 6
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr14){
        auto testCase = stringToTestCase("1 1 1 1 1 1 "
                                         "0 0 1 0 0 0 "
                                         "1 1 1 1 1 0 "
                                         "0 1 0 1 1 1 "
                                         "0 0 0 1 0 1 "
                                         "1 0 0 0 1 1 "
                                         "0 1 0 0 0 1 ",
                                             7,6
        );
        runTestCase(testCase);
    }
    TEST(ColAddition,cmr15){
        auto testCase = stringToTestCase("1",1,1);
        runTestCase(testCase);
    }


    TEST(ColAddition,SixBy10_1){
        runTestCase(seedToTestCase(871520243141229702,6,10));
    }
    TEST(ColAddition,SixBySix2){
        runTestCase(seedToTestCase(28585879706,6,6));
    }
    TEST(ColAddition,SixBySix3){
        runTestCase(seedToTestCase(28662816944,6,6));
    }
    TEST(ColAddition,SixBySix4){
        runTestCase(seedToTestCase(55897751056,6,6));
    }
    TEST(ColAddition,SixBySix5){
        runTestCase(seedToTestCase(50242711234,6,6));
    }
    TEST(ColAddition,SevenByEight1){
        runTestCase(seedToTestCase(48221311036474554,7,8));
    }
    TEST(ColAddition,RandomSixBySix){
        randomlySample(6,6,100'000,42);
    }
    TEST(ColAddition,Random7By8){
        randomlySample(7,8,100'000,3);
    }
    TEST(ColAddition,Random8By7){
        randomlySample(8,7,100'000,16);
    }
    TEST(ColAddition,Random10By6){
        randomlySample(10,6,100'000,20);
    }
    TEST(ColAddition,Random6By10){
        randomlySample(6,10,100'000,38);
    }
}
