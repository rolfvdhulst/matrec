#include "TestHelpers.h"
#include <gtest/gtest.h>
#include <cmr/graphic.h>
#include <matrec/Graphic.h>

namespace GraphicTest {
    CMR_ERROR cmrTestGraphicness(bool &result, const TestCase &testCase) {
        CMR *cmr = NULL;

        CMR_CALL(CMRcreateEnvironment(&cmr));
        CMR_CHRMAT *matrix = NULL;
        std::string filename = "temptestfile.sparse";
        {
            std::ofstream writeStream(filename, std::ofstream::out);
            std::size_t nonzeros = std::accumulate(testCase.matrix.begin(), testCase.matrix.end(),
                                                   0ul, [](std::size_t b,
                                                           const std::vector<MATREC_col> &rowData) -> std::size_t {
                        return b + rowData.size();
                    });
            writeStream << testCase.rows << " " << testCase.cols << " " << nonzeros << "\n";

            for (std::size_t row = 0; row < testCase.rows; ++row) {
                for (const auto &col_idx: testCase.matrix[row]) {
                    writeStream << row + 1 << " " << col_idx + 1 << " 1\n";

                }
            }
        }

        FILE *readFile = fopen(filename.c_str(), "r");
        CMR_CALL(CMRchrmatCreateFromSparseStream(cmr, readFile, &matrix));
        fclose(readFile);
        unlink(filename.c_str());//Destroy the file again.

        CMR_GRAPHIC_STATISTICS stats;
        CMR_CALL(CMRstatsGraphicInit(&stats));
        bool isGraphic;
        CMR_CALL(CMRtestGraphicMatrix(cmr, matrix, &isGraphic, NULL, NULL, NULL, NULL, &stats));

        result = isGraphic;

        CMR_CALL(CMRchrmatFree(cmr, &matrix));
        CMR_CALL(CMRfreeEnvironment(&cmr));
        return CMR_OKAY;
    }

    bool CMRisGraphic(const TestCase &testCase) {
        bool result;
        CMR_ERROR code = cmrTestGraphicness(result, testCase);
        EXPECT_EQ(code, CMR_OKAY);
        return result;
    }

//TODO: should probably write this using a googletest fixture to prevent reallocating all the time
    MATREC_ERROR runDecomposition(const TestCase &testCase, bool &isGood, bool guaranteedGraphic) {
        MATREC *env = NULL;
        MATREC_CALL(MATRECcreateEnvironment(&env));
        MATRECGraphicDecomposition *dec = NULL;
        MATREC_CALL(MATRECGraphicDecompositionCreate(env, &dec, testCase.rows, testCase.cols));
        MATRECGraphicRowAddition *newRow = NULL;
        MATREC_CALL(MATRECcreateGraphicRowAddition(env, &newRow));
        MATRECGraphicColumnAddition *newCol = NULL;
        MATREC_CALL(MATRECcreateGraphicColumnAddition(env, &newCol));

        std::vector<std::vector<MATREC_row>> column_cycles(testCase.cols, std::vector<MATREC_row>());

        {
            std::vector<int> ordering;
            for (int i = 0; i < testCase.rows; ++i) {
                ordering.push_back(i);
            }
            for (int i = 0; i < testCase.cols; ++i) {
                ordering.push_back(-i - 1);
            }

            auto rng = std::mt19937(testCase.seed);
            std::shuffle(ordering.begin(),ordering.end(),rng);

            std::vector<MATREC_row> column_storage(testCase.rows, -1);// A buffer to hold the computed cycles in

            std::vector<MATREC_col> effectiveCols;
            std::vector<MATREC_row> effectiveRows;
            const ColTestCase colTestCase(testCase);

            bool isGraphic = true;
            for (int orderIndex: ordering) {
                if (orderIndex < 0) {
                    int column = -orderIndex - 1;
                    effectiveRows.clear();
                    for (MATREC_row row: colTestCase.matrix[column]) {
                        if (MATRECGraphicDecompositionContainsRow(dec, row)) {
                            effectiveRows.push_back(row);
                        }
                    }
                    MATREC_CALL(MATRECGraphicColumnAdditionCheck(dec, newCol, column, effectiveRows.data(),
                                                             effectiveRows.size()));
                    if (MATRECGraphicColumnAdditionRemainsGraphic(newCol)) {
                        MATREC_CALL(MATRECGraphicColumnAdditionAdd(dec, newCol));
                    } else {
                        isGraphic = false;
                        break;
                    }
                    assert(column_cycles[column].empty());
                    for (MATREC_row row: effectiveRows) {
                        column_cycles[column].push_back(row);
                    }
                } else {
                    int row = orderIndex;
                    effectiveCols.clear();
                    for (MATREC_col col: testCase.matrix[row]) {
                        if (MATRECGraphicDecompositionContainsColumn(dec, col)) {
                            effectiveCols.push_back(col);
                        }
                    }

                    MATREC_CALL(MATRECGraphicRowAdditionCheck(dec, newRow, row, effectiveCols.data(), effectiveCols.size()));
                    if (MATRECGraphicRowAdditionRemainsGraphic(newRow)) {
                        MATREC_CALL(MATRECGraphicRowAdditionAdd(dec, newRow));
                    } else {
                        isGraphic = false;
                        break;
                    }
                    for (MATREC_col col: effectiveCols) {
                        column_cycles[col].push_back(row);
                    }
                }

                bool isMinimal = MATRECGraphicDecompositionIsMinimal(dec);
                EXPECT_TRUE(isMinimal); //Check that there are no series-series and bond-bond connections
                if (!isMinimal) {
                    std::cout << "seed: " << testCase.seed << "\n";
                }


                bool fundamental_cycles_good = true;
                for (std::size_t column = 0; column < column_cycles.size(); ++column) {
                    bool correct = MATRECGraphicDecompositionVerifyCycle(dec, column, column_cycles[column].data(),
                                                                       column_cycles[column].size(),
                                                                       column_storage.data());
                    EXPECT_TRUE(correct);
                    if (!correct) {
                        fundamental_cycles_good = false;
                        std::cout << "seed: " << testCase.seed << "\n";
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
                    std::cout << "seed: " << testCase.seed << "\n";
                }
            }

        }
        MATRECfreeGraphicColumnAddition(env, &newCol);
        MATRECfreeGraphicRowAddition(env, &newRow);
        MATRECGraphicDecompositionFree(&dec);
        MATREC_CALL(MATRECfreeEnvironment(&env));
        return MATREC_OKAY;
    }

    void runTestCase(const TestCase &testCase, bool guaranteedGraphic = false) {
        ASSERT_EQ(testCase.rows, testCase.matrix.size());
        bool isGood = true;

        MATREC_ERROR error = runDecomposition(testCase, isGood, guaranteedGraphic);
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
        for (uint64_t i = 0; i < max; ++i) {
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

    TEST(GraphicInterleaved, AllFourByFour)

    {
        runAllMatrices(4, 4);
    }

    TEST(GraphicInterleaved,FourByFour1){
        TestCase testCase = seedToTestCase(819,4,4);
        runTestCase(testCase);
    }
    TEST(GraphicInterleaved,FourByFour2){
        TestCase testCase = seedToTestCase(51,4,4);
        runTestCase(testCase);
    }
    TEST(GraphicInterleaved,FourByFour3){
        TestCase testCase = seedToTestCase(55,4,4);
        runTestCase(testCase);
    }
    TEST(GraphicInterleaved,FourByFour4){
        TestCase testCase = seedToTestCase(827,4,4);
        runTestCase(testCase);
    }
    TEST(GraphicInterleaved,FourByFour5){
        TestCase testCase = seedToTestCase(1019,4,4);
        runTestCase(testCase);
    }
    TEST(GraphicInterleaved,FourByFour6){
        TestCase testCase = seedToTestCase(15290,4,4);
        runTestCase(testCase);
    }
    TEST(GraphicInterleaved,FourByFour7){
        TestCase testCase = seedToTestCase(28637,4,4);
        runTestCase(testCase);
    }
    TEST(GraphicInterleaved, AllFiveByFive)
    {
        runAllMatrices(5,5);
    }
    TEST(GraphicInterleaved, FiveByFive1)
    {
        TestCase testCase = seedToTestCase(129535,5,5);
        runTestCase(testCase);
    }

    TEST(GraphicInterleaved, FiveByFive2)
    {
        TestCase testCase = seedToTestCase(3325839,5,5);
        runTestCase(testCase);
    }
    TEST(GraphicInterleaved, SevenBySeven)
    {
        randomlySample(7,7,100'000,19);
    }
}