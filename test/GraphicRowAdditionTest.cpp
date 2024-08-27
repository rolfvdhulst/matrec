#include "TestHelpers.h"
#include <gtest/gtest.h>
#include <cmr/graphic.h>
#include <matrec/Graphic.h>

CMR_ERROR cmrTestGraphicness(bool& result, const TestCase& testCase){
    CMR * cmr = NULL;

    CMR_CALL(CMRcreateEnvironment(&cmr));
    CMR_CHRMAT * matrix = NULL;
    std::string filename = "temptestfile.sparse";
    {
        std::ofstream writeStream(filename,std::ofstream::out);
        std::size_t nonzeros = std::accumulate(testCase.matrix.begin(),testCase.matrix.end(),
                                               0ul,[](std::size_t b,const std::vector<MATREC_col>& rowData) -> std::size_t{
                    return b+rowData.size();
                });
        writeStream << testCase.rows <<" "<<testCase.cols <<" "<<nonzeros<<"\n";

        for (std::size_t row = 0; row < testCase.rows; ++row) {
            for(const auto& col_idx : testCase.matrix[row]){
                writeStream << row+1<<" "<<col_idx+1<<" 1\n";

            }
        }
    }

    FILE * readFile = fopen(filename.c_str() ,"r");
    CMR_CALL(CMRchrmatCreateFromSparseStream(cmr, readFile, &matrix));
    fclose(readFile);
    unlink(filename.c_str());//Destroy the file again.

    CMR_GRAPHIC_STATISTICS stats;
    CMR_CALL(CMRgraphicStatsInit(&stats));
    bool isGraphic;
    CMR_CALL(CMRgraphicTestMatrix(cmr, matrix, &isGraphic, NULL, NULL, NULL, NULL, &stats, std::numeric_limits<double>::infinity()));

    result = isGraphic;

    CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    CMR_CALL(CMRfreeEnvironment(&cmr));
    return CMR_OKAY;
}

bool CMRisGraphic(const TestCase& testCase){
    bool result;
    CMR_ERROR code = cmrTestGraphicness(result,testCase);
    EXPECT_EQ(code,CMR_OKAY);
    return result;
}

//TODO: should probably write this using a googletest fixture to prevent reallocating all the time
MATREC_ERROR runDecomposition(const TestCase &testCase, std::vector<std::vector<MATREC_col>> &column_cycles, bool &isGood, bool guaranteedGraphic) {
    MATREC *env = NULL;
    MATREC_CALL(MATRECcreateEnvironment(&env));
    MATRECGraphicDecomposition *dec = NULL;
    MATREC_CALL(MATRECGraphicDecompositionCreate(env, &dec, testCase.rows, testCase.cols));
    MATRECGraphicRowAddition *newRow = NULL;
    MATREC_CALL(MATRECcreateGraphicRowAddition(env, &newRow));
    {
        std::vector<MATREC_row> column_storage(testCase.rows, -1);// A buffer to hold the computed cycles in
        bool isGraphic = true;
        for (std::size_t row = 0; row < testCase.rows; ++row) {
            MATREC_CALL(MATRECGraphicRowAdditionCheck(dec, newRow, row, testCase.matrix[row].data(), testCase.matrix[row].size()));
            if (MATRECGraphicRowAdditionRemainsGraphic(newRow)) {
                MATREC_CALL(MATRECGraphicRowAdditionAdd(dec, newRow));
            } else {
                isGraphic = false;
                break;
            }

            for (const auto &col_index : testCase.matrix[row]) {
                column_cycles[col_index].push_back(row);
            }
            bool isMinimal = MATRECGraphicDecompositionIsMinimal(dec);
            EXPECT_TRUE(isMinimal); //Check that there are no series-series and bond-bond connections
            if(!isMinimal){
                std::cout<<"seed: "<<testCase.seed<<"\n";
            }
            bool fundamental_cycles_good = true;
            for (std::size_t column = 0; column < column_cycles.size(); ++column) {
                bool correct = MATRECGraphicDecompositionVerifyCycle(dec, column, column_cycles[column].data(),
                                                      column_cycles[column].size(), column_storage.data());
                EXPECT_TRUE(correct);
                if (!correct) {
                    fundamental_cycles_good = false;
                    std::cout<<"seed: "<<testCase.seed<<"\n";
                }
            }
            if (!fundamental_cycles_good) {
                isGood = false;
                break;
            }
        }

        if(!isGraphic){
            if(guaranteedGraphic){
                EXPECT_FALSE(true);
            }
            bool cmrSaysGraphic = CMRisGraphic(testCase);
            EXPECT_FALSE(cmrSaysGraphic);
            if(cmrSaysGraphic){
                std::cout<<"seed: "<<testCase.seed<<"\n";
            }
        }
    }
    MATRECfreeGraphicRowAddition(env, &newRow);
    MATRECGraphicDecompositionFree(&dec);
    MATREC_CALL(MATRECfreeEnvironment(&env));
    return MATREC_OKAY;
}

void runTestCase(const TestCase &testCase, bool guaranteedGraphic = false) {
    ASSERT_EQ(testCase.rows, testCase.matrix.size());
    std::vector<std::vector<MATREC_col>> column_cycles(testCase.cols);
    bool isGood = true;

    MATREC_ERROR error = runDecomposition(testCase, column_cycles, isGood, guaranteedGraphic);
    EXPECT_EQ(error, MATREC_OKAY);
    EXPECT_TRUE(isGood);
}
void runAllMatrices(uint64_t rows, uint64_t cols){
    uint64_t numElems= rows*cols;
    if(numElems >= 64ul){
        std::cout<<"cannot run on matrices this large...\n";
        EXPECT_TRUE(false);
        return;
    }

    uint64_t max = 1ul << numElems;
    float minPercent = 0.05;
    for (uint64_t i = 0; i < max; ++i) {
        TestCase testCase = seedToTestCase(i,rows,cols);
        runTestCase(testCase);

        float percentage = float(i) / float(max);
        if(percentage > minPercent){
            std::cout<<std::round(minPercent*100.0)<< "% complete"<<std::endl;
            minPercent+=0.05;
        }
    }
}

void randomlySample(uint64_t rows, uint64_t cols, int numRepetitions, uint64_t firstSeed){
    uint64_t numElems= rows*cols;
    if(numElems >= 64ul){
        std::cout<<"cannot run on matrices this large...\n";
        EXPECT_TRUE(false);
        return;
    }
    std::cout<<"running with first seed: "<<firstSeed<<"\n";
    uint64_t max = 1ul << numElems;

    auto dist = std::uniform_int_distribution<uint64_t>(0,max);
    std::mt19937_64 gen(firstSeed);

    float minPercent = 0.05;
    for (int i = 0; i < numRepetitions; ++i) {
        uint64_t num = dist(gen);
        TestCase testCase = seedToTestCase(num,rows,cols);
        runTestCase(testCase);

        float percentage = float(i) / float(numRepetitions);
        if(percentage > minPercent){
            std::cout<<std::round(minPercent*100.0)<< "% complete"<<std::endl;
            minPercent+=0.05;
        }
    }
}
TEST(RowAddition, SingleRow){
    TestCase testCase({{0,1}}, 1, 3);
    runTestCase(testCase);
}
TEST(RowAddition,PS1){
    TestCase testCase({{0},{0}}, 2, 3);
    runTestCase(testCase);
}
TEST(RowAddition,PS2){
    TestCase testCase({{0,1},{0,1}}, 2, 3);
    runTestCase(testCase);
}
TEST(RowAddition,PS3){
    TestCase testCase({{0,1,2},{0,1}}, 2, 3);
    runTestCase(testCase);
}
TEST(RowAddition,PS4){
    TestCase testCase({{0,1},{0}}, 2, 3);
    runTestCase(testCase);
}
TEST(RowAddition,PS5){
    TestCase testCase({{0,1},{0,1},{0,1}}, 3, 3);
    runTestCase(testCase);
}

TEST(RowAddition,largeButSimple){
    //A large case that only uses single series and parallel expansions
    TestCase testCase({
        {2,4,5,6},
        {2,3,4,5,6},
        {2,4,5,6},
        {4,5},
        {6,7,8},
        {0,1,2},
        {0},
        {2},
        {2,9}},9,10);
    runTestCase(testCase);
}

TEST(RowAddition,simpleMultipleComponents){
    TestCase testCase({{0},{1},{0,1}},3,2);
    runTestCase(testCase);
}

TEST(RowAddition,ZeroFirstRow){
    TestCase testCase({{},
                       {0},
                       {0}}, 3, 3);
    runTestCase(testCase);
}
TEST(RowAddition, Triconnected1){
    TestCase testCase({{0,1},
                       {1, 2},
                       {0,2}}, 3, 3);
    runTestCase(testCase);
}
TEST(RowAddition, Triconnected2){
    TestCase testCase({{0,1},
                       {1, 2},
                       {0,1,2}}, 3, 3);
    runTestCase(testCase);
}

TEST(RowAddition, Triconnected3){
    TestCase testCase({{0,2},
                       {0,1,2},
                       {0,1}}, 3, 3);
    runTestCase(testCase);
}
TEST(RowAddition,Triconnected4){
    TestCase testCase({{1,2,3},
                       {0, 2},
                       {0,1}}, 3, 5);
    runTestCase(testCase);
}
TEST(RowAddition,Triconnected5){
    TestCase testCase({{1,2,3,4},
                       {0, 3},
                       {0,1,2}}, 3, 5);
    runTestCase(testCase);
}
TEST(RowAddition,Triconnected6){
    TestCase testCase({{0,1,3},
                       {0, 1,2},
                       {0,1,2,3}}, 3, 5);
    runTestCase(testCase);
}
TEST(RowAddition, MultipleComponents1){
    TestCase testCase({{1,2},
                       {0},
                       {0,1}}, 3, 3);
    runTestCase(testCase);
}
TEST(RowAddition, PropagatingToOne1) {
    //Simple case where series should be propagated to a parallel
    TestCase testCase({{0,1},
                       {1},
                       {0, 1}}, 3, 3);
    runTestCase(testCase);
}
TEST(RowAddition, PropagatingToOne2) {
    //Simple case where series should be propagated to a parallel
    TestCase testCase({{1},
                       {0, 1},
                       {0, 1}}, 3, 3);
    runTestCase(testCase);
}
TEST(RowAddition, PropagatingToOne3) {
    //Simple case where series should be propagated to a parallel
    TestCase testCase({{0,1},
                       {0},
                       {0, 1}}, 3, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedSplit1){
    TestCase testCase({{1,2},
                       {0, 2},
                       {0, 1},
                       {0}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedSplit2){
    TestCase testCase({{0,2},
                       {0, 1},
                       {1, 2},
                       {0}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedSplit3){
    TestCase testCase({{1,2},
                       {0, 2},
                       {0, 1},
                       {0,1}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedSplit4){
    TestCase testCase({{0,2},
                       {0,1,2},
                       {0, 2},
                       {0,1}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedSplit5){
    TestCase testCase({{0,2},
                       {1, 2},
                       {0, 1},
                       {0,1}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedArticulationNode1){
    TestCase testCase({{1,2},
                       {0, 2},
                       {0, 1},
                       {0,1,2}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedArticulationNode2){
    TestCase testCase({{0,1,2},
                       {0, 2},
                       {0, 1},
                       {0,1,2}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedArticulationNode3){
    TestCase testCase({{0,1,2},
                       {1, 2},
                       {0, 1},
                       {0,1,2}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedArticulationNode4){
    TestCase testCase({{0,1,3},
                       {1, 2},
                       {0, 2},
                       {0,1,}}, 4, 4);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedMerging1){
    TestCase testCase({{1,2,3},
                       {0,3},
                       {0,2},
                       {0,1}}, 4, 4);
    runTestCase(testCase);
}
TEST(RowAddition,TriconnectedMerging2){
    TestCase testCase({{0,2},
                       {0,1},
                       {0,1,2},
                       {0,1,2}}, 4, 3);
    runTestCase(testCase);
}
TEST(RowAddition,RigidRootMerging){
    runTestCase(seedToTestCase(29658,4,4));
}
TEST(RowAddition,NonGraphicArticulationRigid){
    runTestCase(seedToTestCase(31151,4,4));
}
TEST(RowAddition,NotGraphicMerging){
    runTestCase(seedToTestCase(31711,4,4));
}
TEST(RowAddition,TriconnectedMergingArticulation1){
    runTestCase(seedToTestCase(62302,4,4));
}
TEST(RowAddition,TriconnectedArticulationPropagation){
    runTestCase(seedToTestCase(62303,4,4));
}
TEST(RowAddition,TriconnectedMergeIntoOtherArticulation){
    runTestCase(seedToTestCase(62333,4,4));
}
TEST(RowAddition,ArticulationPointRoot){
    runTestCase(seedToTestCase(62423,4,4));
}
TEST(RowAddition,SeriesFailure){
    runTestCase(seedToTestCase(63411,4,4));
}
TEST(RowAddition,TriconnectedMerging3){
    runTestCase(seedToTestCase(4085,4,3));
}
TEST(RowAddition,TriconnectedPropagation1){
    runTestCase(seedToTestCase(29535,4,4));
}
TEST(RowAddition,NotGraphic1){
    TestCase testCase({{0,1,2},
                       {1, 2},
                       {0, 2},
                       {0,1}}, 4, 3);
    runTestCase(testCase);
}

TEST(RowAddition,RigidIntoParallel1){
    runTestCase(seedToTestCase(13790,4,4));
}

TEST(RowAddition,SearchFailure){
    runTestCase(seedToTestCase(13794,4,4));
}
TEST(RowAddition,NotMinimal){
    runTestCase(seedToTestCase(14041,5,3));
}
TEST(RowAddition,CutEdgeNoMemory){
    runTestCase(seedToTestCase(29653,5,3));
}
TEST(RowAddition,reorderingBug){
    runTestCase(seedToTestCase(1652,4,3));
}
TEST(RowAddition,RigidNoCuts){
    runTestCase(seedToTestCase(111290,4,5));
}
TEST(RowAddition,WrongSplit){
    runTestCase(seedToTestCase(239554,4,5));
}
TEST(RowAddition,CMRBug){
    runTestCase(seedToTestCase(241337,4,5));
}
TEST(RowAddition,ArticulationEdgeRigidIntoParent){
    runTestCase(seedToTestCase(509643,4,5));
}
TEST(RowAddition,RigidMergingArticulation){
    runTestCase(seedToTestCase(1019318,4,5));
}
TEST(RowAddition,RigidMergingArticulation2){
    runTestCase(seedToTestCase(1019319,4,5));
}
TEST(RowAddition,FiveByFourNotGraphic1){
    runTestCase(seedToTestCase(206238,5,4));
}
TEST(RowAddition,FiveByFourCase2){
    runTestCase(seedToTestCase(214474,5,4));
}
TEST(RowAddition,FiveByFourIncorrectMerge){
    runTestCase(seedToTestCase(472493,5,4));
}
TEST(RowAddition,FiveByFourBadColoring){
    //case where coloring for articulation node is not bipartite
    runTestCase(seedToTestCase(472508,5,4));
}

TEST(RowAddition,FiveByFourUnknown){
    runTestCase(seedToTestCase(472523,5,4));
}
TEST(RowAddition,ArticulationPointSearchBug1){
    runTestCase(seedToTestCase(996783,5,4));
}
TEST(RowAddition,ArticulationPointSearchBug2){
    runTestCase(seedToTestCase(996797,5,4));
}
TEST(RowAddition,SixBySix1){
    runTestCase(seedToTestCase(47091484421,6,6));
}
TEST(RowAddition,SixBySixThreeArticulationPoints){
    runTestCase(seedToTestCase(16273736387,6,6));
}
TEST(RowAddition,SixByFourThreeArticulationPoints){
    runTestCase(seedToTestCase(7559087,6,4));
}
TEST(RowAddition,FiveByFiveOne){
    runTestCase(seedToTestCase(31804393,5,5));
}
TEST(RowAddition,FiveByFiveParallelSplitting){
    runTestCase(seedToTestCase(33223995,5,5));
}
TEST(RowAddition,FiveByFiveWrong){
    runTestCase(seedToTestCase(13616375,5,5));
}
TEST(RowAddition, Python1) {
    //Simple case where series should be propagated to a parallel
    TestCase testCase({{0},
                       {0, 1},
                       {0, 1}}, 3, 3);
    runTestCase(testCase);
}


TEST(RowAddition, Python2) {
    TestCase testCase({{0, 1, 2},
                       {1, 2},
                       {0, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python3) {
    TestCase testCase({{1, 2},
                       {0, 2},
                       {2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python4) {
    TestCase testCase({{1, 2},
                       {1},
                       {0},
                       {0, 1, 2}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python5) {
    TestCase testCase({{0, 2},
                       {0, 1},
                       {0, 1, 2}}, 3, 3);
    runTestCase(testCase);
}

TEST(RowAddition, Python6) {
    TestCase testCase({{0, 1},
                       {1, 2},
                       {0, 2}}, 3, 3);
    runTestCase(testCase);
}

TEST(RowAddition, Python7) {
    TestCase testCase({{1, 2, 3},
                       {0, 2},
                       {0, 1}}, 3, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python8) {
    TestCase testCase({{0, 1, 3},
                       {0, 1, 2},
                       {0, 1, 2, 3}}, 3, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python9) {
    TestCase testCase({{1, 2},
                       {0, 2},
                       {0, 1},
                       {0}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python10) {
    TestCase testCase({{1, 2},
                       {0, 2},
                       {0, 1},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python11) {
    TestCase testCase({{1, 2},
                       {0, 2},
                       {0, 1},
                       {0}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python12) {
    TestCase testCase({{1, 2},
                       {0, 2},
                       {0},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}


TEST(RowAddition, Python13) {
    TestCase testCase({{0, 1, 2},
                       {1, 2},
                       {0, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python14) {
    TestCase testCase({{1, 2, 3},
                       {0, 3},
                       {0, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python15) {
    TestCase testCase({{2, 3},
                       {0, 1, 3},
                       {0, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python16) {
    TestCase testCase({{0, 1, 3},
                       {0, 1, 2, 3},
                       {0, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python17) {
    TestCase testCase({{0, 1, 2},
                       {0, 2, 3},
                       {0, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python18) {
    TestCase testCase({{2, 3},
                       {0, 3},
                       {0, 1, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python19) {
    TestCase testCase({{1, 2, 3},
                       {0, 3},
                       {0, 1, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python20) {
    TestCase testCase({{0, 2, 3},
                       {1, 3},
                       {0, 1, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python21) {
    TestCase testCase({{0, 1, 2, 3},
                       {1, 3},
                       {0, 1, 2},
                       {0, 1}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python22) {
    TestCase testCase({{1, 2},
                       {0, 2},
                       {0, 1},
                       {0, 1, 2}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python23) {
    TestCase testCase({{0, 1, 2},
                       {0, 2},
                       {0, 1},
                       {0, 1, 2}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python24) {
    TestCase testCase({{1, 2, 3},
                       {0, 3},
                       {0, 1},
                       {0, 1, 2}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python25) {
    TestCase testCase({{1, 3},
                       {1, 2},
                       {0, 3},
                       {0, 1, 2}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python26) {
    TestCase testCase({{1, 3},
                       {0, 2},
                       {0, 1},
                       {0, 1, 2, 3}}, 4, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python27) {
    TestCase testCase({{0, 1, 2},
                       {0, 1, 2},
                       {0, 2},
                       {0, 1},
                       {0, 1, 2}}, 5, 3);
    runTestCase(testCase);
}

TEST(RowAddition, Python28) {
    TestCase testCase({{2, 3},
                       {0, 1, 3},
                       {0, 2},
                       {0},
                       {0, 1}}, 5, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python29) {
    TestCase testCase({{1, 3},
                       {0, 3},
                       {1, 2},
                       {0, 2},
                       {0, 1}}, 5, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python30) {
    TestCase testCase({{2, 3},
                       {0, 3},
                       {1, 2},
                       {0, 2},
                       {0, 1}}, 5, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python31) {
    TestCase testCase({{0, 1, 2, 3},
                       {0, 3},
                       {0, 2},
                       {0, 2, 3},
                       {0, 1}}, 5, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python32) {
    TestCase testCase({{1, 2, 3},
                       {1, 3},
                       {1, 2},
                       {0, 3},
                       {0, 1, 2}}, 5, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python33) {
    TestCase testCase({{1, 3, 4},
                       {0, 2, 4},
                       {2, 3},
                       {0, 1}}, 4, 5);
    runTestCase(testCase);
}

TEST(RowAddition, Python34) {
    TestCase testCase({{0, 1, 3, 4},
                       {0, 1, 2},
                       {0, 4},
                       {0, 1, 2, 3}}, 4, 5);
    runTestCase(testCase);
}

TEST(RowAddition, Python35) {
    TestCase testCase({{2, 4, 5},
                       {2, 5},
                       {3, 5},
                       {0, 2, 4},
                       {0, 1, 3, 5},
                       {1}}, 6, 6);
    runTestCase(testCase);
}

TEST(RowAddition, Python36) {
    TestCase testCase({{2, 3},
                       {0, 2, 5},
                       {1, 2, 5},
                       {4, 5},
                       {0, 3, 4, 5},
                       {1, 3}}, 6, 6);
    runTestCase(testCase);
}

TEST(RowAddition, Python37) {
    TestCase testCase({{4, 5},
                       {0, 1, 3},
                       {0, 1, 2, 4},
                       {0, 2, 4},
                       {2, 4, 5},
                       {0, 1, 4, 5}}, 6, 6);
    runTestCase(testCase);
}

TEST(RowAddition, Python38) {
    TestCase testCase({{0, 2},
                       {0, 2},
                       {0, 1},
                       {0, 1},
                       {1},
                       {0, 2, 3},
                       {1, 2},
                       {0, 2},
                       {0, 1},
                       {0, 3}}, 10, 4);
    runTestCase(testCase);
}

TEST(RowAddition, Python39) {
    TestCase testCase({{4, 5, 6},
                       {4, 6},
                       {0, 5, 6},
                       {1, 4, 5, 6},
                       {0, 1, 2, 3, 4, 5, 6},
                       {0, 3, 4, 5},
                       {3, 6}}, 7, 7);
    runTestCase(testCase);
}

TEST(RowAddition, Python40) {
    TestCase testCase({{0, 2, 3, 5, 6},
                       {0, 1, 4, 6},
                       {1, 6},
                       {3, 5, 6},
                       {0, 1, 4, 5, 6},
                       {2, 4},
                       {3, 4, 5, 6}}, 7, 7);
    runTestCase(testCase);
}

TEST(RowAddition, Python41) {
    TestCase testCase({{0, 1, 4, 5, 6, 10},
                       {1, 2, 4, 6, 7, 8},
                       {7, 10},
                       {0, 2, 5},
                       {8},
                       {3, 4, 5, 6, 10},
                       {4, 5},
                       {3, 6, 10},
                       {1, 4, 5, 6, 10},
                       {1, 4, 6, 7, 8},
                       {6, 7, 9},
                       {4, 8, 9},
                       {},
                       {7, 10},
                       {8, 9},
                       {4, 8, 9},
                       {6, 7},
                       {}}, 18, 11);
    runTestCase(testCase);
}

TEST(RowAddition,AllTwoBySix){
    //This should cover most, if not all, cases with two rows.
    //Note that these matrices are always graphic by definition
    runAllMatrices(2,6);
}
TEST(RowAddition,AllThreeByThree){
    runAllMatrices(3,3);
}
TEST(RowAddition,AllThreeBySix){
    runAllMatrices(3,6);
}
TEST(RowAddition,AllFourByThree){
    runAllMatrices(4,3);
}
TEST(RowAddition,AllFourByFour){
    runAllMatrices(4,4);
}
TEST(RowAddition,AllFiveByThree){
    runAllMatrices(5,3);
}
TEST(RowAddition,AllSixByThree){
    runAllMatrices(6,3);
}

TEST(RowAddition,SixBySix4){
    runTestCase(seedToTestCase(31095735969,6,6));
}
TEST(RowAddition,RandomSixBySix){
    randomlySample(6,6,100'000,42);
}
TEST(RowAddition,RandomSevenBySeven){
    randomlySample(7,7,100'000,3);
}
TEST(RowAddition,Random8By7){
    randomlySample(8,7,100'000,16);
}
TEST(RowAddition,Random10By6){
    randomlySample(10,6,100'000,20);
}
TEST(RowAddition,Random6By10){
    randomlySample(6,10,100'000,38);
}
TEST(RowAddition,SixBySix3){
    runTestCase(seedToTestCase(11258910622,6,6));
}
TEST(RowAddition,SevenBySeven2){
    runTestCase(seedToTestCase(473114690647189,7,7));
}
TEST(RowAddition,SevenBySeven1){
    //Has a relatively complex merge with a few edge cases. Originally produced a cycle of tree edges.
    runTestCase(seedToTestCase(527698575261016,7,7));
}
TEST(RowAddition,TenBySix1){
    runTestCase(seedToTestCase(500425441027923186,10,6));
}
TEST(RowAddition,TenBySix2){
    //Failed to zero out all colors of a previous iteration of articulation points?
    runTestCase(seedToTestCase(1032551395777074495,10,6));
}
TEST(RowAddition,TenBySix3){
    //Failed to zero out all colors of a previous iteration of articulation points?
    runTestCase(seedToTestCase(533753527437188797,10,6));
}
TEST(RowAddition,SevenBySeven3){
    //Contains more than four articulation points, not all on Q, but is still graphic
    runTestCase(seedToTestCase(545165200943522,7,7));
}
TEST(RowAddition,SevenBySeven4){
    //Contains more than four articulation points, not all on Q, but not graphic?
    runTestCase(seedToTestCase(562571741887140,7,7));
}
TEST(RowAddition,AllFourByFive){
    runAllMatrices(4,5);
}
TEST(RowAddition,AllFiveByFour){
    runAllMatrices(5,4);
}

TEST(RowAddition,ErdosRenyi8Density20){
    for(std::size_t seed = 0; seed < 1000; seed++){
        TestCase testCase = createErdosRenyiTestcase(8,0.2,seed);
        runTestCase(testCase);
    }
}
TEST(RowAddition,ErdosRenyi20Density20){
    for(std::size_t seed = 0; seed < 1000; seed++){
        TestCase testCase = createErdosRenyiTestcase(20,0.2,seed);
        runTestCase(testCase);
    }
}

TEST(RowAddition,ErdosRenyi20Density50){
    for(std::size_t seed = 0; seed < 1000; seed++){
        TestCase testCase = createErdosRenyiTestcase(20,0.5,seed);
        runTestCase(testCase);
    }
}
TEST(RowAddition,ErdosRenyi100Density10){
    TestCase testCase = createErdosRenyiTestcase(100,0.1,0);
    runTestCase(testCase);

}

TEST(RowAddition,SixBySixMemoryProblem){
    runTestCase(seedToTestCase(68254337919,6,6));
}

//Below tests take quite long ~ up to an hour together
//TEST(RowAddition,AllSixByFour){
//    runAllMatrices(6,4);
//}
//TEST(RowAddition,AllFiveByFive){
//    runAllMatrices(5,5);
//}
//TEST(RowAddition,ErdosRenyi500Density5){
//    TestCase testCase = createErdosRenyiTestcase(500,0.05,42);
//    runTestCase(testCase);
//
//}
//TEST(RowAddition,ErdosRenyi1000Density1){
//    TestCase testCase = createErdosRenyiTestcase(1000,0.01,1);
//    //999 by 3780 matrix
//    runTestCase(testCase);
//}
