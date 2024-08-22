#include "TestHelpers.h"

#include <utility>
#include <iostream>

TestCase::TestCase(std::vector<std::vector<MATREC_col>> matrix, std::size_t rows, std::size_t cols, std::size_t seed )
        : rows{rows}, cols{cols}, matrix{std::move(matrix)},seed{seed} {};

ColTestCase::ColTestCase(const TestCase& testCase) : rows{testCase.rows}, cols{testCase.cols},
                                                     matrix(testCase.cols,std::vector<MATREC_row>()), seed{testCase.seed}{
    for(std::size_t row = 0; row < testCase.rows; ++row){
        for(std::size_t col : testCase.matrix[row]){
            matrix[col].push_back(row);
        }
    }
}

TestCase stringToTestCase(std::string mat, uint64_t rows, uint64_t cols){
    std::size_t count = 0;
    std::size_t seed = 0;
    std::size_t power = 1;
    for(auto it = mat.begin(); it != mat.end(); ++it){
        if(*it== '1'){
            seed += power;
            power *= 2;
            ++count;
        }else if(*it == '0'){
            power *= 2;
            ++count;
        }
        if(count == rows*cols){
            break;
        }
    }
    return seedToTestCase(seed,rows,cols);
}
TestCase seedToTestCase(uint64_t seed, uint64_t rows, uint64_t cols){
    assert(rows*cols < 64ul);
    assert(seed < (1ul << rows*cols));
    std::vector<std::vector<MATREC_col>> matrix(rows);
    uint64_t seedCopy = seed;
    int nRows = rows;
    int nCols = cols;
    for (MATREC_row row = 0; row < nRows; ++row) {
        for (MATREC_col col = 0; col < nCols; ++col) {
            if(seed & 1){
                matrix[row].push_back(col);
            }
            seed >>= 1;
        }
    }
    return {matrix,rows,cols,seedCopy};
}


bool Edge::operator==(const Edge &other) const {
    return head == other.head && tail == other.tail;
}

int findRepresentative(int entry, std::vector<int>& representative){
    int current = entry;
    int next;

    //traverse down tree to find the root
    while ((next = representative[current]) >= 0) {
        current = next;
        assert(current < representative.size());
    }

    int root = current;
    current = entry;

    while ((next = representative[current]) >= 0) {
        representative[current] = root;
        current = next;
        assert(current < representative.size());
    }
    return root;
}
int makeUnion(std::vector<int>& representative, int first, int second){
    assert(representative[first] < 0 && representative[second] < 0);

    //The rank is stored as a negative number: we decrement it making the negative number larger.
    // We want the new root to be the one with 'largest' rank, so smallest number. If they are equal, we decrement.
    int firstRank = representative[first];
    int secondRank = representative[second];
    if (firstRank > secondRank) {
        std::swap(first,second);
    }

    representative[second] = first;
    if (firstRank == secondRank) {
        --representative[first];
    }
    return first;
}

TestCase createErdosRenyiTestcase(std::size_t nodes, double density, std::size_t seed){
    assert(nodes >= 1);
    if(nodes <= 1){
        return {{},0,0,seed};
    }


    //Create erdos renyi graph:
    std::minstd_rand gen(seed);
    std::bernoulli_distribution dist(density);
    std::vector<Edge> edges;
    for (std::size_t i = 1; i < nodes; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            if(dist(gen)){
                edges.push_back(Edge{.head = i, .tail = j});
            }
        }
    }

    std::vector<int> nodeRepresentative(nodes,-1);
    std::shuffle(edges.begin(),edges.end(),gen);
    std::vector<Edge> treeEdges;
    for(const auto& edge: edges){
        int firstRep = findRepresentative(edge.head,nodeRepresentative);
        int secondRep = findRepresentative(edge.tail,nodeRepresentative);
        if( firstRep != secondRep ){
            //Add edge to tree
            makeUnion(nodeRepresentative,firstRep,secondRep);
            treeEdges.push_back(edge);
            if(treeEdges.size() == (nodes-1)){
                break;
            }
        }
    }
    std::size_t matrixRows = treeEdges.size();
    std::size_t matrixCols = edges.size() - matrixRows;
    if(matrixCols == 0){
        return {{},0,0,seed};
    }
    struct TreeEdge{
        std::size_t otherNode;
        std::size_t edgeNum;
    };
    std::vector<std::vector<TreeEdge>> treeAdjacencyMatrix(nodes);
    for(std::size_t edgeID = 0; edgeID < matrixRows; edgeID++){
        const Edge& edge = treeEdges[edgeID];
        treeAdjacencyMatrix[edge.head].push_back(TreeEdge{.otherNode = edge.tail,.edgeNum = edgeID});
        treeAdjacencyMatrix[edge.tail].push_back(TreeEdge{.otherNode = edge.head,.edgeNum = edgeID});
    }
    std::unordered_set<Edge> treeEdgeSet(treeEdges.begin(),treeEdges.end());
    std::vector<std::vector<MATREC_col>> matrix(matrixRows);
    int colIdx = 0;

    struct CallStackInfo{
        std::size_t node;
        std::size_t adjacencyPos;
    };
    std::vector<CallStackInfo> callStack(nodes);
    for(const auto& edge : edges){
        if(treeEdgeSet.find(edge) != treeEdgeSet.end()) continue;
        std::size_t source = edge.head;
        std::size_t target = edge.tail;

        //Do DFS to determine path from source to target

        callStack[0].node = source;
        callStack[0].adjacencyPos = 0;
        std::size_t callStackSize = 1;
        do{
            CallStackInfo& info = callStack[callStackSize-1];
            if(info.node == target){
                break;
            }
            if(info.adjacencyPos == treeAdjacencyMatrix[info.node].size()){
                --callStackSize;
                continue;
            }
            const TreeEdge& treeEdge = treeAdjacencyMatrix[info.node][info.adjacencyPos];
            ++info.adjacencyPos;
            //Cannot take parent edge
            if(callStackSize <= 1 || treeEdge.otherNode != callStack[callStackSize-2].node){
                callStack[callStackSize].node = treeEdge.otherNode;
                callStack[callStackSize].adjacencyPos = 0;
                ++callStackSize;
            }
        }while(callStackSize > 0);
        for (size_t i = 0; i < callStackSize-1; ++i) {
            const std::vector<TreeEdge>& neighbours = treeAdjacencyMatrix[callStack[i].node];
            const TreeEdge& columnEdge = neighbours[callStack[i].adjacencyPos-1];
            matrix[columnEdge.edgeNum].push_back(colIdx);
        }
        colIdx++;
    }

    return {matrix,matrixRows,edges.size()-matrixRows,seed};
}
// for string delimiter
std::vector<std::string> split(const std::string& s, const std::string& delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        if(pos_end != pos_start){
            res.push_back (token);
        }
        pos_start = pos_end + delim_len;

    }

    res.push_back (s.substr (pos_start));
    return res;
}
DirectedTestCase stringToDirectedTestCase(std::string mat, uint64_t rows, uint64_t cols) {
    auto splitString = split(mat," ");
    std::size_t row = 0;
    std::size_t column = 0;

    std::vector<std::vector<Nonzero>> matrix;
    std::vector<Nonzero> buffer;
    for(const auto& string : splitString){
        double value = std::stod(string);
        if(value != 0.0){
            buffer.push_back(Nonzero{.index = column,.value = value});
        }
        ++column;
        if(column == cols){
            matrix.push_back(buffer);
            buffer.clear();

            ++row;
            column = 0;
            if(row == rows){
                break;
            }
        }
    }
    return {matrix,rows,cols};
}

DirectedTestCase seedToDirectedTestCase(std::size_t seed, uint64_t rows, uint64_t cols){
    std::uniform_int_distribution<int> dist(-1,1);
    std::mt19937 gen(seed);
    std::vector<std::vector<Nonzero>> matrix;

    for(std::size_t i = 0; i < rows; ++i){
        std::vector<Nonzero> buffer;
//        std::cout<<"\"";
        for(std::size_t j = 0 ; j < cols; ++j){
            int x = dist(gen);
//            std::cout<<x<<" ";
            if( x == 0) continue;
            buffer.push_back(Nonzero{.index = j, .value = double(x)});
        }
//        std::cout<<"\"\n";
        matrix.push_back(buffer);
    }
    return {matrix,rows,cols};
}

DirectedTestCase::DirectedTestCase(std::vector<std::vector<Nonzero>> matrix, std::size_t rows, std::size_t cols) :
matrix{std::move(matrix)},rows{rows},cols{cols}{

}

DirectedColTestCase::DirectedColTestCase(const DirectedTestCase &testCase) :  rows{testCase.rows},cols{testCase.cols},
                                                                              matrix(testCase.cols,std::vector<Nonzero>()){
    for(std::size_t row = 0; row < testCase.rows; ++row){
        for(auto nonzero : testCase.matrix[row]){
            matrix[nonzero.index].push_back(Nonzero{.index = row, .value = nonzero.value});
        }
    }
}

DirectedTestCase erdosRenyiDirectedTestCase(std::size_t nodes, double density, std::size_t seed){
    assert(nodes >= 1);
    if(nodes <= 1){
        return DirectedTestCase({},0,0);
    }

    //Create erdos renyi graph:
    std::minstd_rand gen(seed);
    std::bernoulli_distribution dist(density);
    std::vector<Edge> edges;
    for (std::size_t i = 0; i < nodes; ++i) {
        for (std::size_t j = 0; j < nodes; ++j) {
            if(i != j && dist(gen)){
                edges.push_back(Edge{.head = i, .tail = j});
            }
        }
    }

    std::vector<int> nodeRepresentative(nodes,-1);
    std::shuffle(edges.begin(),edges.end(),gen);
    std::vector<Edge> treeEdges;
    for(const auto& edge: edges){
        int firstRep = findRepresentative(edge.head,nodeRepresentative);
        int secondRep = findRepresentative(edge.tail,nodeRepresentative);
        if( firstRep != secondRep ){
            //Add edge to tree
            makeUnion(nodeRepresentative,firstRep,secondRep);
            treeEdges.push_back(edge);
            if(treeEdges.size() == (nodes-1)){
                break;
            }
        }
    }
    std::size_t matrixRows = treeEdges.size();
    std::size_t matrixCols = edges.size() - matrixRows;
    if(matrixCols == 0){
        return {{},0,0};
    }
    struct TreeArc{
        std::size_t otherNode;
        std::size_t edgeNum;
        bool forwards;
    };
    std::vector<std::vector<TreeArc>> treeAdjacencyMatrix(nodes);
    for(std::size_t edgeID = 0; edgeID < matrixRows; edgeID++){
        const Edge& edge = treeEdges[edgeID];
        treeAdjacencyMatrix[edge.head].push_back(TreeArc{.otherNode = edge.tail,.edgeNum = edgeID,
                                                         .forwards = false});
        treeAdjacencyMatrix[edge.tail].push_back(TreeArc{.otherNode = edge.head,.edgeNum = edgeID,
                                                         .forwards = true});
    }
    std::unordered_set<Edge> treeEdgeSet(treeEdges.begin(),treeEdges.end());
    std::vector<std::vector<Nonzero>> matrix(matrixRows);
    int colIdx = 0;

    struct CallStackInfo{
        std::size_t node;
        std::size_t adjacencyPos;
    };
    std::vector<CallStackInfo> callStack(nodes);
    for(const auto& edge : edges){
        if(treeEdgeSet.find(edge) != treeEdgeSet.end()) continue;
        std::size_t source = edge.head;
        std::size_t target = edge.tail;

        //Do DFS to determine path from source to target

        callStack[0].node = source;
        callStack[0].adjacencyPos = 0;
        std::size_t callStackSize = 1;
        do{
            CallStackInfo& info = callStack[callStackSize-1];
            if(info.node == target){
                break;
            }
            if(info.adjacencyPos == treeAdjacencyMatrix[info.node].size()){
                --callStackSize;
                continue;
            }
            const TreeArc& treeEdge = treeAdjacencyMatrix[info.node][info.adjacencyPos];
            ++info.adjacencyPos;
            //Cannot take parent edge
            if(callStackSize <= 1 || treeEdge.otherNode != callStack[callStackSize-2].node){
                callStack[callStackSize].node = treeEdge.otherNode;
                callStack[callStackSize].adjacencyPos = 0;
                ++callStackSize;
            }
        }while(callStackSize > 0);
        for (size_t i = 0; i < callStackSize-1; ++i) {
            const std::vector<TreeArc>& neighbours = treeAdjacencyMatrix[callStack[i].node];
            const TreeArc& columnEdge = neighbours[callStack[i].adjacencyPos-1];
            double value = columnEdge.forwards ? 1.0 : -1.0;
            matrix[columnEdge.edgeNum].push_back(Nonzero{.index = static_cast<MATREC_matrix_size>(colIdx),
                                                         .value = value});
        }
        colIdx++;
    }

    return {matrix,matrixRows,edges.size()-matrixRows};
}