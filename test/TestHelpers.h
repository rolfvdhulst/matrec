#ifndef MATREC_TESTHELPERS_H
#define MATREC_TESTHELPERS_H

#include <matrec/Shared.h>
#include <random>
#include <fstream>
#include <unordered_set>
#include <algorithm>

struct Nonzero {
    MATREC_matrix_size index;
    double value;
};
struct DirectedTestCase {
    DirectedTestCase(std::vector<std::vector<Nonzero>> matrix, std::size_t rows, std::size_t cols);
    std::size_t rows;
    std::size_t cols;
    std::vector<std::vector<Nonzero>> matrix;
};
struct DirectedColTestCase {
    explicit DirectedColTestCase(const DirectedTestCase& testCase);
    std::size_t rows;
    std::size_t cols;
    std::vector<std::vector<Nonzero>> matrix;
};
DirectedTestCase stringToDirectedTestCase(std::string string, uint64_t rows, uint64_t cols);

DirectedTestCase seedToDirectedTestCase(std::size_t seed, uint64_t rows, uint64_t cols);

DirectedTestCase erdosRenyiDirectedTestCase(std::size_t nodes, double density, std::size_t seed);

class TestCase {
public:
    TestCase(std::vector<std::vector<MATREC_col>> matrix, std::size_t rows, std::size_t cols, std::size_t seed = 0);
    std::size_t rows;
    std::size_t cols;
    std::vector<std::vector<MATREC_col>> matrix;
    std::size_t seed;
};

class ColTestCase {
public:
    explicit ColTestCase(const TestCase& testCase);

    std::size_t rows;
    std::size_t cols;
    std::vector<std::vector<MATREC_row>> matrix;
    std::size_t seed;
};

TestCase seedToTestCase(uint64_t seed, uint64_t rows, uint64_t cols);

TestCase stringToTestCase(std::string string, uint64_t rows, uint64_t cols);

struct Edge{
    std::size_t head;
    std::size_t tail;

    bool operator == (const Edge& other) const;
};

int findRepresentative(int entry, std::vector<int>& representative);

int makeUnion(std::vector<int>& representative, int first, int second);

namespace std {

    template <>
    struct hash<Edge>
    {
        std::size_t operator()(const Edge& k) const
        {
            return (k.head << 32) + k.tail;
        }
    };

}

TestCase createErdosRenyiTestcase(std::size_t nodes, double density, std::size_t seed);


#endif //MATREC_TESTHELPERS_H
