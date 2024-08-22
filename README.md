## MATREC: Matrix recognition algorithms
This repository contains algorithms to recognize graphic matrices and network matrices.
Graphic matrices and Network matrices are important objects in optimization and combinatorics.
Network matrices and transposed are a large subclass of totally unimodular matrices,
which makes them particularly interesting for mixed-integer linear programming applications.

Graphic.h contains two algorithms; a column-wise algorithm, which is based on:

Robert E. Bixby, Donald K. Wagner, (1988) An Almost Linear-Time Algorithm for Graph Realization. Mathematics of Operations Research 13(1):99-123.

and a row-wise algorithm, which is based on our recent work:

Rolf van der Hulst, Matthias Walter, (2024) A Row-wise Algorithm for Graph Realization. [Arxiv link]()

In Network.h, we adapted the algorithms from Graphic.h to the network matrix setting.
If you use this software in a publication, please cite our preprint.

### Dependencies
The library has no dependencies. 

If users which to build the tests, then they need to install the following packages:
- [Google Test](https://github.com/google/googletest)
- [CMR]( https://github.com/discopt/cmr )
### Building the software
1. Create a build directory
 
`mkdir build`

2. Navigate to the build directory

`cd build`

3. Configure the build

`cmake ..`

Optionally, users can add `-DBUILD_TESTS=ON` to build the tests. 
Note that for these, dependencies are required.

4. Compile:

`make`

5. (Optional) Install the library

`make install`


