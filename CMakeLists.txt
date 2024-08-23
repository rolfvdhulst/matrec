cmake_minimum_required(VERSION 3.5)

project(matrec
        VERSION 0.1.0
        DESCRIPTION "Algorithms for recognizing graphic and network matrices"
        LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

option(BUILD_TESTS "Build the tests, require GTest and CMR to be installed" OFF)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release"
            CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()


add_library(matrec
src/Graphic.c
src/Incidence.c
src/Matrix.c
src/Network.c
src/Shared.c
        src/SignCheckRowAddition.c
include/matrec/Graphic.h
include/matrec/Incidence.h
include/matrec/Matrix.h
include/matrec/Shared.h
include/matrec/Network.h
        include/matrec/SignCheckRowAddition.h
)

add_library(matrec::matrec ALIAS matrec)

target_include_directories(matrec
PUBLIC include/
)
target_compile_options(matrec PRIVATE
        -Wall
        -Wextra # reasonable and standard
        -Wshadow # warn the user if a variable declaration shadows one from a parent context
        -Wcast-align # warn for potential performance problem casts
        -Wunused # warn on anything being unused
        -Wpedantic # warn if non-standard C++ is used
        -Wconversion # warn on type conversions that may lose data
        -Wsign-conversion # warn on sign conversions
        -Wnull-dereference # warn if a null dereference is detected
        -Wdouble-promotion # warn if float is implicit promoted to double
        -Wformat=2 # warn on security issues around functions that format output (ie printf)
        -Wimplicit-fallthrough # warn on statements that fallthrough without an explicit annotation)
        -Wmisleading-indentation # warn if indentation implies blocks where blocks do not exist
        -Wduplicated-cond # warn if if / else chain has duplicated conditions
        -Wduplicated-branches # warn if if / else branches have duplicated code
        -Wlogical-op # warn about logical operations being used where bitwise were probably wanted
)

set_target_properties(
        matrec
        PROPERTIES VERSION ${PROJECT_VERSION}
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN YES
)


if(BUILD_TESTS)
    find_package(CMR REQUIRED)
    find_package(GTest REQUIRED)
    add_executable(matrec_test
            test/main.cpp
            test/TestHelpers.cpp
            test/TestHelpers.h
            test/GraphicColumnAdditionTest.cpp
            test/GraphicRowAdditionTest.cpp
            test/GraphicTest.cpp
            test/IncidenceTest.cpp #TODO
            test/NetworkTest.cpp
    )

    target_link_libraries(matrec_test
            PUBLIC matrec::matrec
            PUBLIC GTest::GTest
            PUBLIC CMR::cmr
    )
endif()

