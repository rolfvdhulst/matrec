#ifndef MATREC_SHARED_H
#define MATREC_SHARED_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h> //defines bool when c++ is not defined
#include <limits.h>

#ifdef __cplusplus
extern "C"{
#endif

typedef enum
{
    MATREC_OKAY = 0, ///No error
    MATREC_ERROR_MEMORY = 1, ///Error in (re)allocation
    MATREC_ERROR_INPUT = 2, ///Error in input matrix/
} MATREC_ERROR;


#define MATREC_CALL(call)                             \
do                                                  \
{                                                   \
    MATREC_ERROR _spqr_error = call;                  \
    if(_spqr_error){                                \
        switch(_spqr_error){                        \
            case MATREC_ERROR_MEMORY:                 \
            {   printf("Memory allocation failed"); \
                break;}                             \
            case MATREC_ERROR_INPUT:                  \
            {   printf("User input error");         \
                break;}                             \
            default:                                \
                printf("Unrecognized error code "); \
        }                                           \
        printf(" in %s:%d.\n", __FILE__, __LINE__); \
        return _spqr_error;                         \
    }                                               \
} while(false)                                      \

struct MATREC_ENVIRONMENT{
FILE * output;
};

typedef struct MATREC_ENVIRONMENT MATREC;

MATREC_ERROR MATRECcreateEnvironment(MATREC** pSpqr);
MATREC_ERROR MATRECfreeEnvironment(MATREC** pSpqr);

#define MATRECallocBlockArray(spqr, ptr, length) \
    MATRECimplAllocBlockArray(spqr,(void **) (ptr), sizeof(**(ptr)),length)

#define MATRECreallocBlockArray(spqr, ptr, length) \
    MATRECimplReallocBlockArray(spqr,(void **) (ptr), sizeof(**(ptr)),length)

#define MATRECfreeBlockArray(spqr, ptr) \
    MATRECimplFreeBlockArray(spqr,(void **) (ptr))

MATREC_ERROR MATRECimplAllocBlockArray(MATREC * env, void** ptr, size_t size, size_t length);
MATREC_ERROR MATRECimplReallocBlockArray(MATREC* env, void** ptr, size_t size, size_t length);
void MATRECimplFreeBlockArray(MATREC* env, void ** ptr);


#define MATRECallocBlock(spqr, ptr) \
    MATRECimplAllocBlock(spqr,(void **) (ptr), sizeof(**(ptr)))

#define MATRECfreeBlock(spqr, ptr) \
    MATRECimplFreeBlock(spqr,(void **) (ptr))

MATREC_ERROR MATRECimplAllocBlock(MATREC * env, void **ptr, size_t size);

void MATRECimplFreeBlock(MATREC * env, void **ptr);

///Types which define matrix sizes
///Aliased so that switching is much easier if ever desired

typedef size_t MATREC_matrix_size;
typedef MATREC_matrix_size MATREC_row;
typedef MATREC_matrix_size MATREC_col;

#define MATREC_INVALID ULLONG_MAX //TODO: this should maybe be ULLONG_MAX or SSIZE_MAX? Check...
#define MATREC_INVALID_ROW MATREC_INVALID
#define MATREC_INVALID_COL MATREC_INVALID

bool MATRECrowIsInvalid(MATREC_row row);
bool MATRECrowIsValid(MATREC_row row);
bool MATRECcolIsInvalid(MATREC_col col);
bool MATRECcolIsValid(MATREC_col col);




#ifdef __cplusplus
}
#endif

#endif //MATREC_SHARED_H

