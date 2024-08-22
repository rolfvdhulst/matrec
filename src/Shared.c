#include "matrec/Shared.h"

#ifndef NDEBUG
//Only necessary for overflow check assertions
#include <limits.h>
#endif

MATREC_ERROR MATRECcreateEnvironment(MATREC** pSpqr){
    *pSpqr = (MATREC*) malloc(sizeof(MATREC));
    MATREC * env = *pSpqr;
    if(!env){
        return MATREC_ERROR_MEMORY;
    }
    env->output = stdout;
    return MATREC_OKAY;
}
MATREC_ERROR MATRECfreeEnvironment(MATREC** pSpqr){
    if(!pSpqr){
        return MATREC_ERROR_MEMORY;
    }
    MATREC * env = *pSpqr;
    if(!env){
        return MATREC_ERROR_MEMORY;
    }

    free(*pSpqr);
    *pSpqr = NULL;
    return MATREC_OKAY;
}


//TODO: implement other malloc-type functions such as reallocarray and calloc throughout the codebase?

MATREC_ERROR MATRECimplAllocBlockArray(__attribute__((unused)) MATREC * env, void** ptr, size_t size, size_t length){
    assert(env);
    assert(ptr);
    //assert(*ptr == NULL); //TODO: why is this check here, is it necessary?
    assert(!(size > 0 && length > UINT_MAX / size)); //overflow check

    *ptr = malloc(size * length); //TODO check for overflows, possibly use realloc here?

    return *ptr ? MATREC_OKAY : MATREC_ERROR_MEMORY;
}
MATREC_ERROR MATRECimplReallocBlockArray(__attribute__((unused)) MATREC* env, void** ptr, size_t size, size_t length)
{
    assert(env);
    assert(ptr);
    assert(!(size > 0 && length > UINT_MAX / size)); //overflow check
    *ptr = realloc(*ptr, size * length);
    //*ptr = reallocarray(*ptr,length,size); //reallocarray can also work and is a bit safer, but is non-standard...
    return *ptr ? MATREC_OKAY : MATREC_ERROR_MEMORY;
}
void MATRECimplFreeBlockArray(__attribute__((unused)) MATREC* env, void ** ptr){
    assert(env);
    assert(ptr);
    free(*ptr);
    *ptr = NULL;
}

MATREC_ERROR MATRECimplAllocBlock(__attribute__((unused)) MATREC * env, void **ptr, size_t size){
    assert(env);
    assert(ptr);
    *ptr = malloc(size);

    return *ptr ? MATREC_OKAY : MATREC_ERROR_MEMORY;
}

void MATRECimplFreeBlock(__attribute__((unused)) MATREC * env, void **ptr){
    assert(env);
    assert(ptr);
    assert(*ptr);
    free(*ptr);
    *ptr = NULL;
}

