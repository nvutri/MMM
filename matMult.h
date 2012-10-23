/**
 * All MMM Helper methods
 */
#include <papi.h>
#include <papiStdEventDefs.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>


#ifndef _MATMULT
#define _MATMULT
#define X(x, y) X[ (y) * N + (x)]

/**
 * Matrix functions: allocate, index, and destroy
 * N is global variable
 */
extern int N;

double* alloc(int SIZE) {
    return (double*) malloc(sizeof(double) * SIZE * SIZE);
}

void flushCache(double* matrix) {
    /* destroy the matrix */
    free(matrix);
}

void printMatrix(double* X, unsigned N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << X(i, j) << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * Initialize Matrixes
 */
void initialize(double* const X, const double VAL, unsigned N) {
    //pointer pointing to value
    double* pointer;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (VAL == 0)
                X(i, j) = VAL;
            else
                X(i, j) = i + j;
        }
}

/**
 * Papi Handler
 */
void handle_error(int retval) {
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

void init_papi() {
    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT && retval < 0) {
        printf("PAPI library version mismatch!\n");
        exit(1);
    }
    if (retval < 0)
        handle_error(retval);

//    std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(retval)
//            << " MINOR: " << PAPI_VERSION_MINOR(retval) << " REVISION: "
//            << PAPI_VERSION_REVISION(retval) << "\n";
//
}

int begin_papi(int Event) {
    int EventSet = PAPI_NULL;
    int rv;
    /* Create the Event Set */
    if ((rv = PAPI_create_eventset(&EventSet)) != PAPI_OK)
        handle_error(rv);
    if ((rv = PAPI_add_event(EventSet, Event)) != PAPI_OK)
        handle_error(rv);
    /* Start counting events in the Event Set */
    if ((rv = PAPI_start(EventSet)) != PAPI_OK)
        handle_error(rv);
    return EventSet;
}

long_long end_papi(int EventSet) {
    long_long retval;
    int rv;

    /* get the values */
    if ((rv = PAPI_stop(EventSet, &retval)) != PAPI_OK)
        handle_error(rv);

    /* Remove all events in the eventset */
    if ((rv = PAPI_cleanup_eventset(EventSet)) != PAPI_OK)
        handle_error(rv);

    /* Free all memory and data structures, EventSet must be empty. */
    if ((rv = PAPI_destroy_eventset(&EventSet)) != PAPI_OK)
        handle_error(rv);

    return retval;
}

#endif
