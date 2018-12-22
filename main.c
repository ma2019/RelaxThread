//
// Created by moo-ease on 05/11/18.
//

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<pthread.h>
#include<time.h>

double* globalTestArray;

pthread_t* threadArray;

int dimension;
int numberOfThreads;
double precision;

int relaxableElements;
int numberOfElementsPerThread;
int* indicesToRelax;

double precisionReached  = 0;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_barrier_t barrier;

/*
 * The following function is relaxes the array in parallel.
 * Every thread obtains its thread index, and works on a certain
 * range of elements within the array.
 * A copy with the same dimensions as the problem is created
 * which contains values from the previous iteration.
 * This allows the threads to calculate the precisions.
 * A mutex lock and a barrier is used to avoid race conditions.
 */
void *relaxThreaded(void *arg) {

    double achievedPrecision;

    pthread_t selfThread = pthread_self();

    int threadIndex = 0;

    for (int i = 0; i <= numberOfThreads; i++) {
        if (pthread_equal(selfThread, threadArray[i])) {
            threadIndex = i;
        }
    }

    double copyOfArray[dimension*dimension];

    int startingIndex = threadIndex*numberOfElementsPerThread;

    int stoppingIndex;

    if (threadIndex != (numberOfThreads - 1)) {
        stoppingIndex = startingIndex + numberOfElementsPerThread;
    } else {
        stoppingIndex = relaxableElements;
    }

    while (precisionReached == 0) {
        for (int i = startingIndex; i < stoppingIndex; i++) {
            int indexToRelax = indicesToRelax[i];
            copyOfArray[indexToRelax] = globalTestArray[indexToRelax];
            globalTestArray[indexToRelax] = (((globalTestArray[(indexToRelax) - 1]) +
                                        (globalTestArray[(indexToRelax) + 1]) +
                                        (globalTestArray[(indexToRelax) - dimension]) +
                                        (globalTestArray[(indexToRelax) + dimension])) / 4);
            pthread_barrier_wait(&barrier);

            pthread_mutex_lock(&mutex);
            achievedPrecision = copyOfArray[indexToRelax] - globalTestArray[indexToRelax];
            if (achievedPrecision < 0) {
                achievedPrecision = achievedPrecision * (-1);
            }
            if (achievedPrecision < precision) {
                precisionReached = 1;
            } else {
                precisionReached = 0;
            }
            pthread_mutex_unlock(&mutex);
        }
    }
}

/*
 * This function solves the array sequentiall.
 * The array is transformed into a 1 dimensional array,
 * with size dimension*dimension, so that
 * an element in a position (row,column) in a 2D array is
 * positioned at (dimesnion*row+column).
 *
 */
void sequentialRelax(void* arg) {
    double achievedPrecision;

    double copyOfArray[dimension*dimension];

    while (precisionReached == 0) {
        for (int row = 1; row < dimension - 1; row++) {
            for (int column = 1; column < dimension - 1; column++) {
                copyOfArray[dimension * row + column] =
                        globalTestArray[dimension * row + column];
                globalTestArray[dimension * row + column] =
                        (((globalTestArray[(dimension * row + column) - 1])
                        + (globalTestArray[(dimension * row + column) + 1])
                        + (globalTestArray[(dimension * row + column)
                        - dimension]) + (globalTestArray[(dimension * row
                        + column) + dimension])) / 4);
                achievedPrecision = copyOfArray[dimension * row + column]
                        - globalTestArray[dimension * row + column];
                if (achievedPrecision < 0) {
                    achievedPrecision = achievedPrecision * (-1);
                }
                if (achievedPrecision < precision) {
                    precisionReached = 1;
                } else {
                    precisionReached = 0;
                }
            }

        }
    }
}

/*
 * Function to generate random numbers that populate the array.
 */

double generateRandom(double randMin, double randMax) {
    double rangeOfNumbers = randMax - randMin;
    double divide = RAND_MAX/rangeOfNumbers;
    return randMin + (rand()/divide);
}

/*
 * The main function takes in 6 arguments. It calls the
 * sequential algorithm on the problem array, and then creates
 * threads that enter the parallel algorithm to work on the arrays.
 * The clock() function is used to find the time at the beginning
 * and at the end of calling the sequential and parallel functions.
 * The difference between the two times is taken to be the time
 * required for the program to come up with a solution.
 */

int main(int argc, char **argv) {

    clock_t startTime, endTime;
    double seqTimeTaken, parallelTimeTaken;

    if(argc < 6) {
        printf("Too few arguments.\n");
        return 1;
    }

    dimension = atoi(argv[1]);

    numberOfThreads = atoi(argv[2]);

    precision = strtod(argv[3], NULL);

    double randMin = strtod(argv[4], NULL);

    double randMax = strtod(argv[5], NULL);

    if(numberOfThreads < 1 || numberOfThreads > (dimension-2)) {
        numberOfThreads = (dimension-2);
        printf("Invalid number of threads. Switching to default number: %d.\n",
                numberOfThreads);

    }
    double *relaxArray = (double*)malloc(dimension*dimension* sizeof(double));

    for (int i = 0; i< dimension*dimension; i++) {
        relaxArray[i] = generateRandom(randMin, randMax);
    }

    globalTestArray = (double*)malloc(dimension*dimension* sizeof(double));

    for (int i = 0; i < dimension*dimension; i++) {
        globalTestArray[i] = relaxArray[i];
    }

    free(relaxArray);

    startTime = clock();

    sequentialRelax((void*)globalTestArray);

    endTime = clock();

    seqTimeTaken = ((double)(endTime-startTime));

    relaxableElements = pow((dimension-2), 2);

    indicesToRelax = (int*)malloc(relaxableElements * sizeof(int));

    for(int i = 0; i < relaxableElements; i++) {
        indicesToRelax[i] = -1;
    }
    int indexToAddTo = 0;

    for (int row = 1; row < dimension - 1; row++) {
        for (int column = 1; column < dimension - 1; column++) {
            int index = dimension * row + column;
            indicesToRelax[indexToAddTo] = index;
            indexToAddTo++;
        }
    }

    numberOfElementsPerThread = relaxableElements/numberOfThreads;

    pthread_barrier_init(&barrier, NULL, numberOfThreads);

    pthread_t tid[numberOfThreads];

    threadArray = (pthread_t*)malloc(sizeof(pthread_t)*numberOfThreads);

    startTime = clock();

    for (int i = 0; i < numberOfThreads; i++) {
        pthread_create(&tid[i], NULL, relaxThreaded, NULL);
        threadArray[i] = tid[i];
    }

    for (int i = 0; i < numberOfThreads; i++) {
        pthread_join(tid[i], NULL);
    }

    endTime = clock();

    free(globalTestArray);

    parallelTimeTaken = ((double)(endTime - startTime));

    printf("Sequential time for %d dimensions with %f precision and values from "
           "%f to %f: %f\n",
            dimension, precision, randMin, randMax, seqTimeTaken);

    printf("Parallel time for %d dimensions with %d threads, %f precision "
           "and values from %f to %f: %f\n",
            dimension, numberOfThreads, precision, randMin, randMax, parallelTimeTaken);

}