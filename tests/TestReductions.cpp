#include "mpi.h"
#include <stdlib.h>
#include <limits>
#include <iostream>
#include <vector>
#include "csparse.h"
#include "assertions.hpp"
#include "MPIReductions.hpp"
using namespace std;
using namespace Tungsten;


void test_MPIvectorAllMaxloc() {
    const int RANK = MPI::COMM_WORLD.Get_rank();
    const int SIZE = MPI::COMM_WORLD.Get_size();
    const int sendLength = 100*SIZE;
    const int recvLength = sendLength/SIZE;

    int trueMaxLoc = 0;
    float trueMax = 0; 
    vector<float> sendbuf(sendLength);
    vector<float> recvbuf(recvLength);

    if (RANK == 0) {
        srand (time(NULL));
        trueMaxLoc = 0;
        trueMax = -std::numeric_limits<float>::max();
        for (int i = 0; i < sendLength; i++) {
            sendbuf[i] = drand48();
            if (sendbuf[i] > trueMax) {
                trueMax = sendbuf[i];
                trueMaxLoc = i;
            }
        }
    }

    MPI::COMM_WORLD.Scatter(&sendbuf[0], recvLength, MPI_FLOAT,
                            &recvbuf[0], recvLength, MPI_FLOAT, 0);
    triplet t = MPIvectorAllMaxloc(recvbuf);
    if (RANK == 0)
        ASSERT_EQUAL(trueMax, t.value);
}

int test_MPIcsAdd(const std::string& method) {
    const int MASTER = 0;
    const int SIZE = MPI::COMM_WORLD.Get_size();
    const int RANK = MPI::COMM_WORLD.Get_rank();
    const int N = 10; // the size of the matrix
    const int K = 100; // the number of entries on each rank

    // create a bunch of (i,j) pairs. We only use these directly
    // on MASTER
    vector<int> sendII(K*SIZE);
    vector<int> sendJJ(K*SIZE);
    for (int i = 0; i < K*SIZE; i++) {
        sendII[i] = rand() % N;
        sendJJ[i] = rand() % N;
    }
    vector<int> recvII(K*SIZE);
    vector<int> recvJJ(K*SIZE);
    MPI::COMM_WORLD.Scatter(&sendII[0], K, MPI_INT, &recvII[0], K, MPI_FLOAT, MASTER);
    MPI::COMM_WORLD.Scatter(&sendJJ[0], K, MPI_INT, &recvJJ[0], K, MPI_FLOAT, MASTER);

    cs* T = cs_spalloc(N, N, 1, 1, 1);    
    for (int i = 0; i < K; i++)
        cs_entry(T, recvII[i], recvJJ[i], 1.0);   
    // convert to csc
    cs* M = cs_triplet(T);
    cs_dupl(M);
    cs_free(T);

    cs* result;
    if (method.compare("naive") == 0)
        result = MPIcsAdd(M);
    else if (method.compare("efficient") == 0)
        result = MPIcsAdd_efficient(M);

    if (RANK == MASTER) {
        // compute the expected result by adding
        // entries for all of the data that was scattered
        // to separate nodes
        cs* E = cs_spalloc(N, N, 1, 1, 1);    
        for (int i = 0; i < K*SIZE; i++)
            cs_entry(E, sendII[i], sendJJ[i], 1.0);   
        cs* expected  = cs_triplet(E);
        cs_dupl(expected);
        cs_free(E);

        printf("THese should be equal\n");
        cs_print(result, 1);
        cs_print(expected, 1);

        ASSERT_EQUAL(expected->nzmax, result->nzmax);
        ASSERT_EQUAL(expected->m, result->m);
        ASSERT_EQUAL(expected->n, result->n);
        for (int i = 0; i < result->nzmax; i++) {
            ASSERT_EQUAL((expected->i)[i], (result->i)[i]);
            ASSERT_EQUAL((expected->x)[i], (result->x)[i]);
        }
    }
}


int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    const int RANK = MPI::COMM_WORLD.Get_rank();
    const int SIZE = MPI::COMM_WORLD.Get_size();
    int status = 0;

    try {
        test_MPIvectorAllMaxloc();
        test_MPIcsAdd("naive");
        test_MPIcsAdd("efficient");
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        status = 1;
    }

    vector<int> statuses(SIZE);
    MPI::COMM_WORLD.Allgather(&status, 1, MPI_INT, &statuses[0], 1, MPI_INT);
    for (int i = 0; i < SIZE; i++)
        if (statuses[i] != 0)
            status = statuses[i];
    if (status == 0 && RANK == 0)
        cout << "Done" << endl;
    MPI::Finalize();
    return status;
}
