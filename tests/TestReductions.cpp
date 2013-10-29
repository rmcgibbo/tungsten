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
    //float* a = (float*) malloc(100);
    //a[0] = 1;

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


int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    const int RANK = MPI::COMM_WORLD.Get_rank();
    const int SIZE = MPI::COMM_WORLD.Get_size();
    int status = 0;

    try {
        test_MPIvectorAllMaxloc();
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
