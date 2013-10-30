#include <stdlib.h> // rand
#include <iostream>
#include <vector>
#include "assertions.hpp"
#include "csparse.h"
using namespace std;
using namespace Tungsten;

void test_csget() {
    const int N = 3;
    const int K = 5;
    cs* T = cs_spalloc(N, N, 1, 1, 1);
    vector<vector<double> > dense(N, vector<double>(N));
   
    for (int i = 0; i < K; i++) {
        int row = rand() % N;
        int col = rand() % N;
        double entry = drand48();
        cs_entry(T, row, col, entry);
        dense[row][col] += entry;
    }
    cs* M = cs_triplet(T);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double X1 = cs_get(T, i, j);
            double X2 = cs_get(M, i, j);
            
            ASSERT_EQUAL(dense[i][j], X1);
            ASSERT_EQUAL(dense[i][j], X2);
        }
    }
    cs_free(T);
    cs_free(M);
}


int main(int argc, char* argv[]) {
    try {
        test_csget();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
