//////////////////////////////////////////////////////////////////////////
// This file is part of Tungsten
//
// Tungsten is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 2.1 of the License, or
// (at your option) any later version.
//
// Tungsten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Tungsten. If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013 Stanford University
// Author: Robert McGibbon
// Contributors:

#include "mpi.h"
#include "stdlib.h"
#include <algorithm>
#include <vector>
#include <utility>
#include <numeric>
#include "config.h"
#include "csparse.h"
#include "typedefs.hpp"
#include "utilities.hpp"
#include "ParallelMSM.hpp"
#if defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#include <sys/time.h>
#else
#include <time.h>
#endif

namespace Tungsten {

using std::vector;
using std::pair;
using std::map;

static const int MASTER = 0;

// Utilities
template <class T1, class T2> static bool pairComparatorSecond(const pair<T1, T2>& lhs, const pair<T1, T2>& rhs) {
    return lhs.second < rhs.second;
}

#ifdef __GNUC__
  /* Test for GCC > 2.95 */
  #if __GNUC__ > 2 || (__GNUC__ == 2 && (__GNUC_MINOR__ > 95))
    #define likely(x)   __builtin_expect(!!(x), 1)
    #define unlikely(x) __builtin_expect(!!(x), 0)
  #else /* __GNUC__ > 2 ... */
    #define likely(x)   (x)
    #define unlikely(x) (x)
  #endif /* __GNUC__ > 2 ... */
#else /* __GNUC__ */
  #define likely(x)   (x)
  #define unlikely(x) (x)
#endif /* __GNUC__ */


// Class implementation

ParallelMSM::ParallelMSM(const vector<gindex>& assignments,
                         const vector<gindex>& centers,
                         const vector<float>& time)
    : assignments_(assignments)
    , numStates_(centers.size())
    , centers_(centers)
    , time_(time)
    , rank_(MPI::COMM_WORLD.Get_rank())
    , size_(MPI::COMM_WORLD.Get_size())
    , countsMatrix_(NULL)
{
    // build MSM
    computeStateLabels();
    computeTransitionCounts();
}


void ParallelMSM::computeStateLabels() {
    invLabelMap_.resize(centers_.size());
    for (int i = 0; i < centers_.size(); i++) {
        labelMap_[centers_[i]] = i;
        invLabelMap_[i] = centers_[i];
    }

    stateLabels_.resize(assignments_.size());
    for (int i = 0; i < assignments_.size(); i++)
        stateLabels_[i] = labelMap_[assignments_[i]];
}



gindex ParallelMSM::scatterMinCountStates() {
    // this array is computed on master and then scattered
    // to the ranks -- totally unused on other nodes
    vector<int> minCountLabels(size_);

    if (rank_ == MASTER) {
        // compute the number of inbound transition into
        // each state. these are the row sums of the counts
        // matrix, which has the "to" state as the first index
        // and the "from" index as the second index.
        vector<double> ones(numStates_, 1.0);
        vector<double> rowSums(numStates_, 0.0);
        cs_gaxpy(countsMatrix_, &ones[0], &rowSums[0]);

        //for (int i = 0; i < numStates_; i++)
        //printf("rowSums[%d]=%f\n", i, rowSums[i]);

        vector<pair<int, double> > rowSumsWithIndex(numStates_);
        for (int i = 0; i < numStates_; i++) {
            rowSumsWithIndex[i] = pair<int, double>(i, rowSums[i]);
        }
        std::sort(rowSumsWithIndex.begin(), rowSumsWithIndex.end(), pairComparatorSecond<int, double>);

        // for (int i = 0; i < numStates_; i++)
        // printf("minCountLabels[%d]=%d\n", i, rowSumsWithIndex[i].first);

        for (int i = 0; i < size_; i++)
            minCountLabels[i] = rowSumsWithIndex[i % numStates_].first;
    }

    int myNewLabel;
    MPI::COMM_WORLD.Scatter(&minCountLabels[0], 1, MPI_INT, &myNewLabel, 1, MPI_INT, MASTER);
    return invLabelMap_[myNewLabel];
}


void ParallelMSM::computeTransitionCounts() {
#if defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
    struct timeval startTime;
    gettimeofday(&startTime, NULL);
#else
    time_t startTime = (NULL);
#endif

    if (stateLabels_.size() == 0)
        exitWithMessage("NO ASSIGNMENTS!");

    cs* T = cs_spalloc(numStates_, numStates_, 1, 1, 1);
    // Add the pairs
    for (int i = 0; i < stateLabels_.size()-1; i++)
        // in tungsten, each of the trajectory files contain
        // disjoint segments, since each adaptive sampling
        // round is written to the same trajectory file. we dont
        // want to spurriously mark counts between the end of one
        // trajectory and the beginning of the other
        if (likely(time_[i+1] > time_[i]))
            cs_entry(T, stateLabels_[i+1], stateLabels_[i], 1.0);

    // Convert to CSC form and sum duplicace entries
    countsMatrix_ = cs_triplet(T);
    cs_dupl(countsMatrix_);
    cs_free(T);

    // Gather nzmax on root, the maximum number of entries
    // each each rank's countsMatrix_
    vector<int> rootNzmax(size_);
    MPI::COMM_WORLD.Gather(&countsMatrix_->nzmax, 1, MPI_INT, &rootNzmax[0], 1, MPI_INT, MASTER);

    cs* newCounts;
    if (rank_ != MASTER) {
        // All of the slave nodes send their buffers to to MASTER
        // for accumulation
        MPI::COMM_WORLD.Isend(countsMatrix_->p, countsMatrix_->n+1, MPI_INT, MASTER, 0);
        MPI::COMM_WORLD.Isend(countsMatrix_->i, countsMatrix_->nzmax, MPI_INT, MASTER, 1);
        MPI::COMM_WORLD.Isend(countsMatrix_->x, countsMatrix_->nzmax, MPI_DOUBLE, MASTER, 2);
    } else {
        for (int j = 1; j < size_; j++) {
            // The master node receives these entries and uses them to
            // reconstruct a sparse matrix, using cs_add to then
            // add it to its own.
            vector<int> p(numStates_+1);
            vector<int> i(rootNzmax[j]);
            vector<double> x(rootNzmax[j]);
            MPI::Request rP, rI, rX;
            rP = MPI::COMM_WORLD.Irecv(&p[0], numStates_+1, MPI_INT, j, 0);
            rI = MPI::COMM_WORLD.Irecv(&i[0], rootNzmax[j], MPI_INT, j, 1);
            rX = MPI::COMM_WORLD.Irecv(&x[0], rootNzmax[j], MPI_DOUBLE, j, 2);
            rP.Wait();
            rI.Wait();
            rX.Wait();

            // place this data in a struct
            cs M;
            M.nzmax = rootNzmax[j];
            M.m = numStates_;
            M.n = numStates_;
            M.p = &p[0];
            M.i = &i[0];
            M.x = &x[0];
            M.nz = -1;

            newCounts = cs_add(countsMatrix_, &M, 1.0, 1.0);
            cs_free(countsMatrix_);
            countsMatrix_ = newCounts;
        }
    }

#if defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
    struct timeval endTime;
    gettimeofday(&endTime, NULL);
    long long elapsedLL = (endTime.tv_sec - startTime.tv_sec)*1000000LL + (endTime.tv_usec-startTime.tv_usec);
    double elapsed = (double) elapsed / 1000000.0;
#else
    time_t endTime = time(NULL);
    double elapsed = difftime(endTime, startTime);
#endif
    printfM("Transition Counts Matrix Built (t=%.2f s)\n\n", elapsed / 1000000.0);

}

ParallelMSM::~ParallelMSM() {
    if (countsMatrix_ != NULL) {
        cs_free(countsMatrix_);
    }
}

}  // namespace
