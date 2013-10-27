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
// along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013 Stanford University
// Author: Robert McGibbon
// Contributors:

#include "mpi.h"
#include "omp.h"
#include <math.h>
#include <cstdlib>
#include <limits>
#include <algorithm>    // std::copy
#include <vector>
#include "config.h"
#include "typedefs.hpp"
#include "utilities.hpp"
#include "aligned_allocator.hpp"
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"
#include "theobald_rmsd.h"
#if defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#include <sys/time.h>
#else
#include <time.h>
#endif

namespace Tungsten {

using std::vector;
using std::pair;
using std::copy;

static const int MASTER = 0;
static const int MAX_KCENTERS_LINES = 10;

typedef struct {
    int rank;
    int frame;
    float value;
} triplet;

/* *
 * MPI Parallel Reduction. Each rank provides input data, and the return
 * value, on each node, is a triplet containing the rank, index and value of
 * the maximum entry. It's a global argmax.
 */
static triplet maxLocAllReduce(const vector<float>& input) {
    static const int rank = MPI::COMM_WORLD.Get_rank();
    struct {
        float value;
        int   index;
    } localMaxLoc, globalMaxLoc;
    int trajLength = input.size();
    int maxTrajLength;

    // Share the longest trajectory length with all of the nodes so that
    // each can properly compute a unique index of
    // localIndex + rank*maxLengthOfAnyLocalIndex
    MPI::COMM_WORLD.Allreduce(&trajLength, &maxTrajLength, 1, MPI_INT, MPI_MAX);

    // local maxloc
    localMaxLoc.value = input[0];
    localMaxLoc.index = 0;
    for (int i = 1; i < trajLength; i++)
        if (localMaxLoc.value < input[i]) {
            localMaxLoc.value = input[i];
            localMaxLoc.index = i;
        }

    // give the local maxloc a globaly-resolvably index
    localMaxLoc.index = localMaxLoc.index + rank * maxTrajLength;

    // global maxloc
    MPI::COMM_WORLD.Allreduce(&localMaxLoc, &globalMaxLoc, 1, MPI_FLOAT_INT, MPI_MAXLOC);
    int outRank = globalMaxLoc.index / maxTrajLength;
    int outFrame = globalMaxLoc.index % maxTrajLength;
    float outValue = globalMaxLoc.value;

    triplet t = {outRank, outFrame, outValue};
    return t;
}


ParallelKCenters::ParallelKCenters(const NetCDFTrajectoryFile& ncTraj,
                                   int stride, const vector<int>& atomIndices)
    : stride_(stride)
    , rank_(MPI::COMM_WORLD.Get_rank())
    , size_(MPI::COMM_WORLD.Get_size()) {
    atomIndices_ = vector<int>(atomIndices);  // copy
    // If atomIndices is empty, we use ALL of the atoms in the trajectory
    if (atomIndices_.size() == 0)
        for (int i = 0; i < ncTraj.getNumAtoms(); i++)
            atomIndices_.push_back(i);
    numPaddedAtoms_ = ncTraj.getNumPaddedAtoms<4>(atomIndices);
    numAtoms_ = atomIndices_.size();
    numFrames_ = ncTraj.getNumFrames();
    numCoordinates_ = numFrames_ * numPaddedAtoms_ * 3;

    traces_.resize(numFrames_);
    coordinates_ = ncTraj.loadAllAxisMajorPositions<4>(stride, atomIndices);

    centerCoordinates();
    computeTraces();
}

void ParallelKCenters::cluster(double rmsdCutoff, int seedRank, int seedIndex) {
    gindex newCenter = {seedRank, seedIndex};
    vector<float> distances(numFrames_);
    assignments_.resize(numFrames_);
    centers_.resize(0);
    fill(distances.begin(), distances.end(), std::numeric_limits<float>::max());

    printfM("\nParallel KCenters Clustering\n");
    printfM("----------------------------\n");
#if defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
    struct timeval startTime;
    gettimeofday(&startTime, NULL);
#else
    time_t startTime = (NULL);
#endif

    for (int i = 0; true; i++) {
        triplet max = maxLocAllReduce(distances);

        if (i > 0 && i < MAX_KCENTERS_LINES)
            // don't print when the
            printfM("Finishing when %.4f nm < %.4f nm   ", max.value, rmsdCutoff);
        if (max.value < rmsdCutoff)
            break;

        if (i > 0 && i < MAX_KCENTERS_LINES)
            printfM("Found new center (%d, %d)\n", max.rank, max.frame);
        if (i == MAX_KCENTERS_LINES)
            printfM("... [truncating output] ...\n\n");


        gindex newCenter = {max.rank,  max.frame};
        vector<float> newDistances = getRmsdsFrom(newCenter.rank, newCenter.frame);
        for (int j = 0; j < numFrames_; j++)
            if (newDistances[j] < distances[j]) {
                distances[j] = newDistances[j];
                assignments_[j] = newCenter;
            }
        centers_.push_back(newCenter);
        MPI::COMM_WORLD.Barrier();

        #ifdef VERBOSE
        MPI::COMM_WORLD.Barrier();
        printfM("RMSDs to (%d, %d)\n", newCenter.rank, newCenter.frame);
        printMPIVector(distances);
        fflush(stdout);
        #endif
    }
#if defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
    struct timeval endTime;
    gettimeofday(&endTime, NULL);
    long long elapsedLL = (endTime.tv_sec - startTime.tv_sec)*1000000LL + (endTime.tv_usec-startTime.tv_usec);
    double elapsed = elapsed / 1000000.0;
#else
    time_t endTime = time(NULL);
    double elapsed = difftime(endTime, startTime);
#endif
    printfM("LOCATED k=%lu states (t=%.2f s)\n", centers_.size(), elapsed / 1000000.0);

}


void ParallelKCenters::centerCoordinates() {
    for (int i = 0; i < numFrames_; i++) {
        double center[] = {0, 0, 0};
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < numAtoms_; k++)
                center[j] += coordinates_[i*3*numPaddedAtoms_ + j*numPaddedAtoms_ + k];

        for (int j = 0; j < 3; j++)
            center[j] /= numAtoms_;
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < numAtoms_; k++)
                coordinates_[i*3*numPaddedAtoms_ + j*numPaddedAtoms_ + k] -= center[j];
    }
}

void ParallelKCenters::computeTraces() {
    for (int i = 0; i < numFrames_; i++) {
        double trace = 0;
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < numAtoms_; k++)
                trace += coordinates_[i*3*numPaddedAtoms_ + j*numPaddedAtoms_ + k] *
                         coordinates_[i*3*numPaddedAtoms_ + j*numPaddedAtoms_ + k];
        traces_[i] = trace;
    }
}


vector<float> ParallelKCenters::getRmsdsFrom(int targetRank, int targetIndex) const {
    if (targetRank >= size_)
        exitWithMessage("IndexError: No such rank.");
    if (rank_ == targetRank && targetIndex >= numFrames_)
        exitWithMessage("IndexError: No such rank.");

    int err = 0;
    float g = 0;
    // we just want to use the same allocator as coordinates_, but you
    // can't call get_allocator() because the template needs to be
    // compile time constaint.
    vector<float, aligned_allocator<float, 4*sizeof(float)> > frame(numPaddedAtoms_*3);

    if (rank_ == targetRank) {
        // copy appropriate coordinate into the send buffer
        std::copy(&coordinates_[(targetIndex)*numPaddedAtoms_*3],
                  &coordinates_[(targetIndex+1)*numPaddedAtoms_*3],
                  frame.begin());
        g = traces_[targetIndex];
    }

    MPI::COMM_WORLD.Bcast(&frame[0], numPaddedAtoms_*3, MPI_FLOAT, targetRank);
    MPI::COMM_WORLD.Bcast(&g, 1, MPI_FLOAT, targetRank);

    vector<float> result(numFrames_);
    #pragma omp for
    for (int i = 0; i < numFrames_; i++)
        result[i] = sqrtf(msd_axis_major(numAtoms_, numPaddedAtoms_, numPaddedAtoms_,
                                         &frame[0], &coordinates_[i*numPaddedAtoms_*3],
                                         g, traces_[i]));

    return result;
}

}  // namespace Tungsten
