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

#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <vector>
#include "utilities.hpp"
#include "typedefs.hpp"
#include "aligned_allocator.hpp"
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"
#include "theobald_rmsd.h"

namespace Tungsten {

using std::vector;
using std::pair;

static const int MASTER = 0;
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
    atomIndices_ = vector<int>(atomIndices); // copy
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
    fill(distances.begin(), distances.end(), FLT_MAX);

    printfM("\nParallel KCenters Clustering\n");
    printfM("----------------------------\n");

    for (int i = 0; true; i++) {
        triplet max = maxLocAllReduce(distances);
        printfM("Finishing when %.4f < %.4f.    ", max.value, rmsdCutoff);
        if (max.value < rmsdCutoff)
            break;

        gindex newCenter = {max.rank,  max.frame};
        printfM("Found new center (%d, %d)\n", newCenter.rank, newCenter.frame);

        vector<float> newDistances = getRmsdsFrom(newCenter.rank, newCenter.frame);
        for (int j = 0; j < numFrames_; j++)
            if (newDistances[j] < distances[j]) {
                distances[j] = newDistances[j];
                assignments_[j] = newCenter;
            }
        centers_.push_back(newCenter);

        MPI::COMM_WORLD.Barrier();
        printfM("RMSDs to (%d, %d)\n", newCenter.rank, newCenter.frame);
        printMPIVector(distances);
        fflush(stdout);
    }

    printfM("Located k=%lu clusters\n\n", centers_.size());
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


vector<float> ParallelKCenters::getRmsdsFrom(int rank, int index) const {
    if (rank >= size_)
        exitWithMessage("IndexError: No such rank.");
    if (index >= numFrames_)
        exitWithMessage("IndexError: No such frame.");
    int err = 0;
    float g = traces_[index];
    float* frame;

    // Broadcast the frame of interest to all of the nodes
    if (rank_ != rank) {
        err = posix_memalign(reinterpret_cast<void**>(&frame), 16,
                             numPaddedAtoms_*3*sizeof(float));
        if (err != 0) exitWithMessage("Malloc error");
    } else
        frame = const_cast<float*>(&coordinates_[index*numPaddedAtoms_*3]);


    MPI::COMM_WORLD.Bcast(frame, numPaddedAtoms_*3, MPI_FLOAT, rank);
    MPI::COMM_WORLD.Bcast(&g, 1, MPI_FLOAT, rank);

    vector<float> result(numFrames_);
    #pragma omp for
    for (int i = 0; i < numFrames_; i++)
        result[i] = sqrtf(msd_axis_major(numAtoms_, numPaddedAtoms_, numPaddedAtoms_,
                                         frame, &coordinates_[i*numPaddedAtoms_*3],
                                         g, traces_[i]));

    if (rank_ != rank)
        free(frame);
    return result;
}

}
