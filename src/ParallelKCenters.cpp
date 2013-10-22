// Copyright 2013 Robert McGibbon
#include <mpi.h>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <vector>
#include <utility>
#include <cstring>

#include "utilities.hpp"
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"
#include "theobald_rmsd.h"
namespace Tungsten {

using std::vector;
using std::pair;

static const int MASTER = 0;
typedef pair<pair<int, int>, float> triplet;

triplet maxLocAllReduce(const vector<float>& input) {
  /* MPI Parallel Reduction. Each rank provides input data
   * (it must be the same length on each rank. this is a
   * constraint), and the return value, on each node, is
   * a triplet containing the rank, index and value of the
   * maximum entry. It's a global argmax.
   */
  static const int rank = MPI::COMM_WORLD.Get_rank();
  struct {
    float value;
    int   index;
  } local, global;
  int count = input.size();

  // local maxloc
  local.value = input[0];
  local.index = 0;
  for (int i = 1; i < count; i++)
    if (local.value < input[i]) {
      local.value = input[i];
      local.index = i;
    }

  // global maxloc
  local.index += rank * count;
  MPI::COMM_WORLD.Allreduce(&local, &global, 1, MPI_FLOAT_INT, MPI_MAXLOC);
  int outRank = global.index / count;
  int outIndex = global.index % count;
  float outValue = global.value;

  triplet t(pair<int, int>(outRank, outIndex), outValue);
  return t;
}


ParallelKCenters::ParallelKCenters(const NetCDFTrajectoryFile& ncTraj,
				   int stride, const vector<int>& atomIndices) :
  stride_(stride),
  rank_(MPI::COMM_WORLD.Get_rank()),
  size_(MPI::COMM_WORLD.Get_size())
{
  atomIndices_ = vector<int>(atomIndices); // copy
  // If atomIndices is empty, we use ALL of the atoms in the trajectory
  if (atomIndices_.size() == 0)
    for (int i = 0; i < ncTraj.getNumAtoms(); i++)
      atomIndices_.push_back(i);
  numAtoms_ = atomIndices_.size();
  numPaddedAtoms_ = ((numAtoms_ + 3) / 4) * 4;
  numFrames_ = ncTraj.getNumFrames();
  numCoordinates_ = numFrames_ * numPaddedAtoms_ * 3;

  traces_.resize(numFrames_);
  ncTraj.readAxisMajorPositions(1, atomIndices_, 4, coordinates_);

  if (rank_ == MASTER) {
    printf("Coordinates\n");
    for (int j = 0; j < 3; j++) {
      printf("[");
      for (int k = 0; k < numPaddedAtoms_; k++)
        printf("%f   ", coordinates_[j*numPaddedAtoms_ + k]);
      printf("]\n");
    }
  }

  centerCoordinates();
  computeTraces();
}

void ParallelKCenters::cluster(float rmsdCutoff, const pair<int, int>& seed) {
  pair<int, int> newCenter = seed;
  vector<float> distances(numFrames_);
  assignments_.resize(numFrames_);
  centers_.resize(0);
  fill(distances.begin(), distances.end(), FLT_MAX);

  if (rank_ == MASTER) {
    printf("\nParallel KCenters Clustering\n");
    printf("============================\n");
  }

  for (int i = 0; true; i++) {
    triplet max = maxLocAllReduce(distances);
    if (rank_ == MASTER)
      printf("Finishing when %.4f < %.4f.    ", max.second, rmsdCutoff);
    if (max.second < rmsdCutoff)
      break;

    pair<int, int> newCenter = max.first;
    if (rank_ == MASTER)
      printf("Found new center (%d, %d)\n", newCenter.first, newCenter.second);

    vector<float> newDistances = getRmsdsFrom(newCenter);
    for (int j = 0; j < numFrames_; j++)
      if (newDistances[j] < distances[j]) {
        distances[j] = newDistances[j];
        assignments_[j] = max.first;
      }
    centers_.push_back(max.first);

    printVector(newDistances);
  }

  if (rank_ == MASTER) {
    printf("\n=============================\n");
    printf("Located k=%lu clusters\n", centers_.size());
  }
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


vector<float> ParallelKCenters::getRmsdsFrom(const pair<int, int> &ref) const {
  if (ref.first >= size_)
    exitWithMessage("IndexError: No such rank.");
  if (ref.second >= numFrames_)
    exitWithMessage("IndexError: No such frame.");
  int err = 0;
  float g = traces_[ref.second];
  float* frame;

  // Broadcast the frame of interest to all of the nodes
  if (rank_ != ref.first) {
    err = posix_memalign((void**) &frame, 16, numPaddedAtoms_*3*sizeof(float));
    if (err != 0) exitWithMessage("Malloc error");
  } else
    frame = const_cast<float*>(&coordinates_[ref.second*numPaddedAtoms_*3]);


  MPI::COMM_WORLD.Bcast(frame, numPaddedAtoms_*3, MPI_FLOAT, ref.first);
  MPI::COMM_WORLD.Bcast(&g, 1, MPI_FLOAT, ref.first);

  vector<float> result(numFrames_);
  #pragma omp for
  for (int i = 0; i < numFrames_; i++)
    result[i] = sqrtf(msd_axis_major(numAtoms_, numPaddedAtoms_, numPaddedAtoms_,
                                     frame, &coordinates_[i*numPaddedAtoms_*3],
                                     g, traces_[i]));

  if (rank_ != ref.first)
    free(frame);
  return result;
}

}
