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

using std::vector;
using std::pair;

#define MASTER 0
typedef pair<pair<int, size_t>, float> triplet;

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

  triplet t(pair<int, size_t>(outRank, outIndex), outValue);
  return t;
}


ParallelKCenters::ParallelKCenters(const NetCDFTrajectoryFile& ncTraj,
         int stride, const vector<int>& atomIndices_) : stride(stride) {
  static const int rank = MPI::COMM_WORLD.Get_rank();
  atomIndices = vector<int>(atomIndices_); // copy
  // If atomIndices is empty, we use ALL of the atoms in the trajectory
  if (atomIndices.size() == 0)
    for (int i = 0; i < ncTraj.getNumAtoms(); i++)
      atomIndices.push_back(i);
  numAtoms = atomIndices.size();
  numPaddedAtoms = ((numAtoms + 3) / 4) * 4;
  numFrames = ncTraj.getNumFrames();
  numCoordinates = numFrames * numPaddedAtoms * 3;

  int err1 = posix_memalign((void**) &coordinates, 16, numCoordinates*sizeof(float));
  int err2 = posix_memalign((void**) &traces, 16, numFrames*sizeof(float));
  if (err1 != 0 || err2 != 0)
    exitWithMessage("Malloc error");

  memset(coordinates, 0, numCoordinates*sizeof(float));
  ncTraj.readAxisMajorPositions(1, atomIndices, 4, coordinates);

  if (rank == MASTER) {
    printf("Coordinates\n");
    for (int j = 0; j < 3; j++) {
      printf("[");
      for (int k = 0; k < numPaddedAtoms; k++)
        printf("%f   ", coordinates[j*numPaddedAtoms + k]);
      printf("]\n");
    }
  }

  center();
  computeTraces();
}

void ParallelKCenters::cluster(float rmsdCutoff, const pair<int, size_t>& seed) {
  static const int rank = MPI::COMM_WORLD.Get_rank();
  pair<int, size_t> newCenter = seed;
  vector<float> distances(numFrames);
  vector< pair<int, size_t> > assignments(numFrames);
  vector< pair<int, size_t> > centers;
  fill(distances.begin(), distances.end(), FLT_MAX);

  if (rank == MASTER) {
    printf("\nParallel KCenters Clustering\n");
    printf("============================\n");
  }

  for (int i = 0; true; i++) {
    triplet max = maxLocAllReduce(distances);
    if (rank == MASTER)
      printf("Finishing when %.4f < %.4f.    ", max.second, rmsdCutoff);
    if (max.second < rmsdCutoff)
      break;

    pair<int, size_t> newCenter = max.first;
    if (rank == MASTER)
      printf("Found new center (%d, %lu)\n", newCenter.first, newCenter.second);

    vector<float> newDistances = getRmsdsFrom(newCenter);
    for (int j = 0; j < numFrames; j++)
      if (newDistances[j] < distances[j]) {
        distances[j] = newDistances[j];
        assignments[j] = max.first;
      }
    centers.push_back(max.first);

    printVector(newDistances);
  }

  if (rank == MASTER) {
    printf("\n=============================\n");
    printf("Located k=%lu clusters\n", centers.size());
  }
}


void ParallelKCenters::center() {
  for (int i = 0; i < numFrames; i++) {
    double center[] = {0, 0, 0};
    float* frame = coordinates + i*3*numPaddedAtoms;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < numAtoms; k++)
        center[j] += frame[j*numPaddedAtoms + k];

    for (int j = 0; j < 3; j++)
      center[j] /= numAtoms;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < numAtoms; k++)
        frame[j*numPaddedAtoms + k] -= center[j];
  }
}

void ParallelKCenters::computeTraces() {
  for (int i = 0; i < numFrames; i++) {
    double trace = 0;
    float* frame = coordinates + i*3*numPaddedAtoms;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < numAtoms; k++)
        trace += frame[j*numPaddedAtoms+k] * frame[j*numPaddedAtoms+k];
    traces[i] = trace;
  }
}


vector<float> ParallelKCenters::getRmsdsFrom(const pair<int, size_t> &ref) {
  static const int rank = MPI::COMM_WORLD.Get_rank();
  static const int size = MPI::COMM_WORLD.Get_size();
  if (ref.first >= size)
    exitWithMessage("IndexError: No such rank.");
  if (ref.second >= numFrames)
    exitWithMessage("IndexError: No such frame.");
  int err = 0;
  float g = traces[ref.second];
  float* frame = coordinates + ref.second*numPaddedAtoms*3;
  vector<float> result(numFrames);


  // Broadcast the frame of interest to all of the nodes
  if (rank != ref.first) {
    err = posix_memalign((void**) &frame, 16, numPaddedAtoms*3*sizeof(float));
    if (err != 0) exitWithMessage("Malloc error");
  }
  MPI::COMM_WORLD.Bcast(frame, numPaddedAtoms*3, MPI_FLOAT, ref.first);
  MPI::COMM_WORLD.Bcast(&g, 1, MPI_FLOAT, ref.first);

  //  #pragma omp for
  for (size_t i = 0; i < numFrames; i++)
    result[i] = sqrtf(msd_axis_major(numAtoms, numPaddedAtoms, numPaddedAtoms,
                                     frame, &coordinates[i*numPaddedAtoms*3],
                                     g, traces[i]));

  if (rank != ref.first)
    free(frame);
  return result;
}
