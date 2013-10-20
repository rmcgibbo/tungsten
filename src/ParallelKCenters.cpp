#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <limits>
#include <utility>
#include <mpi.h>
#include <omp.h>

#include "utilities.hpp"
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"
#include "theobald_rmsd.h"

#define MASTER 0

ParallelKCenters::ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride) : stride(stride) {
  int rank = MPI::COMM_WORLD.Get_rank();
  numFrames = ncTraj.getNumFrames();
  numAtoms = ncTraj.getNumAtoms();
  numPaddedAtoms = ((numAtoms + 3) / 4) * 4;
  numCoordinates = numFrames * numPaddedAtoms * 3;

  int err1 = posix_memalign((void**) &coordinates, 16, numCoordinates*sizeof(float));
  int err2 = posix_memalign((void**) &traces, 16, numFrames*sizeof(float));
  if (err1 != 0 || err2 != 0)
    exitWithMessage("Malloc error");

  memset(coordinates, 0, numCoordinates*sizeof(float));
  ncTraj.readPositions(1, coordinates);
  center();
  computeTraces();
}


void ParallelKCenters::center() {
  for (size_t i = 0; i < numFrames; i++) {
    double center[] = {0, 0, 0};
    float* frame = coordinates + i*numPaddedAtoms*3;
    for (size_t j = 0; j < numAtoms; j++)
      for (size_t k = 0; k < 3; k++)
	center[k] += frame[j*3 + k]; 
    
    for (size_t k = 0; k < 3; k++)
      center[k] /= numAtoms;
    for (size_t j = 0; j < numAtoms; j++)
      for (size_t k = 0; k < 3; k++)
	frame[j*3+k] -= center[k];
  }
}

void ParallelKCenters::computeTraces() {
  for (size_t i = 0; i < numFrames; i++) {
    double trace = 0;
    float* frame = coordinates + i*numPaddedAtoms*3;
    for (size_t j = 0; j < numAtoms; j++)
      for (size_t k = 0; k < 3; k++)
	trace += frame[j*3+k] * frame[j*3+k];
    traces[i] = trace;
  }
}

std::pair<std::pair<int, size_t>, float> ParallelKCenters::collectFarthestFrom(std::pair<int, size_t> &ref) {
  static const int size = MPI::COMM_WORLD.Get_size();
  static const int rank = MPI::COMM_WORLD.Get_rank();
  float* rootRmsds;
  int* rootFarthestFrame;
  float maxRmsd = 0;
  int farthestNode;
  size_t farthestFrame = -1;

  std::vector<float> rmsds = getRmsdsFrom(ref);

  // Compute the farthest point from the ref on each MPI rank
  for (size_t i = 0; i < rmsds.size(); i++)
    if (rmsds[i] > maxRmsd) {
      maxRmsd = rmsds[i];
      farthestFrame = i;
    }

  if (rank == MASTER) {
    int err1 = posix_memalign((void**) &rootRmsds, 16, size*sizeof(float));
    int err2 = posix_memalign((void**) &rootFarthestFrame, 16, size*sizeof(int));
    if (err1 != 0 || err2 != 0) exitWithMessage("Malloc error");
  }

  // Gather the best points on each node on the root.
  MPI::COMM_WORLD.Gather(&maxRmsd, 1, MPI_FLOAT, rootRmsds, size, MPI_FLOAT, MASTER);
  MPI::COMM_WORLD.Gather(&farthestFrame, 1, MPI_INT, rootFarthestFrame, size, MPI_INT, MASTER);

  maxRmsd = 0;
  if (rank == MASTER)
    // Compute which node had the farthest away point
    for (size_t i = 0; i < size; i++)
      if (rootRmsds[i] > maxRmsd) {
	maxRmsd = rootRmsds[i];
	farthestNode = i;
	farthestFrame = rootFarthestFrame[i];
      }

  // Broacast the result to all of the nodes
  MPI::COMM_WORLD.Bcast(&farthestNode, 1, MPI_FLOAT, MASTER);
  MPI::COMM_WORLD.Bcast(&farthestFrame, 1, MPI_FLOAT, MASTER);

  std::pair<int, size_t> farthest(farthestNode, farthestFrame);
  std::pair<std::pair<int, size_t>, float> triplet(farthest, maxRmsd);

  if (rank == MASTER) {
    free(rootRmsds);
    free(rootFarthestFrame);
  }

  return triplet;
}

std::vector<float> ParallelKCenters::getRmsdsFrom(std::pair<int, size_t> &ref) {
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  int err = 0;
  float g = traces[ref.second];
  float* frame = coordinates + ref.first*numPaddedAtoms*3;
  std::vector<float> result(numFrames);

  if (ref.first >= numFrames)  
    exitWithMessage("IndexError: No such frame.");
  if (ref.second >= size)
    exitWithMessage("IndexError: No such rank.");

  // Broadcast the frame of interest to all of the nodes
  if (rank != ref.first) {
    err = posix_memalign((void**) &frame, 16, numPaddedAtoms*3*sizeof(float));
    if (err != 0) exitWithMessage("Malloc error");
  }
  MPI::COMM_WORLD.Bcast(frame, numPaddedAtoms*3, MPI_FLOAT, ref.first);
  MPI::COMM_WORLD.Bcast(&g, 1, MPI_FLOAT, ref.first);

  #pragma omp for
  for (size_t i = 0; i < numFrames; i++)
      result[i] = sqrtf(msd_atom_major(numAtoms, numPaddedAtoms, frame, &coordinates[i*numPaddedAtoms*3], g, traces[i]));
  
  if (rank != ref.first)
    free(frame);
  return result;
}
