#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <utility>
#include <mpi.h>

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

  if (rank == 0) {
    for (int i = 0; i < numFrames; i++)
      printf("g[%d]=%f\n", i, traces[i]);

    printf("coordinates[0,0]=%f %f %f\n",
	   coordinates[0],
	   coordinates[1],
	   coordinates[2]);
    printf("coordinates[1,0]=%f %f %f\n",
	   coordinates[numPaddedAtoms*3 + 0],
	   coordinates[numPaddedAtoms*3 + 1],
	   coordinates[numPaddedAtoms*3 + 2]);
    

  }

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

std::vector<float> ParallelKCenters::getRmsdsTo(std::pair<int, int> &ref) {
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


  for (size_t i = 0; i < numFrames; i++)
      result[i] = sqrtf(msd_atom_major(numAtoms, numPaddedAtoms, frame, &coordinates[i*numPaddedAtoms*3], g, traces[i]));
  
  if (rank != ref.first)
    //free(frame);
  return result;
  }
