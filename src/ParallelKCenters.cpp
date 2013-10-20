#include <cstdio>
#include <vector>
#include <mpi.h>
#include <malloc.h>
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"

ParallelKCenters::ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride) : stride(stride) {
  int rank = MPI::COMM_WORLD.Get_rank();

  numCoordinates = ncTraj.getNumFrames()*ncTraj.getNumAtoms()*3;
  coordinate = (float*) malloc(numCoordinates*sizeof(float));
  if (coordinates == NULL) {
    fprintf(stderr, "MEMORY ERROR");
  }
  ncTraj.readPositions(1, coordinates);


}

