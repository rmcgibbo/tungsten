#include <cstdio>
#include <vector>
#include <mpi.h>
#include <malloc.h>

#include "utilities.hpp"
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"

ParallelKCenters::ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride) : stride(stride) {
  int rank = MPI::COMM_WORLD.Get_rank();

  numCoordinates = ncTraj.getNumFrames()*ncTraj.getNumAtoms()*3;
  coordinates = (float*) malloc(numCoordinates*sizeof(float));
  if (coordinates == NULL) {
    exitWithMessage("Malloc error");
  }
  ncTraj.readPositions(1, coordinates);


}

