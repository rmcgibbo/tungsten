#include <cstdio>
#include <vector>
#include <mpi.h>
#include <malloc.h>
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"

ParallelKCenters::ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride) : stride(stride) {
  int rank = MPI::COMM_WORLD.Get_rank();
  int numCoordinates = ncTraj.getNumFrames()*ncTraj.getNumAtoms()*3;
  float *positions, *rpositions;
  positions = (float*) malloc(numCoordinates*sizeof(float));
  if (positions == NULL) {
    fprintf(stderr, "MEMORY ERROR");
  }
  ncTraj.readPositions(1, positions);



  //  MPI::COMM_WORLD.Gather(positions, numCoordinates, MPI_FLOAT, 

  free(positions);
}

