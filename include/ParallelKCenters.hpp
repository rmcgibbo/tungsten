#ifndef TUNGSTEN_PARALLELKCENTERS_H
#define TUNGSTEN_PARALLELKCENTERS_H
#include <malloc.h>
#include "NetCDFTrajectoryFile.hpp"


class ParallelKCenters {
public:
  ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride);
  ~ParallelKCenters() {
    free(coordinates);
  }

 private:
  int stride;
  float* coordinates;
  size_t numCoordinates;
};


#endif
