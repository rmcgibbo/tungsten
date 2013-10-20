#ifndef TUNGSTEN_PARALLELKCENTERS_H
#define TUNGSTEN_PARALLELKCENTERS_H
#include "NetCDFTrajectoryFile.hpp"


class ParallelKCenters {
public:
  ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride);

 private:
   int stride;
};


#endif
