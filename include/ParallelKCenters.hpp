#ifndef TUNGSTEN_PARALLELKCENTERS_H
#define TUNGSTEN_PARALLELKCENTERS_H
#include <vector>
#include <utility>
#include <malloc.h>
#include "NetCDFTrajectoryFile.hpp"


class ParallelKCenters {
public:
  ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride);
  ~ParallelKCenters() {
    free(coordinates);
  }
  std::vector<float> getRmsdsTo(std::pair<int, int> &ref);

 private:
  void center();
  void computeTraces();

  int stride;
  float* coordinates;
  float* traces;
  size_t numCoordinates;
  size_t numAtoms;
  size_t numPaddedAtoms;
  size_t numFrames;

};


#endif
