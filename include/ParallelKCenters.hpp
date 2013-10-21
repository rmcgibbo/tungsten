#ifndef TUNGSTEN_PARALLELKCENTERS_H
#define TUNGSTEN_PARALLELKCENTERS_H
#include <vector>
#include <utility>
#include <malloc.h>
#include "NetCDFTrajectoryFile.hpp"


class ParallelKCenters {
public:
  ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride,
		   const std::vector<int>& atomIndices);
  ~ParallelKCenters() {
    free(coordinates);
  }
  std::vector<float> getRmsdsFrom(const std::pair<int, size_t> &ref);
  void cluster(float rmsdCutoff, const std::pair<int, size_t>& seed);


 private:
  void center();
  void computeTraces();
  std::vector<int> atomIndices;

  int stride;
  float* coordinates;
  float* traces;
  size_t numCoordinates;
  size_t numAtoms;
  size_t numPaddedAtoms;
  size_t numFrames;

};


#endif
