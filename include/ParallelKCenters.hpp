#ifndef TUNGSTEN_PARALLELKCENTERS_H
#define TUNGSTEN_PARALLELKCENTERS_H
#include <vector>
#include <utility>
#include <malloc.h>
#include "NetCDFTrajectoryFile.hpp"
namespace Tungsten {


class ParallelKCenters {
public:
  ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride,
		   const std::vector<int>& atomIndices);
  ~ParallelKCenters() {}
  std::vector<float> getRmsdsFrom(const std::pair<int, int> &ref) const;
  void cluster(float rmsdCutoff, const std::pair<int, int>& seed);


 private:
  void center();
  void computeTraces();
  std::vector<int> atomIndices;

  int stride;
  fvector16 coordinates;
  fvector16 traces;

  const int rank;
  const int size;
  size_t numCoordinates;
  size_t numAtoms;
  size_t numPaddedAtoms;
  size_t numFrames;

};

}
#endif
