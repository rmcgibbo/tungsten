#ifndef TUNGSTEN_PARALLELMSM_H
#define TUNGSTEN_PARALLELMSM_H
#include <vector>
#include <utility>
#include "csparse.h"
namespace Tungsten {

class ParallelMSM {
public:
  ParallelMSM(const std::vector<std::pair<int, int> >& assignments,
	      const std::vector<std::pair<int, int> >& centers);
  ~ParallelMSM();
  std::pair<int, int> scatterMinCountStateLabels();


private:
  void computeStateLabels();
  void computeTransitionCounts();

  const int rank_;
  const int size_;
  const int numStates_;
  const std::vector<std::pair<int, int > > assignments_;
  const std::vector<std::pair<int, int> > centers_;


  // This class deals externally with assignments and states
  // as std::pair<int, int> giving the node on which the conformation
  // resides and its offset in that trajectory. But internally,
  // it deals with states indexed as a single integer from 0
  // to numStates-1. The maps going between these two index schemes,
  // both forward and backwards, are here
  std::map<std::pair<int, int>, int> labelMap_;
  std::vector<std::pair<int, int> > invLabelMap_;

  std::vector<int> stateLabels_;
  cs* countsMatrix_;
};

}
#endif
