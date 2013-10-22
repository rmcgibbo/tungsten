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

 private:
  void computeStateLabels();
  void computeTransitionCounts();

  const int rank_;
  const int size_;
  const int numStates_;
  const std::vector<std::pair<int, int > > assignments_;
  const std::vector<std::pair<int, int> > centers_;
  std::vector<int> stateLabels_;
  cs* countsMatrix_;
};

}
#endif
