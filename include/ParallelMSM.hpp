#ifndef TUNGSTEN_PARALLELMSM_H
#define TUNGSTEN_PARALLELMSM_H
#include <vector>
#include "csparse.h"

class ParallelMSM {
public:
  ParallelMSM(const std::vector<int>& assignments, int numStates);
  ~ParallelMSM();

 private:
  const std::vector<int> assignments;
  const int numStates;
  cs* counts;
};


#endif
