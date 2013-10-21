#include "mpi.h"
#include <vector>
#include "csparse.h"
#include "ParallelMSM.hpp"

ParallelMSM::ParallelMSM(const std::vector<int>& assignments, int numStates=0):
  assignments(assignments), 
  numStates(numStates) {

  int length = assignments.size();
  cs* T = cs_spalloc(0, 0, 1, 1, 1);
  // Add the pairs
  for (int i = 0; i < length-1; i++)
    cs_entry(T, assignments[i], assignments[i+1], 1.0);
  // Convert to CSC form and sum duplicace entries
  counts = cs_triplet(T);
  cs_dupl(counts);
  cs_free(T);

 
  if (MPI::COMM_WORLD.Get_rank() == 0) {
    cs_print(counts, 0);
    printf("Largest column norm: %f\n", cs_norm(counts));
  }

}

ParallelMSM::~ParallelMSM() {
  cs_free(counts);
}
