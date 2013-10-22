#include "mpi.h"
#include "stdlib.h"
#include <vector>
#include <numeric>
#include "csparse.h"
#include "utilities.hpp"
#include "ParallelMSM.hpp"
namespace Tungsten {

static const int MASTER = 0;

using std::vector;
using std::pair;
using std::map;


ParallelMSM::ParallelMSM(const std::vector<pair<int, int> >& assignments,
			 const std::vector<pair<int, int> >& centers):
  assignments_(assignments), numStates_(centers.size()), centers_(centers),
  rank_(MPI::COMM_WORLD.Get_rank()), size_(MPI::COMM_WORLD.Get_size())
{
  // build MSM
  computeStateLabels();
}

void ParallelMSM::computeStateLabels() {
  map<pair<int, int>, int> labelMap;
  for (int i = 0; i < centers_.size(); i++)
    labelMap[centers_[i]] = i;

  stateLabels_.resize(assignments_.size());
  for (int i = 0; i < assignments_.size(); i++)
    stateLabels_[i] = labelMap[assignments_[i]];
}

void ParallelMSM::computeTransitionCounts() {
    // this will probably fail if assignments is empty
  cs* T = cs_spalloc(numStates_, numStates_, 1, 1, 1);
  // Add the pairs
  for (int i = 0; i < stateLabels_.size()-1; i++)
    cs_entry(T, stateLabels_[i], stateLabels_[i+1], 1.0);

  // Convert to CSC form and sum duplicace entries
  countsMatrix_ = cs_triplet(T);
  cs_dupl(countsMatrix_);
  cs_free(T);

  // Gather nzmax on root, the maximum number of entries
  // each each rank's countsMatrix_
  vector<int> rootNzmax(size_);
  MPI::COMM_WORLD.Gather(&countsMatrix_->nzmax, 1, MPI_INT, &rootNzmax[0], 1, MPI_INT, MASTER);

  cs* newCounts;
  if (rank_ != MASTER) {
    // All of the slave nodes send their buffers to to MASTER
    // for accumulation
    MPI::COMM_WORLD.Isend(countsMatrix_->p, countsMatrix_->n+1, MPI_INT, MASTER, 0);
    MPI::COMM_WORLD.Isend(countsMatrix_->i, countsMatrix_->nzmax, MPI_INT, MASTER, 1);
    MPI::COMM_WORLD.Isend(countsMatrix_->x, countsMatrix_->nzmax, MPI_DOUBLE, MASTER, 2);
  } else {
    for (int j = 1; j < size_; j++) {
      // The master node receives these entries and uses them to
      // reconstruct a sparse matrix, using cs_add to then
      // add it to its own.
      vector<int> p(numStates_+1);
      vector<int> i(rootNzmax[j]);
      vector<double> x(rootNzmax[j]);
      MPI::Request rP, rI, rX;
      rP = MPI::COMM_WORLD.Irecv(&p[0], numStates_+1, MPI_INT, j, 0);
      rI = MPI::COMM_WORLD.Irecv(&i[0], rootNzmax[j], MPI_INT, j, 1);
      rX = MPI::COMM_WORLD.Irecv(&x[0], rootNzmax[j], MPI_DOUBLE, j, 2);
      rP.Wait();
      rI.Wait();
      rX.Wait();

      // place this data in a struct
      cs M;
      M.nzmax = rootNzmax[j];
      M.m = numStates_;
      M.n = numStates_;
      M.p = &p[0];
      M.i = &i[0];
      M.x = &x[0];
      M.nz = -1;

      newCounts = cs_add(countsMatrix_, &M, 1.0, 1.0);
      cs_free(countsMatrix_);
      countsMatrix_ = newCounts;
    }
  }

  if (rank_ == MASTER) {
    printf("\nFinal Data\n");
    cs_print(countsMatrix_, 0);
  }

}

ParallelMSM::~ParallelMSM() {
  //cs_free(countsMatrix_);
}

}
