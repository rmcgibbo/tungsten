//////////////////////////////////////////////////////////////////////////
// This file is part of Tungsten
//
// Tungsten is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 2.1 of the License, or
// (at your option) any later version.
//
// Tungsten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013 Stanford University
// Author: Robert McGibbon
// Contributors:

#ifndef TUNGSTEN_PARALLELMSM_H
#define TUNGSTEN_PARALLELMSM_H
#include <vector>
#include "typedefs.hpp"
#include "csparse.h"
namespace Tungsten {

class ParallelMSM {
public:

    /**
     * Construct a Markov state model in parallel over MPI. This method
     * should be called simultaniously by all of the MPI ranks.
     *
     * @param assignments   Each rank is associated with a single trajectory
     *                      and passes in the assignments of that trajctory,
     *                      which give the <rank,frame> pair of the cluster
     *                      center to which each frame is assigned.
     * @param centers       Centers, which should be the SAME on all of the
     *                      MPI ranks, is a list of all of the <rank,index>
     *                      pairs that are cluster centers.
     * @param time          The simulation time of each frame in the assignments.
     *                      Currently, this does NOT support variable timestep
     *                      integrator methods, but it does support trajectories
     *                      whose time index "restarts", indicating disjoint
     *                      segments of a trajectory
     */
    ParallelMSM(const std::vector<gindex>& assignments,
                const std::vector<gindex>& centers,
                const std::vector<float>& time);

    /**
     * Scatter the <rank,index> pair of the N states which have the fewest
     * incomming transitions. Each MPI rank will receive a single index,
     * corresponding to the new state to which it is tasked with (presumably)
     * simulating further.
     */
    gindex scatterMinCountStates();

    /**
     * Free this class.
     */
    ~ParallelMSM();


private:
    /**
     * The input argumnets for the constructor give the assignments in terms
     * of the global <rank,index> coordinates of the generator to which each
     * frame is assigned. But internally when we build a sparse matrix, we
     * want to use contiguous [0,...,numStates-1] indexing. This method
     * computes the mapping and inverse mapping between these two index systems
     * and populates the private variable `stateLabels_`, which gives the
     * "mapped" version of the constructor argument `assignments`, in
     * [0,...,numStates-1] indexing.
     */
    void computeStateLabels();

    /**
     * Compute the transition counts matrix, countsMatrix_, which
     * are stored in a csparse sparse matrix, with
     * countsMatrix_[j, i] as the number of transitions from state i
     * to state j. The "to" is the first index, and the "from" is
     * the second index.
     *
     * NOTE: ONLY THE ROOT NODE GETS THE "FULL" COUNTS MATRIX
     * (it does a reduction, getting the data from the other
     * ranks). The other nodes only have their *local* counts.
     */
    void computeTransitionCounts();

    const int rank_;  // my MPI rank
    const int size_;  // total MPI communicator size
    const int numStates_;  // number of states in the MSM
    const std::vector<gindex> assignments_;
    const std::vector<gindex> centers_;
    const std::vector<float> time_;

    // This class deals externally with assignments and states
    // as std::pair<int, int> giving the node on which the conformation
    // resides and its offset in that trajectory. But internally,
    // it deals with states indexed as a single integer from 0
    // to numStates-1. The maps going between these two index schemes,
    // both forward and backwards, are here
    std::map<gindex, int> labelMap_;
    std::vector<gindex> invLabelMap_;

    std::vector<int> stateLabels_;
    cs* countsMatrix_;
};

}
#endif
