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

#ifndef TUNGSTEN_PARALLELKCENTERS_H
#define TUNGSTEN_PARALLELKCENTERS_H
#include <vector>
#include "typedefs.hpp"
#include "aligned_allocator.hpp"
#include "NetCDFTrajectoryFile.hpp"
namespace Tungsten {


class ParallelKCenters {
public:
    /**
     * Initialize k-centers RMSD clustering in parallel over MPI. Each MPI
     * rank is associated with a single trajectory.
     *
     * @param ncTraj        The trajectory that that this rank loads from
     * @param stride        Load only every stride-th frame, to save memory
     * @param atomIndices   The indices of the atoms to load and use for the
     *                      RMSD calculation.
     */
    ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride,
                     const std::vector<int>& atomIndices);

    /**
     * Compute the RMSDs from one conformation that may be located on ANY
     * MPI rank to all of the conformations on THIS mpi rank. If any of
     * the positions contain NANs in them, the distance is reported as zero.
     * (This is good for kcenters because it ensures that these conformations
     * are never chosen as a center.)
     *
     * @param rank    The index of the node (MPI rank) on which the reference
     *                conformation resides
     * @param index   The index of the reference conformation on its node.
     *
     * @return        A vector of distances, of length equal to the number of
     *                frames owned by THIS node, giving the distance from the
     *                reference conformation to all of THIS node's conformations.
     */
    std::vector<float> getRmsdsFrom(int rank, int index) const;

    /**
     * Run k-centers clustering
     *
     * @param rmsdCutoff    New state will be added until the distance
     *                      from every conformation to its assigned cluster
     *                      center is less than this cutoff.
     * @param seedRank      The zeroth cluster must be picked arbirarily for
     *                      k-centers. seedRank and seedIndex specify it.
     * @param seedIndex     The zeroth cluster must be picked arbirarily for
     *                      k-centers. seedRank and seedIndex specify it.
     */
    void cluster(double rmsdCutoff, int seedRank, int seedIndex);

    /**
     * Get the assignments for each conformation, the global index of the
     * conformation which serves as the center of the cluster to which each
     * of the frames owned by this node are assigned
     *
     */
    std::vector<gindex> getAssignments() {
        return assignments_;
    }

    /**
     * Get the global index of the conformations which are serving as cluster
     * centers. These are the unique elements of getAssignments(), which the
     * cavait that the uniqueness is over all of the MPI ranks. A given rank
     * might not have any assignments to one of the clusters, but that
     * cluster will still be listed in the output of getCenters(), which
     * is identical on every MPI rank.
     */
    std::vector<gindex> getCenters() {
        return centers_;
    }


private:
    void centerCoordinates();
    void computeTraces();
    std::vector<int> atomIndices_;

    int stride_;
    std::vector<float, aligned_allocator<float, 4*sizeof(float)> > coordinates_;
    std::vector<float> traces_;

    std::vector<gindex> assignments_;
    std::vector<gindex> centers_;

    const int rank_;
    const int size_;
    size_t numCoordinates_;
    size_t numAtoms_;
    size_t numPaddedAtoms_;
    size_t numFrames_;

};

}
#endif
