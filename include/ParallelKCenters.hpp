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
     * Run k-centers RMSD clustering in parallel over MPI. Each MPI rank is
     * associated with a single trajectory. EACH TRAJECTORY MUST BE THE
     * SAME LENGTH. This is currently assumed, but unchecked. It could
     * be circumvented with a more clever parallel reductio step.
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
     * this is good for kcenters because it ensures that these conformations
     * are never chosen as a center.
     *
     * @param rank
     * @param index
     *
     * @return
     */
    std::vector<float> getRmsdsFrom(int rank, int index) const;

    /**
     *
     *
     */
    void cluster(double rmsdCutoff, int seedRank, int seedIndex);

    std::vector<gindex> getAssignments() {
        return assignments_;
    }
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
