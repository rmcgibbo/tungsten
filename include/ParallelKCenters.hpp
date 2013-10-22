#ifndef TUNGSTEN_PARALLELKCENTERS_H
#define TUNGSTEN_PARALLELKCENTERS_H
#include <vector>
#include <utility>
#include "aligned_allocator.hpp"
#include "NetCDFTrajectoryFile.hpp"
namespace Tungsten {


class ParallelKCenters {
public:
    ParallelKCenters(const NetCDFTrajectoryFile& ncTraj, int stride,
                     const std::vector<int>& atomIndices);
    ~ParallelKCenters() {}
    std::vector<float> getRmsdsFrom(const std::pair<int, int> &ref) const;
    void cluster(double rmsdCutoff, const std::pair<int, int>& seed);
    std::vector< std::pair<int, int> > getAssignments() {
        return assignments_;
    }
    std::vector< std::pair<int, int> > getCenters() {
        return centers_;
    }


private:
    void centerCoordinates();
    void computeTraces();
    std::vector<int> atomIndices_;

    int stride_;
    std::vector<float, aligned_allocator<float, 4*sizeof(float)> > coordinates_;
    std::vector<float> traces_;

    std::vector< std::pair<int, int> > assignments_;
    std::vector< std::pair<int, int> > centers_;

    const int rank_;
    const int size_;
    size_t numCoordinates_;
    size_t numAtoms_;
    size_t numPaddedAtoms_;
    size_t numFrames_;

};

}
#endif
