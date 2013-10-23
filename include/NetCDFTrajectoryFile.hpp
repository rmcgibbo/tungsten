#ifndef TUNGSTEN_NETCDFTRAJECTORY_FILE_H_
#define TUNGSTEN_NETCDFTRAJECTORY_FILE_H_
#include <stddef.h>
#include <vector>
#include <string>
#include <netcdfcpp.h>
#include "OpenMM.h"
#include "typedefs.hpp"
#include "aligned_allocator.hpp"
namespace Tungsten {

class NetCDFTrajectoryFile {
public:
    NetCDFTrajectoryFile(const std::string& filename, const std::string& mode, int numAtoms);
    ~NetCDFTrajectoryFile(void) {
        delete handle_;
    }

    /*
     * Write an OpenMM state to disk. If any nonfinite numbers are detected
     * in the positions, the positions are written with all NANs.
     */
    int write(OpenMM::State state);


    /*
     * Fetch the `index`-th frame of positions from the trajectory
     * residing on the `rank`th MPI rank.
     *
     * @param rank    the MPI rank of the process to load from
     * @param index   the index of the frame to load
     */
    PositionsAndPeriodicBox loadNonlocalStateMPI(int rank, int index);

    /*
     * Read an entire trajectory from disk, returning it into an aligned vector.
     *
     * @tparam N            Template parameter for the alignment you request. To
     *                      get 4-float (16 byte) aligned memory, use N=4.
     * @param stride        Load only every `stride`-th frame from the tajectory,
     *                      skipping the others
     * @param atomIndices   Only these (zero-based) indices are loaded from disk.
     *
     * @returns             An aligned vector of size numFrames*3*numPaddedAtoms,
     *                      which is implicitly a 3-dimensional array of shape
     *                      (numFrames, 3, numPaddedAtoms). numPaddedAtoms is
     *                      equal to the size of atomIndices rounded UP to the
     *                      the nearest N, so that each frame is propertly
     *                      aligned. The coordinates of the "padding" atoms are
     *                      set to (0, 0, 0).
     */
    template<std::size_t N> std::vector<float, aligned_allocator<float, N*sizeof(float)> >
    loadAllAxisMajorPositions(int stride, const std::vector<int>& atomIndices) const {
        // total number of frames in the dataset
        int numTotalFrames = handle_->get_dim("frame")->size();
        // number of frames that we're actually going to read
        int numFrames = (numTotalFrames + stride - 1) / stride;
        // number of atoms in the output dimension
        int numPaddedAtoms = getNumPaddedAtoms<N>(atomIndices);
        NcVar* coord = handle_->get_var("coordinates");
        std::vector<float> frame(numAtoms_*3);
        std::vector<float, aligned_allocator<float, N*sizeof(float)> > out(numFrames*3*numPaddedAtoms);

        int ii = 0;
        for (int i = 0; i < numTotalFrames; i += stride, ii++) {
            // i  is the index of the frame to read off the disk,
            // ii is the index int `out` where we want to put it
            coord->set_cur(i, 0, 0);
            coord->get(&frame[0], 1, numAtoms_, 3);
            for (int jj = 0; jj < atomIndices.size(); jj++) {
                int j = atomIndices[jj];
                // j is the index of the atom on disk
                // jj is the index of the atom in atomIndices

                for (int k = 0; k < 3; k++) {
                    float v = frame[j*3 + k] / 10.0;  // angstroms to nm
                    out[ii*numPaddedAtoms*3 + k*numPaddedAtoms + jj] = v;
                }
            }
        }
        return out;
    }

    /*
     * Compute the number of atoms, including padding, which will be in the
     * results of loadAllAxisMajorPositions
     *
     * @param N             Template parameter for the alignment you request. To
     *                      get 4-float (16 byte) aligned memory, use N=4.
     * @param atomIndices   The indices of the atoms you're loading.
     */
    template<std::size_t N> inline
    int getNumPaddedAtoms(const std::vector<int>& atomIndices) const {
        return ((atomIndices.size() + N - 1) / N) * N;
    }

    /* 
     * Get the number of atoms stored in this trajectory
     */ 
    int getNumAtoms() const {
        return numAtoms_;
    }
    
    /*
     * Get the number of frames/snapshots stored in this trajectory
     */ 
    int getNumFrames() const {
        if (handle_ == NULL)
            return 0;
        return handle_->get_dim("frame")->size();
    }
    
    /*
     * Flush the trajectory, writing any contents in internal buffers to disk
     */
    void flush(void) {
        if (handle_ != NULL)
            handle_->sync();
    }

private:
    const int rank_;
    const int size_;
    int numAtoms_;
    NcFile* handle_;
    const std::string mode_;

    int initializeHeaders(void);
};

}
#endif
