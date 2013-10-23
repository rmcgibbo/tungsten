// Copyright 2013 Robert McGibbon
#ifndef TUNGSTEN_NETCDFTRAJECTORY_FILE_H_
#define TUNGSTEN_NETCDFTRAJECTORY_FILE_H_
#include <stddef.h>
#include "netcdf.h"
#include <vector>
#include <string>
#include "OpenMM.h"
#include "typedefs.hpp"
#include "aligned_allocator.hpp"
namespace Tungsten {
static const int NC_INVALID = -1;

#define NC_ERR(e) {printf("Error: %s at %s, line %d\n", nc_strerror(e), __FILE__, __LINE__); exitWithMessage("");}

class NetCDFTrajectoryFile {
public:
    NetCDFTrajectoryFile(const std::string& filename, const std::string& mode, int numAtoms);
    ~NetCDFTrajectoryFile(void) {
        if (int r = nc_close(ncid_)) NC_ERR(r);
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
        size_t numTotalFrames = getNumFrames();
        int numAtoms = getNumAtoms();
        // number of frames that we're actually going to read
        int numFrames = (numTotalFrames + stride - 1) / stride;
        // number of atoms in the output dimension
        int numPaddedAtoms = getNumPaddedAtoms<N>(atomIndices);
        std::vector<float> frame(numAtoms*3);
        std::vector<float, aligned_allocator<float, N*sizeof(float)> > out(numFrames*3*numPaddedAtoms);

        int ii = 0;
        for (int i = 0; i < numTotalFrames; i += stride, ii++) {
            // i  is the index of the frame to read off the disk,
            // ii is the index int `out` where we want to put it
            //coord->set_cur(i, 0, 0);
            //coord->get(&frame[0], 1, numAtoms_, 3);
            size_t start[] = {i, 0, 0};
            size_t count[] = {1, numAtoms, 3};
            if (int r = nc_get_vara_float(ncid_, coordVar_, start, count, &frame[0])) NC_ERR(r);
            
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
    size_t getNumAtoms() const {
        if (!isvalid()) {
            printf("invalid");
        }
        size_t numAtoms;
        if (int r = nc_inq_dimlen(ncid_, atomDim_, &numAtoms)) NC_ERR(r);
        return numAtoms;
    }
    
    /*
     * Get the number of frames/snapshots stored in this trajectory
     */ 
    size_t getNumFrames() const {
        // total number of frames in the dataset
        size_t numTotalFrames;
        if (int r = nc_inq_dimlen(ncid_, frameDim_, &numTotalFrames)) NC_ERR(r);
        return numTotalFrames;
    }
    
    /**
     * Flush the trajectory, writing any contents in internal buffers to disk
     */
    void flush(void) {
        if (int r = nc_sync(ncid_)) NC_ERR(r);
    }
    
    /**
     * Check whether the netcdf file handle is valid.
     */
    bool isvalid(void) const {
        return (ncid_ != NC_INVALID);
    }

private:
    int initializeHeaders(int numAtoms);
    int loadHeaders();
    
    const int rank_;
    const int size_;
    const std::string mode_;

    // the file handle
    int ncid_; 
    // dimensions
    int frameDim_;
    int spatialDim_;
    int atomDim_;
    int cellSpatialDim_;
    int cellAngularDim_;
    int labelDim_;
    // variables
    int cellSpatialVar_;
    int cellAngularVar_;
    int cellLengthsVar_;
    int cellAnglesVar_;
    int timeVar_;
    int coordVar_;
};

}
#endif
