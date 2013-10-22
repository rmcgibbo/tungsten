#ifndef TUNGSTEN_NETCDFTRAJECTORY_FILE_H_
#define TUNGSTEN_NETCDFTRAJECTORY_FILE_H_
#include <vector>
#include <string>
#include <netcdfcpp.h>
#include "OpenMM.h"
#include "typedefs.hpp"
namespace Tungsten {

class NetCDFTrajectoryFile {
public:
  NetCDFTrajectoryFile(const std::string& filename, const char* mode, int n_atom);
  ~NetCDFTrajectoryFile(void) {
    delete handle_;
  }

  /*
   * Write an OpenMM state to disk.
   *
   */
  int write(OpenMM::State state);


  /*
   * Fetch the `index`-th frame of positions from the trajectory
   * residing on the `rank`th MPI rank.
   */
  std::vector<OpenMM::Vec3> loadNonlocalPositionsMPI(int rank, int index);

  /*
   * Read an entire trajectory from disk, into aligned memory.
   *
   */
  void loadAllAxisMajorPositions(int stride, const std::vector<int>& atomIndices,
			      int atomAlignment, fvector16& out) const;


  int getNumAtoms() const {
    return numAtoms_;
  }
  int getNumFrames() const {
    if (handle_ == NULL)
      return 0;
    return handle_->get_dim("frame")->size();
  }
  void flush(void) {
    if (handle_ != NULL)
      handle_->sync();
  }

private:
  const int rank_;
  const int size_;
  int numAtoms_;
  NcFile* handle_;
  const char* mode_;

  int initializeHeaders(void);
};

}
#endif
