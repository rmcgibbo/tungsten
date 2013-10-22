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
    delete handle;
  }
  int write(OpenMM::State state);
  void readPositions(int stride, float* out) const;
  void readAxisMajorPositions(int stride, const std::vector<int>& atomIndices,
			      int atomAlignment, fvector16& out) const;

  int getNumAtoms() const {
    return n_atoms;
  }
  int getNumFrames() const {
    if (handle == NULL)
      return 0;
    return handle->get_dim("frame")->size();
  }
  void flush(void) {
    handle->sync();
  }

private:
  NcFile* handle;
  const char* mode;
  int n_atoms;
  int initializeHeaders(void);
};

}
#endif
