// Copyright 2013 Robert McGibbon
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include "utilities.hpp"
#include "typedefs.hpp"
#include "NetCDFTrajectoryFile.hpp"
namespace Tungsten {
static const int MASTER = 0;

using std::string;
using std::vector;

NetCDFTrajectoryFile::NetCDFTrajectoryFile(const string& filename,
       const char* mode, int numAtoms = 0):
  mode_(mode), numAtoms_(numAtoms), rank_(MPI::COMM_WORLD.Get_rank()),
  size_(MPI::COMM_WORLD.Get_size())
{
  NcFile::FileMode ncMode;
  if (strcmp(mode_, "r") == 0) {
    handle_ = new NcFile(filename.c_str(), NcFile::ReadOnly);
    numAtoms_ = handle_->get_dim("atom")->size();
  } else if (strcmp(mode_, "w") == 0) {
    ncMode = NcFile::Replace;
    if (numAtoms_ <= 0) {
      exitWithMessage("Invalid number of atoms");
    }
    handle_ = new NcFile(filename.c_str(), NcFile::Replace, NULL, 0,
                        NcFile::Offset64Bits);
  } else {
    exitWithMessage("NetCDFTrajectoryFile: Bad Mode");
  }

  if (!handle_->is_valid()) {
    printf("------------------\n");
    printf("ERROR OPENING FILE\n");
    printf("Mode=%s\n", mode);
    exit(1);
  }

  if (strcmp(mode_, "w") == 0) {
    initializeHeaders();
  }
}


int NetCDFTrajectoryFile::initializeHeaders() {
  handle_->add_att("title", "");
  handle_->add_att("application", "OpenMM");
  handle_->add_att("program", "Tungsten");
  handle_->add_att("programVersion", "0.1");
  handle_->add_att("Conventions", "AMBER");
  handle_->add_att("ConventionVersion", "1.0");

  NcDim* frameDim = handle_->add_dim("frame", 0);
  NcDim* spatialDim = handle_->add_dim("spatial", 3);
  NcDim* atomDim = handle_->add_dim("atom", numAtoms_);
  NcDim* cellSpatialDim = handle_->add_dim("cell_spatial", 3);
  NcDim* cellAngularDim = handle_->add_dim("cell_angular", 3);
  NcDim* labelDim = handle_->add_dim("label", 5);

  NcVar* cellSpatialVar = handle_->add_var("cell_spatial", ncChar, spatialDim);
  NcVar* cellAngularVar = handle_->add_var("cell_angular", ncChar, spatialDim, labelDim);
  NcVar* cellLengthsVar = handle_->add_var("cell_lengths", ncDouble, frameDim, cellSpatialDim);
  NcVar* cellAnglesVar = handle_->add_var("cell_angles", ncDouble, frameDim, cellAngularDim);
  cellAnglesVar->add_att("units", "degree");
  cellLengthsVar->add_att("units", "angstrom");


  cellSpatialVar->put("XYZ", 3);
  cellAngularVar->put("alpha", 1, 5);
  cellAngularVar->set_cur(1, 0);
  cellAngularVar->put("beta", 1, 4);
  cellAngularVar->set_cur(2, 0);
  cellAngularVar->put("gamma", 1, 5);

  NcVar* timeVar = handle_->add_var("time", ncFloat, frameDim);
  timeVar->add_att("units", "picosecond");
  NcVar* coordVar = handle_->add_var("coordinates", ncFloat, frameDim, atomDim, spatialDim);
  coordVar->add_att("units", "angstrom");

  return 1;
}


int NetCDFTrajectoryFile::write(OpenMM::State state) {
  if (strcmp(mode_, "w") != 0)
    exitWithMessage("Writing is not allowed in this mode");
  int frame = handle_->get_dim("frame")->size();

  OpenMM::Vec3 a;
  OpenMM::Vec3 b;
  OpenMM::Vec3 c;
  double time = state.getTime();
  state.getPeriodicBoxVectors(a, b, c);
  double cellLengths[] = {a[0]*10.0, b[1]*10.0, c[2]*10.0};
  double cellAngles[] = {90.0, 90.0, 90.0};
  vector<OpenMM::Vec3> positions = state.getPositions();
  for (size_t i = 0; i < positions.size(); i++)
    positions[i] *= 10;

  handle_->get_var("time")->put_rec(&time, frame);
  handle_->get_var("cell_lengths")->put_rec(cellLengths, frame);
  handle_->get_var("cell_angles")->put_rec(cellAngles, frame);
  handle_->get_var("coordinates")->put_rec(&positions[0][0], frame);

  return 1;
}


vector<OpenMM::Vec3> NetCDFTrajectoryFile::loadNonlocalPositionsMPI(int rank, int index) {

  // This could be done more efficiently by only exchanging the necessary pairs, but
  // its complex to get right without deadlocking because many nodes may request 
  // data that all needs to come from one node. So lets just share *all* of the indices.
  // with Allgather.
  vector<int> gatheredRank(size_);
  vector<int> gatheredIndex(size_);
  MPI::COMM_WORLD.Allgather(&rank, 1, MPI_INT, &gatheredRank[0], 1, MPI_INT);
  MPI::COMM_WORLD.Allgather(&index, 1, MPI_INT, &gatheredIndex[0], 1, MPI_INT);
  NcVar* coord = handle_->get_var("coordinates");

  // Now, each node needs to look through the gatheredRank/gatheredSize and make 
  // the appropriate send calls to give the data to the other nodes. It also
  // needs to obviously make its own receive call.
  const int tag = 123;
  vector<double> sendBuffer(numAtoms_*3);
  vector<double> recvBuffer(numAtoms_*3, -12345);
  // each node only gets one receive, but it may need to
  // make 0 or more sends
  MPI::Request recvRequest;
  vector<MPI::Request> sendRequests;

  for (int i = 0; i < size_; i++) {
    int j = gatheredRank[i];
    // j needs to send data to node i.
    // lets do the pair of send/recv.

    if (j == rank_) {
      // load the data into the send buffer
      coord->set_cur(gatheredIndex[i], 0, 0);
      coord->get(&sendBuffer[0], 1, numAtoms_, 3);
      sendRequests.push_back(MPI::COMM_WORLD.Isend(&sendBuffer[0], numAtoms_*3, MPI_DOUBLE, i, tag));
    }

    if (i == rank_) {
      // receive the data
      recvRequest = MPI::COMM_WORLD.Irecv(&recvBuffer[0], numAtoms_*3, MPI_DOUBLE, j, tag);
    }
  }
  
  // make sure the nonblocking requests went through
  recvRequest.Wait();
  for (int i = 0; i < sendRequests.size(); i++)
    sendRequests[i].Wait();

  // reformat our recvBuffer into a vector of OpenMM::Vec3
  vector<OpenMM::Vec3> retval;
  for (int i = 0; i < numAtoms_; i++) {
    OpenMM::Vec3 v(recvBuffer[i*3 + 0], recvBuffer[i*3 + 1], recvBuffer[i*3 + 2]);
    retval.push_back(v);
  }

  return retval;
}

void NetCDFTrajectoryFile::loadAllAxisMajorPositions(
   int stride, const vector<int>& atomIndices, int atomAlignment, fvector16& out) const
{
  int numTotalFrames = handle_->get_dim("frame")->size();
  int numFrames = (numTotalFrames+stride-1)/stride;
  int numPaddedAtoms = ((atomIndices.size() + 3) / atomAlignment) * atomAlignment;
  NcVar* coord = handle_->get_var("coordinates");
  out.resize(numFrames*3*numPaddedAtoms);
  vector<float> frame(numAtoms_*3);
  
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
}

} //namespace
