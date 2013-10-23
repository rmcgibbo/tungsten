// Copyright 2013 Robert McGibbon
#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <limits>  // for portable isinf/isnan
#include <vector>
#include <string>
#include <iostream>
#include "utilities.hpp"
#include "aligned_allocator.hpp"
#include "typedefs.hpp"
#include "NetCDFTrajectoryFile.hpp"
namespace Tungsten {
static const int MASTER = 0;

using std::string;
using std::vector;
using std::numeric_limits;

#include <iostream>
#include "math.h"
using namespace std;

// C++ doesn't seem to contain a portable isfinite() function unless you're
// using C++11. I got this idea from stack overflow. I'm not sure that it
// covers all of the corner cases, but it seems to work for common infs/nans.
template <typename T>
static inline bool isfinite(T x)  {
    return (x <= std::numeric_limits<T>::max() && x >= -std::numeric_limits<T>::max());
}

NetCDFTrajectoryFile::NetCDFTrajectoryFile(const string& filename,
        const string& mode, int numAtoms = 0)
    : mode_(mode)
    , numAtoms_(numAtoms)
    , rank_(MPI::COMM_WORLD.Get_rank())
    , size_(MPI::COMM_WORLD.Get_size()) {
    NcFile::FileMode ncMode;
    if (mode_.compare("r") == 0) {
        handle_ = new NcFile(filename.c_str(), NcFile::ReadOnly);
        numAtoms_ = handle_->get_dim("atom")->size();
    } else if (mode_.compare("w") == 0) {
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
        fprintf(stderr, "ERROR OPENING FILE Mode=%s\n", mode.c_str());
        exit(1);
    }

    if (mode_.compare("w") == 0)
        initializeHeaders();
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
    if (mode_.compare("w") != 0)
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

    bool hasInfinite = false;
    for (int i = 0; i < positions.size(); i++) {
        positions[i] *= 10;
        if (!isfinite(positions[i][0]) || !isfinite(positions[i][1]) || !isfinite(positions[i][2]))
            hasInfinite = true;
    }

    if (hasInfinite)
        for (int i = 0; i < positions.size(); i++)
            positions[i] = OpenMM::Vec3(NAN, NAN, NAN);

    handle_->get_var("time")->put_rec(&time, frame);
    handle_->get_var("cell_lengths")->put_rec(cellLengths, frame);
    handle_->get_var("cell_angles")->put_rec(cellAngles, frame);
    handle_->get_var("coordinates")->put_rec(&positions[0][0], frame);

    return 1;
}


PositionsAndPeriodicBox NetCDFTrajectoryFile::loadNonlocalStateMPI(int rank, int index) {

    // This could be done more efficiently by only exchanging the necessary pairs, but
    // its complex to get right without deadlocking because many nodes may request
    // data that all needs to come from one node. So lets just share *all* of the indices.
    // with Allgather.
    vector<int> gatheredRank(size_);
    vector<int> gatheredIndex(size_);
    MPI::COMM_WORLD.Allgather(&rank, 1, MPI_INT, &gatheredRank[0], 1, MPI_INT);
    MPI::COMM_WORLD.Allgather(&index, 1, MPI_INT, &gatheredIndex[0], 1, MPI_INT);
    NcVar* coord = handle_->get_var("coordinates");
    NcVar* cellLengths = handle_->get_var("cell_lengths");

    // Now, each node needs to look through the gatheredRank/gatheredSize and make
    // the appropriate send calls to give the data to the other nodes. It also
    // needs to obviously make its own receive call.
    const int tag1 = 123;
    const int tag2 = 234;
    vector<float> sendPositionsBuffer(numAtoms_*3);
    vector<float> recvPositionsBuffer(numAtoms_*3, -12345);
    vector<double> sendCellLengthsBuffer(3, -23456);
    vector<double> recvCellLengthsBuffer(3);
    // each node only gets one receive, but it may need to
    // make 0 or more sends
    MPI::Request recvRequest2, recvRequest1;
    vector<MPI::Request> sendRequests1, sendRequests2;


    for (int i = 0; i < size_; i++) {
        int j = gatheredRank[i];
        // j needs to send data to node i.
        // lets do the pair of send/recv.
        if (j == rank_) {
            // load the data into the send buffer
            coord->set_cur(gatheredIndex[i], 0, 0);
            coord->get(&sendPositionsBuffer[0], 1, numAtoms_, 3);
            coord->set_cur(gatheredIndex[i], 0, 0);
            cellLengths->get(&sendCellLengthsBuffer[0], 1, 3);
            sendRequests1.push_back(MPI::COMM_WORLD.Isend(&sendPositionsBuffer[0],
                                    numAtoms_*3, MPI_FLOAT, i, tag1));
            sendRequests2.push_back(MPI::COMM_WORLD.Isend(&sendCellLengthsBuffer[0],
                                    3, MPI_DOUBLE, i, tag2));
        }

        if (i == rank_) {
            // receive the data
            recvRequest1 = MPI::COMM_WORLD.Irecv(&recvPositionsBuffer[0],
                                                 numAtoms_*3, MPI_FLOAT, j, tag1);
            recvRequest2 = MPI::COMM_WORLD.Irecv(&recvCellLengthsBuffer[0], 3,
                                                 MPI_DOUBLE, j, tag2);
        }
    }
    // make sure the nonblocking requests went through
    recvRequest1.Wait();
    for (int i = 0; i < sendRequests1.size(); i++)
        sendRequests1[i].Wait();
    recvRequest2.Wait();
    for (int i = 0; i < sendRequests2.size(); i++)
        sendRequests2[i].Wait();

    // reformat our recvBuffer into a vector of OpenMM::Vec3. Also convert
    // from loaded data in A to OpenMM's nm unit system.
    vector<OpenMM::Vec3> newPositions;
    for (int i = 0; i < numAtoms_; i++) {
        OpenMM::Vec3 v(recvPositionsBuffer[i*3 + 0] / 10.0,
                       recvPositionsBuffer[i*3 + 1] / 10.0,
                       recvPositionsBuffer[i*3 + 2] / 10.0);
        newPositions.push_back(v);
    }

    PositionsAndPeriodicBox retval;
    retval.positions = newPositions;
    retval.boxA = OpenMM::Vec3(recvCellLengthsBuffer[0] / 10.0, 0, 0);
    retval.boxB = OpenMM::Vec3(0, recvCellLengthsBuffer[1] / 10.0, 0);
    retval.boxC = OpenMM::Vec3(0, 0, recvCellLengthsBuffer[2] / 10.0);
    return retval;
}


}  // namespace
