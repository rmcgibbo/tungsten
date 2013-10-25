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

#include "mpi.h"
#include "netcdf.h"
#include "sys/utsname.h"  // for uname info in title
#include "time.h"         // for creation info in title
#include "math.h"
#include "sys/stat.h"
#include "stdio.h"
#include <cstring>         // strlen
#include <cstdlib>
#include <limits>         // for portable isinf/isnan
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
using namespace std;

// C++ doesn't seem to contain a portable isfinite() function unless you're
// using C++11. I got this idea from stack overflow. I'm not sure that it
// covers all of the corner cases, but it seems to work for common infs/nans.
template <typename T>
static inline bool isfinite(T x)  {
    return (x <= std::numeric_limits<T>::max() && x >= -std::numeric_limits<T>::max());
}

inline bool file_exists (const string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}


NetCDFTrajectoryFile::NetCDFTrajectoryFile(const string& filename,
        const string& mode, int numAtoms = 0)
    : mode_(mode)
    , rank_(MPI::COMM_WORLD.Get_rank())
    , size_(MPI::COMM_WORLD.Get_size())
    , ncid_(NC_INVALID) {
    int r;

    if (mode_.compare("r") == 0) {
        if ((r = nc_open(filename.c_str(), NC_NOWRITE, &ncid_))) NC_ERR(r);
    } else if (mode_.compare("w") == 0) {
        if (file_exists(filename)) {
            if ((r = nc_open(filename.c_str(), NC_WRITE, &ncid_))) NC_ERR(r);
        } else
            if ((r = nc_create(filename.c_str(), NC_64BIT_OFFSET, &ncid_))) NC_ERR(r);
    } else
        exitWithMessage("NetCDFTrajectoryFile: Bad Mode");

    if (!isvalid())
        exitWithMessage("ERROR OPENING FILE Mode=%s\n", mode.c_str());

    int numDims;
    if ((r = nc_inq_ndims(ncid_, &numDims))) NC_ERR(r);
    if (numDims == 6) {
        loadHeaders();
    } else if (mode_.compare("w") == 0 && numDims == 0)
        initializeHeaders(numAtoms);
    else
        exitWithMessage("Malformed AMBER NetCDF trajector file");

    printfM("\nLoading/Creating output trajectories\n");
    printfM("------------------------------------\n");
    printfAOrd("Rank %d: Initializing NetCDF trj %s; currently contains %d frames\n", rank_, filename.c_str(), getNumFrames());


}

int NetCDFTrajectoryFile::loadHeaders() {
    int r;
    if ((r = nc_inq_dimid(ncid_, "frame", &frameDim_))) NC_ERR(r);
    if ((r = nc_inq_dimid(ncid_, "spatial", &spatialDim_))) NC_ERR(r);
    if ((r = nc_inq_dimid(ncid_, "atom", &atomDim_))) NC_ERR(r);
    if ((r = nc_inq_dimid(ncid_, "cell_spatial", &cellSpatialDim_))) NC_ERR(r);
    if ((r = nc_inq_dimid(ncid_, "cell_angular", &cellAngularDim_))) NC_ERR(r);
    if ((r = nc_inq_dimid(ncid_, "label", &labelDim_))) NC_ERR(r);

    if ((r = nc_inq_varid(ncid_, "cell_spatial", &cellSpatialVar_))) NC_ERR(r);
    if ((r = nc_inq_varid(ncid_, "cell_angular", &cellAngularVar_))) NC_ERR(r);
    if ((r = nc_inq_varid(ncid_, "cell_lengths", &cellLengthsVar_))) NC_ERR(r);
    if ((r = nc_inq_varid(ncid_, "cell_angles", &cellAnglesVar_))) NC_ERR(r);
    if ((r = nc_inq_varid(ncid_, "time", &timeVar_))) NC_ERR(r);
    if ((r = nc_inq_varid(ncid_, "coordinates", &coordVar_))) NC_ERR(r);

    return 1;
}

vector<float> NetCDFTrajectoryFile::loadTime() {
    vector<float> time(getNumFrames());
    int r;
    if ((r = nc_get_var_float(ncid_, timeVar_, &time[0]))) NC_ERR(r);
    return time;
}


int NetCDFTrajectoryFile::initializeHeaders(int numAtoms) {
    int r;
    if (numAtoms <= 0) {
        exitWithMessage("Invalid number of atoms");
    }

    // produce the title to use.
    // NOTE. I tried implementing this using a stringstream and a couple
    // of the entries were garbled. When printing the string inside this code,
    // it looked fine, by ncdump or python to load the resulting file, you could
    // tell. There seems to be no issue with the sprintf solution though, despite
    // being really ugly.
    char title[1024];  // go a little overkill here, I would prefer not to hit a buffer overflow
    time_t curtime = time (NULL);
    struct utsname sysinfo;
    uname(&sysinfo);
    strftime(title, 1024, "Tungsten Adaptive Sampling\nCreated: %Y-%m-%d %H:%M:%S\n", localtime(&curtime));
    sprintf(&title[56], "System: %s %s %s %s\nOpenMMVersion: %s", sysinfo.nodename, sysinfo.sysname, sysinfo.release, sysinfo.machine, OpenMM::Platform::getOpenMMVersion().c_str());

    // declare the global attribues
    const char* application = "OpenMM";
    const char* program = "Tungsten";
    const char* programVersion = "0.1";
    const char* Conventions = "AMBER";
    const char* ConventionVersion = "1.0";

    // add the global attributes
    if ((r = nc_put_att_text(ncid_, NC_GLOBAL, "title", strlen(title)+1, title))) NC_ERR(r);
    if ((r = nc_put_att_text(ncid_, NC_GLOBAL, "application", strlen(application)+1, application))) NC_ERR(r);
    if ((r = nc_put_att_text(ncid_, NC_GLOBAL, "program", strlen(program)+1, program))) NC_ERR(r);
    if ((r = nc_put_att_text(ncid_, NC_GLOBAL, "programVersion", strlen(programVersion)+1, programVersion))) NC_ERR(r);
    if ((r = nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", strlen(Conventions)+1, Conventions))) NC_ERR(r);
    if ((r = nc_put_att_text(ncid_, NC_GLOBAL, "ConventionVersion", strlen(ConventionVersion)+1, ConventionVersion))) NC_ERR(r);

    if ((r = nc_def_dim(ncid_, "frame", NC_UNLIMITED, &frameDim_))) NC_ERR(r);
    if ((r = nc_def_dim(ncid_, "spatial", 3, &spatialDim_))) NC_ERR(r);
    if ((r = nc_def_dim(ncid_, "atom", numAtoms, &atomDim_))) NC_ERR(r);
    if ((r = nc_def_dim(ncid_, "cell_spatial", 3, &cellSpatialDim_))) NC_ERR(r);
    if ((r = nc_def_dim(ncid_, "cell_angular", 3, &cellAngularDim_))) NC_ERR(r);
    if ((r = nc_def_dim(ncid_, "label", 5, &labelDim_))) NC_ERR(r);

    int cellSpatialDimids[] = {spatialDim_};
    int cellAngularDimids[] = {spatialDim_, labelDim_};
    int cellLengthsDimids[] = {frameDim_, cellSpatialDim_};
    int cellAnglesDimids[]  = {frameDim_, cellAngularDim_};
    int timeDimids[]        = {frameDim_};
    int coordinatesDimids[] = {frameDim_, atomDim_, spatialDim_};
    if ((r = nc_def_var(ncid_, "cell_spatial", NC_CHAR, 1, cellSpatialDimids, &cellSpatialVar_))) NC_ERR(r);
    if ((r = nc_def_var(ncid_, "cell_angular", NC_CHAR, 2, cellAngularDimids, &cellAngularVar_))) NC_ERR(r);
    if ((r = nc_def_var(ncid_, "cell_lengths", NC_DOUBLE, 2, cellLengthsDimids, &cellLengthsVar_))) NC_ERR(r);
    if ((r = nc_def_var(ncid_, "cell_angles", NC_DOUBLE, 2, cellAnglesDimids, &cellAnglesVar_))) NC_ERR(r);
    if ((r = nc_def_var(ncid_, "time", NC_FLOAT, 1, timeDimids, &timeVar_))) NC_ERR(r);
    if ((r = nc_def_var(ncid_, "coordinates", NC_FLOAT, 3, coordinatesDimids, &coordVar_))) NC_ERR(r);

    const char* degree = "degree";
    const char* angstrom = "angstrom";
    const char* picoseconds = "picoseconds";
    const char* xyz = "XYZ";
    const char* alphabetagamma = "alphabeta\0gamma";
    if ((r = nc_put_att_text(ncid_, cellAnglesVar_, "units", strlen(degree), degree))) NC_ERR(r);
    if ((r = nc_put_att_text(ncid_, cellLengthsVar_, "units", strlen(angstrom), angstrom))) NC_ERR(r);
    if ((r = nc_put_att_text(ncid_, timeVar_, "units", strlen(picoseconds), picoseconds))) NC_ERR(r);
    if ((r = nc_put_att_text(ncid_, coordVar_, "units", strlen(angstrom), angstrom))) NC_ERR(r);

    if ((r = nc_enddef(ncid_))) NC_ERR(r);

    if ((r = nc_put_var_text(ncid_, cellSpatialVar_, xyz))) NC_ERR(r);
    if ((r = nc_put_var_text(ncid_, cellAngularVar_, alphabetagamma))) NC_ERR(r);
    return 1;
}


int NetCDFTrajectoryFile::write(OpenMM::State state) {
    if (mode_.compare("w") != 0)
        exitWithMessage("Writing is not allowed in this mode");
    int frame = getNumFrames();
    int numAtoms = getNumAtoms();

    OpenMM::Vec3 a;
    OpenMM::Vec3 b;
    OpenMM::Vec3 c;
    float time = state.getTime();
    state.getPeriodicBoxVectors(a, b, c);
    double cellLengths[] = {a[0]*10.0, b[1]*10.0, c[2]*10.0};
    double cellAngles[] = {90.0, 90.0, 90.0};
    vector<OpenMM::Vec3> mmPositions = state.getPositions();
    vector<float> positions(numAtoms*3);

    bool hasInfinite = false;
    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < 3; j++) {
            if (!isfinite(mmPositions[i][0]) || !isfinite(mmPositions[i][1]) || !isfinite(mmPositions[i][2]))
                hasInfinite = true;
            positions[i*3 + j] = mmPositions[i][j] * 10.0;  // convert to angstroms
        }

    if (hasInfinite)
        for (int i = 0; i < positions.size(); i++)
            positions[i] = NAN;

    int r;  // errors
    const size_t index1[] = {frame};
    if ((r = nc_put_var1_float(ncid_, timeVar_, index1, &time))) NC_ERR(r);

    const size_t index2[] = {frame, 0};
    const size_t count2[] = {1, 3};
    if ((r = nc_put_vara_double(ncid_, cellLengthsVar_, index2, count2, cellLengths))) NC_ERR(r);
    if ((r = nc_put_vara_double(ncid_, cellAnglesVar_, index2, count2, cellAngles))) NC_ERR(r);

    const size_t index3[] = {frame, 0, 0};
    const size_t count3[] = {1, getNumAtoms(), 3};
    if ((r = nc_put_vara_float(ncid_, coordVar_, index3, count3, &positions[0]))) NC_ERR(r);

    return 1;
}


PositionsAndPeriodicBox NetCDFTrajectoryFile::loadNonlocalStateMPI(int rank, int index) {
    const int magic1 = -12345;
    const int magic2 = -23456;
    const int tag1 = 123;
    const int tag2 = 234;
    int r;

    // This could be done more efficiently by only exchanging the necessary pairs, but
    // its complex to get right without deadlocking because many nodes may request
    // data that all needs to come from one node. So lets just share *all* of the indices.
    // with Allgather.
    size_t numAtoms = getNumAtoms();
    vector<int> gatheredRank(size_);
    vector<int> gatheredIndex(size_);
    MPI::COMM_WORLD.Allgather(&rank, 1, MPI_INT, &gatheredRank[0], 1, MPI_INT);
    MPI::COMM_WORLD.Allgather(&index, 1, MPI_INT, &gatheredIndex[0], 1, MPI_INT);

    // Now, each node needs to look through the gatheredRank/gatheredSize and make
    // the appropriate send calls to give the data to the other nodes. It also
    // needs to obviously make its own receive call.
    vector<float> sendPositionsBuffer(numAtoms*3, magic1);
    vector<float> recvPositionsBuffer(numAtoms*3, magic2);
    vector<double> sendCellLengthsBuffer(3, magic1);
    vector<double> recvCellLengthsBuffer(3, magic2);
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
            size_t start3[] = {gatheredIndex[i], 0, 0};
            size_t count3[] = {1, numAtoms, 3};
            if ((r = nc_get_vara_float(ncid_, coordVar_, start3, count3, &sendPositionsBuffer[0]))) NC_ERR(r);
            size_t start2[] = {gatheredIndex[i], 0};
            size_t count2[] = {1, 3};
            if ((r = nc_get_vara_double(ncid_, cellLengthsVar_, start2, count2, &sendCellLengthsBuffer[0]))) NC_ERR(r);
            sendRequests1.push_back(MPI::COMM_WORLD.Isend(&sendPositionsBuffer[0],
                                    numAtoms*3, MPI_FLOAT, i, tag1));
            sendRequests2.push_back(MPI::COMM_WORLD.Isend(&sendCellLengthsBuffer[0],
                                    3, MPI_DOUBLE, i, tag2));
        }

        if (i == rank_) {
            // receive the data
            recvRequest1 = MPI::COMM_WORLD.Irecv(&recvPositionsBuffer[0],
                                                 numAtoms*3, MPI_FLOAT, j, tag1);
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
    for (int i = 0; i < numAtoms; i++) {
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
