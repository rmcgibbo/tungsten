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

#include <limits.h> /* PATH_MAX */
#include <sys/utsname.h>
#include <mpi.h>
#include <stdarg.h>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "openmm/serialization/XmlSerializer.h"
#include "utilities.hpp"
#include "INIReader.h"

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

namespace Tungsten {
static const int MASTER = 0;
using std::string;
using std::ifstream;
using std::vector;
using std::stringstream;

using OpenMM::Context;
using OpenMM::System;
using OpenMM::Platform;
using OpenMM::Integrator;
using OpenMM::XmlSerializer;


Context* createContext(ifstream& systemXml, ifstream& integratorXml, const string& platformName) {
    const int rank = MPI::COMM_WORLD.Get_rank();
    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());

    Platform* platform;
    System* system = XmlSerializer::deserialize<OpenMM::System>(systemXml);
    Integrator* integrator =XmlSerializer::deserialize<OpenMM::Integrator>(integratorXml);
    resetRandomNumberSeed(system, integrator);
    try {
     platform = &Platform::getPlatformByName(platformName);
    } catch (OpenMM::OpenMMException e) {
        printf("An exception occured on rank=%d: %s\n", rank, e.what());
        exitWithMessage("Exit Failure");
    }
    Context* context = new Context(*system, *integrator, *platform);

    if (rank == MASTER) {
        printf("\nContext created on OpenMM platform: %s\n", (*context).getPlatform().getName().c_str());
        vector<string> names = (*context).getPlatform().getPropertyNames();
        for (int i = 0; i < names.size(); i++) {
            string value = (*context).getPlatform().getPropertyValue(*context, names[i]);
            printf("  %s: %s\n", names[i].c_str(), value.c_str());
        }
    }

    return context;
}


void exitWithMessage(const string& format, ...) {
    const int rank = MPI::COMM_WORLD.Get_rank();
    va_list args;
    va_start(args, format);
    int r;
    r = vfprintf(stderr, format.c_str(), args);
    va_end(args);
    fflush(stderr);
    MPI::COMM_WORLD.Abort(EXIT_FAILURE);
}

int printfM(const string& format, ...) {
    const int rank = MPI::COMM_WORLD.Get_rank();
    va_list args;
    va_start(args, format);
    int r;
    if (rank == MASTER)
        r = vprintf(format.c_str(), args);
    va_end(args);
    return r;
}

int printfAOrd(const string& format, ...) {
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();
    va_list args;
    va_start(args, format);
    int r;
    fflush(stdout);
    MPI::COMM_WORLD.Barrier();

    for (int i = 0; i < size; i++) {
        MPI::COMM_WORLD.Barrier();
        if (i == rank)
            r = vprintf(format.c_str(), args);
    }
    fflush(stdout);
    MPI::COMM_WORLD.Barrier();

    va_end(args);
    return r;
}


void printPerformance(double mdTime, time_t endWallTime, time_t startWallTime) {
    const int size = MPI::COMM_WORLD.Get_size();
    const int rank = MPI::COMM_WORLD.Get_rank();
    static const double SEC_PER_DAY = 86400.0;
    // number of simulated ns per day of wall time
    double nsPerDay = (mdTime/1000.0) / (difftime(endWallTime, startWallTime)/SEC_PER_DAY);
    vector<double> recvBuffer(size);
    MPI::COMM_WORLD.Gather(&nsPerDay, 1, MPI_DOUBLE, &recvBuffer[0], 1, MPI_DOUBLE, MASTER);
    if (rank == MASTER) {
        double sum;
        for (int i = 0; i < size; i++)
            sum += recvBuffer[i];
        double mean = sum / size;
        double sumSquare = 0;
        for (int i = 0; i < size; i++)
            sumSquare += (recvBuffer[i] - mean)*(recvBuffer[i] - mean);
        double stdDev = sqrt(sumSquare / size);

        printf("%.2f +/- %.2f ns/day/node\n", mean, stdDev);
        printf("Aggregate: %.2f ns/day\n", sum);
    }
}


ConfigOpts parseConfigFile(const char* configFileName) {
    INIReader reader(configFileName);
    if (reader.ParseError() < 0) {
        exitWithMessage("Could not find config file");
    }

    ConfigOpts out;
    out.numRounds = reader.GetInteger("", "numRounds", -1);
    out.numStepsPerRound = reader.GetInteger("", "numStepsPerRound", -1);
    out.numStepsPerWrite = reader.GetInteger("", "numStepsPerWrite", -1);
    out.outputRootPath = reader.Get("", "outputRootPath", ".");
    out.openmmPlatform = reader.Get("", "openmmPlatform", "Reference");
    out.kcentersRmsdCutoff = reader.GetReal("", "kcentersRmsdCutoff", 1.0);

    const string& atom_indices_file = reader.Get("", "kcentersRmsdIndicesFile", "");
    if (atom_indices_file.size() > 0) {
        char buf[PATH_MAX + 1];
        realpath(atom_indices_file.c_str(), buf);
        ifstream is(buf);
        if (!is) {
            stringstream ss;
            ss << "cannot access " << buf << ": No such file or directory";
            exitWithMessage(ss.str());
        }
        std::istream_iterator<int> start(is), end;
        vector<int> kcentersRmsdIndices(start, end);
        out.kcentersRmsdIndices = kcentersRmsdIndices;
    } else {
        out.kcentersRmsdIndices = vector<int>();
    }

    if (out.numStepsPerRound <= 0)
        exitWithMessage(">numStepsPerRound must be given and greater than 0");
    if (out.numRounds <= 0)
        exitWithMessage("numRounds must be given and greater than 0");
    if (out.numStepsPerWrite <= 0)
        exitWithMessage("numStepsPerWrite must be given and greater than 0");

    return out;
}


bool hasPeriodicBoundaries(const System& system) {
    const int numForces = system.getNumForces();
    bool isPeriodic = false;
    for (int i = 0; i < numForces; i++) {
        const OpenMM::NonbondedForce* force  = dynamic_cast<const OpenMM::NonbondedForce *>(&system.getForce(i));
        if (force != NULL) {
            int nbMethod = force->getNonbondedMethod();
            if (nbMethod == OpenMM::NonbondedForce::NoCutoff) {
                isPeriodic = false;
                break;
            } else if (nbMethod == OpenMM::NonbondedForce::CutoffNonPeriodic) {
                isPeriodic = false;
                break;
            } else if (nbMethod == OpenMM::NonbondedForce::CutoffPeriodic) {
                isPeriodic = true;
                break;
            } else if (nbMethod == OpenMM::NonbondedForce::Ewald) {
                isPeriodic = true;
                break;
            } else if (nbMethod == OpenMM::NonbondedForce::PME) {
                isPeriodic = true;
            }
        }
    }
    return isPeriodic;
}


void printUname(void) {
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();
    struct utsname sysinfo;
    uname(&sysinfo);

    for (int i = 0; i < size; i++) {
        MPI::COMM_WORLD.Barrier();
        if (i == rank) {
            printf("Rank %d Information\n", rank);
            printf("  System: %s %s %s\n", sysinfo.sysname, sysinfo.release, sysinfo.machine);
            printf("  Host Name: %s\n", sysinfo.nodename);
            printf("  OMP_NUM_THREADS: %s\n", getenv("OMP_NUM_THREADS"));
            printf("  OpenMM Version: %s\n", Platform::getOpenMMVersion().c_str());
        }
    }
}


void resetRandomNumberSeed(System* system, Integrator* integrator)  {
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int randomSeed = (rand() % 4294967296) * (rank+1);

    for (int i = 0; i < system->getNumForces(); i++) {
        OpenMM::Force &force = system->getForce(i);
        try {
            OpenMM::AndersenThermostat &ATForce = dynamic_cast<OpenMM::AndersenThermostat &>(force);
            ATForce.setRandomNumberSeed(randomSeed);
            continue;
        } catch(const std::bad_cast &bc) {}
        try {
            OpenMM::MonteCarloBarostat &MCBForce = dynamic_cast<OpenMM::MonteCarloBarostat &>(force);
            MCBForce.setRandomNumberSeed(randomSeed);
            continue;
        } catch(const std::bad_cast &bc) {}
        try {
            OpenMM::MonteCarloAnisotropicBarostat &MCBForce = dynamic_cast<OpenMM::MonteCarloAnisotropicBarostat &>(force);
            MCBForce.setRandomNumberSeed(randomSeed);
            continue;
        } catch(const std::bad_cast &bc) {}
    }

    try {
        OpenMM::LangevinIntegrator &li = dynamic_cast<OpenMM::LangevinIntegrator &>(*integrator);
        li.setRandomNumberSeed(randomSeed);
    } catch(const std::bad_cast &bc) {}
}


double getTemperature(const System& system, const Integrator& integrator)  {
    double temperature = -1;

    try {
        const OpenMM::LangevinIntegrator &li = dynamic_cast<const OpenMM::LangevinIntegrator &>(integrator);
        temperature = li.getTemperature();
    } catch(const std::bad_cast &bc) {}

    try {
        const OpenMM::VariableLangevinIntegrator &li = dynamic_cast<const OpenMM::VariableLangevinIntegrator &>(integrator);
        temperature = li.getTemperature();
    } catch(const std::bad_cast &bc) {}

    try {
        const OpenMM::BrownianIntegrator &li = dynamic_cast<const OpenMM::BrownianIntegrator &>(integrator);
        temperature = li.getTemperature();
    } catch(const std::bad_cast &bc) {}

    for (int i=0; i < system.getNumForces(); i++) {
        const OpenMM::Force &force = system.getForce(i);
        try {
            const OpenMM::AndersenThermostat &atforce = dynamic_cast<const OpenMM::AndersenThermostat &>(force);
            temperature = atforce.getDefaultTemperature();
            continue;
        } catch(const std::bad_cast &bc) {}
    }
    return temperature;
}

}  //  namespace
