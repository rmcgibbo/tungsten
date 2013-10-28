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
// along with Tungsten. If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013 Stanford University
// Author: Robert McGibbon
// Contributors:

#ifndef TUNGSTEN_UTILITIES_H
#define TUNGSTEN_UTILITIES_H
#include "mpi.h"
#include "time.h"
#include <iostream>
#include <fstream>
#include <string>
#include "OpenMM.h"
namespace Tungsten {

typedef struct {
    int numRounds;
    int numStepsPerRound;
    int numStepsPerWrite;
    std::string outputRootPath;
    std::vector<int> kcentersRmsdIndices;
    double kcentersRmsdCutoff;
    std::string openmmPlatform;
} ConfigOpts;


/**
 * Boot up an OpenMM context from a serialized System and Integrator. Reports
 * some output to stdout, and resets the random number seeds to ensure that
 * multiple trajectories initialized on different MPI ranks diverge from
 * each other.
 */
OpenMM::Context* createContext(std::ifstream& systemXml, std::ifstream& integratorXml, const std::string& platformName);

/**
 * Convenience method to exit the process, printing an error message on the
 * MPI root node.
 */
void exitWithMessage(const std::string& fmt, ...);

/**
 * printf only on the master MPI rank.
 */
int printfM(const std::string& fmt, ...);

/**
 * printf on all the ranks, with the output coming
 * to the terminal ordered (first rank 0, then 1, etc)
 */
int printfAOrd(const std::string& format, ...);

/**
 * Calculate and print the performance of an MD simulation
 * from the elapsed simulation time on each node (mdTime)
 * and the elapsed wall time
 */
void printPerformance(double elapsedMDTime, double elapsedWallTime);


/**
 * Parse the config file
 */
ConfigOpts parseConfigFile(const char* configFileName);

/**
 * Does an OpenMM System contain periodic boundary conditions?
 */
bool hasPeriodicBoundaries(const OpenMM::System& system);

/**
 * Print a description of this system, like a linux `uname`, to stdout
 */
void printUname(void);

/**
 * Reset the random seed for any stochastic elements (mc barostat,
 * langevin integrator, etc) in the simulation.
 */
void resetRandomNumberSeed(OpenMM::System* system, OpenMM::Integrator* integrator);

/**
 * Print a each MPI rank's copy of a vector
 */
template <typename T> void printMPIVector(std::vector<T> const & d, bool breif=true) {
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();

    std::cout.flush();
    for (int i = 0; i < size; i++) {
        MPI::COMM_WORLD.Barrier();
        if (rank == i) {
            std::cout << "Rank " << rank << ": [";
            for (int j = 0; j < d.size(); j++) {
                std::cout << d[j];
                if (j < d.size()-1)
                    std::cout << ",  ";
                if (breif && i > 20) {
                    std::cout << "...";
                    break;
                }
            }
            std::cout << "]" << std::endl;
        }
    }
}

/**
 * Get the target temperature of an OpenMM simulation. A value of -1 is
 * returned if the simulation is not under temperature control.
 */
double getTemperature(const OpenMM::System& system, const OpenMM::Integrator& integrator);

} // namespace
#endif
