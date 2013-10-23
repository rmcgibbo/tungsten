// Copyright 2013 Robert McGibbon
#ifndef TUNGSTEN_UTILITIES_H
#define TUNGSTEN_UTILITIES_H
#include <fstream>
#include <string>
#include "mpi.h"
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
