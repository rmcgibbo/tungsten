#ifndef TUNGSTEN_UTILITIES_H
#define TUNGSTEN_UTILITIES_H
#include <fstream>
#include <string>
#include "mpi.h"
#include "OpenMM.h"
namespace Tungsten {

typedef struct{
  int numRounds;
  int numStepsPerRound;
  int numStepsPerWrite;
  std::string outputRootPath;
  std::vector<int> kcentersRmsdIndices;
  double kcentersRmsdCutoff;
  std::string openmmPlatform;     
} ConfigOpts;


/**
 * Main entry point
 */
int main(int argc, char* argv[]);

/**
 * Boot up a context from a serialized system and integrator
 */
OpenMM::Context* createContext(std::ifstream& systemXml, std::ifstream& integratorXml, const std::string& platformName);

/**
 * Convenience method to exit from MPI
 */
void exitWithMessage(const char* message);

/**
 * Parse the config file
 */
void parseConfigFile(const char* configFileName, ConfigOpts* out);

/**
 * Does an OpenMM System contain periodic boundary conditions?
 */
bool hasPeriodicBoundaries(const OpenMM::System& system);


/**
 * Log the uname to stdout
 */
void printUname(void);

/**
 * Reset the random seed for any stochastic elements (mc barostat, langevin integrator, etc)
 */
void resetRandomNumberSeed(OpenMM::System* system, OpenMM::Integrator* integrator);


template <typename T> void printMPIVector(std::vector<T> const & d) {
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();

  for (int i = 0; i < size; i++) {
    MPI::COMM_WORLD.Barrier();
    if (rank == i) {
      std::cout << "Rank " << rank << ": [";
      for (int j = 0; j < d.size(); j++)
	std::cout << d[j] << ",  ";
      std::cout << "]" << std::endl;
    }
  }
}

}
#endif
