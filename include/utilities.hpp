#ifndef TUNGSTEN_UTILITIES_H
#define TUNGSTEN_UTILITIES_H

#include <fstream>
#include <string>
#include "OpenMM.h"

 typedef struct{
     int numRounds;
     int numStepsPerRound;
     int numStepsPerWrite;
     std::string outputRootPath;
     std::vector<int> rmsdAtomIndices;
     float kcentersRmsdCutoff;
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

void printVector(const std::vector<float>& d);

#endif
