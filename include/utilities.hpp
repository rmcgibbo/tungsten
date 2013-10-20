#ifndef TUNGSTEN_UTILITIES_H
#define TUNGSTEN_UTILITIES_H

#include <cstdlib>
#include <fstream>
#include <string>
#include <mpi.h>
#include "OpenMM.h"
#include "openmm/serialization/XmlSerializer.h"

#include "INIReader.h"

#define MASTER 0

 typedef struct{
     int n_rounds;
     int n_steps_per_round;
     int save_frequency;
     std::string output_root_path;
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















OpenMM::Context* createContext(std::ifstream& systemXml, std::ifstream& integratorXml, const std::string& platformName) {
  int rank = MPI::COMM_WORLD.Get_rank();
  OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

  OpenMM::System* system = OpenMM::XmlSerializer::deserialize<OpenMM::System>(systemXml);
  OpenMM::Integrator* integrator = OpenMM::XmlSerializer::deserialize<OpenMM::Integrator>(integratorXml);
  OpenMM::Platform* platform = &OpenMM::Platform::getPlatformByName(platformName);


  OpenMM::Context* context = new OpenMM::Context(*system, *integrator, *platform);
  if (rank == MASTER) {
    printf("Context created on OpenMM platform: %s\n", (*context).getPlatform().getName().c_str());
    std::vector<std::string> names = (*context).getPlatform().getPropertyNames();
    for (int i = 0; i < names.size(); i++) {
      std::string value = (*context).getPlatform().getPropertyValue(*context, names[i]);
      printf("  %s: %s\n", names[i].c_str(), value.c_str());
    }
  }

  return context;
}


void exitWithMessage(const char* message) {
  int rank = MPI::COMM_WORLD.Get_rank();
  if (rank == MASTER)
      fprintf(stderr, "%s\n", message);
  MPI::Finalize();
  exit(EXIT_FAILURE);
}


void parseConfigFile(const char* configFileName, ConfigOpts* out) {
  INIReader reader(configFileName);
  if (reader.ParseError() < 0) {
      exitWithMessage("Could not find config file");
  }
  
  out->n_steps_per_round = reader.GetInteger("", "n_steps_per_round", -1);
  out->n_rounds = reader.GetInteger("", "n_rounds", -1);
  out->save_frequency = reader.GetInteger("", "save_frequency", -1);
  out->output_root_path = reader.Get("", "output_root_path", ".");

  
  if (out->n_steps_per_round <= 0)
      exitWithMessage("n_steps_per_round must be given and greater than 0");
  if (out->n_rounds <= 0)
    exitWithMessage("n_rounds must be given and greater than 0");
  if (out->save_frequency <= 0)
    exitWithMessage("save_frequency must be given and greater than 0");
}


bool hasPeriodicBoundaries(const OpenMM::System& system) {
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


#endif
