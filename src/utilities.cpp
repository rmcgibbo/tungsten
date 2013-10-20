#include <cstdlib>
#include <ctime>
#include <sys/utsname.h>
#include <mpi.h>
#include "openmm/serialization/XmlSerializer.h"
#include "utilities.hpp"
#include "INIReader.h"
#define MASTER 0

OpenMM::Context* createContext(std::ifstream& systemXml, std::ifstream& integratorXml, const std::string& platformName) {
  int rank = MPI::COMM_WORLD.Get_rank();
  OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

  OpenMM::System* system = OpenMM::XmlSerializer::deserialize<OpenMM::System>(systemXml);
  OpenMM::Integrator* integrator = OpenMM::XmlSerializer::deserialize<OpenMM::Integrator>(integratorXml);
  resetRandomNumberSeed(integrator);
  OpenMM::Platform* platform = &OpenMM::Platform::getPlatformByName(platformName);
  OpenMM::Context* context = new OpenMM::Context(*system, *integrator, *platform);

  if (rank == MASTER) {
    printf("\nContext created on OpenMM platform: %s\n", (*context).getPlatform().getName().c_str());
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


void printUname(void) {
  static const int rank = MPI::COMM_WORLD.Get_rank();
  static const int size = MPI::COMM_WORLD.Get_size();
  struct utsname sysinfo;
  uname(&sysinfo);
  
  for (int i = 0; i < size; i++) {
    MPI::COMM_WORLD.Barrier();
    if (i == rank) {
      printf("Rank %d Information\n", rank);
      printf("  System: %s %s %s\n", sysinfo.sysname, sysinfo.release, sysinfo.machine);
      printf("  Host Name: %s\n", sysinfo.nodename);
      printf("  OMP_NUM_THREADS: %s\n", getenv("OMP_NUM_THREADS"));
      printf("  OpenMM Version: %s\n", OpenMM::Platform::getOpenMMVersion().c_str());
    }
  }
}


void resetRandomNumberSeed(OpenMM::Integrator* integrator) {
  srand(time(NULL));
  int seed = rand() % 4294967296;

  OpenMM::LangevinIntegrator* lIntegrator  = dynamic_cast<OpenMM::LangevinIntegrator*>(integrator);
  if (lIntegrator != NULL) {
    printf("LangevinIntegrator\n");
    lIntegrator->setRandomNumberSeed(seed);
  }
}
