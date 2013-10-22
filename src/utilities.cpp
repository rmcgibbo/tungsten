// Copyright 2013 Robert McGibbon
#include <limits.h> /* PATH_MAX */
#include <sys/utsname.h>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "openmm/serialization/XmlSerializer.h"
#include "utilities.hpp"
#include "INIReader.h"

#define MASTER 0
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

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
  static const int rank = MPI::COMM_WORLD.Get_rank();
  Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());

  System* system = XmlSerializer::deserialize<OpenMM::System>(systemXml);
  Integrator* integrator =XmlSerializer::deserialize<OpenMM::Integrator>(integratorXml);
  resetRandomNumberSeed(system, integrator);
  Platform* platform = &Platform::getPlatformByName(platformName);
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


void exitWithMessage(const char* message) {
  int rank = MPI::COMM_WORLD.Get_rank();
  if (rank == MASTER)
      fprintf(stderr, "%s\n", message);
  MPI::Finalize();
  exit(EXIT_FAILURE);
}
void exitWithMessage(const string& message) {
  exitWithMessage(message.c_str());
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
  out->openmm_platform = reader.Get("", "openmm_platform", "Reference");
  out->kcenters_rmsd_cutoff = reader.GetFloat("", "kcenters_rmsd_cutoff", 1.0);


  const string& atom_indices_file = reader.Get("", "kcenters_rmsd_indices_file", "");
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
    vector<int> atomIndices(start, end);
    out->atomIndices = atomIndices;
  } else {
    out->atomIndices = vector<int>();
  }

  if (out->n_steps_per_round <= 0)
      exitWithMessage("n_steps_per_round must be given and greater than 0");
  if (out->n_rounds <= 0)
    exitWithMessage("n_rounds must be given and greater than 0");
  if (out->save_frequency <= 0)
    exitWithMessage("save_frequency must be given and greater than 0");
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
      printf("  OpenMM Version: %s\n", Platform::getOpenMMVersion().c_str());
    }
  }
}


void resetRandomNumberSeed(System* system, Integrator* integrator)  {
  const int rank = MPI::COMM_WORLD.Get_rank(); 
  const int randomSeed = (rand() % 4294967296) * (rank+1);

  for(int i=0; i< system->getNumForces(); i++) {
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

void printVector(const vector<float>& d) {
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();

  for (int i = 0; i < size; i++) {
    MPI::COMM_WORLD.Barrier();
    if (rank == i) {
      printf("Rank %d: [", rank);
      for (int j = 0; j < d.size(); j++)
        printf("%f,  ", d[j]);
      printf("]\n");
    }
  }
}
