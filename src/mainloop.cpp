// Copyright 2013 Robert McGibbon
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <utility>
#include <iomanip>
#include "OpenMM.h"
#include "openmm/serialization/XmlSerializer.h"

#include "utilities.hpp"
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"
#include "ParallelMSM.hpp"

#define MASTER 0
#define WIDTH 5    // used to create the output format
#define PLATFORM_NAME "OpenCL"

using std::string;
using std::vector;
using std::pair;
using std::ifstream;
using std::fstream;
using std::stringstream;
using namespace Tungsten;


int main(int argc, char* argv[]) {
  MPI::Init(argc, argv);
  static const int size = MPI::COMM_WORLD.Get_size();
  static const int rank = MPI::COMM_WORLD.Get_rank();

  // Parse the command line
  if (argc != 5) {
    fprintf(stderr, "usage: %s <system.xml> <integrator.xml> <state.xml> <config.ini>\n", argv[0]);
    exitWithMessage("");
  }

  printUname();
  // Create the context from the input files
  ifstream systemXml(argv[1]);
  ifstream integratorXml(argv[2]);
  OpenMM::Context* context =  createContext(systemXml, integratorXml, PLATFORM_NAME);
  OpenMM::Integrator& integrator = context->getIntegrator();
  const OpenMM::System& system = context->getSystem();
  int numAtoms = system.getNumParticles();
  bool isPeriodic = hasPeriodicBoundaries(system);

  // Set the initial state
  fstream stateXml(argv[3]);
  OpenMM::State* state = OpenMM::XmlSerializer::deserialize<OpenMM::State>(stateXml);
  context->setState(*state);

  // Get the config file
  ConfigOpts opts;
  parseConfigFile(argv[4], &opts);

  // Create the trajectory for our work on this one
  stringstream s;
  s <<  opts.output_root_path << "/trj-" << std::setw(WIDTH) << std::setfill('0') << rank << ".nc";
  string fileName = s.str();
  NetCDFTrajectoryFile file(fileName, "w", numAtoms);

  for (int round = 0; round < opts.n_rounds; round++) {
    file.write(context->getState(OpenMM::State::Positions, isPeriodic));
    for (int step = 0; step < opts.n_steps_per_round; step += opts.save_frequency) {
      integrator.step(opts.save_frequency);
      file.write(context->getState(OpenMM::State::Positions, isPeriodic));
    }
    file.flush();


    ParallelKCenters clusterer(file, 1, opts.atomIndices);
    clusterer.cluster(0.05, pair<int, int>(0, 0));
    
    vector<pair<int, int> > assignments;
    assignments.push_back(pair<int, int>(0,0));
    assignments.push_back(pair<int, int>(0,0));
    assignments.push_back(pair<int, int>(0,0));
    assignments.push_back(pair<int, int>(0,1));
    assignments.push_back(pair<int, int>(0,1));
    assignments.push_back(pair<int, int>(0,1));
    assignments.push_back(pair<int, int>(0,1));
    vector<pair<int, int> > centers;
    centers.push_back(pair<int, int>(0,0));
    centers.push_back(pair<int, int>(0,1));
    ParallelMSM(assignments, centers);



  }
  delete context;
  MPI::Finalize();
}
