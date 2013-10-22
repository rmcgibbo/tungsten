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

using std::string;
using std::vector;
using std::pair;
using std::ifstream;
using std::fstream;
using std::stringstream;
using namespace Tungsten;


int main(int argc, char* argv[]) {
  MPI::Init(argc, argv);
  const int size = MPI::COMM_WORLD.Get_size();
  const int rank = MPI::COMM_WORLD.Get_rank();

  // Parse the command line
  if (argc != 5) {
    fprintf(stderr, "usage: %s <system.xml> <integrator.xml> <state.xml> <config.ini>\n", argv[0]);
    exitWithMessage("");
  }

  printUname();

  // Get the config file
  ConfigOpts opts;
  parseConfigFile(argv[4], &opts);

  // Create the context from the input files
  ifstream systemXml(argv[1]);
  ifstream integratorXml(argv[2]);
  OpenMM::Context* context =  createContext(systemXml, integratorXml, opts.openmmPlatform);
  OpenMM::Integrator& integrator = context->getIntegrator();
  const OpenMM::System& system = context->getSystem();
  int numAtoms = system.getNumParticles();
  bool isPeriodic = hasPeriodicBoundaries(system);

  // Set the initial state
  fstream stateXml(argv[3]);
  OpenMM::State* state = OpenMM::XmlSerializer::deserialize<OpenMM::State>(stateXml);
  context->setState(*state);

  // Create the trajectory for our work on this one
  stringstream s;
  s <<  opts.outputRootPath << "/trj-" << std::setw(WIDTH) << std::setfill('0') << rank << ".nc";
  string fileName = s.str();
  NetCDFTrajectoryFile file(fileName, "w", numAtoms);

  for (int round = 0; round < opts.numRounds; round++) {
    if (rank == MASTER)
      printf("\nBeginning Round %d\n===================\n", round);

    file.write(context->getState(OpenMM::State::Positions, isPeriodic));
    for (int step = 0; step < opts.numStepsPerRound; step += opts.numStepsPerWrite) {
      if (rank == MASTER)
	printf("#");
      integrator.step(opts.numStepsPerWrite);
      file.write(context->getState(OpenMM::State::Positions, isPeriodic));
    }
    if (rank == MASTER)
      printf("\n");
    file.flush();

    ParallelKCenters clusterer(file, 1, opts.kcentersRmsdIndices);
    clusterer.cluster(opts.kcentersRmsdCutoff, pair<int, int>(0, 0));
    ParallelMSM markovModel(clusterer.getAssignments(), clusterer.getCenters());
    pair<int, int> newCoordinateIndex = markovModel.scatterMinCountStateLabels();

    MPI::COMM_WORLD.Barrier();
    fflush(stdout);
    printf("New Coordinates on Rank=%d: (%d, %d)\n", rank, newCoordinateIndex.first, newCoordinateIndex.second);
    

    file.loadPositionsMPI(newCoordinateIndex.first, newCoordinateIndex.second);

    

  }
  delete context;
  MPI::Finalize();
}
