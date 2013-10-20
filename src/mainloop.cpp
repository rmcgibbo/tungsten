#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <utility>
#include <iomanip>
#include <mpi.h>
#include "OpenMM.h"
#include "openmm/serialization/XmlSerializer.h"

#include "utilities.hpp"
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"


 #define MASTER 0
 #define WIDTH 5    //used to create the output format
 #define PLATFORM_NAME "OpenCL"


 int main(int argc, char* argv[]) {
     MPI::Init(argc, argv);
     int size = MPI::COMM_WORLD.Get_size();
     int rank = MPI::COMM_WORLD.Get_rank();

     // Parse the command line
     if (argc != 5) {
	 fprintf(stderr, "usage: %s <system.xml> <integrator.xml> <state.xml> <config.ini>\n", argv[0]);
	 exitWithMessage("");
     }

     // Create the context from the input files
     std::ifstream systemXml(argv[1]);
     std::ifstream integratorXml(argv[2]);
     OpenMM::Context* context =  createContext(systemXml, integratorXml, PLATFORM_NAME);
     OpenMM::Integrator& integrator = context->getIntegrator();
     const OpenMM::System& system = context->getSystem();
     int numAtoms = system.getNumParticles();
     bool isPeriodic = hasPeriodicBoundaries(system);

     // Set the initial state
     std::ifstream stateXml(argv[3]);
     OpenMM::State* state = OpenMM::XmlSerializer::deserialize<OpenMM::State>(stateXml);
     context->setState(*state);

     // Get the config file
     ConfigOpts opts;
     parseConfigFile(argv[4], &opts);

     // Create the trajectory for our work on this one
     std::stringstream s;
     s <<  opts.output_root_path << "/trj-" << std::setw(WIDTH) << std::setfill('0') << rank << ".nc";
     std::string fileName = s.str();
     NetCDFTrajectoryFile file(fileName, "w", numAtoms);

     for (int round = 0; round < opts.n_rounds; round++) {
       file.write(context->getState(OpenMM::State::Positions, isPeriodic));
       for (int step = 0; step < opts.n_steps_per_round; step += opts.save_frequency) {
	 integrator.step(opts.save_frequency);
	 file.write(context->getState(OpenMM::State::Positions, isPeriodic));
       }
       file.flush();


       ParallelKCenters clusterer(file, 1);
       std::pair<int, int>pair (0, 0);

       std::vector<float> distance = clusterer.getRmsdsTo(pair);
       MPI::COMM_WORLD.Barrier();
       for (int i = 0; i < size; i++) {
	 MPI::COMM_WORLD.Barrier();
	 if (i == rank) {
	   printf("\n");
	   for (int j = 0; j < distance.size(); j++)
	     printf("RANK %d:: distance[%d] = %f\n", rank, j, distance[j]);
	 }
       }
     }

    delete context;
    MPI::Finalize();
}
