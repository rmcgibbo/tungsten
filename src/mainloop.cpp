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
// along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013 Stanford University
// Author: Robert McGibbon
// Contributors:

#include "mpi.h"   // parallelism
#include "time.h"  // for timing the speed
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>  // set fill with for building filename
#include "OpenMM.h"
#include "openmm/serialization/XmlSerializer.h"
#include "utilities.hpp"
#include "typedefs.hpp"
#include "NetCDFTrajectoryFile.hpp"
#include "ParallelKCenters.hpp"
#include "ParallelMSM.hpp"


using std::string;
using std::vector;
using std::ifstream;
using std::fstream;
using std::stringstream;
using OpenMM::Context;
using OpenMM::State;
using OpenMM::Integrator;
using OpenMM::System;
using OpenMM::XmlSerializer;
using namespace Tungsten;

static const int MASTER = 0;  // mpi master node, for all terminal IO
// the output trajectories will be trj-00001.nc, with a width given here
static const int FILENAME_NUMBER_WIDTH = 5;


int run(int argc, char* argv[], double* totalMDTime) {
    const int size = MPI::COMM_WORLD.Get_size();
    const int rank = MPI::COMM_WORLD.Get_rank();
    // Parse the command line
    if (argc != 5) {
        printf("usage: [mpirun] %s <system.xml> <integrator.xml> <state.xml> <config.ini>\n", argv[0]);
        MPI::Finalize();
        exit(0);
    }

    printUname();

    // Get the config file
    ConfigOpts opts = parseConfigFile(argv[4]);

    // Create the context from the input files
    ifstream systemXml(argv[1]);
    ifstream integratorXml(argv[2]);
    Context* context =  createContext(systemXml, integratorXml, opts.openmmPlatform);
    Integrator& integrator = context->getIntegrator();
    const System& system = context->getSystem();
    int numAtoms = system.getNumParticles();
    bool isPeriodic = hasPeriodicBoundaries(system);
    const double temperature = getTemperature(system, integrator);
    if (temperature <= 0)
        exitWithMessage("tungsten is only compatible with systems under temperature control");

    // Set the initial state
    fstream stateXml(argv[3]);
    State* state = XmlSerializer::deserialize<State>(stateXml);
    context->setState(*state);

    // Create the trajectory for our work on this one
    stringstream s;
    s <<  opts.outputRootPath << "/trj-" << std::setw(FILENAME_NUMBER_WIDTH) << std::setfill('0') << rank << ".nc";
    string fileName = s.str();
    NetCDFTrajectoryFile file(fileName, "w", numAtoms);

    for (int round = 0; round < opts.numRounds; round++) {
        printfM("\n==================================\n");
        printfM("Running Adaptive Sampling Round %d\n", round+1);
        printfM("==================================\n\n");
        file.write(context->getState(State::Positions, isPeriodic));

        time_t roundStartWallTime = time(NULL);
        for (int step = 0; step < opts.numStepsPerRound; step += opts.numStepsPerWrite) {
            integrator.step(opts.numStepsPerWrite);
            file.write(context->getState(State::Positions, isPeriodic));
            file.flush();
        }
        time_t roundEndWallTime = time(NULL);
        // note, using getTime() here relies on the fact that every round, we
        // reset the simulation clock to zero.
        double mdTime = context->getState(0).getTime();
        printfM("MD Performance: ");
        printPerformance(mdTime, roundEndWallTime, roundStartWallTime);
        (*totalMDTime) += mdTime;

        // set up for the next round
        // Run clustering with a strude of 1

        ParallelKCenters clusterer(file, 1, opts.kcentersRmsdIndices);
        clusterer.cluster(opts.kcentersRmsdCutoff, 0, 0);
        printfM("Scattering new starting min-counts starting confs to each rank\n");
        printfM("--------------------------------------------------------------\n");
        ParallelMSM markovModel(clusterer.getAssignments(), clusterer.getCenters());
        gindex newConformation = markovModel.scatterMinCountStates();
        printfAOrd("Rank %d: received frame %d from rank %d\n", rank, newConformation.frame, newConformation.rank);

        PositionsAndPeriodicBox s = file.loadNonlocalStateMPI(newConformation.rank, newConformation.frame);
        context->setPositions(s.positions);
        context->setPeriodicBoxVectors(s.boxA, s.boxB, s.boxC);
        context->setVelocitiesToTemperature(temperature);
        context->setTime(0.0);
    }
    delete context;
}



int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    const int rank = MPI::COMM_WORLD.Get_rank();
    if (rank == MASTER) {
        printf("Tungsten: Parallel Markov State Model Acceleraed Molecular Dynamics\n\n");
        printf("Copyright (C) 2013 Stanford University. This program\n");
        printf("comes with ABSOLUTELY NO WARRANTY. Tungsten is free\n");
        printf("software, and you are welcome to redistribute it under\n");
        printf("certain conditions; see LICENSE.txt for details\n\n");
    }
    double totalMDTime = 0;
    double aggegateMDTime = 0;
    time_t startWallTime = time(NULL);

    try {
        run(argc, argv, &totalMDTime);
    } catch (std::exception e) {
        fflush(stdout);
        fflush(stderr);
        fprintf(stderr, "============================\n");
        fprintf(stderr, "An exception was triggered on RANK=%d\n", rank);
        std::cerr << e.what() << std::endl;
        fprintf(stderr, "============================\n");
        fflush(stderr);
        MPI::COMM_WORLD.Abort(EXIT_FAILURE);
    }

    time_t endWallTime = time(NULL);
    printfM("\nOverall Performance\n");
    printfM("===================\n");
    printPerformance(totalMDTime, endWallTime, startWallTime);
    MPI::COMM_WORLD.Reduce(&totalMDTime, &aggegateMDTime, 1, MPI_DOUBLE, MPI_SUM, MASTER);
    printfM("Total Sampling: %.3f ns\n\n", aggegateMDTime/1000.0);
    printfM("Starfleet Out.\n");
    MPI::Finalize();
}
