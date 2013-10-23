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

#define MASTER 0
#define WIDTH 5    // used to create the output format

using std::string;
using std::vector;
using std::ifstream;
using std::fstream;
using std::stringstream;
using namespace Tungsten;


int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    const int size = MPI::COMM_WORLD.Get_size();
    const int rank = MPI::COMM_WORLD.Get_rank();

    // Parse the command line
    if (argc != 5)
        exitWithMessage("usage: %s <system.xml> <integrator.xml> <state.xml> <config.ini>\n", argv[0]);

    printUname();

    // Get the config file
    ConfigOpts opts = parseConfigFile(argv[4]);

    // Create the context from the input files
    ifstream systemXml(argv[1]);
    ifstream integratorXml(argv[2]);
    OpenMM::Context* context =  createContext(systemXml, integratorXml, opts.openmmPlatform);
    OpenMM::Integrator& integrator = context->getIntegrator();
    const OpenMM::System& system = context->getSystem();
    int numAtoms = system.getNumParticles();
    bool isPeriodic = hasPeriodicBoundaries(system);
    const double temperature = getTemperature(system, integrator);
    if (temperature <= 0)
        exitWithMessage("tungsten is only compatible with systems under temperature control");

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
        printfM("\nBeginning Round %d\n===================\n", round);

        file.write(context->getState(OpenMM::State::Positions, isPeriodic));
        for (int step = 0; step < opts.numStepsPerRound; step += opts.numStepsPerWrite) {
            if (rank == MASTER)
                printf("#");
            try {
                integrator.step(opts.numStepsPerWrite);
            } catch (OpenMM::OpenMMException e) {
                printf("An exception occured on rank=%d: %s\n", rank, e.what());
                exitWithMessage("Exit Failure");
            }
            file.write(context->getState(OpenMM::State::Positions, isPeriodic));
        }
        printfM("\n");
        file.flush();

        // Run clustering with a strude of 1
        ParallelKCenters clusterer(file, 1, opts.kcentersRmsdIndices);
        clusterer.cluster(opts.kcentersRmsdCutoff, 0, 0);
        ParallelMSM markovModel(clusterer.getAssignments(), clusterer.getCenters());
        gindex newConformation = markovModel.scatterMinCountStates();

        MPI::COMM_WORLD.Barrier();
        fflush(stdout);
        printf("New Coordinates on Rank=%d: (%d, %d)\n", rank,
               newConformation.rank, newConformation.frame);

        PositionsAndPeriodicBox s = file.loadNonlocalStateMPI(
                                        newConformation.rank, newConformation.frame);
        context->setPositions(s.positions);
        context->setPeriodicBoxVectors(s.boxA, s.boxB, s.boxC);
        printfM("Randomizing velocities to T=%.2f\n", temperature);
        context->setVelocitiesToTemperature(temperature);
        printfM("Setting simulation clock back to 0");
        context->setTime(0.0);


    }
    delete context;
    MPI::Finalize();
}
