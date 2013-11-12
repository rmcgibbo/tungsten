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
// along with Tungsten.  If not, see <http://www.gnu.org/licenses/>.
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

int MPI_allFilesExist(const string& str) {
    int outputFileExists = fileExists(str);
    int allOutputFilesExists;
    bool startWithClustering;
    MPI::COMM_WORLD.Allreduce(&outputFileExists, &allOutputFilesExists, 1, MPI_INT, MPI_LAND);
    if (allOutputFilesExists)
        return true;
    return false;
}

void run(int argc, char* argv[], double* totalMDTime, double* totalWallTimeInMD) {
    const int SIZE = MPI::COMM_WORLD.Get_size();
    const int RANK = MPI::COMM_WORLD.Get_rank();
    // Parse the command line
    if (argc != 5) {
        printfM("usage: [mpirun] %s <system.xml> <integrator.xml> <state.xml|inp.nc> <config.ini>\n", argv[0]);
        MPI::Finalize();
        exit(0);
    }

    printUname();
    fflush(stdout);
    MPI::COMM_WORLD.Barrier();

    // Get the config file
    ConfigOpts opts = parseConfigFile(argv[4]);

    // Create the context from the input files
    ifstream systemXml(argv[1]);
    if (systemXml == NULL) exitWithMessage("No such file or directory: '%s'\n", argv[1]);
    ifstream integratorXml(argv[2]);
    if (integratorXml == NULL) exitWithMessage("No such file or directory: '%s'\n", argv[2]);
    Context* context =  createContext(systemXml, integratorXml, opts.openmmPlatform);
    Integrator& integrator = context->getIntegrator();
    const System& system = context->getSystem();
    int numAtoms = system.getNumParticles();
    bool isPeriodic = hasPeriodicBoundaries(system);
    const double temperature = getTemperature(system, integrator);
    if (temperature <= 0)
        exitWithMessage("tungsten is only compatible with systems under temperature control\n");

    // Set the initial state
    if (endswith(argv[3], ".nc")) {
        printfM("\n(netCDF) Loading starting coordinates from: %s\n", argv[3]);
        NetCDFTrajectoryFile nc(argv[3], "r", numAtoms);
        int numFrames = nc.getNumFrames();
        PositionsAndPeriodicBox s = nc.loadState(RANK % numFrames);
        context->setPositions(s.positions);
        context->setPeriodicBoxVectors(s.boxA, s.boxB, s.boxC);
        context->setVelocitiesToTemperature(temperature);
    } else {
        printfM("\n(OpenMM xml) Loading starting coordinates from: %s\n", argv[3]);
        fstream stateXml(argv[3]);
        if (stateXml == NULL) exitWithMessage("No such file or directory: '%s'\n", argv[1]);
        State* state = XmlSerializer::deserialize<State>(stateXml);
        context->setState(*state);
    }
    // Create per-node output trajectory
    stringstream s;
    s <<  opts.outputRootPath << "/trj-" << std::setw(FILENAME_NUMBER_WIDTH) << std::setfill('0') << RANK << ".nc";
    NetCDFTrajectoryFile file(s.str().c_str(), "w", numAtoms);
    int localStartWithClustering = (file.getNumFrames() > 1);
    int startWithClustering = 0;
    MPI::COMM_WORLD.Allreduce(&localStartWithClustering, &startWithClustering, 1, MPI_INT, MPI_LAND);

    for (int round = 0; round < opts.numRounds; round++) {
        if (startWithClustering) {
            printfM("\n====================================\n");
            printfM("Skipping simulation round %d and\n", round+1);
            printfM("proceeding straight to clustering\n");
            printfM("====================================\n\n");
            startWithClustering = 0;
        } else {
            printfM("\n====================================\n");
            printfM("Running Adaptive Sampling Round %3d\n", round+1);
            printfM("====================================\n\n");
            context->setTime(0.0);
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
            double elapsedWallTime = difftime(roundEndWallTime, roundStartWallTime);
            printfM("MD Performance: ");
            printPerformance(mdTime, elapsedWallTime);
            (*totalMDTime) += mdTime;
            (*totalWallTimeInMD) += elapsedWallTime;
        }

        // set up for the next round
        // Run clustering with a strude of 1
        int numClusters = -1;
        if (opts.kcentersNumClusterMultiplier > 0) {
            int numFrames = file.getNumFrames();
            int totalNumFrames = 0;
            MPI::COMM_WORLD.Allreduce(&numFrames, &totalNumFrames, 1, MPI_INT, MPI_SUM);
            numClusters = (int) (totalNumFrames * opts.kcentersNumClusterMultiplier);
            if (numClusters <= 0)
                numClusters = 1;
            printfM("\nTotal numFrames=%d. numClusters=%d\n", totalNumFrames, numClusters);
        }

        ParallelKCenters clusterer(file, 1, opts.kcentersRmsdIndices);
        clusterer.cluster(opts.kcentersRmsdCutoff, numClusters, 0, 0);
        ParallelMSM markovModel(clusterer.getAssignments(), clusterer.getCenters(), file.loadTime());
        printfM("Scattering new starting min-counts starting confs to each rank\n");
        printfM("--------------------------------------------------------------\n");
        gindex newConformation = markovModel.scatterMinCountStates();

        PositionsAndPeriodicBox s = file.loadNonlocalStateMPI(newConformation.rank, newConformation.frame);
        context->setPositions(s.positions);
        context->setPeriodicBoxVectors(s.boxA, s.boxB, s.boxC);
        context->setVelocitiesToTemperature(temperature);
    }
    delete context;
}



int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    const int RANK = MPI::COMM_WORLD.Get_rank();
    printfM("Tungsten: Parallel Markov State Model Acceleraed Molecular Dynamics\n\n");
    printfM("Copyright (C) 2013 Stanford University. This program\n");
    printfM("comes with ABSOLUTELY NO WARRANTY. Tungsten is free\n");
    printfM("software, and you are welcome to redistribute it under\n");
    printfM("certain conditions; see LICENSE.txt for details\n\n");
    double totalMDTime = 0;
    double totalTimeInMD = 0;
    time_t startWallTime = time(NULL);

    try {
        run(argc, argv, &totalMDTime, &totalTimeInMD);
    } catch (std::exception e) {
        fflush(stdout);
        fflush(stderr);
        fprintf(stderr, "============================\n");
        fprintf(stderr, "An exception was triggered on RANK=%d\n", RANK);
        std::cerr << e.what() << std::endl;
        fprintf(stderr, "============================\n");
        fflush(stderr);
        MPI::COMM_WORLD.Abort(EXIT_FAILURE);
    }

    time_t endWallTime = time(NULL);
    double elapsedWallTime = difftime(endWallTime, startWallTime);

    double aggegateTimeInMD = 0;  // core*seconds in which we were *actually* running MD
    double aggregateWallTime = 0; // core*seconds which the process consumed overall
    MPI::COMM_WORLD.Reduce(&totalTimeInMD, &aggegateTimeInMD, 1, MPI_DOUBLE, MPI_SUM, MASTER);
    MPI::COMM_WORLD.Reduce(&elapsedWallTime, &aggregateWallTime, 1, MPI_DOUBLE, MPI_SUM, MASTER);
    printfM("\nOverall Performance\n");
    printfM("===================\n");
    printPerformance(totalMDTime, elapsedWallTime);
    printfM("Wall-Clock Efficiency (core*hours doing MD / core*hours overall): %.2f%%\n\n",
            100*(aggegateTimeInMD / aggregateWallTime));

    printfM("Starfleet out.\n");
    MPI::Finalize();
}
