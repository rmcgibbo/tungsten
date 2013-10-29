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
// along with Tungsten. If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013 Stanford University
// Author: Robert McGibbon
// Contributors:

#include "mpi.h"
#include "stdint.h"
#include <vector>
#include "MPIReductions.hpp"

namespace Tungsten {

using std::vector;
static int fastlog2(uint32_t v);
static const int MASTER = 0;

//////////////////////////////////////////////////////////////////////////////
////                        Global Argmax                                 ////
//////////////////////////////////////////////////////////////////////////////

/*
 * MaxLoc Reduction. Each rank provides input data, and the return
 * value, on each node, is a triplet containing the rank, index and value of
 * the maximum entry. It's a global argmax.
 */
triplet MPIvectorAllMaxloc(const vector<float>& input) {
    static const int rank = MPI::COMM_WORLD.Get_rank();
    struct {
        float value;
        int   index;
    } localMaxLoc, globalMaxLoc;
    int trajLength = input.size();
    int maxTrajLength;

    // Share the longest trajectory length with all of the nodes so that
    // each can properly compute a unique index of
    // localIndex + rank*maxLengthOfAnyLocalIndex
    MPI::COMM_WORLD.Allreduce(&trajLength, &maxTrajLength, 1, MPI_INT, MPI_MAX);

    // local maxloc
    localMaxLoc.value = input[0];
    localMaxLoc.index = 0;
    for (int i = 1; i < trajLength; i++)
        if (localMaxLoc.value < input[i]) {
            localMaxLoc.value = input[i];
            localMaxLoc.index = i;
        }

    // give the local maxloc a globaly-resolvably index
    localMaxLoc.index = localMaxLoc.index + rank * maxTrajLength;

    // global maxloc
    MPI::COMM_WORLD.Allreduce(&localMaxLoc, &globalMaxLoc, 1, MPI_FLOAT_INT, MPI_MAXLOC);
    int outRank = globalMaxLoc.index / maxTrajLength;
    int outFrame = globalMaxLoc.index % maxTrajLength;
    float outValue = globalMaxLoc.value;

    triplet t = {outRank, outFrame, outValue};
    return t;
}


//////////////////////////////////////////////////////////////////////////////
////                     SUMMING sparse matrices                         ////
//////////////////////////////////////////////////////////////////////////////


cs* MPIcsAdd_efficient(cs* m) {
    const int SIZE = MPI::COMM_WORLD.Get_size();
    const int RANK = MPI::COMM_WORLD.Get_rank();
    const int numStates = m->n;
    const int lastpower = 1 << fastlog2(SIZE);
    MPI::Request rP, rI, rX;

    vector<int> Nzmax(SIZE);
    MPI::COMM_WORLD.Allgather(&m->nzmax, 1, MPI_INT, &Nzmax[0], 1, MPI_INT);

    // each of the ranks greater than the last power of 2 less than size
    // need to downshift their data, since the binary tree reduction below
    // only works when N is a power of two.
    for (int i = lastpower; i < SIZE; i++)
        if (RANK == i) {
            MPI::COMM_WORLD.Isend(m->p, m->n+1,   MPI_INT,    i-lastpower, 0);
            MPI::COMM_WORLD.Isend(m->i, m->nzmax, MPI_INT,    i-lastpower, 1);
            MPI::COMM_WORLD.Isend(m->x, m->nzmax, MPI_DOUBLE, i-lastpower, 2);
        }
    for (int i = 0; i < SIZE-lastpower; i++)
        if (RANK == i) {
            vector<int>    p(numStates+1);
            vector<int>   ii(Nzmax[i+lastpower]);
            vector<double> x(Nzmax[i+lastpower]);
            rP = MPI::COMM_WORLD.Irecv(&p[0],  numStates+1,        MPI_INT,    i+lastpower, 0);
            rI = MPI::COMM_WORLD.Irecv(&ii[0], Nzmax[i+lastpower], MPI_INT,    i+lastpower, 1);
            rX = MPI::COMM_WORLD.Irecv(&x[0],  Nzmax[i+lastpower], MPI_DOUBLE, i+lastpower, 2);
            rP.Wait();
            rI.Wait();
            rX.Wait();

            cs M;
            M.nzmax = Nzmax[i+lastpower];
            M.m = numStates;
            M.n = numStates;
            M.p = &p[0];
            M.i = &ii[0];
            M.x = &x[0];
            M.nz = -1;
            cs* newCounts = cs_add(m, &M, 1.0, 1.0);
            cs_free(m);
            m = newCounts;
        }

    // Now to the binary tree reduction
    for (int d = 0; d < fastlog2(lastpower); d++) {
        for (int k = 0; k < lastpower; k += 1 << (d + 1)) {
            const int receiver = k;
            const int sender = k + (1 << d);

            if (RANK == sender) {
                MPI::COMM_WORLD.Isend(m->p, m->n+1,   MPI_INT,    receiver, 0);
                MPI::COMM_WORLD.Isend(m->i, m->nzmax, MPI_INT,    receiver, 1);
                MPI::COMM_WORLD.Isend(m->x, m->nzmax, MPI_DOUBLE, receiver, 2);
            } else if (RANK == receiver) {
                vector<int>    p(numStates+1);
                vector<int>   ii(Nzmax[sender]);
                vector<double> x(Nzmax[sender]);
                rP = MPI::COMM_WORLD.Irecv(&p[0],  numStates+1,   MPI_INT,    sender, 0);
                rI = MPI::COMM_WORLD.Irecv(&ii[0], Nzmax[sender], MPI_INT,    sender, 1);
                rX = MPI::COMM_WORLD.Irecv(&x[0],  Nzmax[sender], MPI_DOUBLE, sender, 2);
                rP.Wait();
                rI.Wait();
                rX.Wait();

                cs M;
                M.nzmax = Nzmax[sender];
                M.m = numStates;
                M.n = numStates;
                M.p = &p[0];
                M.i = &ii[0];
                M.x = &x[0];
                M.nz = -1;
                cs* newCounts = cs_add(m, &M, 1.0, 1.0);
                cs_free(m);
                m = newCounts;
            }
        }
    }
    return m;
}


cs* MPIcsAdd(cs* m) {
    const int SIZE = MPI::COMM_WORLD.Get_size();
    const int RANK = MPI::COMM_WORLD.Get_rank();

    // The matrices on each rank must have the same dimensions, and must be
    // squre
    const int numStates = m->n;

    // Gather nzmax on root, the maximum number of entries
    // each each rank's countsMatrix_
    vector<int> rootNzmax(SIZE);
    MPI::COMM_WORLD.Gather(&m->nzmax, 1, MPI_INT, &rootNzmax[0], 1, MPI_INT, MASTER);

    cs* newCounts;
    if (RANK != MASTER) {
        // All of the slave nodes send their buffers to to MASTER
        // for accumulation
        MPI::COMM_WORLD.Isend(m->p, m->n+1, MPI_INT, MASTER, 0);
        MPI::COMM_WORLD.Isend(m->i, m->nzmax, MPI_INT, MASTER, 1);
        MPI::COMM_WORLD.Isend(m->x, m->nzmax, MPI_DOUBLE, MASTER, 2);
    } else {
        for (int j = 1; j < SIZE; j++) {
            // The master node receives these entries and uses them to
            // reconstruct a sparse matrix, using cs_add to then
            // add it to its own.
            vector<int> p(numStates+1);
            vector<int> i(rootNzmax[j]);
            vector<double> x(rootNzmax[j]);
            MPI::Request rP, rI, rX;
            rP = MPI::COMM_WORLD.Irecv(&p[0], numStates+1, MPI_INT, j, 0);
            rI = MPI::COMM_WORLD.Irecv(&i[0], rootNzmax[j], MPI_INT, j, 1);
            rX = MPI::COMM_WORLD.Irecv(&x[0], rootNzmax[j], MPI_DOUBLE, j, 2);
            rP.Wait();
            rI.Wait();
            rX.Wait();

            // place this data in a struct
            cs M;
            M.nzmax = rootNzmax[j];
            M.m = numStates;
            M.n = numStates;
            M.p = &p[0];
            M.i = &i[0];
            M.x = &x[0];
            M.nz = -1;

            newCounts = cs_add(m, &M, 1.0, 1.0);
            cs_free(m);
            m = newCounts;
        }
    }
    return m;
}


static int fastlog2(uint32_t v) {
    // http://graphics.stanford.edu/~seander/bithacks.html
    int r;
    static const int MultiplyDeBruijnBitPosition[32] = {
        0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };

    v |= v >> 1; // first round down to one less than a power of 2
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;

    r = MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
    return r;
}

}  // namespace Tungsten
