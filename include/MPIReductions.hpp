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

#ifndef TUNGSTEN_MPIREDUCTIONS_H
#define TUNGSTEN_MPIREDUCTIONS_H
#include <vector>
#include "csparse.h"
namespace Tungsten {


typedef struct {
    int rank;
    int index;
    float value;
} triplet;


/**
 * Find the argmax of a distributed set of vectors.
 *
 * @param input       Each MPI rank calls this function with one vector
 * @return triplet    The return value gives the global index (rank, index)
 *                    of the max value, and that value itself.
 */
triplet MPI_vectorAllMaxloc(const std::vector<float>& input);

/**
 * Add a distributed collection of sparse matrices. All of the
 * matrices _must_ have the same dimensions.
 *
 * Both MPI_csAdd and MPI_csAdd_efficient do the same thing, but have
 * different implementations in terms of the order of the communication
 * between the MPI processes.
 *
 * @param m          Each MPI rank calls this function with one sparse matrix
 * @return sum       The result is the sum of the matrices
 */
cs* MPI_csAdd_efficient(cs* m);

/**
 * Add a distributed collection of sparse matrices. All of the
 * matrices _must_ have the same dimensions.
 *
 * Both MPI_csAdd and MPI_csAdd_efficient do the same thing, but have
 * different implementations in terms of the order of the communication
 * between the MPI processes.
 *
 * @param m          Each MPI rank calls this function with one sparse matrix
 * @return sum       The result is the sum of the matrices
 */
cs* MPI_csAdd(cs* m);

} // namespace Tungsten
#endif
