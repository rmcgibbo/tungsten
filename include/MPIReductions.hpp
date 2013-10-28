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
    int frame;
    float value;
} triplet;

triplet MPIvectorAllMaxloc(const std::vector<float>& input);

cs* MPIcsAdd_efficient(cs* m);

cs* MPIcsAdd(cs* m);

} // namespace Tungsten
#endif