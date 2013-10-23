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

#ifndef TUNGSTEN_TYPEDEFS_H_
#define TUNGSTEN_TYPEDEFS_H_
#include <vector>
#include "OpenMM.h"
namespace Tungsten {


struct gindex {
    /**
     * Small struct to hold the globally (MPI) addressable coordinates
     * of an individual simulation frame, refering to both the MPI rank
     * it is owned by and the frame number.
     */
    int rank;
    int frame;

    /**
     * This is required for using the struct as a std::map key.
     */
    bool operator<( const gindex & that ) const {
        if (this->rank < that.rank)
            return true;
        if (this->rank > that.rank)
            return false;
        return this->frame < that.frame;
    }
};

struct PositionsAndPeriodicBox {
    /**
     * The positions and periodic box vectors of an individual simulation
     * frame.
     */
    std::vector<OpenMM::Vec3> positions;
    OpenMM::Vec3 boxA;
    OpenMM::Vec3 boxB;
    OpenMM::Vec3 boxC;
};

} //namespace
#endif
