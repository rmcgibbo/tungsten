// Copyright 2013 Robert McGibbon
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

