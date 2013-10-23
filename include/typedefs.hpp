#ifndef TUNGSTEN_TYPEDEFS_H_
#define TUNGSTEN_TYPEDEFS_H_
#include <vector>
#include "OpenMM.h"
namespace Tungsten {


struct gindex {
    int rank;
    int frame;

    bool operator<( const gindex & that ) const {
        if (this->rank < that.rank)
            return true;
        if (this->rank > that.rank)
            return false;
        return this->frame < that.frame;
    }
};

struct PositionsAndPeriodicBox {
    std::vector<OpenMM::Vec3> positions;
    OpenMM::Vec3 boxA;
    OpenMM::Vec3 boxB;
    OpenMM::Vec3 boxC;
};

} //namespace
#endif
    
    