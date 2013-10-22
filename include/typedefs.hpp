#ifndef _TUNGSTEN_TYPEDEFS_H_
#define _TUNGSTEN_TYPEDEFS_H_
#include "aligned_allocator.hpp"
namespace Tungsten {

typedef std::vector<float, aligned_allocator<float, 16> > fvector16;

}
#endif
