/**
 * \file:
 */
#ifndef BLOCK_LATTICE_2D_CUH
#define BLOCK_LATTICE_2D_CUH

#if !defined(BLOCK_LATTICE_2D_H)
#error Do not include this file directly. Include 'blockLattice2D.h'.
#endif

#include "core/plbCuda.h"
#include "core/geometry2D.h"
#include <core/cell.h>
namespace plb
{
void device_collide ( Cell<T,Descriptor> **grid, Box2D domain )
{

};
}
#endif
