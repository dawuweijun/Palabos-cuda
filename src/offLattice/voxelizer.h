/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: HuiJie Zhang */

#ifndef VOXELIZER_H
#define VOXELIZER_H
#include "core/globalDefs.h"
#include "core/plbDebug.h"
namespace plb
{
namespace voxelFlag
{
/// It is a requirement that "undetermined" equals zero, because the
///   default initialization value for scalar-fields is zero.
static const int undetermined = 0;
static const int outside      = 1;
/// Cells which are outside, but which have
///   inside neighbors.
static const int outerBorder  = 2;
static const int inside       = 3;
/// Cells which are inside, but which have
///   outside neighbors.
static const int innerBorder  = 4;

inline int invert ( int arg )
{
    switch ( arg )
    {
    case inside:
        return outside;
    case outside:
        return inside;
    case innerBorder:
        return outerBorder;
    case outerBorder:
        return innerBorder;
    case undetermined:
        return undetermined;
    default:
        PLB_ASSERT ( false );
    }
    return undetermined;
}
inline int bulkFlag ( int arg )
{
    if ( arg==innerBorder || arg==inside )
    {
        return inside;
    }
    else if ( arg==outerBorder || arg==outside )
    {
        return outside;
    }
    else
    {
        return undetermined;
    }
}
inline int borderFlag ( int arg )
{
    if ( arg==inside || arg==innerBorder )
    {
        return innerBorder;
    }
    else if ( arg==outside || arg==outerBorder )
    {
        return outerBorder;
    }
    else
    {
        return undetermined;
    }
}
inline bool insideFlag ( int arg )
{
    return arg==inside || arg==innerBorder;
}
inline bool outsideFlag ( int arg )
{
    return arg==outside || arg==outerBorder;
}

}//namespace voxelFlag
}//namespace plb
#endif
