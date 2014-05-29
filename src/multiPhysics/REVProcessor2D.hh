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

#ifndef REV_PROCESSORS_2D_HH
#define REV_PROCESSORS_2D_HH
#include "multiPhysics/REVProcessor2D.h"
#include <atomicBlock/dataField2D.h>

namespace plb
{

template<typename T1, template<typename U> class Descriptor, typename T2>
void BrinkmanProcessor2D<T1,Descriptor,T2>::process ( Box2D domain, BlockLattice2D<T1,Descriptor>& lattice,
        TensorField2D<T2,4>& negNiuInvsK )
{
    enum
    {
        forceOffset = Descriptor<T1>::ExternalField::forceBeginsAt
    };
    Dot2D offset = computeRelativeDisplacement ( lattice, negNiuInvsK );
    Array< T1, 2  > vel;
    T1 *force;
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            Array< T2, 4  > &negNiuInvsK_=negNiuInvsK.get ( offset.x+iX,offset.y+iY );

            lattice.get ( iX,iY ).computeVelocity ( vel );
            force = lattice.get ( iX,iY ).getExternal ( forceOffset );

            force[0]= negNiuInvsK_[0]*vel[0]+negNiuInvsK_[1]*vel[1];
            force[1]= negNiuInvsK_[2]*vel[0]+negNiuInvsK_[3]*vel[1];
        }
    }
}

template<typename T1, template<typename U> class Descriptor, typename T2>
BrinkmanProcessor2D<T1,Descriptor,T2> *
BrinkmanProcessor2D<T1,Descriptor,T2>:: clone() const
{
    return new BrinkmanProcessor2D<T1,Descriptor,T2>();
}
}// namespace plb
#endif
