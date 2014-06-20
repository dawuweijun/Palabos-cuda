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

#ifndef SPMRevProcessor3D_HH
#define SPMRevProcessor3D_HH
#include "multiPhysics/REVProcessor3D.h"
#include <atomicBlock/dataField3D.h>
namespace plb
{
template<typename T1, template<typename U> class Descriptor, typename T2>
void BrinkmanProcessor3D<T1,Descriptor,T2>::process ( Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
        TensorField3D<T2,9>& negNiuInvsK )
{
    enum
    {
        forceOffset = Descriptor<T1>::ExternalField::forceBeginsAt
    };
    Dot3D offset = computeRelativeDisplacement ( lattice, negNiuInvsK );
    Array< T1, 3  > vel;
    T1 *force;
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            for ( plint iZ=domain.z0; iZ<=domain.z1; ++iZ )
            {

                Array< T2, 9 > &negNiuInvsK_= negNiuInvsK.get ( offset.x+iX,offset.y+iY,offset.z+iZ );
                lattice.get ( iX,iY,iZ ).computeVelocity ( vel );

                force = lattice.get ( iX,iY,iZ ).getExternal ( forceOffset );

                force[0]=negNiuInvsK_[0]*vel[0]+negNiuInvsK_[1]*vel[1]+negNiuInvsK_[2]*vel[2];
                force[1]=negNiuInvsK_[3]*vel[0]+negNiuInvsK_[4]*vel[1]+negNiuInvsK_[5]*vel[2];
                force[2]=negNiuInvsK_[6]*vel[0]+negNiuInvsK_[7]*vel[1]+negNiuInvsK_[8]*vel[2];
            }
        }
    }
}

template<typename T1, template<typename U> class Descriptor, typename T2>
BrinkmanProcessor3D<T1,Descriptor,T2> *
BrinkmanProcessor3D<T1,Descriptor,T2>:: clone() const
{
    return new BrinkmanProcessor3D<T1,Descriptor,T2>();
}

/****************************************************************************/
template<typename T, template<typename U> class Descriptor>
void BrinkmanProcessor3DL<T,Descriptor>::process ( Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    enum
    {
        forceOffset = Descriptor<T>::ExternalField::forceBeginsAt
    };
    Array< T, 3  > vel;
    T *force;
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            for ( plint iZ=domain.z0; iZ<=domain.z1; ++iZ )
            {
                lattice.get ( iX,iY,iZ ).computeVelocity ( vel );
                force = lattice.get ( iX,iY,iZ ).getExternal ( forceOffset );

                force[0]=negNiuInvsK[0]*vel[0]+negNiuInvsK[1]*vel[1]+negNiuInvsK[2]*vel[2];
                force[1]=negNiuInvsK[3]*vel[0]+negNiuInvsK[4]*vel[1]+negNiuInvsK[5]*vel[2];
                force[2]=negNiuInvsK[6]*vel[0]+negNiuInvsK[7]*vel[1]+negNiuInvsK[8]*vel[2];
            }
        }
    }
}

}// namespace plb
#endif
