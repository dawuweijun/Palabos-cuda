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

            Array< T2, 4  > &localP=negNiuInvsK.get ( offset.x+iX,offset.y+iY );

            lattice.get ( iX,iY ).computeVelocity ( vel );
            force = lattice.get ( iX,iY ).getExternal ( forceOffset );

            force[0]= localP[0]*vel[0]+localP[1]*vel[1];
            force[1]= localP[2]*vel[0]+localP[3]*vel[1];
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
