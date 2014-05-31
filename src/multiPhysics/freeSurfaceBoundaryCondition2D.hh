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

#ifndef FREE_SURFACE_BOUNDARY_CONDITION_2D_HH
#define FREE_SURFACE_BOUNDARY_CONDITION_2D_HH

#include "multiPhysics/freeSurfaceBoundaryCondition2D.h"
#include "multiPhysics/freeSurfaceUtil2D.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include <cmath>
#include <iostream>

namespace plb
{

template<typename T, template<typename U> class Descriptor>
FreeSurfaceFadingArea2D<T,Descriptor>::FreeSurfaceFadingArea2D ( T factor_ )
    : factor ( factor_ )
{ }

template<typename T, template<typename U> class Descriptor>
void FreeSurfaceFadingArea2D<T,Descriptor>::process ( Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    std::vector<T> decomposedVariables;

    enum
    {
        forceOffset          = Descriptor<T>::ExternalField::forceBeginsAt,
        momentumStoredOffset = Descriptor<T>::ExternalField::momentumBeginsAt,
        densityStoredOffset  = Descriptor<T>::ExternalField::densityBeginsAt,
    };

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Cell<T,Descriptor>& cell = lattice.get ( iX,iY );
            plint order = 0;
            cell.getDynamics().decompose ( cell, decomposedVariables, order );

            T density = Descriptor<T>::fullRho ( decomposedVariables[0] );
            if ( density > T ( 0 ) ) density *= factor;
            decomposedVariables[0] = Descriptor<T>::rhoBar ( density );
            cell.getDynamics().recompose ( cell, decomposedVariables, order );

            *cell.getExternal ( densityStoredOffset ) = density;

            Array<T,Descriptor<T>::d> j;
            j.resetToZero();
            T rhoBar;
            momentTemplates<T,Descriptor>::get_rhoBar_j ( cell, rhoBar, j );

            // TODO: What about mass, volumeFraction, flagStatus?
            j.to_cArray ( cell.getExternal ( momentumStoredOffset ) );
        }
    }
}


template<typename T, template<typename U> class Descriptor>
FreeSurfaceFadingArea2D<T,Descriptor>* FreeSurfaceFadingArea2D<T,Descriptor>::clone() const
{
    return new FreeSurfaceFadingArea2D<T,Descriptor> ( *this );
}


template<typename T, template<typename U> class Descriptor>
void RemoveMass2D<T,Descriptor>::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            //param.attributeDynamics(iX,iY, new NoDynamics<T,Descriptor>());
            param.setDensity ( iX,iY, ( T ) 1. );
            param.setMomentum ( iX,iY,Array<T,3> ( 0.,0.,0. ) );
            param.mass ( iX,iY ) = ( T ) 0;
            param.volumeFraction ( iX,iY ) = ( T ) 0;
            //param.flag(iX,iY) = empty;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
RemoveMass2D<T,Descriptor>* RemoveMass2D<T,Descriptor>::clone() const
{
    return new RemoveMass2D<T,Descriptor> ( *this );
}


template<typename T, template<typename U> class Descriptor>
void PouringLiquid2D<T,Descriptor>::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            T iniRho = T ( 1 );
            param.attributeDynamics (
                iX,iY, dynamicsTemplate->clone() );
            iniCellAtEquilibrium ( param.cell ( iX,iY ), iniRho, injectionVelocity );
            param.setDensity ( iX,iY, iniRho );
            param.setMomentum ( iX,iY, iniRho*injectionVelocity );
            param.mass ( iX,iY ) = iniRho;
            param.volumeFraction ( iX,iY ) = ( T ) 1;
            param.flag ( iX,iY ) = fluid;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
PouringLiquid2D<T,Descriptor>* PouringLiquid2D<T,Descriptor>::clone() const
{
    return new PouringLiquid2D<T,Descriptor> ( *this );
}


template<typename T, template<typename U> class Descriptor>
void ShortenBounceBack2D<T,Descriptor>::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    typedef Descriptor<T> D;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    Box2D extDomain = domain.enlarge ( 1 );

    for ( plint iX=extDomain.x0; iX<=extDomain.x1; ++iX )
    {
        bool xBoundary = iX==extDomain.x0 || iX==extDomain.x1;
        for ( plint iY=extDomain.y0; iY<=extDomain.y1; ++iY )
        {
            bool yBoundary = xBoundary || iY==extDomain.y0 || iY==extDomain.y1;
            if ( param.flag ( iX,iY ) ==wall )
            {
                for ( plint iNeighbor=1; iNeighbor<D::q; ++iNeighbor )
                {
                    plint nextX = iX+D::c[iNeighbor][0];
                    plint nextY = iY+D::c[iNeighbor][1];
                    if ( contained ( nextX,nextY, domain ) )
                    {
                        if ( isWet ( param.flag ( nextX,nextY ) ) )
                        {
                            plint opp = indexTemplates::opposite<D> ( iNeighbor );
                            param.cell ( nextX,nextY ) [iNeighbor] = param.cell ( iX,iY ) [opp];
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ShortenBounceBack2D<T,Descriptor>* ShortenBounceBack2D<T,Descriptor>::clone() const
{
    return new ShortenBounceBack2D<T,Descriptor> ( *this );
}

}  // namespace plb

#endif  // FREE_SURFACE_BOUNDARY_CONDITION_2D_HH

