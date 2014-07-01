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

/** \file
 * Neumann boundary conditions -- generic implementation.
 */
#ifndef NEUMANN_CONDITION_2D_HH
#define NEUMANN_CONDITION_2D_HH

#include "boundaryCondition/neumannCondition2D.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void CopyUnknownPopulationsFunctional2D<T,Descriptor,direction,orientation>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    std::vector<plint> const& unknownIndices
        = indexTemplates::subIndex<Descriptor<T>, direction, -orientation>();
    enum {
        normalX = direction==0 ? orientation : 0,
        normalY = direction==1 ? orientation : 0
    };
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (pluint fIndex=0; fIndex<unknownIndices.size(); ++fIndex) {
                plint iPop = unknownIndices[fIndex];
                lattice.get(iX,iY)[iPop] = lattice.get(iX-normalX, iY-normalY)[iPop];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
CopyUnknownPopulationsFunctional2D<T,Descriptor,direction,orientation>*
    CopyUnknownPopulationsFunctional2D<T,Descriptor,direction,orientation>::clone() const
{
    return new CopyUnknownPopulationsFunctional2D<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void CopyUnknownPopulationsFunctional2D<T,Descriptor,direction,orientation>::
         getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
BlockDomain::DomainT CopyUnknownPopulationsFunctional2D<T,Descriptor,direction,orientation>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyAllPopulationsFunctional2D<T,Descriptor,normalX,normalY>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                lattice.get(iX,iY)[iPop] = lattice.get(iX-normalX, iY-normalY)[iPop];
            }
        }
    }
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
CopyAllPopulationsFunctional2D<T,Descriptor,normalX,normalY>*
    CopyAllPopulationsFunctional2D<T,Descriptor,normalX,normalY>::clone() const
{
    return new CopyAllPopulationsFunctional2D<T,Descriptor,normalX,normalY>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyAllPopulationsFunctional2D<T,Descriptor, normalX,normalY>::
         getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
BlockDomain::DomainT CopyAllPopulationsFunctional2D<T,Descriptor,normalX,normalY>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyVelocityFunctional2D<T,Descriptor,normalX,normalY>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> u;
            lattice.get(iX-normalX, iY-normalY).computeVelocity(u);
            lattice.get(iX,iY).defineVelocity(u);
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
CopyVelocityFunctional2D<T,Descriptor,normalX,normalY>*
    CopyVelocityFunctional2D<T,Descriptor,normalX,normalY>::clone() const
{
    return new CopyVelocityFunctional2D<T,Descriptor,normalX,normalY>(*this);
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyVelocityFunctional2D<T,Descriptor, normalX,normalY>::
         getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::dynamicVariables;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
BlockDomain::DomainT CopyVelocityFunctional2D<T,Descriptor,normalX,normalY>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyTangentialVelocityFunctional2D<T,Descriptor,normalX,normalY>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> u;
            lattice.get(iX-normalX, iY-normalY).computeVelocity(u);
            if (normalX!=0) {
                u[0] = T();
            }
            if (normalY!=0) {
                u[1] = T();
            }
            lattice.get(iX,iY).defineVelocity(u);
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
CopyTangentialVelocityFunctional2D<T,Descriptor,normalX,normalY>*
    CopyTangentialVelocityFunctional2D<T,Descriptor,normalX,normalY>::clone() const
{
    return new CopyTangentialVelocityFunctional2D<T,Descriptor,normalX,normalY>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
BlockDomain::DomainT CopyTangentialVelocityFunctional2D<T,Descriptor,normalX,normalY>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyTangentialVelocityFunctional2D<T,Descriptor, normalX,normalY>::
         getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::dynamicVariables;
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyNormalVelocityFunctional2D<T,Descriptor,normalX,normalY>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> u;
            lattice.get(iX-normalX, iY-normalY).computeVelocity(u);
            if (normalX==0) {
                u[0] = T();
            }
            if (normalY==0) {
                u[1] = T();
            }
            lattice.get(iX,iY).defineVelocity(u);
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
CopyNormalVelocityFunctional2D<T,Descriptor,normalX,normalY>*
    CopyNormalVelocityFunctional2D<T,Descriptor,normalX,normalY>::clone() const
{
    return new CopyNormalVelocityFunctional2D<T,Descriptor,normalX,normalY>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyNormalVelocityFunctional2D<T,Descriptor, normalX,normalY>::
         getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::dynamicVariables;
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
BlockDomain::DomainT CopyNormalVelocityFunctional2D<T,Descriptor,normalX,normalY>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyDensityFunctional2D<T,Descriptor,normalX,normalY>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).defineDensity (
                    lattice.get(iX-normalX, iY-normalY).computeDensity() );
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
CopyDensityFunctional2D<T,Descriptor,normalX,normalY>*
    CopyDensityFunctional2D<T,Descriptor,normalX,normalY>::clone() const
{
    return new CopyDensityFunctional2D<T,Descriptor,normalX,normalY>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
void CopyDensityFunctional2D<T,Descriptor, normalX,normalY>::
         getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::dynamicVariables;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY> 
BlockDomain::DomainT CopyDensityFunctional2D<T,Descriptor,normalX,normalY>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
VirtualOutlet2D<T,Descriptor>::VirtualOutlet2D(T outsideDensity_, Box2D globalDomain_, int type_)
        : outsideDensity(outsideDensity_),
          globalDomain(globalDomain_),
          type(type_)
{
    PLB_ASSERT(type == 0 || type == 1);
}

template<typename T, template<typename U> class Descriptor>
void VirtualOutlet2D<T,Descriptor>::processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice2D<T,Descriptor> *lattice = dynamic_cast<BlockLattice2D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(lattice);
    ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[1]);
    PLB_ASSERT(rhoBar);
    TensorField2D<T,2> *j = dynamic_cast<TensorField2D<T,2>*>(blocks[2]);
    PLB_ASSERT(j);

    Dot2D absOfs = lattice->getLocation();

    Dot2D ofsRB = computeRelativeDisplacement(*lattice, *rhoBar);
    Dot2D ofsJ = computeRelativeDisplacement(*lattice, *j);

    if (type == 0) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
{
                    Cell<T,Descriptor>& cell = lattice->get(iX, iY);
                    for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
                        plint prevX = iX - Descriptor<T>::c[iPop][0];
                        plint prevY = iY - Descriptor<T>::c[iPop][1];
                        
                        if (!contained(prevX + absOfs.x, prevY + absOfs.y, globalDomain)) {
                            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
                            T savedPop = lattice->get(prevX, prevY)[opp];

                            // Velocity is simply taken from the previous time step.
                            Array<T,2> J = j->get(iX + ofsJ.x, iY + ofsJ.y);
                            T Jsqr = dot<T,2>(J, J);
                            // Density is prescribed as a boundary condition.
                            T outsideRhoBar = Descriptor<T>::rhoBar(outsideDensity);
                            T feq_i = cell.computeEquilibrium(iPop, outsideRhoBar, J, Jsqr);
                            T feq_opp_i = cell.computeEquilibrium(opp, outsideRhoBar, J, Jsqr);
                            cell[iPop] = feq_i + feq_opp_i - savedPop;
                        }
                    }
                }
            }
        }
    } else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
{
                    Cell<T,Descriptor>& cell = lattice->get(iX, iY);
                    for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
                        plint prevX = iX - Descriptor<T>::c[iPop][0];
                        plint prevY = iY - Descriptor<T>::c[iPop][1];
                        
                        if (!contained(prevX + absOfs.x, prevY + absOfs.y, globalDomain)) {
                            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
                            T savedPop = lattice->get(prevX, prevY)[opp];

                            // Velocity is simply taken from the previous time step.
                            Array<T,2> J = j->get(iX + ofsJ.x, iY + ofsJ.y);
                            T Jsqr = dot<T,2>(J, J);

                            // For the density there are several choices:
                            // First choice: density from the previous time step.
                            //T RhoBar = rhoBar->get(iX + ofsRB.x, iY + ofsRB.y, iZ +ofsRB.z);

                            // Second choice: Laplacian smoothed value from the previous time step.
                            T RhoBar = 0.0;
                            int n = 0;
                            for (plint dx = -1; dx <= 1; dx++) {
                                plint i = iX + dx;
                                for (plint dy = -1; dy <= 1; dy++) {
                                    plint j = iY + dy;
                                        if (!(dx == 0 && dy == 0)) {
                                            if (contained(i + absOfs.x, j + absOfs.y, globalDomain)) {
                                                RhoBar += rhoBar->get(i + ofsRB.x, j + ofsRB.y);
                                                n++;
                                            }
                                        }
                                    }
                                }
                            if (n != 0) {
                                RhoBar /= n;
                            } else {
                                RhoBar = Descriptor<T>::rhoBar(outsideDensity);
                            }

                            // Third choice: Mean value of the Laplacian smoothed value and the local
                            //               value from the previous time step. This is believed
                            //               to kill checkerboard modes.
                            RhoBar = 0.5 * (RhoBar + rhoBar->get(iX + ofsRB.x, iY + ofsRB.y));

                            T feq_i = cell.computeEquilibrium(iPop, RhoBar, J, Jsqr);
                            T feq_opp_i = cell.computeEquilibrium(opp, RhoBar, J, Jsqr);
                            cell[iPop] = feq_i + feq_opp_i - savedPop;
                        }
                    }
                }
            }
        }
    }
}

}  // namespace plb

#endif  // NEUMANN_CONDITION_2D_HH
