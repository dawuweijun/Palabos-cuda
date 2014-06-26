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

#ifndef FILIPPOVA_HAENEL_2D_HH
#define FILIPPOVA_HAENEL_2D_HH

#include "offLattice/filippovaHaenel2D.h"
#include "offLattice/nextNeighbors2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "core/dynamics.h"
#include <algorithm>
#include <vector>
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
FilippovaHaenelModel2D<T,Descriptor>::FilippovaHaenelModel2D (
        BoundaryShape2D<T,Array<T,2> >* shape_, int flowType_, bool useAllDirections_ )
    : OffLatticeModel2D<T,Array<T,2> >(shape_, flowType_),
      computeStat(true)
{ }

template<typename T, template<typename U> class Descriptor>
FilippovaHaenelModel2D<T,Descriptor>* FilippovaHaenelModel2D<T,Descriptor>::clone() const {
    return new FilippovaHaenelModel2D(*this);
}

template<typename T, template<typename U> class Descriptor>
plint FilippovaHaenelModel2D<T,Descriptor>::getNumNeighbors() const {
    return 1;
}

template<typename T, template<typename U> class Descriptor>
void FilippovaHaenelModel2D<T,Descriptor>::prepareCell (
        Dot2D const& cellLocation,
        AtomicContainerBlock2D& container )
{
    typedef Descriptor<T> D;
    Dot2D offset = container.getLocation();
    OffLatticeInfo2D* info = dynamic_cast<OffLatticeInfo2D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<int> liquidNeighbors;
    std::vector<plint> ids;
    if (!this->isFluid(cellLocation+offset)) {
        for (int iPop=0; iPop<D::q; ++iPop) {
            Dot2D neighbor(cellLocation.x+D::c[iPop][0], cellLocation.y+D::c[iPop][1]);
            // If the non-fluid node has a fluid neighbor ...
            if (this->isFluid(neighbor+offset)) {
                // ... check how many fluid nodes it has ahead of it ...
                plint iTriangle=-1;
                global::timer("intersect").start();
                Array<T,2> locatedPoint;
                T distance;
                Array<T,2> wallNormal;
                Array<T,2> surfaceData;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface (
                            cellLocation+offset, Dot2D(D::c[iPop][0],D::c[iPop][1]), locatedPoint, distance,
                            wallNormal, surfaceData, bdType, iTriangle );
                // In the following, the importance of directions is sorted wrt. how well they
                //   are aligned with the wall normal. It is better to take the continuous normal,
                //   because it is not sensitive to the choice of the triangle when we shoot at
                //   an edge.
                //wallNormal = this->computeContinuousNormal(locatedPoint, iTriangle);
                global::timer("intersect").stop();
                PLB_ASSERT( ok );
                // ... then add this node to the list.
                liquidNeighbors.push_back(iPop);
                ids.push_back(iTriangle);
            }
        }
        if (!liquidNeighbors.empty()) {
            info->getDryNodes().push_back(cellLocation);
            info->getDryNodeFluidDirections().push_back(liquidNeighbors);
            info->getDryNodeIds().push_back(ids);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ContainerBlockData*
    FilippovaHaenelModel2D<T,Descriptor>::generateOffLatticeInfo() const
{
    return new OffLatticeInfo2D;
}

template<typename T, template<typename U> class Descriptor>
Array<T,2> FilippovaHaenelModel2D<T,Descriptor>::getLocalForce (
                AtomicContainerBlock2D& container ) const
{
    OffLatticeInfo2D* info =
        dynamic_cast<OffLatticeInfo2D*>(container.getData());
    PLB_ASSERT( info );
    return info->getLocalForce();
}

template<typename T, template<typename U> class Descriptor>
void FilippovaHaenelModel2D<T,Descriptor>::boundaryCompletion (
        AtomicBlock2D& nonTypeLattice,
        AtomicContainerBlock2D& container,
        std::vector<AtomicBlock2D const*> const& args )
{
    BlockLattice2D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice2D<T,Descriptor>&> (nonTypeLattice);
    OffLatticeInfo2D* info =
        dynamic_cast<OffLatticeInfo2D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<Dot2D> const&
        dryNodes = info->getDryNodes();
    std::vector<std::vector<int > > const&
        dryNodeFluidDirections = info->getDryNodeFluidDirections();
    std::vector<std::vector<plint> > const&
        dryNodeIds = info->getDryNodeIds();
    PLB_ASSERT( dryNodes.size() == dryNodeFluidDirections.size() );

    Dot2D absoluteOffset = lattice.getLocation();

    Array<T,2>& localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iDry=0; iDry<dryNodes.size(); ++iDry) {
        cellCompletion (
            lattice, dryNodes[iDry], dryNodeFluidDirections[iDry],
            dryNodeIds[iDry], absoluteOffset, localForce, args );
    }
}

template<typename T, template<typename U> class Descriptor>
void FilippovaHaenelModel2D<T,Descriptor>::cellCompletion (
        BlockLattice2D<T,Descriptor>& lattice, Dot2D const& guoNode,
        std::vector<int> const& dryNodeFluidDirections,
        std::vector<plint> const& dryNodeIds, Dot2D const& absoluteOffset,
        Array<T,2>& localForce, std::vector<AtomicBlock2D const*> const& args )
{
    typedef Descriptor<T> D;
    Cell<T,Descriptor>& s_cell =
        lattice.get( guoNode.x, guoNode.y);
    int noDynId = NoDynamics<T,Descriptor>().getId();
    PLB_ASSERT( s_cell.getDynamics().getId() == noDynId );
    for (plint iDirection=0; iDirection<(plint)dryNodeFluidDirections.size(); ++iDirection)
    {
        int iOpp = dryNodeFluidDirections[iDirection];
        int iPop = indexTemplates::opposite<Descriptor<T> >(iOpp);
        Dot2D fluidDirection(D::c[iOpp][0],D::c[iOpp][1]);
        plint dryNodeId = dryNodeIds[iDirection];

        Array<T,2> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;

        Cell<T,Descriptor> const& f_cell =
            lattice.get( guoNode.x+fluidDirection.x,
                         guoNode.y+fluidDirection.y);

        Cell<T,Descriptor> const& ff_cell =
            lattice.get( guoNode.x+2*fluidDirection.x,
                         guoNode.y+2*fluidDirection.y);

        Cell<T,Descriptor> collidedCell(f_cell);
        BlockStatistics statsCopy(lattice.getInternalStatistics());
        collidedCell.collide(statsCopy);

        T f_rhoBar, ff_rhoBar;
        Array<T,2> f_j, ff_j;
        Array<T,2> wallNormal;
        f_cell.getDynamics().computeRhoBarJ(f_cell, f_rhoBar, f_j);
        ff_cell.getDynamics().computeRhoBarJ(ff_cell, ff_rhoBar, ff_j);
        T f_rho = D::fullRho(f_rhoBar);
        T f_jSqr = normSqr(f_j);

#ifdef PLB_DEBUG
        bool ok =
#endif
        this->pointOnSurface( guoNode+absoluteOffset, fluidDirection,
                              wallNode, wallDistance, wallNormal,
                              wall_vel, bdType, dryNodeId );
        PLB_ASSERT( ok );

        Array<T,2> w_j = wall_vel*f_rho;
        T d = sqrt(D::cNormSqr[iOpp]);
        PLB_ASSERT( wallDistance <= d );
        T delta = 1.0-wallDistance / d;

        T kappa = 0.;
        Array<T,2> wf_j; wf_j.resetToZero();
        T omega = f_cell.getDynamics().getOmega();


        if (delta<0.5) {
            //wf_j = f_j;
            //kappa = (omega*(2.0*delta-1.0))/(1.0-omega);

            wf_j = ff_j;
            kappa = (omega*(2*delta-1.0))/(1.0-2.0*omega);
        }
        else {
            //wf_j = f_j * ((delta-1.0)/delta) + w_j / delta;
            //kappa = omega*(2.0*delta-1.0)

            wf_j = (1.0-3.0/(2.0*delta))*f_j+3.0/(2.0*delta)*w_j;
            kappa = (2.0*omega*(2.0*delta-1.0))/(2.0+omega);
        }

        T c_i_wf_j_f_j = D::c[iPop][0]*(wf_j[0]-f_j[0]) +
                         D::c[iPop][1]*(wf_j[1]-f_j[1]) ;

        T c_i_w_j     = D::c[iPop][0]*w_j[0] + D::c[iPop][1]*w_j[1];


        T f_ieq = f_cell.getDynamics().computeEquilibrium(iPop, f_rhoBar, f_j, f_jSqr) +
                      D::t[iPop]*D::invCs2*c_i_wf_j_f_j;
        s_cell[iOpp] = (1.0-kappa)*collidedCell[iPop]+kappa*f_ieq+2.0*D::t[iPop]*D::invCs2*c_i_w_j;

        localForce[0] += D::c[iPop][0]*(s_cell[iPop]+s_cell[iOpp]);
        localForce[1] += D::c[iPop][1]*(s_cell[iPop]+s_cell[iOpp]);
    }
}

}  // namespace plb

#endif  // FILIPPOVA_HAENEL_2D_HH

