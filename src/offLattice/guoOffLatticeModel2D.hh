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

#ifndef GUO_OFF_LATTICE_MODEL_2D_HH
#define GUO_OFF_LATTICE_MODEL_2D_HH

#include "offLattice/guoOffLatticeModel2D.h"
#include "offLattice/nextNeighbors2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include <algorithm>
#include <vector>
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
GuoOffLatticeModel2D<T,Descriptor>::LiquidNeighbor::LiquidNeighbor
            (plint iNeighbor_, plint depth_, plint iSegment_, Array<T,2> wallNormal)
        : iNeighbor(iNeighbor_),
          depth(depth_),
          iSegment(iSegment_)
{
    int const* c = NextNeighbor<T>::c[iNeighbor];
    Array<T,2> neighborVect(c[0],c[1]);
    cosAngle = fabs(dot(neighborVect,wallNormal))*NextNeighbor<T>::invD[iNeighbor];
}

template<typename T, template<typename U> class Descriptor>
bool GuoOffLatticeModel2D<T,Descriptor>::LiquidNeighbor::
         operator<(LiquidNeighbor const& rhs) const
{
    return cosAngle < rhs.cosAngle;
}

template<typename T, template<typename U> class Descriptor>
GuoOffLatticeModel2D<T,Descriptor>::GuoOffLatticeModel2D (
        BoundaryShape2D<T,Array<T,2> >* shape_, int flowType_, bool useAllDirections_ )
    : OffLatticeModel2D<T,Array<T,2> >(shape_, flowType_),
      useAllDirections(useAllDirections_),
      regularizedModel(true),
      secondOrderFlag(true),
      computeStat(true)
{ }

template<typename T, template<typename U> class Descriptor>
GuoOffLatticeModel2D<T,Descriptor>* GuoOffLatticeModel2D<T,Descriptor>::clone() const {
    return new GuoOffLatticeModel2D(*this);
}

template<typename T, template<typename U> class Descriptor>
plint GuoOffLatticeModel2D<T,Descriptor>::getNumNeighbors() const {
    return 2;
}

template<typename T, template<typename U> class Descriptor>
void GuoOffLatticeModel2D<T,Descriptor>::prepareCell (
        Dot2D const& cellLocation,
        AtomicContainerBlock2D& container )
{
    Dot2D offset = container.getLocation();
    GuoOffLatticeInfo2D* info =
        dynamic_cast<GuoOffLatticeInfo2D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<LiquidNeighbor> liquidNeighbors;
    if (!this->isFluid(cellLocation+offset)) {
        for (int iNeighbor=0; iNeighbor<NextNeighbor<T>::numNeighbors; ++iNeighbor) {
            int const* c = NextNeighbor<T>::c[iNeighbor];
            Dot2D neighbor(cellLocation.x+c[0], cellLocation.y+c[1]);
            // If the non-fluid node has a fluid neighbor ...
            if (this->isFluid(neighbor+offset)) {
                // ... check how many fluid nodes it has ahead of it ...
                int depth = 1;
                for (int iDepth=2; iDepth<=getNumNeighbors(); ++iDepth) {
                    Dot2D nextNeighbor(cellLocation.x+iDepth*c[0],
                                       cellLocation.y+iDepth*c[1]);
                    if (this->isFluid(nextNeighbor+offset)) {
                        depth = iDepth;
                    }
                    else {
                        break;
                    }
                }
                plint iSegment=-1;
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
                            cellLocation+offset, Dot2D(c[0],c[1]), locatedPoint, distance,
                            wallNormal, surfaceData, bdType, iSegment );
                // In the following, the importance of directions is sorted wrt. how well they
                //   are aligned with the wall normal. It is better to take the continuous normal,
                //   because it is not sensitive to the choice of the triangle when we shoot at
                //   an edge.
                //wallNormal = this->computeContinuousNormal(locatedPoint, iSegment);
                global::timer("intersect").stop();
                PLB_ASSERT( ok );
                // ... then add this node to the list.
                liquidNeighbors.push_back(LiquidNeighbor(iNeighbor, depth, iSegment, wallNormal));
            }
        }
        if (!liquidNeighbors.empty()) {
            info->getDryNodes().push_back(cellLocation);
            std::sort(liquidNeighbors.begin(), liquidNeighbors.end());
            std::vector<std::pair<int,int> > neighborDepthPairs;
            std::vector<plint> ids;
            if (useAllDirections) {
                for (pluint i=0; i<liquidNeighbors.size(); ++i) {
                    neighborDepthPairs.push_back(std::make_pair(liquidNeighbors[i].iNeighbor, liquidNeighbors[i].depth));
                    ids.push_back(liquidNeighbors[i].iSegment);
                }
            }
            else {
                plint i = liquidNeighbors.size()-1;
                neighborDepthPairs.push_back(std::make_pair(liquidNeighbors[i].iNeighbor, liquidNeighbors[i].depth));
                ids.push_back(liquidNeighbors[i].iSegment);
            }
            info->getDryNodeFluidDirections().push_back(neighborDepthPairs);
            info->getDryNodeIds().push_back(ids);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ContainerBlockData*
    GuoOffLatticeModel2D<T,Descriptor>::generateOffLatticeInfo() const
{
    return new GuoOffLatticeInfo2D;
}

template<typename T, template<typename U> class Descriptor>
Array<T,2> GuoOffLatticeModel2D<T,Descriptor>::getLocalForce (
                AtomicContainerBlock2D& container ) const
{
    GuoOffLatticeInfo2D* info =
        dynamic_cast<GuoOffLatticeInfo2D*>(container.getData());
    PLB_ASSERT( info );
    return info->getLocalForce();
}

template<typename T, template<typename U> class Descriptor>
void GuoOffLatticeModel2D<T,Descriptor>::boundaryCompletion (
        AtomicBlock2D& nonTypeLattice,
        AtomicContainerBlock2D& container,
        std::vector<AtomicBlock2D const*> const& args )
{
    BlockLattice2D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice2D<T,Descriptor>&> (nonTypeLattice);
    GuoOffLatticeInfo2D* info =
        dynamic_cast<GuoOffLatticeInfo2D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<Dot2D> const&
        dryNodes = info->getDryNodes();
    std::vector<std::vector<std::pair<int,int> > > const&
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
class GuoAlgorithm2D {
public:
    typedef Descriptor<T> D;
    GuoAlgorithm2D (
        OffLatticeModel2D<T,Array<T,2> >& model_,
        BlockLattice2D<T,Descriptor>& lattice_,
        Dot2D const& guoNode_,
        std::vector<std::pair<int,int> > const& dryNodeFluidDirections_,
        std::vector<plint> const& dryNodeIds_, Dot2D const& absoluteOffset_,
        Array<T,2>& localForce_, std::vector<AtomicBlock2D const*> const& args_,
        bool computeStat_, bool secondOrder_);
    virtual ~GuoAlgorithm2D() { }
    bool computeNeighborData();
    void finalize();

    virtual void extrapolateVariables (
              Dot2D const& fluidDirection, int depth, Array<T,2> const& wallNode, T delta,
              Array<T,2> const& wall_vel, OffBoundary::Type bdType,
              Array<T,2> const& wallNormal, plint triangleId, plint iDirection ) =0;
    virtual void reduceVariables(T sumWeights) =0;
    virtual void complete() =0;

protected:
    OffLatticeModel2D<T,Array<T,2> >& model;
    BlockLattice2D<T,Descriptor>& lattice;
    Dot2D const& guoNode;
    Cell<T,Descriptor>& cell;
    std::vector<std::pair<int,int> > const& dryNodeFluidDirections;
    std::vector<plint> const& dryNodeIds;
    Dot2D absoluteOffset;
    Array<T,2>& localForce;
    std::vector<AtomicBlock2D const*> const& args;

    plint numDirections;
    std::vector<T> weights;
    std::vector<T> rhoBarVect;
    std::vector<Array<T,Descriptor<T>::d> > jVect;

    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    bool computeStat, secondOrder;
};

template<typename T, template<typename U> class Descriptor>
GuoAlgorithm2D<T,Descriptor>::GuoAlgorithm2D (
            OffLatticeModel2D<T,Array<T,2> >& model_,
            BlockLattice2D<T,Descriptor>& lattice_,
            Dot2D const& guoNode_,
            std::vector<std::pair<int,int> > const& dryNodeFluidDirections_,
            std::vector<plint> const& dryNodeIds_, Dot2D const& absoluteOffset_,
            Array<T,2>& localForce_, std::vector<AtomicBlock2D const*> const& args_,
            bool computeStat_, bool secondOrder_ )
    : model(model_),
      lattice(lattice_),
      guoNode(guoNode_),
      cell(lattice.get(guoNode.x, guoNode.y)),
      dryNodeFluidDirections(dryNodeFluidDirections_),
      dryNodeIds(dryNodeIds_),
      absoluteOffset(absoluteOffset_),
      localForce(localForce_),
      args(args_),
      computeStat(computeStat_),
      secondOrder(secondOrder_)
{
    numDirections = (plint)dryNodeFluidDirections.size();
    weights.resize(numDirections);
    rhoBarVect.resize(numDirections);
    jVect.resize(numDirections);
}

template<typename T, template<typename U> class Descriptor>
bool GuoAlgorithm2D<T,Descriptor>::computeNeighborData()
{
    T sumWeights = T();
    Array<T,2> wallNormal;
    for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
        int iNeighbor = dryNodeFluidDirections[iDirection].first;
        int const* c = NextNeighbor<T>::c[iNeighbor];
        Dot2D fluidDirection(c[0],c[1]);
        plint dryNodeId = dryNodeIds[iDirection];
        int depth = dryNodeFluidDirections[iDirection].second;

        Array<T,2> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;
#ifdef PLB_DEBUG
        bool ok =
#endif
        this->model.pointOnSurface( guoNode+absoluteOffset, fluidDirection,
                                    wallNode, wallDistance, wallNormal,
                                    wall_vel, bdType, dryNodeId );
        PLB_ASSERT( ok );
        if (! ( bdType==OffBoundary::dirichlet || bdType==OffBoundary::neumann ||
                bdType==OffBoundary::freeSlip || bdType==OffBoundary::constRhoInlet || bdType==OffBoundary::densityNeumann) )
        {
            return false;
        }
        if (bdType==OffBoundary::dirichlet) {
            for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                // Use the formula uLB = uP - 1/2 g. If there is no external force,
                //   the force term automatically evaluates to zero.
                wall_vel[iD] -= (T)0.5*getExternalForceComponent(cell,iD);
            }
        }
        T invDistanceToNeighbor = NextNeighbor<T>::invD[iNeighbor];
        PLB_ASSERT( wallDistance <= NextNeighbor<T>::d[iNeighbor] );
        T delta = (T)1. - wallDistance * invDistanceToNeighbor;
        Array<T,2> normalFluidDirection((T)fluidDirection.x, (T)fluidDirection.y);
        normalFluidDirection *= invDistanceToNeighbor;
        weights[iDirection] = fabs(dot(normalFluidDirection, wallNormal));
        sumWeights += weights[iDirection];
        this->extrapolateVariables (
                fluidDirection, depth, wallNode, delta, wall_vel,
                bdType, wallNormal, dryNodeId, iDirection );
    }
    this->reduceVariables(sumWeights);

    return true;
}


template<typename T, template<typename U> class Descriptor>
void GuoAlgorithm2D<T,Descriptor>::finalize()
{
    Array<T,D::d> deltaJ;
    deltaJ.resetToZero();
    if (computeStat) {
        for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
            int iNeighbor = dryNodeFluidDirections[iDirection].first;
            int iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
            if (iPop>=0) {
                plint oppPop = indexTemplates::opposite<D>(iPop);
                deltaJ[0] += D::c[oppPop][0]*cell[oppPop];
                deltaJ[1] += D::c[oppPop][1]*cell[oppPop];
            }
        }
    }

    this->complete();

    if (computeStat) {
        Cell<T,Descriptor> collidedCell(cell);
        BlockStatistics statsCopy(lattice.getInternalStatistics());
        collidedCell.collide(statsCopy);

        for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
            int iNeighbor = dryNodeFluidDirections[iDirection].first;
            plint iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
            if (iPop>=0) {
                deltaJ[0] -= D::c[iPop][0]*collidedCell[iPop];
                deltaJ[1] -= D::c[iPop][1]*collidedCell[iPop];
            }
        }
        // Don't divide by rho. Here we just divide by rho0=1. Remember that,
        //   while j is conserved along a channel, rho is not due to the
        //   pressure drop. Dividing by rho instead of rho0 would yield a
        //   larger result upstream than downstream, independent of whether
        //   the velocity is u or j.
        localForce += deltaJ;
    }
}

template<typename T, template<typename U> class Descriptor>
class GuoPiNeqAlgorithm2D : public GuoAlgorithm2D<T,Descriptor>
{
public:
    typedef Descriptor<T> D;
    GuoPiNeqAlgorithm2D (
        OffLatticeModel2D<T,Array<T,2> >& model_,
        BlockLattice2D<T,Descriptor>& lattice_,
        Dot2D const& guoNode_,
        std::vector<std::pair<int,int> > const& dryNodeFluidDirections_,
        std::vector<plint> const& dryNodeIds_, Dot2D const& absoluteOffset_,
        Array<T,2>& localForce_, std::vector<AtomicBlock2D const*> const& args_,
        bool computeStat_, bool secondOrder_ );
    virtual void extrapolateVariables (
              Dot2D const& fluidDirection, int depth, Array<T,2> const& wallNode, T delta,
              Array<T,2> const& wall_vel, OffBoundary::Type bdType,
              Array<T,2> const& wallNormal, plint triangleId, plint iDirection );
    virtual void reduceVariables(T sumWeights);
    virtual void complete();
private:
    std::vector<Array<T,SymmetricTensor<T,Descriptor>::n> > PiNeqVect;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
};

template<typename T, template<typename U> class Descriptor>
GuoPiNeqAlgorithm2D<T,Descriptor>::GuoPiNeqAlgorithm2D (
            OffLatticeModel2D<T,Array<T,2> >& model_,
            BlockLattice2D<T,Descriptor>& lattice_,
            Dot2D const& guoNode_,
            std::vector<std::pair<int,int> > const& dryNodeFluidDirections_,
            std::vector<plint> const& dryNodeIds_, Dot2D const& absoluteOffset_,
            Array<T,2>& localForce_, std::vector<AtomicBlock2D const*> const& args_,
            bool computeStat_, bool secondOrder_ )
    : GuoAlgorithm2D<T,Descriptor> (
            model_, lattice_, guoNode_, dryNodeFluidDirections_, dryNodeIds_,
            absoluteOffset_, localForce_, args_, computeStat_, secondOrder_ )
{
    PiNeqVect.resize(this->numDirections);
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void GuoPiNeqAlgorithm2D<T,Descriptor>::extrapolateVariables (
          Dot2D const& fluidDirection, int depth, Array<T,2> const& wallNode, T delta,
          Array<T,2> const& wall_vel, OffBoundary::Type bdType,
          Array<T,2> const& wallNormal, plint triangleId, plint iDirection )
{
    if (!this->secondOrder) {
        depth=1;
    }
    T rhoBar1;
    Array<T,Descriptor<T>::d> j1, j2;
    Cell<T,Descriptor> const& cell1 =
        this->lattice.get( this->guoNode.x+fluidDirection.x,
                           this->guoNode.y+fluidDirection.y);
    cell1.getDynamics().computeRhoBarJPiNeq(cell1, rhoBar1, j1, this->PiNeqVect[iDirection]);
    if (this->args.empty()) {
        Cell<T,Descriptor> const& cell2 =
            this->lattice.get( this->guoNode.x+2*fluidDirection.x,
                               this->guoNode.y+2*fluidDirection.y);

        T tmpRhoBar;
        cell2.getDynamics().computeRhoBarJ(cell2, tmpRhoBar, j2);
    }
    else {
        PLB_ASSERT( (plint)this->args.size()==1 );
        NTensorField2D<T> const* macroField =
            dynamic_cast<NTensorField2D<T> const*>( this->args[0] );
        PLB_ASSERT( macroField );
        // 1 Variable for rhoBar, 3 variables for j.
        PLB_ASSERT( macroField->getNdim()==4 );
        Dot2D offset = computeRelativeDisplacement(this->lattice, *macroField);
        T const* macroscopic = macroField->get (
                this->guoNode.x+2*fluidDirection.x+offset.x,
                this->guoNode.y+2*fluidDirection.y+offset.y);
        j2.from_cArray(macroscopic+1);
    }

    if ( bdType==OffBoundary::constRhoInlet ||
         bdType==OffBoundary::densityNeumann )
    {
        this->rhoBarVect[iDirection]=Descriptor<T>::rhoBar(wall_vel[0]);
    }
    else {
        this->rhoBarVect[iDirection] = rhoBar1;
    }
    Array<T,2> wall_j (
            this->model.velIsJ() ? wall_vel :
                             (Descriptor<T>::fullRho(this->rhoBarVect[iDirection])*wall_vel) );
    if (depth < 2) {
        if (delta < (T)0.25) {
            this->jVect[iDirection] = wall_j;
        }
        else {
            //j = (T)1./delta * (wall_j+(delta-(T)1.)*j1);
            this->jVect[iDirection] = wall_j; // Temporary fix for complex geometries.
        }
    }
    else {  // depth >= 2           d=1   d=0.5   d=0
        if (delta < (T)0.75) {   // x---|=========o-------------o
            this->jVect[iDirection] = wall_j + (delta-(T)1.)*j1 +
                ((T)1.-delta)/((T)1.+delta)*((T)2.*wall_j+(delta-(T)1.)*j2);
        }  //       d=1   d=0.5   d=0
        else {   // x===|---------o-------------o
            this->jVect[iDirection] = (T)1./delta * (wall_j+(delta-(T)1.)*j1);
        }
    }
    if ( bdType==OffBoundary::neumann )
    {
        this->jVect[iDirection] = j1;
    }
    else if ( bdType==OffBoundary::densityNeumann ) {
        this->jVect[iDirection] = dot(j1,wallNormal)*wallNormal;
    }
    else if ( bdType==OffBoundary::freeSlip ) {
        Array<T,2> continuousNormal =
            this->model.computeContinuousNormal(wallNode, triangleId);
        this->jVect[iDirection] = j1-dot(j1,continuousNormal)*continuousNormal;
    }
}

template<typename T, template<typename U> class Descriptor>
void GuoPiNeqAlgorithm2D<T,Descriptor>::reduceVariables(T sumWeights) {
    this->rhoBar = T();
    this->j.resetToZero();
    PiNeq.resetToZero();
    for (plint iDirection=0; iDirection<this->numDirections; ++iDirection) {
        this->rhoBar += this->rhoBarVect[iDirection] * this->weights[iDirection];
        this->j += this->jVect[iDirection] * this->weights[iDirection];
        PiNeq += PiNeqVect[iDirection] * this->weights[iDirection];
    }
    this->rhoBar /= sumWeights;
    this->j /= sumWeights;
    PiNeq /= sumWeights;
}

template<typename T, template<typename U> class Descriptor>
void GuoPiNeqAlgorithm2D<T,Descriptor>::complete() {
    Dynamics<T,Descriptor> const& dynamics = this->cell.getDynamics();
    T jSqr = normSqr(this->j);
    if (this->model.getPartialReplace()) {
        Cell<T,Descriptor> saveCell(this->cell);
        dynamics.regularize(this->cell, this->rhoBar, this->j, jSqr, PiNeq);
        for (plint iDirection=0; iDirection<this->numDirections; ++iDirection) {
            int iNeighbor = this->dryNodeFluidDirections[iDirection].first;
            plint iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
            plint oppPop = indexTemplates::opposite<D>(iPop);
            this->cell[oppPop] = saveCell[oppPop];
        }
    }
    else {
        dynamics.regularize(this->cell, this->rhoBar, this->j, jSqr, PiNeq);
    }
}


template<typename T, template<typename U> class Descriptor>
class GuoOffPopAlgorithm2D : public GuoAlgorithm2D<T,Descriptor>
{
public:
    typedef Descriptor<T> D;
    GuoOffPopAlgorithm2D (
        OffLatticeModel2D<T,Array<T,2> >& model_,
        BlockLattice2D<T,Descriptor>& lattice_,
        Dot2D const& guoNode_,
        std::vector<std::pair<int,int> > const& dryNodeFluidDirections_,
        std::vector<plint> const& dryNodeIds_, Dot2D const& absoluteOffset_,
        Array<T,2>& localForce_, std::vector<AtomicBlock2D const*> const& args_, bool computeStat_, bool secondOrder_ );
    virtual void extrapolateVariables (
              Dot2D const& fluidDirection, int depth, Array<T,2> const& wallNode, T delta,
              Array<T,2> const& wall_vel, OffBoundary::Type bdType,
              Array<T,2> const& wallNormal, plint triangleId, plint iDirection );
    virtual void reduceVariables(T sumWeights);
    virtual void complete();
private:
    std::vector<Array<T,Descriptor<T>::q> > fNeqVect;
    Array<T,Descriptor<T>::q> fNeq;
};

template<typename T, template<typename U> class Descriptor>
GuoOffPopAlgorithm2D<T,Descriptor>::GuoOffPopAlgorithm2D (
            OffLatticeModel2D<T,Array<T,2> >& model_,
            BlockLattice2D<T,Descriptor>& lattice_,
            Dot2D const& guoNode_,
            std::vector<std::pair<int,int> > const& dryNodeFluidDirections_,
            std::vector<plint> const& dryNodeIds_, Dot2D const& absoluteOffset_,
            Array<T,2>& localForce_, std::vector<AtomicBlock2D const*> const& args_,
            bool computeStat_, bool secondOrder_ )
    : GuoAlgorithm2D<T,Descriptor> (
            model_, lattice_, guoNode_, dryNodeFluidDirections_, dryNodeIds_,
            absoluteOffset_, localForce_, args_, computeStat_, secondOrder_ )
{
    fNeqVect.resize(this->numDirections);
    fNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void GuoOffPopAlgorithm2D<T,Descriptor>::extrapolateVariables (
          Dot2D const& fluidDirection, int depth, Array<T,2> const& wallNode, T delta,
          Array<T,2> const& wall_vel, OffBoundary::Type bdType,
          Array<T,2> const& wallNormal, plint triangleId, plint iDirection )
{
    if (!this->secondOrder) {
        depth=1;
    }
    T rhoBar1, rhoBar2;
    Array<T,Descriptor<T>::d> j1, j2;
    Array<T,Descriptor<T>::q> fNeq1, fNeq2;
    Cell<T,Descriptor> const& cell1 =
        this->lattice.get( this->guoNode.x+fluidDirection.x,
                           this->guoNode.y+fluidDirection.y);
    Cell<T,Descriptor> const& cell2 =
        this->lattice.get( this->guoNode.x+2*fluidDirection.x,
                           this->guoNode.y+2*fluidDirection.y);
    cell1.getDynamics().computeRhoBarJ(cell1, rhoBar1, j1);
    cell2.getDynamics().computeRhoBarJ(cell2, rhoBar2, j2);
    T j1sqr = normSqr(j1);
    T j2sqr = normSqr(j2);
    for(plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        fNeq1[iPop] = cell1[iPop]-cell1.getDynamics().computeEquilibrium(iPop, rhoBar1, j1, j1sqr);
        fNeq2[iPop] = cell2[iPop]-cell2.getDynamics().computeEquilibrium(iPop, rhoBar2, j2, j2sqr);
    }

    if ( bdType==OffBoundary::constRhoInlet ||
         bdType==OffBoundary::densityNeumann )
    {
        this->rhoBarVect[iDirection]=Descriptor<T>::rhoBar(wall_vel[0]);
    }
    else {
        this->rhoBarVect[iDirection] = rhoBar1;
    }
    Array<T,2> wall_j (
            this->model.velIsJ() ? wall_vel :
                             (Descriptor<T>::fullRho(this->rhoBarVect[iDirection])*wall_vel) );
    Array<T,2> jNeumann;
    this->fNeqVect[iDirection] = fNeq1;
    jNeumann = j1;
    if (depth < 2) {
        if (delta < (T)0.25) {
            this->jVect[iDirection] = wall_j;
        }
        else {
            this->jVect[iDirection] = (T)1./delta * (wall_j+(delta-(T)1.)*j1);
        }
    }
    else {  // depth >= 2
        if (delta < (T)0.75) {
            this->jVect[iDirection] = wall_j + (delta-(T)1.)*j1 +
                ((T)1.-delta)/((T)1.+delta)*((T)2.*wall_j+(delta-(T)1.)*j2);
        }
        else {
            this->jVect[iDirection] = (T)1./delta * (wall_j+(delta-(T)1.)*j1);
        }
    }
    if ( bdType==OffBoundary::neumann )
    {
        this->jVect[iDirection] = jNeumann;
    }
    else if ( bdType==OffBoundary::densityNeumann ) {
        this->jVect[iDirection] = dot(jNeumann,wallNormal)*wallNormal;
    }
    else if ( bdType==OffBoundary::freeSlip ) {
        Array<T,2> continuousNormal =
            this->model.computeContinuousNormal(wallNode, triangleId);
        this->jVect[iDirection] = jNeumann-dot(jNeumann,continuousNormal)*continuousNormal;
    }
}


template<typename T, template<typename U> class Descriptor>
void GuoOffPopAlgorithm2D<T,Descriptor>::reduceVariables(T sumWeights) {
    this->rhoBar = T();
    this->j.resetToZero();
    fNeq.resetToZero();
    for (plint iDirection=0; iDirection<this->numDirections; ++iDirection) {
        this->rhoBar += this->rhoBarVect[iDirection] * this->weights[iDirection];
        this->j += this->jVect[iDirection] * this->weights[iDirection];
        fNeq += fNeqVect[iDirection] * this->weights[iDirection];
    }
    this->rhoBar /= sumWeights;
    this->j /= sumWeights;
    fNeq /= sumWeights;
}

template<typename T, template<typename U> class Descriptor>
void GuoOffPopAlgorithm2D<T,Descriptor>::complete() {
    T jSqr = normSqr(this->j);
    if (this->model.getPartialReplace()) {
        for (plint iDirection=0; iDirection<this->numDirections; ++iDirection) {
            int iNeighbor = this->dryNodeFluidDirections[iDirection].first;
            plint iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
            this->cell[iPop] = this->cell.computeEquilibrium(iPop, this->rhoBar, this->j, jSqr)+fNeq[iPop];
        }
    }
    else {
        for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
            this->cell[iPop] = this->cell.computeEquilibrium(iPop, this->rhoBar, this->j, jSqr)+fNeq[iPop];
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void GuoOffLatticeModel2D<T,Descriptor>::cellCompletion (
        BlockLattice2D<T,Descriptor>& lattice,
        Dot2D const& guoNode,
        std::vector<std::pair<int,int> > const& dryNodeFluidDirections,
        std::vector<plint> const& dryNodeIds, Dot2D const& absoluteOffset,
        Array<T,2>& localForce, std::vector<AtomicBlock2D const*> const& args )
{
    GuoAlgorithm2D<T,Descriptor>* algorithm=0;
    if (this->regularizedModel) {
        algorithm = new GuoPiNeqAlgorithm2D<T,Descriptor> (
                *this, lattice, guoNode, dryNodeFluidDirections,
                dryNodeIds, absoluteOffset, localForce, args, computesStat(), usesSecondOrder() );
    }
    else {
        algorithm = new GuoOffPopAlgorithm2D<T,Descriptor> (
                *this, lattice, guoNode, dryNodeFluidDirections,
                dryNodeIds, absoluteOffset, localForce, args, computesStat(), usesSecondOrder() );
    }
    bool ok = algorithm -> computeNeighborData();
    PLB_ASSERT( ok );
    algorithm->finalize();
    delete algorithm;
}

}  // namespace plb

#endif  // GUO_OFF_LATTICE_MODEL_2D_HH

