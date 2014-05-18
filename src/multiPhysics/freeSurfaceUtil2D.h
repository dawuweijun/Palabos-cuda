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

#ifndef FREE_SURFACE_UTIL_2D_H
#define FREE_SURFACE_UTIL_2D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiPhysics/freeSurfaceUtil.h"
#include <vector>
#include <set>
#include <string>

namespace plb
{


/// Create a parameter-list for most free-surface data processors.
template< typename T,template<typename U> class Descriptor>
std::vector<MultiBlock2D*> aggregateFreeSurfaceParams2D (
    MultiBlockLattice2D<T,Descriptor>& fluid, MultiScalarField2D<T>& rhoBar,
    MultiTensorField2D<T,2>& j, MultiScalarField2D<T>& mass,
    MultiScalarField2D<T>& volumeFraction, MultiScalarField2D<int>& flag,
    MultiTensorField2D<T,2>& normal,
    MultiContainerBlock2D& interfaceLists, MultiScalarField2D<T>& curvature,
    MultiScalarField2D<T>& outsideDensity )
{
    std::vector<MultiBlock2D*> aggregation;

    aggregation.push_back ( &fluid );
    aggregation.push_back ( &rhoBar );
    aggregation.push_back ( &j );
    aggregation.push_back ( &mass );
    aggregation.push_back ( &volumeFraction );
    aggregation.push_back ( &flag );
    aggregation.push_back ( &normal );
    aggregation.push_back ( &interfaceLists );
    aggregation.push_back ( &curvature );
    aggregation.push_back ( &outsideDensity );

    return aggregation;
}

/// A wrapper offering convenient access to the free-surface data provided to
/// data processors. Avoids verbous casting, asserting, etc.
template<typename T,template<typename U> class Descriptor>
class FreeSurfaceProcessorParam2D
{
public:
    typedef typename InterfaceLists<T,Descriptor>::Node Node;
    FreeSurfaceProcessorParam2D ( std::vector<AtomicBlock2D*>& atomicBlocks )
    {
        PLB_ASSERT ( atomicBlocks.size() >= 10 );

        fluid_ = dynamic_cast<BlockLattice2D<T,Descriptor>*> ( atomicBlocks[0] );
        PLB_ASSERT ( fluid_ );

        rhoBar_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[1] );
        PLB_ASSERT ( rhoBar_ );

        j_ = dynamic_cast<TensorField2D<T,2>*> ( atomicBlocks[2] );
        PLB_ASSERT ( j_ );

        mass_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[3] );
        PLB_ASSERT ( mass_ );

        volumeFraction_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[4] );
        PLB_ASSERT ( volumeFraction_ );

        flag_ = dynamic_cast<ScalarField2D<int>*> ( atomicBlocks[5] );
        PLB_ASSERT ( flag_ );

        normal_ = dynamic_cast<TensorField2D<T,2>*> ( atomicBlocks[6] );
        PLB_ASSERT ( normal_ );

        containerInterfaceLists_ = dynamic_cast<AtomicContainerBlock2D*> ( atomicBlocks[7] );
        PLB_ASSERT ( containerInterfaceLists_ );

        interfaceLists_ = dynamic_cast<InterfaceLists<T,Descriptor>*> ( containerInterfaceLists_->getData() );
        //PLB_ASSERT(interfaceLists_);
        //Put the assertion at the usage of interfaceLists, so we can still work with both freeSurfaceProcessorParam and twoPhaseProcessorParam.


        curvature_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[8] );
        PLB_ASSERT ( curvature_ );

        outsideDensity_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[9] );
        PLB_ASSERT ( outsideDensity_ );

        absoluteOffset       = fluid_->getLocation();
        relativeOffsetRhoBar = computeRelativeDisplacement ( *fluid_,*rhoBar_ );
        relativeOffsetJ      = computeRelativeDisplacement ( *fluid_,*j_ );
        relativeOffsetMass   = computeRelativeDisplacement ( *fluid_,*mass_ );
        relativeOffsetVF     = computeRelativeDisplacement ( *fluid_,*volumeFraction_ );
        relativeOffsetFS     = computeRelativeDisplacement ( *fluid_,*flag_ );
        relativeOffsetNormal = computeRelativeDisplacement ( *fluid_,*normal_ );
        relativeOffsetC      = computeRelativeDisplacement ( *fluid_,*curvature_ );
        relativeOffsetOD     = computeRelativeDisplacement ( *fluid_,*outsideDensity_ );
    }
    Cell<T,Descriptor>& cell ( plint iX, plint iY )
    {
        return fluid_->get ( iX,iY );
    }
    T& mass ( plint iX, plint iY )
    {
        return mass_->get ( iX+relativeOffsetMass.x,iY+relativeOffsetMass.y );
    }
    T& volumeFraction ( plint iX, plint iY )
    {
        return volumeFraction_->get ( iX+relativeOffsetVF.x,iY+relativeOffsetVF.y );
    }
    T& curvature ( plint iX, plint iY )
    {
        return curvature_->get ( iX+relativeOffsetC.x,iY+relativeOffsetC.y );
    }
    T& outsideDensity ( plint iX, plint iY )
    {
        return outsideDensity_->get ( iX+relativeOffsetOD.x,iY+relativeOffsetOD.y );
    }
    int& flag ( plint iX, plint iY )
    {
        return flag_->get ( iX+relativeOffsetFS.x,iY+relativeOffsetFS.y );
    }
    void setForce ( plint iX, plint iY, Array<T,2> const& force )
    {
        force.to_cArray ( cell ( iX,iY ).getExternal ( forceOffset ) );
    }
    Array<T,2> getForce ( plint iX, plint iY )
    {
        Array<T,2> force;
        force.from_cArray ( cell ( iX,iY ).getExternal ( forceOffset ) );
        return force;
    }
    void setMomentum ( plint iX, plint iY, Array<T,2> const& momentum )
    {
        j_->get ( iX+relativeOffsetJ.x,iY+relativeOffsetJ.y ) = momentum;
    }
    Array<T,2> getMomentum ( plint iX, plint iY )
    {
        return j_->get ( iX+relativeOffsetJ.x,iY+relativeOffsetJ.y );
    }
    T getDensity ( plint iX, plint iY )
    {
        return Descriptor<T>::fullRho (
                   rhoBar_->get ( iX+relativeOffsetRhoBar.x, iY+relativeOffsetRhoBar.y ) );
    }
    void setDensity ( plint iX, plint iY, T rho )
    {
        rhoBar_->get ( iX+relativeOffsetRhoBar.x, iY+relativeOffsetRhoBar.y )
            = Descriptor<T>::rhoBar ( rho );
    }
    void setNormal ( plint iX, plint iY, Array<T,2> const& normal )
    {
        normal_->get ( iX+relativeOffsetNormal.x,iY+relativeOffsetNormal.y ) = normal;
    }
    Array<T,2> getNormal ( plint iX, plint iY )
    {
        return normal_->get ( iX+relativeOffsetNormal.x,iY+relativeOffsetNormal.y );
    }

    void attributeDynamics ( plint iX, plint iY, Dynamics<T,Descriptor>* dynamics )
    {
        fluid_->attributeDynamics ( iX,iY, dynamics );
    }

    bool isBoundary ( plint iX, plint iY )
    {
        return cell ( iX, iY ).getDynamics().isBoundary();
    }

    void addToTotalMass ( T addedTotalMass )
    {
        fluid_->getInternalStatistics().gatherSum ( 0, addedTotalMass );
    }
    void addToLostMass ( T addedLostMass )
    {
        fluid_->getInternalStatistics().gatherSum ( 1, addedLostMass );
    }
    void addToInterfaceCells ( plint addedInterfaceCells )
    {
        fluid_->getInternalStatistics().gatherIntSum ( 0, addedInterfaceCells );
    }
    T getSumMassMatrix() const
    {
        return fluid_->getInternalStatistics().getSum ( 0 );
    }
    T getSumLostMass() const
    {
        return fluid_->getInternalStatistics().getSum ( 1 );
    }
    T getTotalMass() const
    {
        return getSumMassMatrix() + getSumLostMass();
    }
    plint getNumInterfaceCells() const
    {
        return fluid_->getInternalStatistics().getIntSum ( 0 );
    }

    T smoothVolumeFraction ( plint iX, plint iY )
    {
        using namespace twoPhaseFlag;

        if ( flag_->get ( iX+relativeOffsetFS.x,iY+relativeOffsetFS.y ) == wall )
        {
            return volumeFraction_->get ( iX+relativeOffsetVF.x,iY+relativeOffsetVF.y );
        }

        T val = 0.0;
        int n = 0;
        for ( int i = -1; i < 2; i++ )
        {
            plint nextX = iX + i;
            for ( int j = -1; j < 2; j++ )
            {
                plint nextY = iY + j;
                if ( ! ( i == 0 && j == 0 )
                        &&flag_->get ( nextX+relativeOffsetFS.x,nextY+relativeOffsetFS.y ) != wall )
                {
                    n++;
                    val += volumeFraction_->get ( nextX+relativeOffsetVF.x,nextY+relativeOffsetVF.y );
                }
            }
        }
        if ( n != 0 )
        {
            val /= ( T ) n;
        }
        else
        {
            val = volumeFraction_->get ( iX+relativeOffsetVF.x,iY+relativeOffsetVF.y );
        }

        return val;
    }

    std::map<Node,T>& massExcess()
    {
        PLB_ASSERT ( interfaceLists_ );
        return interfaceLists_ -> massExcess;
    }
    std::set<Node>& interfaceToFluid()
    {
        PLB_ASSERT ( interfaceLists_ );
        return interfaceLists_ -> interfaceToFluid;
    }
    std::set<Node>& interfaceToEmpty()
    {
        PLB_ASSERT ( interfaceLists_ );
        return interfaceLists_ -> interfaceToEmpty;
    }
    std::set<Node>& emptyToInterface()
    {
        PLB_ASSERT ( interfaceLists_ );
        return interfaceLists_ -> emptyToInterface;
    }

    Dot2D const& absOffset() const
    {
        return absoluteOffset;
    }
    Box2D getBoundingBox() const
    {
        return volumeFraction_->getBoundingBox();
    }
private:
    BlockLattice2D<T,Descriptor>* fluid_;
    ScalarField2D<T>* rhoBar_;
    TensorField2D<T,2>* j_;
    ScalarField2D<T>* mass_;
    ScalarField2D<T>* volumeFraction_;
    ScalarField2D<int>* flag_;
    TensorField2D<T,2>* normal_;
    AtomicContainerBlock2D* containerInterfaceLists_;
    InterfaceLists<T,Descriptor>* interfaceLists_;
    ScalarField2D<T>* curvature_;
    ScalarField2D<T>* outsideDensity_;

    Dot2D absoluteOffset, relativeOffsetRhoBar, relativeOffsetJ, relativeOffsetMass,
          relativeOffsetVF, relativeOffsetFS, relativeOffsetNormal, relativeOffsetC,
          relativeOffsetOD;

    static const int forceOffset = Descriptor<T>::ExternalField::forceBeginsAt;
};

}  // namespace plb

#endif  // FREE_SURFACE_UTIL_2D_H

