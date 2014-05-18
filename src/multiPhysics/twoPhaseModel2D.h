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

#ifndef TWO_PHASE_MODEL_2D_H
#define TWO_PHASE_MODEL_2D_H

#include "core/globalDefs.h"
#include "multiPhysics/freeSurfaceModel2D.h"
#include "multiPhysics/twoPhaseModel.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include <memory>

namespace plb
{

template<typename T, int nDim>
class IniFilteredDensity2D :
    public BoxProcessingFunctional2D_ST<T,T,nDim>
{
public:
    virtual void process ( Box2D domain, ScalarField2D<T>& scalarField,
                           TensorField2D<T,nDim>& tensorField );
    virtual IniFilteredDensity2D<T,nDim>* clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const;
};

template<typename T, int nDim>
void IniFilteredDensity2D<T,nDim>::process (
    Box2D domain, ScalarField2D<T>& scalarField,
    TensorField2D<T,nDim>& tensorField )
{
    Dot2D offset = computeRelativeDisplacement ( scalarField, tensorField );
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            for ( int iDim=0; iDim<nDim; ++iDim )
            {
                tensorField.get ( iX+offset.x,iY+offset.y ) [iDim] =
                    scalarField.get ( iX,iY );
            }
        }
    }
}

template<typename T, int nDim>
IniFilteredDensity2D<T,nDim>* IniFilteredDensity2D<T,nDim>::clone() const
{
    return new IniFilteredDensity2D<T,nDim> ( *this );
}

template<typename T, int nDim>
void IniFilteredDensity2D<T,nDim>::getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
void iniFilteredDensity ( MultiTensorField2D<T,nDim>& densities, MultiScalarField2D<T>& density, Box2D domain )
{
    applyProcessingFunctional (
        new IniFilteredDensity2D<T,nDim>(), domain, density, densities );
}


template<typename T, int nDim>
class AddFilteredDensity2D :
    public BoxProcessingFunctional2D_ST<T,T,nDim>
{
public:
    virtual void process ( Box2D domain, ScalarField2D<T>& scalarField,
                           TensorField2D<T,nDim>& tensorField );
    virtual AddFilteredDensity2D<T,nDim>* clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const;
};

template<typename T, int nDim>
void AddFilteredDensity2D<T,nDim>::process (
    Box2D domain, ScalarField2D<T>& scalarField,
    TensorField2D<T,nDim>& tensorField )
{
    Dot2D offset = computeRelativeDisplacement ( scalarField, tensorField );
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            for ( int iDim=0; iDim<nDim-1; ++iDim )
            {
                tensorField.get ( iX+offset.x,iY+offset.y ) [iDim] =
                    tensorField.get ( iX+offset.x,iY+offset.y ) [iDim+1];
            }
            tensorField.get ( iX+offset.x,iY+offset.y ) [nDim-1] =
                scalarField.get ( iX,iY );
        }
    }
}

template<typename T, int nDim>
AddFilteredDensity2D<T,nDim>* AddFilteredDensity2D<T,nDim>::clone() const
{
    return new AddFilteredDensity2D<T,nDim> ( *this );
}

template<typename T, int nDim>
void AddFilteredDensity2D<T,nDim>::getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
void addFilteredDensity ( MultiTensorField2D<T,nDim>& densities, MultiScalarField2D<T>& density, Box2D domain )
{
    applyProcessingFunctional (
        new AddFilteredDensity2D<T,nDim>(), domain, density, densities );
}


template<typename T, int nDim>
class GetFilteredDensity2D :
    public BoxProcessingFunctional2D_ST<T,T,nDim>
{
public:
    virtual void process ( Box2D domain, ScalarField2D<T>& scalarField,
                           TensorField2D<T,nDim>& tensorField );
    virtual GetFilteredDensity2D<T,nDim>* clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const;
};

template<typename T, int nDim>
void GetFilteredDensity2D<T,nDim>::process (
    Box2D domain, ScalarField2D<T>& scalarField,
    TensorField2D<T,nDim>& tensorField )
{
    Dot2D offset = computeRelativeDisplacement ( scalarField, tensorField );
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            T averageDensity=T();
            for ( int iDim=0; iDim<nDim; ++iDim )
            {
                averageDensity +=
                    tensorField.get ( iX+offset.x,iY+offset.y ) [iDim];
            }
            scalarField.get ( iX,iY ) = averageDensity;
        }
    }
}

template<typename T, int nDim>
GetFilteredDensity2D<T,nDim>* GetFilteredDensity2D<T,nDim>::clone() const
{
    return new GetFilteredDensity2D<T,nDim> ( *this );
}

template<typename T, int nDim>
void GetFilteredDensity2D<T,nDim>::getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, int nDim>
void getFilteredDensity ( MultiTensorField2D<T,nDim>& densities, MultiScalarField2D<T>& density, Box2D domain )
{
    applyProcessingFunctional (
        new GetFilteredDensity2D<T,nDim>(), domain, density, densities );
}


template<typename T, template<typename U> class Descriptor>
class TwoPhaseAveragePressure2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    TwoPhaseAveragePressure2D ( T densityRatio_, T rhoDefault_, TwoPhaseModel model_ );
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseAveragePressure2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseAveragePressure2D<T,Descriptor> ( *this );
    }
    T getAveragePressure() const;
    T getAveragePressure2() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::nothing;         // Fluid 2.
            modified[11] = modif::nothing;         // rhoBar2.
            modified[12] = modif::nothing;         // j2.
            modified[13] = modif::staticVariables; // mass2.
        }
    }
private:
    plint sumRho1_ID, sumRho2_ID, weight1_ID, weight2_ID;
    T densityRatio, rhoDefault;
    TwoPhaseModel model;
};

template<typename T, template<typename U> class Descriptor>
class TwoPhaseAverageVelocity2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    TwoPhaseAverageVelocity2D ( TwoPhaseModel model_ );
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseAverageVelocity2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseAverageVelocity2D<T,Descriptor> ( *this );
    }
    Array<T,2> getAverageVelocity() const;
    Array<T,2> getAverageVelocity2() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::nothing;         // Fluid 2.
            modified[11] = modif::nothing;         // rhoBar2.
            modified[12] = modif::nothing;         // j2.
            modified[13] = modif::staticVariables; // mass2.
        }
    }
private:
    Array<plint,2> sumVel1_ID, sumVel2_ID;
    plint weight1_ID, weight2_ID;
    TwoPhaseModel model;
};

template<typename T,template<typename U> class Descriptor>
class TwoPhaseComputePressure2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseComputePressure2D ( T densityRatio_, T rhoDefault_, TwoPhaseModel model_ )
        : densityRatio ( densityRatio_ ),
          rhoDefault ( rhoDefault_ ),
          model ( model_ ),
          computeFluid1 ( true ),
          computeFluid2 ( true )
    { }
    TwoPhaseComputePressure2D ( T densityRatio_, T rhoDefault_, TwoPhaseModel model_, bool computeFluid1_, bool computeFluid2_ )
        : densityRatio ( densityRatio_ ),
          rhoDefault ( rhoDefault_ ),
          model ( model_ ),
          computeFluid1 ( computeFluid1_ ),
          computeFluid2 ( computeFluid2_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseComputePressure2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseComputePressure2D ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::nothing;         // Fluid 2.
            modified[11] = modif::nothing;         // rhoBar2.
            modified[12] = modif::nothing;         // j2.
            modified[13] = modif::nothing;         // mass2.
            modified[14] = modif::staticVariables; // resulting density
        }
        else
        {
            modified[10] = modif::staticVariables; // resulting density
        }
    }
private:
    T densityRatio, rhoDefault;
    TwoPhaseModel model;
    bool computeFluid1;
    bool computeFluid2;
};

template<typename T,template<typename U> class Descriptor>
class TwoPhaseComputeVelocity2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseComputeVelocity2D ( T densityRatio_, bool useFreeSurfaceLimit_ )
        : densityRatio ( densityRatio_ ),
          computeFluid1 ( true ),
          computeFluid2 ( true ),
          useFreeSurfaceLimit ( useFreeSurfaceLimit_ )
    { }
    TwoPhaseComputeVelocity2D ( T densityRatio_, bool computeFluid1_, bool computeFluid2_, bool useFreeSurfaceLimit_ )
        : densityRatio ( densityRatio_ ),
          computeFluid1 ( computeFluid1_ ),
          computeFluid2 ( computeFluid2_ ),
          useFreeSurfaceLimit ( useFreeSurfaceLimit_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseComputeVelocity2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseComputeVelocity2D ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if ( !useFreeSurfaceLimit )
        {
            modified[10] = modif::nothing;         // Fluid 2.
            modified[11] = modif::nothing;         // rhoBar2.
            modified[12] = modif::nothing;         // j2.
            modified[13] = modif::nothing;         // mass2.
            modified[14] = modif::staticVariables; // resulting velocity
        }
        else
        {
            modified[10] = modif::staticVariables; // resulting velocity
        }
    }
private:
    T densityRatio;
    bool computeFluid1;
    bool computeFluid2;
    bool useFreeSurfaceLimit;
};


template<typename T, template<typename U> class Descriptor>
class TwoPhasePunchSphere2D : public BoxProcessingFunctional2D
{
public:
    TwoPhasePunchSphere2D (
        Array<T,2> const& center_, T radius_, T rhoEmpty_, T referenceDensity_,
        T densityRatio_, Dynamics<T,Descriptor> *dynamicsTemplate2_, TwoPhaseModel model_ )
        : center ( center_ ),
          radius ( radius_ ),
          rhoEmpty ( rhoEmpty_ ),
          referenceDensity ( referenceDensity_ ),
          densityRatio ( densityRatio_ ),
          dynamicsTemplate2 ( dynamicsTemplate2_ ),
          model ( model_ )
    { }
    TwoPhasePunchSphere2D ( TwoPhasePunchSphere2D<T,Descriptor> const& rhs )
        : center ( rhs.center ),
          radius ( rhs.radius ),
          rhoEmpty ( rhs.rhoEmpty ),
          referenceDensity ( rhs.referenceDensity ),
          densityRatio ( rhs.densityRatio ),
          dynamicsTemplate2 ( rhs.dynamicsTemplate2->clone() ),
          model ( rhs.model )
    { }
    TwoPhasePunchSphere2D<T,Descriptor>& operator= ( TwoPhasePunchSphere2D<T,Descriptor> const& rhs )
    {
        TwoPhasePunchSphere2D<T,Descriptor> ( rhs ).swap ( *this );
        return *this;
    }
    void swap ( TwoPhasePunchSphere2D<T,Descriptor>& rhs )
    {
        std::swap ( center, rhs.center );
        std::swap ( radius, rhs.radius );
        std::swap ( rhoEmpty, rhs.rhoEmpty );
        std::swap ( referenceDensity, rhs.referenceDensity );
        std::swap ( densityRatio, rhs.densityRatio );
        std::swap ( dynamicsTemplate2, rhs.dynamicsTemplate2 );
        std::swap ( model, rhs.model );
    }
    ~TwoPhasePunchSphere2D()
    {
        delete dynamicsTemplate2;
    }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhasePunchSphere2D<T,Descriptor>* clone() const
    {
        return new TwoPhasePunchSphere2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::staticVariables; // Fluid.
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::staticVariables; // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::staticVariables; // Volume-fraction.
        modified[5] = modif::staticVariables; // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::staticVariables; // Curvature.
        modified[9] = modif::staticVariables; // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::dataStructure;  // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }
private:
    Array<T,2> center;
    T radius;
    T rhoEmpty, referenceDensity, densityRatio;
    Dynamics<T,Descriptor> *dynamicsTemplate2;
    TwoPhaseModel model;
};

template<typename T, template<typename U> class Descriptor>
class TwoPhasePunchRectangle2D : public BoxProcessingFunctional2D
{
public:
    TwoPhasePunchRectangle2D (
        Box2D rectangle_, T rhoEmpty_, T referenceDensity_,
        T densityRatio_, Dynamics<T,Descriptor> *dynamicsTemplate2_, TwoPhaseModel model_ )
        : rectangle ( rectangle_ ),
          rhoEmpty ( rhoEmpty_ ),
          referenceDensity ( referenceDensity_ ),
          densityRatio ( densityRatio_ ),
          dynamicsTemplate2 ( dynamicsTemplate2_ ),
          model ( model_ )
    { }
    TwoPhasePunchRectangle2D ( TwoPhasePunchRectangle2D<T,Descriptor> const& rhs )
        : rectangle ( rhs.rectangle ),
          rhoEmpty ( rhs.rhoEmpty ),
          referenceDensity ( rhs.referenceDensity ),
          densityRatio ( rhs.densityRatio ),
          dynamicsTemplate2 ( rhs.dynamicsTemplate2->clone() ),
          model ( rhs.model )
    { }
    TwoPhasePunchRectangle2D<T,Descriptor>& operator= ( TwoPhasePunchRectangle2D<T,Descriptor> const& rhs )
    {
        TwoPhasePunchRectangle2D<T,Descriptor> ( rhs ).swap ( *this );
        return *this;
    }
    void swap ( TwoPhasePunchRectangle2D<T,Descriptor>& rhs )
    {
        std::swap ( rectangle, rhs.rectangle );
        std::swap ( rhoEmpty, rhs.rhoEmpty );
        std::swap ( referenceDensity, rhs.referenceDensity );
        std::swap ( densityRatio, rhs.densityRatio );
        std::swap ( dynamicsTemplate2, rhs.dynamicsTemplate2 );
        std::swap ( model, rhs.model );
    }
    ~TwoPhasePunchRectangle2D()
    {
        delete dynamicsTemplate2;
    }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhasePunchRectangle2D<T,Descriptor>* clone() const
    {
        return new TwoPhasePunchRectangle2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::staticVariables; // Fluid.
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::staticVariables; // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::staticVariables; // Volume-fraction.
        modified[5] = modif::staticVariables; // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::staticVariables; // Curvature.
        modified[9] = modif::staticVariables; // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::dataStructure;  // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }
private:
    Box2D rectangle;
    T rhoEmpty, referenceDensity, densityRatio;
    Dynamics<T,Descriptor> *dynamicsTemplate2;
    TwoPhaseModel model;
};

template<typename T, template<typename U> class Descriptor>
void twoPhasePunchSphere ( std::vector<MultiBlock2D*> const& twoPhaseArgs, Array<T,2> const& center, T radius,
                           T rhoEmpty, T referenceDensity, T densityRatio, Dynamics<T,Descriptor>& dynamics,
                           TwoPhaseModel model, Box2D domain );

template<typename T, template<typename U> class Descriptor>
void twoPhasePunchRectangle ( std::vector<MultiBlock2D*> const& twoPhaseArgs, Box2D rectangle,
                              T rhoEmpty, T referenceDensity, T densityRatio, Dynamics<T,Descriptor>& dynamics,
                              TwoPhaseModel model, Box2D domain );

template<typename T, template<typename U> class Descriptor>
T computeAverageSphereDensity ( std::vector<MultiBlock2D*> const& twoPhaseArgs,
                                Array<T,2> const& center, T radius, Box2D domain );


/// A wrapper offering convenient access to the free-surface data provided to
/// data processors. Avoids verbous casting, asserting, etc.
template<typename T,template<typename U> class Descriptor>
class TwoPhaseProcessorParam2D
{
public:
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::ExtrapolInfo ExtrapolInfo;
    TwoPhaseProcessorParam2D ( std::vector<AtomicBlock2D*>& atomicBlocks )
    {
        if ( atomicBlocks.size() >=14 )
        {
            useFreeSurfaceLimit = false;
        }
        else
        {
            useFreeSurfaceLimit = true;
        }

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

        interfaceLists_ = dynamic_cast<TwoPhaseInterfaceLists<T,Descriptor>*> ( containerInterfaceLists_->getData() );
        //PLB_ASSERT(interfaceLists_);
        //Put the assertion at the usage of interfaceLists, so we can still work with both freeSurfaceProcessorParam and twoPhaseProcessorParam.

        curvature_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[8] );
        PLB_ASSERT ( curvature_ );

        outsideDensity_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[9] );
        PLB_ASSERT ( outsideDensity_ );

        if ( !useFreeSurfaceLimit )
        {
            fluid2_ = dynamic_cast<BlockLattice2D<T,Descriptor>*> ( atomicBlocks[10] );
            PLB_ASSERT ( fluid2_ );

            rhoBar2_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[11] );
            PLB_ASSERT ( rhoBar2_ );

            j2_ = dynamic_cast<TensorField2D<T,2>*> ( atomicBlocks[12] );
            PLB_ASSERT ( j2_ );

            mass2_ = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[13] );
            PLB_ASSERT ( mass2_ );
        }
        else
        {
            fluid2_ = 0;
            rhoBar2_ = 0;
            j2_ = 0;
            mass2_ = 0;
        }

        absoluteOffset       = fluid_->getLocation();
        relativeOffsetRhoBar = computeRelativeDisplacement ( *fluid_,*rhoBar_ );
        relativeOffsetJ      = computeRelativeDisplacement ( *fluid_,*j_ );
        relativeOffsetMass   = computeRelativeDisplacement ( *fluid_,*mass_ );
        relativeOffsetVF     = computeRelativeDisplacement ( *fluid_,*volumeFraction_ );
        relativeOffsetFS     = computeRelativeDisplacement ( *fluid_,*flag_ );
        relativeOffsetNormal = computeRelativeDisplacement ( *fluid_,*normal_ );
        relativeOffsetC      = computeRelativeDisplacement ( *fluid_,*curvature_ );
        relativeOffsetOD     = computeRelativeDisplacement ( *fluid_,*outsideDensity_ );

        if ( !useFreeSurfaceLimit )
        {
            relativeOffsetFluid2 = computeRelativeDisplacement ( *fluid_,*fluid2_ );
            relativeOffsetRhoBar2= computeRelativeDisplacement ( *fluid_,*rhoBar2_ );
            relativeOffsetJ2     = computeRelativeDisplacement ( *fluid_,*j2_ );
            relativeOffsetMass2  = computeRelativeDisplacement ( *fluid_,*mass2_ );
        }
    }
    Cell<T,Descriptor>& cell ( plint iX, plint iY )
    {
        return fluid_->get ( iX,iY );
    }
    Cell<T,Descriptor>& cell2 ( plint iX, plint iY )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        return fluid2_->get ( iX+relativeOffsetFluid2.x,iY+relativeOffsetFluid2.y );
    }
    T& mass ( plint iX, plint iY )
    {
        return mass_->get ( iX+relativeOffsetMass.x,iY+relativeOffsetMass.y );
    }
    T& mass2 ( plint iX, plint iY )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        return mass2_->get ( iX+relativeOffsetMass2.x,iY+relativeOffsetMass2.y );
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
    void setForce2 ( plint iX, plint iY, Array<T,2> const& force )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        force.to_cArray ( cell2 ( iX,iY ).getExternal ( forceOffset ) );
    }
    Array<T,2> getForce2 ( plint iX, plint iY )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        Array<T,2> force;
        force.from_cArray ( cell2 ( iX,iY ).getExternal ( forceOffset ) );
        return force;
    }
    void setMomentum ( plint iX, plint iY, Array<T,2> const& momentum )
    {
        j_->get ( iX+relativeOffsetJ.x,iY+relativeOffsetJ.y ) = momentum;
    }
    void setMomentum2 ( plint iX, plint iY, Array<T,2> const& momentum )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        j2_->get ( iX+relativeOffsetJ2.x,iY+relativeOffsetJ2.y ) = momentum;
    }
    Array<T,2> getMomentum ( plint iX, plint iY )
    {
        return j_->get ( iX+relativeOffsetJ.x,iY+relativeOffsetJ.y );
    }
    Array<T,2> getMomentum2 ( plint iX, plint iY )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        return j2_->get ( iX+relativeOffsetJ2.x,iY+relativeOffsetJ2.y );
    }
    T getDensity ( plint iX, plint iY )
    {
        return Descriptor<T>::fullRho (
                   rhoBar_->get ( iX+relativeOffsetRhoBar.x, iY+relativeOffsetRhoBar.y ) );
    }
    T getDensity2 ( plint iX, plint iY )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        return Descriptor<T>::fullRho (
                   rhoBar2_->get ( iX+relativeOffsetRhoBar2.x, iY+relativeOffsetRhoBar2.y ) );
    }
    void setDensity ( plint iX, plint iY, T rho )
    {
        rhoBar_->get ( iX+relativeOffsetRhoBar.x, iY+relativeOffsetRhoBar.y )
            = Descriptor<T>::rhoBar ( rho );
    }
    void setDensity2 ( plint iX, plint iY, T rho )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        rhoBar2_->get ( iX+relativeOffsetRhoBar2.x, iY+relativeOffsetRhoBar2.y )
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
        fluid_->attributeDynamics ( iX,iY,dynamics );
    }

    void attributeDynamics2 ( plint iX, plint iY, Dynamics<T,Descriptor>* dynamics )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        fluid2_->attributeDynamics ( iX,iY,dynamics );
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
    void addToTotalMass2 ( T addedTotalMass2 )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        fluid2_->getInternalStatistics().gatherSum ( 0, addedTotalMass2 );
    }
    void addToLostMass2 ( T addedLostMass2 )
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        fluid2_->getInternalStatistics().gatherSum ( 1, addedLostMass2 );
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
    T getSumMassMatrix2() const
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        return fluid2_->getInternalStatistics().getSum ( 0 );
    }
    T getSumLostMass2() const
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        return fluid2_->getInternalStatistics().getSum ( 1 );
    }
    T getTotalMass2() const
    {
        PLB_ASSERT ( !useFreeSurfaceLimit );
        return getSumMassMatrix2() + getSumLostMass2();
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
                if ( ! ( i == 0 && j == 0 ) &&
                        flag_->get ( nextX+relativeOffsetFS.x,nextY+relativeOffsetFS.y ) != wall )
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
    std::map<Node,T>& massExcess2()
    {
        PLB_ASSERT ( interfaceLists_ );
        return interfaceLists_ -> massExcess2;
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
    std::map<Node,ExtrapolInfo>& emptyToInterface()
    {
        PLB_ASSERT ( interfaceLists_ );
        return interfaceLists_ -> emptyToInterface;
    }
    std::map<Node,ExtrapolInfo>& fluidToInterface()
    {
        PLB_ASSERT ( interfaceLists_ );
        return interfaceLists_ -> fluidToInterface;
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
    bool useFreeSurfaceLimit;
    BlockLattice2D<T,Descriptor> *fluid_, *fluid2_;
    ScalarField2D<T> *rhoBar_, *rhoBar2_;
    TensorField2D<T,2> *j_, *j2_;
    ScalarField2D<T> *mass_, *mass2_;
    ScalarField2D<T>* volumeFraction_;
    ScalarField2D<int>* flag_;
    TensorField2D<T,2>* normal_;
    AtomicContainerBlock2D* containerInterfaceLists_;
    TwoPhaseInterfaceLists<T,Descriptor>* interfaceLists_;
    ScalarField2D<T>* curvature_;
    ScalarField2D<T>* outsideDensity_;

    Dot2D absoluteOffset, relativeOffsetFluid2, relativeOffsetRhoBar,
          relativeOffsetRhoBar2, relativeOffsetJ,
          relativeOffsetJ2, relativeOffsetMass, relativeOffsetMass2,
          relativeOffsetVF, relativeOffsetFS, relativeOffsetNormal, relativeOffsetC,
          relativeOffsetOD;

    static const int forceOffset = Descriptor<T>::ExternalField::forceBeginsAt;
};

/// Create a parameter-list for most free-surface data processors.
template< typename T,template<typename U> class Descriptor>
std::vector<MultiBlock2D*> aggregateTwoPhaseParams (
    MultiBlockLattice2D<T,Descriptor>& fluid,
    MultiBlockLattice2D<T,Descriptor>* fluid2,
    MultiScalarField2D<T>& rhoBar,
    MultiScalarField2D<T>* rhoBar2,
    MultiTensorField2D<T,2>& j,
    MultiTensorField2D<T,2>* j2,
    MultiScalarField2D<T>& mass,
    MultiScalarField2D<T>* mass2,
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
    if ( fluid2 )
    {
        PLB_ASSERT ( rhoBar2 );
        PLB_ASSERT ( j2 );
        PLB_ASSERT ( mass2 );
        aggregation.push_back ( fluid2 );
        aggregation.push_back ( rhoBar2 );
        aggregation.push_back ( j2 );
        aggregation.push_back ( mass2 );
    }

    return aggregation;
}

template< typename T,template<typename U> class Descriptor>
class DefaultInitializeTwoPhase2D : public BoxProcessingFunctional2D
{
public:
    DefaultInitializeTwoPhase2D ( Dynamics<T,Descriptor> *dynamicsTemplate_,
                                  Dynamics<T,Descriptor> *dynamicsTemplate2_,
                                  Array<T,2> g_, Array<T,2> g2_, T rhoIni_, bool useFreeSurfaceLimit_ )
        : dynamicsTemplate ( dynamicsTemplate_ ),
          dynamicsTemplate2 ( dynamicsTemplate2_ ),
          g ( g_ ), g2 ( g2_ ), rhoIni ( rhoIni_ ),
          useFreeSurfaceLimit ( useFreeSurfaceLimit_ )
    { }
    DefaultInitializeTwoPhase2D ( DefaultInitializeTwoPhase2D<T,Descriptor> const& rhs )
        : dynamicsTemplate ( rhs.dynamicsTemplate->clone() ),
          dynamicsTemplate2 ( rhs.dynamicsTemplate2->clone() ),
          g ( rhs.g ), g2 ( rhs.g2 ), rhoIni ( rhs.rhoIni ),
          useFreeSurfaceLimit ( rhs.useFreeSurfaceLimit )
    { }
    DefaultInitializeTwoPhase2D<T,Descriptor>& operator= ( DefaultInitializeTwoPhase2D<T,Descriptor> const& rhs )
    {
        DefaultInitializeTwoPhase2D<T,Descriptor> ( rhs ).swap ( *this );
        return *this;
    }
    void swap ( DefaultInitializeTwoPhase2D<T,Descriptor>& rhs )
    {
        std::swap ( dynamicsTemplate, rhs.dynamicsTemplate );
        std::swap ( dynamicsTemplate2, rhs.dynamicsTemplate2 );
        std::swap ( g, rhs.g );
        std::swap ( g2, rhs.g2 );
        std::swap ( rhoIni, rhs.rhoIni );
        std::swap ( useFreeSurfaceLimit, rhs.useFreeSurfaceLimit );
    }
    virtual ~DefaultInitializeTwoPhase2D()
    {
        delete dynamicsTemplate;
        delete dynamicsTemplate2;
    }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual DefaultInitializeTwoPhase2D<T,Descriptor>* clone() const
    {
        return new DefaultInitializeTwoPhase2D ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Volume-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if ( !useFreeSurfaceLimit )
        {
            modified[10] = modif::dataStructure;   // Fluid 2.
            modified[11] = modif::staticVariables;   // rhoBar2.
            modified[12] = modif::staticVariables;   // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }
private:
    Dynamics<T,Descriptor> *dynamicsTemplate, *dynamicsTemplate2;
    Array<T,2> g, g2;
    T rhoIni;
    bool useFreeSurfaceLimit;
};


template< typename T,template<typename U> class Descriptor>
class PartiallyDefaultInitializeTwoPhase2D : public BoxProcessingFunctional2D
{
public:
    PartiallyDefaultInitializeTwoPhase2D ( Dynamics<T,Descriptor> *dynamicsTemplate_,
                                           Dynamics<T,Descriptor> *dynamicsTemplate2_,
                                           Array<T,2> g_, Array<T,2> g2_, T rhoIni_, bool useFreeSurfaceLimit_ )
        : dynamicsTemplate ( dynamicsTemplate_ ),
          dynamicsTemplate2 ( dynamicsTemplate2_ ),
          g ( g_ ), g2 ( g2_ ), rhoIni ( rhoIni_ ),
          useFreeSurfaceLimit ( useFreeSurfaceLimit_ )
    { }
    PartiallyDefaultInitializeTwoPhase2D ( PartiallyDefaultInitializeTwoPhase2D<T,Descriptor> const& rhs )
        : dynamicsTemplate ( rhs.dynamicsTemplate->clone() ),
          dynamicsTemplate2 ( rhs.dynamicsTemplate2->clone() ),
          g ( rhs.g ), g2 ( rhs.g2 ), rhoIni ( rhs.rhoIni ),
          useFreeSurfaceLimit ( rhs.useFreeSurfaceLimit )
    { }
    PartiallyDefaultInitializeTwoPhase2D<T,Descriptor>& operator= ( PartiallyDefaultInitializeTwoPhase2D<T,Descriptor> const& rhs )
    {
        PartiallyDefaultInitializeTwoPhase2D<T,Descriptor> ( rhs ).swap ( *this );
        return *this;
    }
    void swap ( PartiallyDefaultInitializeTwoPhase2D<T,Descriptor>& rhs )
    {
        std::swap ( dynamicsTemplate, rhs.dynamicsTemplate );
        std::swap ( dynamicsTemplate2, rhs.dynamicsTemplate2 );
        std::swap ( g, rhs.g );
        std::swap ( g2, rhs.g2 );
        std::swap ( rhoIni, rhs.rhoIni );
        std::swap ( useFreeSurfaceLimit, rhs.useFreeSurfaceLimit );
    }
    virtual ~PartiallyDefaultInitializeTwoPhase2D()
    {
        delete dynamicsTemplate;
        delete dynamicsTemplate2;
    }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual PartiallyDefaultInitializeTwoPhase2D<T,Descriptor>* clone() const
    {
        return new PartiallyDefaultInitializeTwoPhase2D ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Volume-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if ( !useFreeSurfaceLimit )
        {
            modified[10] = modif::dataStructure;    // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }
private:
    Dynamics<T,Descriptor> *dynamicsTemplate, *dynamicsTemplate2;
    Array<T,2> g, g2;
    T rhoIni;
    bool useFreeSurfaceLimit;
};


// Functional to impose an initial constant velocity to one or to both fluids.
// CAUTION: This data processor must be called only after proper initialization of the TwoPhaseFields2D data structure.
template<typename T,template<typename U> class Descriptor, class Function>
class ConstantIniVelocityTwoPhase2D : public BoxProcessingFunctional2D
{
public:
    ConstantIniVelocityTwoPhase2D ( Array<T,2> velocity_, Function f_, bool imposeToFluid1_, bool imposeToFluid2_, bool useFreeSurfaceLimit_ )
        : velocity ( velocity_ ),
          f ( f_ ),
          imposeToFluid1 ( imposeToFluid1_ ),
          imposeToFluid2 ( imposeToFluid2_ ),
          useFreeSurfaceLimit ( useFreeSurfaceLimit_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual ConstantIniVelocityTwoPhase2D<T,Descriptor,Function>* clone() const
    {
        return new ConstantIniVelocityTwoPhase2D<T,Descriptor,Function> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0]  = modif::staticVariables; // Fluid.
        modified[1]  = modif::nothing;         // rhoBar.
        modified[2]  = modif::staticVariables; // j.
        modified[3]  = modif::nothing;         // Mass.
        modified[4]  = modif::nothing;         // Volume fraction.
        modified[5]  = modif::nothing;         // Flag-status.
        modified[6]  = modif::nothing;         // Normal.
        modified[7]  = modif::nothing;         // Interface-lists.
        modified[8]  = modif::nothing;         // Curvature.
        modified[9]  = modif::nothing;         // Outside density.
        if ( !useFreeSurfaceLimit )
        {
            modified[10] = modif::staticVariables; // Fluid 2.
            modified[11] = modif::nothing;         // rhoBar2.
            modified[12] = modif::staticVariables; // j2.
            modified[13] = modif::nothing;         // mass2.
        }
    }
private:
    Array<T,2> velocity;
    Function f;
    bool imposeToFluid1, imposeToFluid2;
    bool useFreeSurfaceLimit;
};


/// Compute the mass balance on every node in the domain, and store in mass matrix.
/** Input:
  *   - Flag-status:   needed in bulk+1
  *   - Mass:          needed in bulk
  *   - Volume fraction: needed in bulk
  *   - Populations:   needed in bulk+1
  * Output:
  *   - mass.
  **/
template< typename T,template<typename U> class Descriptor>
class TwoPhaseMassChange2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseMassChange2D ( bool useFreeSurfaceLimit_ )
        : useFreeSurfaceLimit ( useFreeSurfaceLimit_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseMassChange2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseMassChange2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;         // Fluid.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::nothing;         // Volume fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
        if ( !useFreeSurfaceLimit )
        {
            modified[10] = modif::nothing;        // Fluid 2.
            modified[11] = modif::nothing;   // rhoBar2.
            modified[12] = modif::nothing;   // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }
private:
    bool useFreeSurfaceLimit;
};

/// Completion scheme on the post-collide populations on interface cells.
/** Input:
  *   - Flag-status:   needed in bulk+1
  *   - Volume fraction: needed in bulk+1
  *   - Populations:   needed in bulk+1
  *   - Momentum:      needed in bulk+1
  *   - Density:       needed in bulk+1
  * Output:
  *   - Populations.
  **/
// ASK: This data processor loops over the whole volume. Is this really
//      necessary, or could one of the lists be used instead?
template<typename T, template<typename U> class Descriptor>
class TwoPhaseCompletion2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseCompletion2D ( bool useFreeSurfaceLimit_ )
        : useFreeSurfaceLimit ( useFreeSurfaceLimit_ )
    { }
    virtual TwoPhaseCompletion2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseCompletion2D<T,Descriptor> ( *this );
    }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;         // Fluid. Should be: staticVariables.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::nothing;         // Mass.
        modified[4] = modif::nothing;         // Volume fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
        if ( !useFreeSurfaceLimit )
        {
            modified[10] = modif::nothing;        // Fluid 2. Should be: staticVariables.
            modified[11] = modif::nothing;        // rhoBar2.
            modified[12] = modif::nothing;        // j2.
            modified[13] = modif::nothing;        // mass2.
        }
    }
private:
    bool useFreeSurfaceLimit;
};

/// Compute and store mass-fraction and macroscopic variables.
/** Input:
  *   - Flag-status:   needed in bulk
  *   - Mass:          needed in bulk
  *   - Populations:   needed in bulk
  * Output:
  *   - mass-fraction, density, momentum, flag (because setting bounce-back).
  **/
template<typename T, template<typename U> class Descriptor>
class TwoPhaseMacroscopic2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseMacroscopic2D ( T rhoDefault_, T densityRatio_, T surfaceTension_,
                            TwoPhaseModel model_ )
        : rhoDefault ( rhoDefault_ ),
          densityRatio ( densityRatio_ ),
          surfaceTension ( surfaceTension_ ),
          model ( model_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseMacroscopic2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseMacroscopic2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::staticVariables;  // Fluid. Should be: staticVariables.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass. Should be: staticVariables.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::staticVariables;  // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::staticVariables; // Fluid 2. Should be: staticVariables.
            modified[11] = modif::staticVariables; // rhoBar2.
            modified[12] = modif::staticVariables; // j2.
            modified[13] = modif::staticVariables; // mass2.
        }
    }
private:
    T rhoDefault;
    T densityRatio;
    T surfaceTension;
    TwoPhaseModel model;
};

template<typename T, template<typename U> class Descriptor>
class VerifyTwoPhase2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual VerifyTwoPhase2D<T,Descriptor>* clone() const
    {
        return new VerifyTwoPhase2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;          // Fluid. Should be: nothing.
        modified[1] = modif::nothing;  // rhoBar.
        modified[2] = modif::nothing;  // j.
        modified[3] = modif::nothing;  // Mass. Should be: nothing.
        modified[4] = modif::nothing;  // Volume fraction.
        modified[5] = modif::nothing;  // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        modified[10] = modif::nothing;         // Fluid 2. Should be: nothing.
        modified[11] = modif::nothing; // rhoBar2.
        modified[12] = modif::nothing; // j2.
        modified[13] = modif::nothing; // mass2.
    }
};

template<typename T, template<typename U> class Descriptor>
class TwoPhaseInterfaceFilter2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseInterfaceFilter2D ( T rhoDefault_, T densityRatio_, TwoPhaseModel model_ )
        : rhoDefault ( rhoDefault_ ),
          densityRatio ( densityRatio_ ),
          model ( model_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseInterfaceFilter2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseInterfaceFilter2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::staticVariables;  // Fluid. Should be: staticVariables.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass. Should be: staticVariables.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::staticVariables;  // Flag-status.
        modified[6] = modif::staticVariables;  // Normal.
        modified[7] = modif::staticVariables;  // Interface-lists.
        modified[8] = modif::staticVariables;  // Curvature.
        modified[9] = modif::staticVariables;  // Outside density.
        modified[10] = modif::dataStructure;   // Fluid 2. Should be: staticVariables.
        modified[11] = modif::staticVariables; // rhoBar2.
        modified[12] = modif::staticVariables; // j2.
        modified[13] = modif::staticVariables; // mass2.
    }
private:
    T rhoDefault;
    T densityRatio;
    TwoPhaseModel model;
};

/** Input:
  *   - interface-to-fluid list: needed in bulk+1
  *   - interface-to-empty list: needed in bulk+1
  *   - density: needed in bulk+1
  *   - mass:    needed in bulk+1
  *   - flag:    needed in bulk+1
  * Output:
  *   - flag, dynamics, mass, volumeFraction, density, force, momentum
  *   - mass-excess-list: defined in bulk+1
  **/
template<typename T,template<typename U> class Descriptor>
class TwoPhaseIniInterfaceToAnyNodes2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseIniInterfaceToAnyNodes2D (
        T rhoDefault_,
        Dynamics<T,Descriptor> *dynamicsTemplate_,
        Dynamics<T,Descriptor> *dynamicsTemplate2_,
        Array<T,2> const& force_, Array<T,2> const& force2_,
        T densityRatio_, T surfaceTension_, TwoPhaseModel model_ )
        : rhoDefault ( rhoDefault_ ),
          dynamicsTemplate ( dynamicsTemplate_ ),
          dynamicsTemplate2 ( dynamicsTemplate2_ ),
          force ( force_ ), force2 ( force2_ ), densityRatio ( densityRatio_ ),
          surfaceTension ( surfaceTension_ ), model ( model_ )
    { }
    TwoPhaseIniInterfaceToAnyNodes2D ( TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor> const& rhs )
        : rhoDefault ( rhs.rhoDefault ),
          dynamicsTemplate ( rhs.dynamicsTemplate->clone() ),
          dynamicsTemplate2 ( rhs.dynamicsTemplate2->clone() ),
          force ( rhs.force ), force2 ( rhs.force2 ),
          densityRatio ( rhs.densityRatio ),
          surfaceTension ( rhs.surfaceTension ),
          model ( rhs.model )
    { }
    TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor>& operator= ( TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor> const& rhs )
    {
        TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor> ( rhs ).swap ( *this );
        return *this;
    }
    void swap ( TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor>& rhs )
    {
        std::swap ( rhoDefault, rhs.rhoDefault );
        std::swap ( dynamicsTemplate, rhs.dynamicsTemplate );
        std::swap ( dynamicsTemplate2, rhs.dynamicsTemplate2 );
        std::swap ( force, rhs.force );
        std::swap ( force2, rhs.force2 );
        std::swap ( densityRatio, rhs.densityRatio );
        std::swap ( surfaceTension, rhs.surfaceTension );
        std::swap ( model, rhs.model );
    }
    virtual ~TwoPhaseIniInterfaceToAnyNodes2D()
    {
        delete dynamicsTemplate;
        delete dynamicsTemplate2;
    }
    TwoPhaseIniInterfaceToAnyNodes2D ( T rhoDefault_ );
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );

    virtual TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;           // Fluid. Gets assigned new dynamics. Should be: dataStructure
        modified[1] = modif::staticVariables;   // rhoBar.
        modified[2] = modif::nothing;           // j. Should be: staticVariables.
        modified[3] = modif::staticVariables;   // Mass. Is redistributed and initialized from neighborying density.
        modified[4] = modif::nothing;           // Volume fraction. Is default-initialized. Should be: staticVariables.
        modified[5] = modif::staticVariables;   // Flag-status. Is adapted according to cell-change lists.
        modified[6] = modif::nothing;           // Normal.
        modified[7] = modif::nothing;           // Interface-lists. Read-only.
        modified[8] = modif::nothing;           // Curvature.
        modified[9] = modif::nothing;           // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::nothing;          // Fluid 2. Should be: dataStructure
            modified[11] = modif::staticVariables;   // rhoBar2.
            modified[12] = modif::nothing;   // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }
private:
    T rhoDefault;
    Dynamics<T,Descriptor> *dynamicsTemplate, *dynamicsTemplate2;
    Array<T,2> force, force2;
    T densityRatio, surfaceTension;
    TwoPhaseModel model;
};

/// Based on the previously computed empty->interface list, initialize flow variables for
///   new interface cells.
/** Input:
  *   - Populations: needed in bulk+0
  *   - Momentum:    needed in bulk+1
  *   - Density:     needed in bulk+1
  *   - Flag-status: needed in bulk+0
  * Output:
  *   - flag-status:   initialized to "interface" on corresponding cells.
  *   - lattice:       initialized from neighbor averages on new interface cells.
  *   - mass:          initialized to zero on new interface cells.
  *   - mass-fraction: initialized to zero on new interface cells.
  *   - momentum
  **/
template<typename T,template<typename U> class Descriptor>
class TwoPhaseIniEmptyToInterfaceNodes2D: public BoxProcessingFunctional2D
{
public:
    TwoPhaseIniEmptyToInterfaceNodes2D ( Dynamics<T,Descriptor>* dynamicsTemplate_,
                                         Array<T,Descriptor<T>::d> force_, T densityRatio_, TwoPhaseModel model_ )
        : dynamicsTemplate ( dynamicsTemplate_ ), force ( force_ ), densityRatio ( densityRatio_ ), model ( model_ )
    { }
    TwoPhaseIniEmptyToInterfaceNodes2D ( TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor> const& rhs )
        : dynamicsTemplate ( rhs.dynamicsTemplate->clone() ),
          force ( rhs.force ),
          densityRatio ( rhs.densityRatio ),
          model ( rhs.model )
    { }
    TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor>& operator= (
        TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor> const& rhs )
    {
        TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor> ( rhs ).swap ( *this );
        return *this;
    }
    void swap ( TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor>& rhs )
    {
        std::swap ( dynamicsTemplate, rhs.dynamicsTemplate );
        std::swap ( force, rhs.force );
        std::swap ( densityRatio, rhs.densityRatio );
        std::swap ( model, rhs.model );
    }
    ~TwoPhaseIniEmptyToInterfaceNodes2D()
    {
        delete dynamicsTemplate;
    }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;         // Fluid. Should be: dataStructure
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::nothing;         // j. Should be: staticVariables.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::nothing;         // Volume fraction, read-only. Should be: staticVariables
        modified[5] = modif::staticVariables; // Flag-status, read-only.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists. Read access to gasCellToInitializeData.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::dataStructure;        // Fluid 2.
            modified[11] = modif::staticVariables;   // rhoBar2.
            modified[12] = modif::nothing;   // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }
private:
    Dynamics<T,Descriptor>* dynamicsTemplate;
    Array<T,Descriptor<T>::d> force; // Body force, for initialization of the new interface cell.
    T densityRatio;
    TwoPhaseModel model;
};

/// Isolated cells cannot be part of the interface. This data processor spots and
/// removes them.
/** Input:
  *   - Flag-status: needed in bulk+2
  *   - mass:        needed in bulk+1
  *   - density:     needed in bulk+1
  * Output:
  *   - interfaceToFluidNodes:   initialized in bulk+1
  *   - interfaceToEmptyNodes:   initialized in bulk+1
  *   - massExcess list:         initialized in bulk+1
  *   - mass, density, mass-fraction, dynamics, force, momentum, flag: in bulk+1
  **/
template<typename T,template<typename U> class Descriptor>
class TwoPhaseRemoveFalseInterfaceCells2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseRemoveFalseInterfaceCells2D ( T rhoDefault_, TwoPhaseModel model_ )
        : rhoDefault ( rhoDefault_ ),
          model ( model_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseRemoveFalseInterfaceCells2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseRemoveFalseInterfaceCells2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0]  = modif::nothing;          // Fluid: Gets NoDynamics when node changes to empty. Should be: dataStructure.
        modified[1]  = modif::staticVariables;  // rhoBar.
        modified[2]  = modif::nothing;          // j. Should be: staticVariables.
        modified[3]  = modif::staticVariables;  // Mass.
        modified[4]  = modif::nothing;          // Volume fraction. Should be: staticVariables.
        modified[5]  = modif::staticVariables;  // Flag-status.
        modified[6]  = modif::nothing;          // Normal.
        modified[7]  = modif::nothing;          // Interface lists.
        modified[8]  = modif::nothing;          // Curvature.
        modified[9]  = modif::nothing;          // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::nothing;    // Fluid 2.
            modified[11] = modif::staticVariables;   // rhoBar2.
            modified[12] = modif::nothing;   // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }
private:
    T rhoDefault;
    TwoPhaseModel model;
};

template<typename T,template<typename U> class Descriptor>
class TwoPhaseInitializeInterfaceLists2D : public BoxProcessingFunctional2D
{
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
    {
        PLB_ASSERT ( atomicBlocks.size() ==1 );

        AtomicContainerBlock2D* containerInterfaceLists = dynamic_cast<AtomicContainerBlock2D*> ( atomicBlocks[0] );
        PLB_ASSERT ( containerInterfaceLists );
        TwoPhaseInterfaceLists<T,Descriptor>* interfaceLists = new TwoPhaseInterfaceLists<T,Descriptor>;
        containerInterfaceLists->setData ( interfaceLists );

    }
    virtual TwoPhaseInitializeInterfaceLists2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseInitializeInterfaceLists2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        // Default-assign potential other parameters present in a multi-fluid system.
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::staticVariables;
    }
};

/// Wrapper for execution of TwoPhaseInitializeInterfaceLists2D.
template<typename T,template<typename U> class Descriptor>
void twoPhaseInitializeInterfaceLists2D ( MultiContainerBlock2D& interfaceListBlock )
{
    std::vector<MultiBlock2D*> arg;
    arg.push_back ( &interfaceListBlock );
    applyProcessingFunctional (
        new TwoPhaseInitializeInterfaceLists2D<T,Descriptor>,
        interfaceListBlock.getBoundingBox(), arg );
}

template< typename T,template<typename U> class Descriptor>
class TwoPhaseComputeInterfaceLists2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseComputeInterfaceLists2D ( TwoPhaseModel model_, T rhoDefault_, T densityRatio_, T surfaceTension_ )
        : model ( model_ ),
          rhoDefault ( rhoDefault_ ),
          densityRatio ( densityRatio_ ),
          surfaceTension ( surfaceTension_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseComputeInterfaceLists2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseComputeInterfaceLists2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::nothing;  // Fluid.
        modified[1] = modif::nothing;  // rhoBar.
        modified[2] = modif::nothing;  // j.
        modified[3] = modif::nothing;  // Mass.
        modified[4] = modif::nothing;  // Volume fraction.
        modified[5] = modif::nothing;  // Flag-status.
        modified[6] = modif::nothing;  // Normal.
        modified[7] = modif::nothing;  // Interface-lists.
        modified[8] = modif::nothing;  // Curvature.
        modified[9] = modif::nothing;  // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::nothing; // Fluid 2.
            modified[11] = modif::nothing; // rhoBar2.
            modified[12] = modif::nothing; // j2.
            modified[13] = modif::nothing; // mass2.
        }
    }
private:
    static T kappa; // Safety threshold for state-change, to prevent back-and-forth oscillations.
    TwoPhaseModel model;
    T rhoDefault;
    T densityRatio;
    T surfaceTension;
};

/// Enforce exact mass balance when interface cells become fluid or empty.
/** Input:
  *   - mass-excess list: needed in bulk+1
  *   - Flag-status: needed in bulk+2
  *   - mass:        needed in bulk+2
  *   - density:     needed in bulk+2
  * Output:
  *   - mass, mass-fraction
  **/
template<typename T,template<typename U> class Descriptor>
class TwoPhaseEqualMassExcessReDistribution2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseEqualMassExcessReDistribution2D ( TwoPhaseModel model_ )
        : model ( model_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseEqualMassExcessReDistribution2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseEqualMassExcessReDistribution2D ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::dataStructure;    // Fluid.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if ( model!=freeSurface )
        {
            modified[10] = modif::dataStructure;         // Fluid 2.
            modified[11] = modif::staticVariables;         // rhoBar2.
            modified[12] = modif::staticVariables;         // j2.
            modified[13] = modif::staticVariables; // mass2.
        }
    }
private:
    TwoPhaseModel model;
};

template<typename T, template<typename U> class Descriptor>
struct TwoPhaseFields2D
{
    //static const int envelopeWidth = 2;
    static const int envelopeWidth = 3; // Necessary when we use height functions to compute the curvature,
    // or when double smoothing is used at the data processor that computes the
    // normals from the volume fraction.
    //static const int envelopeWidth = 4; // Necessary when we use height functions to compute the curvature and
    //   use the old contact angle algorithm.
    static const int smallEnvelopeWidth = 2;

    static const int envelopeWidthForImmersedWalls = 4;

    TwoPhaseFields2D ( SparseBlockStructure2D const& blockStructure,
                       Dynamics<T,Descriptor>* dynamics_, Dynamics<T,Descriptor>* dynamics2_,
                       T rhoDefault_, T densityRatio_, T surfaceTension_, T contactAngle_, Array<T,2> force_,
                       Array<T,2> force2_, TwoPhaseModel model_, bool useImmersedWalls = false )
    // model=true --> model adapted for air bubbles, which are prevented from collapsing.
    // model=false --> model adapted for water droplets, where shearing stresses are properly treated.
        : dynamics ( dynamics_ ),
          dynamics2 ( dynamics2_ ),
          rhoDefault ( rhoDefault_ ),
          densityRatio ( densityRatio_ ),
          surfaceTension ( surfaceTension_ ), contactAngle ( contactAngle_ ),
          force ( force_ ), force2 ( force2_ ),
          model ( model_ ),
          lattice (
              MultiBlockManagement2D (
                  blockStructure, defaultMultiBlockPolicy2D().getThreadAttribution(),
                  smallEnvelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>(), dynamics->clone() ),
          lattice2 ( 0 ),
          helperLists ( lattice ),
          mass ( lattice ),
          mass2 ( 0 ),
          flag (
              MultiBlockManagement2D (
                  blockStructure,
                  defaultMultiBlockPolicy2D().getThreadAttribution(),
                  useImmersedWalls ? envelopeWidthForImmersedWalls : envelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiScalarAccess<int>() ),
          volumeFraction ( ( MultiBlock2D& ) flag ),
          curvature (
              MultiBlockManagement2D (
                  blockStructure,
                  defaultMultiBlockPolicy2D().getThreadAttribution(),
                  envelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiScalarAccess<T>() ),
          outsideDensity ( ( MultiBlock2D& ) curvature ),
          rhoBar (
              MultiBlockManagement2D (
                  blockStructure,
                  defaultMultiBlockPolicy2D().getThreadAttribution(),
                  useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiScalarAccess<T>() ),
          rhoBar2 ( 0 ),
          j (
              MultiBlockManagement2D (
                  blockStructure,
                  defaultMultiBlockPolicy2D().getThreadAttribution(),
                  useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiTensorAccess<T,2>() ),
          j2 ( 0 ),
          normal ( ( MultiBlock2D& ) curvature )
    {
        initialization ( blockStructure, useImmersedWalls );
    }
    // Constructor for free-surface case
    TwoPhaseFields2D ( SparseBlockStructure2D const& blockStructure,
                       Dynamics<T,Descriptor>* dynamics_, T rhoDefault_, T surfaceTension_,
                       T contactAngle_, Array<T,2> force_, bool useImmersedWalls = false )
    // model=true --> model adapted for air bubbles, which are prevented from collapsing.
    // model=false --> model adapted for water droplets, where shearing stresses are properly treated.
        : dynamics ( dynamics_->clone() ),
          dynamics2 ( dynamics_ ),
          rhoDefault ( rhoDefault_ ),
          densityRatio ( T() ),
          surfaceTension ( surfaceTension_ ), contactAngle ( contactAngle_ ),
          force ( force_ ), force2 ( force_ ),
          model ( freeSurface ),
          lattice (
              MultiBlockManagement2D (
                  blockStructure, defaultMultiBlockPolicy2D().getThreadAttribution(),
                  smallEnvelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>(), dynamics->clone() ),
          lattice2 ( 0 ),
          helperLists ( lattice ),
          mass ( lattice ),
          mass2 ( 0 ),
          flag (
              MultiBlockManagement2D (
                  blockStructure,
                  defaultMultiBlockPolicy2D().getThreadAttribution(),
                  useImmersedWalls ? envelopeWidthForImmersedWalls : envelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiScalarAccess<int>() ),
          volumeFraction ( ( MultiBlock2D& ) flag ),
          curvature (
              MultiBlockManagement2D (
                  blockStructure,
                  defaultMultiBlockPolicy2D().getThreadAttribution(),
                  envelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiScalarAccess<T>() ),
          outsideDensity ( ( MultiBlock2D& ) curvature ),
          rhoBar (
              MultiBlockManagement2D (
                  blockStructure,
                  defaultMultiBlockPolicy2D().getThreadAttribution(),
                  useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiScalarAccess<T>() ),
          rhoBar2 ( 0 ),
          j (
              MultiBlockManagement2D (
                  blockStructure,
                  defaultMultiBlockPolicy2D().getThreadAttribution(),
                  useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
              defaultMultiBlockPolicy2D().getBlockCommunicator(),
              defaultMultiBlockPolicy2D().getCombinedStatistics(),
              defaultMultiBlockPolicy2D().getMultiTensorAccess<T,2>() ),
          j2 ( 0 ),
          normal ( ( MultiBlock2D& ) curvature )
    {
        initialization ( blockStructure, useImmersedWalls );
    }
    void initialization ( SparseBlockStructure2D const& blockStructure, bool useImmersedWalls )
    {
        Precision precision;
        if ( sizeof ( T ) == sizeof ( float ) )
            precision = FLT;
        else if ( sizeof ( T ) == sizeof ( double ) )
            precision = DBL;
        else if ( sizeof ( T ) == sizeof ( long double ) )
            precision = LDBL;
        else
            PLB_ASSERT ( false );

        T eps = getEpsilon<T> ( precision );
        // The contact angle must take values between 0 and 180 degrees. If it is negative,
        // this means that contact angle effects will not be modeled.
        PLB_ASSERT ( contactAngle < ( T ) 180.0 || fabs ( contactAngle - ( T ) 180.0 ) <= eps );

        if ( fabs ( surfaceTension ) <= eps )
        {
            useSurfaceTension = false;
        }
        else
        {
            useSurfaceTension = true;
        }

        if ( model!=freeSurface )
        {
            lattice2 = new MultiBlockLattice2D<T,Descriptor> (
                MultiBlockManagement2D (
                    blockStructure, defaultMultiBlockPolicy2D().getThreadAttribution(),
                    smallEnvelopeWidth ),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>(), dynamics2->clone() );
            mass2 = new MultiScalarField2D<T> ( lattice );
            rhoBar2 = new MultiScalarField2D<T> (
                MultiBlockManagement2D (
                    blockStructure,
                    defaultMultiBlockPolicy2D().getThreadAttribution(),
                    useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiScalarAccess<T>() );
            j2 = new MultiTensorField2D<T,2> (
                MultiBlockManagement2D (
                    blockStructure,
                    defaultMultiBlockPolicy2D().getThreadAttribution(),
                    useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiTensorAccess<T,2>() );
            lattice2->periodicity().toggleAll ( true );
            mass2->periodicity().toggleAll ( true );
            rhoBar2->periodicity().toggleAll ( true );
            j2->periodicity().toggleAll ( true );
        }

        twoPhaseArgs = aggregateTwoPhaseParams ( lattice, lattice2, rhoBar, rhoBar2, j, j2,
                       mass, mass2, volumeFraction,
                       flag, normal, helperLists, curvature, outsideDensity );

        twoPhaseInitializeInterfaceLists2D<T,Descriptor> ( helperLists );
        lattice.periodicity().toggleAll ( true );
        mass.periodicity().toggleAll ( true );
        flag.periodicity().toggleAll ( true );
        volumeFraction.periodicity().toggleAll ( true );
        curvature.periodicity().toggleAll ( true );
        outsideDensity.periodicity().toggleAll ( true );
        rhoBar.periodicity().toggleAll ( true );
        j.periodicity().toggleAll ( true );
        normal.periodicity().toggleAll ( true );
        setToConstant ( flag, flag.getBoundingBox(), ( int ) twoPhaseFlag::empty );
        setToConstant ( outsideDensity, outsideDensity.getBoundingBox(), rhoDefault );
        rhoBarJparam.push_back ( &lattice );
        rhoBarJparam.push_back ( &rhoBar );
        rhoBarJparam.push_back ( &j );

        if ( model!=freeSurface )
        {
            rhoBarJparam2.push_back ( lattice2 );
            rhoBarJparam2.push_back ( rhoBar2 );
            rhoBarJparam2.push_back ( j2 );
        }

        lattice.internalStatSubscription().subscribeSum();     // Total mass.
        lattice.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice.internalStatSubscription().subscribeIntSum();  // Num interface cells.
        if ( model!=freeSurface )
        {
            lattice2->internalStatSubscription().subscribeSum();     // Total mass 2.
            lattice2->internalStatSubscription().subscribeSum();     // Lost mass 2.
        }

        freeSurfaceDataProcessors();
    }

    TwoPhaseFields2D ( TwoPhaseFields2D<T,Descriptor> const& rhs )
        : dynamics ( rhs.dynamics->clone() ),
          dynamics2 ( rhs.dynamics2->clone() ),
          rhoDefault ( rhs.rhoDefault ),
          densityRatio ( rhs.densityRatio ),
          surfaceTension ( rhs.surfaceTension ),
          contactAngle ( rhs.contactAngle ),
          useSurfaceTension ( rhs.useSurfaceTension ),
          force ( rhs.force ),
          force2 ( rhs.force2 ),
          model ( rhs.model ),
          lattice ( rhs.lattice ),
          lattice2 ( rhs.lattice2 ? new MultiBlockLattice2D<T,Descriptor> ( *rhs.lattice2 ) : 0 ),
          helperLists ( rhs.helperLists ),
          mass ( rhs.mass ),
          mass2 ( rhs.mass2 ? new MultiScalarField2D<T> ( *rhs.mass2 ) : 0 ),
          flag ( rhs.flag ),
          volumeFraction ( rhs.volumeFraction ),
          curvature ( rhs.curvature ),
          outsideDensity ( rhs.outsideDensity ),
          rhoBar ( rhs.rhoBar ),
          rhoBar2 ( rhs.rhoBar2 ? new MultiScalarField2D<T> ( *rhs.rhoBar2 ) : 0 ),
          j ( rhs.j ),
          j2 ( rhs.j2 ? new MultiTensorField2D<T,2> ( *rhs.j2 ) : 0 ),
          normal ( rhs.normal ),
          rhoBarJparam(),
          rhoBarJparam2(),
          twoPhaseArgs()
    {
        twoPhaseArgs = aggregateTwoPhaseParams ( lattice, lattice2, rhoBar, rhoBar2, j, j2,
                       mass, mass2, volumeFraction,
                       flag, normal, helperLists, curvature, outsideDensity );

        rhoBarJparam.push_back ( &lattice );
        rhoBarJparam.push_back ( &rhoBar );
        rhoBarJparam.push_back ( &j );

        if ( model!=freeSurface )
        {
            rhoBarJparam2.push_back ( &lattice2 );
            rhoBarJparam2.push_back ( &rhoBar2 );
            rhoBarJparam2.push_back ( &j2 );
        }
    }

    void swap ( TwoPhaseFields2D<T,Descriptor>& rhs )
    {
        std::swap ( dynamics, rhs.dynamics );
        std::swap ( dynamics2, rhs.dynamics2 );
        std::swap ( rhoDefault, rhs.rhoDefault );
        std::swap ( densityRatio, rhs.densityRatio );
        std::swap ( surfaceTension, rhs.surfaceTension );
        std::swap ( contactAngle, rhs.contactAngle );
        std::swap ( useSurfaceTension, rhs.useSurfaceTension );
        std::swap ( force, rhs.force );
        std::swap ( force2, rhs.force2 );
        std::swap ( model, rhs.model );
        std::swap ( lattice, rhs.lattice );
        std::swap ( lattice2, rhs.lattice2 );
        std::swap ( helperLists, rhs.helperLists );
        std::swap ( mass, rhs.mass );
        std::swap ( mass2, rhs.mass2 );
        std::swap ( flag, rhs.flag );
        std::swap ( volumeFraction, rhs.volumeFraction );
        std::swap ( curvature, rhs.curvature );
        std::swap ( outsideDensity, rhs.outsideDensity );
        std::swap ( rhoBar, rhs.rhoBar );
        std::swap ( rhoBar2, rhs.rhoBar2 );
        std::swap ( j2, rhs.j2 );
        std::swap ( normal, rhs.normal );
        std::swap ( rhoBarJparam, rhs.rhoBarJparam );
        std::swap ( rhoBarJparam2, rhs.rhoBarJparam2 );
        std::swap ( twoPhaseArgs, rhs.twoPhaseArgs );
    }

    TwoPhaseFields2D<T,Descriptor>& operator= ( TwoPhaseFields2D<T,Descriptor> const& rhs )
    {
        TwoPhaseFields2D<T,Descriptor> ( rhs ).swap ( *this );
        return *this;
    }

    TwoPhaseFields2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseFields2D<T,Descriptor> ( *this );
    }

    ~TwoPhaseFields2D()
    {
        delete dynamics;
        delete dynamics2;
        delete lattice2;
        delete rhoBar2;
        delete j2;
        delete mass2;
    }

    void periodicityToggle ( plint direction, bool periodic )
    {
        PLB_ASSERT ( direction == 0 || direction == 1 || direction == 2 );

        lattice.periodicity().toggle ( direction, periodic );
        mass.periodicity().toggle ( direction, periodic );
        flag.periodicity().toggle ( direction, periodic );
        volumeFraction.periodicity().toggle ( direction, periodic );
        curvature.periodicity().toggle ( direction, periodic );
        outsideDensity.periodicity().toggle ( direction, periodic );
        rhoBar.periodicity().toggle ( direction, periodic );
        j.periodicity().toggle ( direction, periodic );
        normal.periodicity().toggle ( direction, periodic );
        if ( model!=freeSurface )
        {
            lattice2->periodicity().toggle ( direction, periodic );
            mass2->periodicity().toggle ( direction, periodic );
            rhoBar2->periodicity().toggle ( direction, periodic );
            j2->periodicity().toggle ( direction, periodic );
        }
    }

    void periodicityToggleAll ( bool periodic )
    {
        lattice.periodicity().toggleAll ( periodic );
        mass.periodicity().toggleAll ( periodic );
        flag.periodicity().toggleAll ( periodic );
        volumeFraction.periodicity().toggleAll ( periodic );
        curvature.periodicity().toggleAll ( periodic );
        outsideDensity.periodicity().toggleAll ( periodic );
        rhoBar.periodicity().toggleAll ( periodic );
        j.periodicity().toggleAll ( periodic );
        normal.periodicity().toggleAll ( periodic );
        if ( model!=freeSurface )
        {
            lattice2->periodicity().toggleAll ( periodic );
            mass2->periodicity().toggleAll ( periodic );
            rhoBar2->periodicity().toggleAll ( periodic );
            j2->periodicity().toggleAll ( periodic );
        }
    }

    void defaultInitialize()
    {
        applyProcessingFunctional (
            new DefaultInitializeTwoPhase2D<T,Descriptor> (
                dynamics->clone(), dynamics2->clone(), force, force2, rhoDefault, model==freeSurface ),
            lattice.getBoundingBox(), twoPhaseArgs );
    }

    void partiallyDefaultInitialize()
    {
        applyProcessingFunctional (
            new PartiallyDefaultInitializeTwoPhase2D<T,Descriptor> (
                dynamics->clone(), dynamics2->clone(), force, force2, rhoDefault, model==freeSurface ),
            lattice.getBoundingBox(), twoPhaseArgs );
    }
    std::auto_ptr<MultiScalarField2D<T> > computePressure ( Box2D domain, bool computeFluid1, bool computeFluid2 )
    {
        std::auto_ptr<MultiScalarField2D<T> > pressure =
            generateMultiScalarField<T> ( lattice, domain );
        std::vector<MultiBlock2D*> args ( twoPhaseArgs );
        args.push_back ( pressure.get() );
        applyProcessingFunctional (
            new TwoPhaseComputePressure2D<T,Descriptor> ( densityRatio, rhoDefault, model, computeFluid1, computeFluid2 ),
            domain, args );
        return pressure;
    }
    std::auto_ptr<MultiScalarField2D<T> > computePressure ( Box2D domain )
    {
        return computePressure ( domain, true, true );
    }
    std::auto_ptr<MultiScalarField2D<T> > computePressure()
    {
        return computePressure ( lattice.getBoundingBox(), true, true );
    }
    void computePressureAverage ( Box2D domain, T& p1, T& p2 )
    {
        TwoPhaseAveragePressure2D<T,Descriptor> functional ( densityRatio, rhoDefault, model );
        applyProcessingFunctional ( functional, domain, twoPhaseArgs );
        p1 = functional.getAveragePressure();
        p2 = functional.getAveragePressure2();
    }
    void computeVelocityAverage ( Box2D domain, Array<T,2>& v1, Array<T,2>& v2 )
    {
        TwoPhaseAverageVelocity2D<T,Descriptor> functional ( model );
        applyProcessingFunctional ( functional, domain, twoPhaseArgs );
        v1 = functional.getAverageVelocity();
        v2 = functional.getAverageVelocity2();
    }
    std::auto_ptr<MultiTensorField2D<T,2> > computeVelocity ( Box2D domain, bool computeFluid1, bool computeFluid2 )
    {
        std::auto_ptr<MultiTensorField2D<T,2> > velocity =
            generateMultiTensorField<T,2> ( lattice, domain );
        std::vector<MultiBlock2D*> args ( twoPhaseArgs );
        args.push_back ( velocity.get() );
        applyProcessingFunctional (
            new TwoPhaseComputeVelocity2D<T,Descriptor> ( densityRatio, computeFluid1, computeFluid2, model==freeSurface ),
            domain, args );
        return velocity;
    }
    std::auto_ptr<MultiTensorField2D<T,2> > computeVelocity ( Box2D domain )
    {
        return computeVelocity ( domain, true, true );
    }
    std::auto_ptr<MultiTensorField2D<T,2> > computeVelocity()
    {
        return computeVelocity ( lattice.getBoundingBox(), true, true );
    }

    void freeSurfaceDataProcessors()
    {
        plint pl; // Processor level.

        /***** Initial level ******/
        pl = 0;

        integrateProcessingFunctional (
            new ExternalRhoJcollideAndStream2D<T,Descriptor>,
            lattice.getBoundingBox(), rhoBarJparam, pl );

        // twophase
        if ( model!=freeSurface )
        {
            integrateProcessingFunctional (
                new ExternalRhoJcollideAndStream2D<T,Descriptor>,
                lattice.getBoundingBox(), rhoBarJparam2, pl );
        }

        integrateProcessingFunctional (
            new TwoPhaseComputeNormals2D<T,Descriptor>,
            lattice.getBoundingBox(), twoPhaseArgs, pl );


        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new TwoPhaseMassChange2D<T,Descriptor> ( model==freeSurface ), lattice.getBoundingBox(),
            twoPhaseArgs, pl );

        integrateProcessingFunctional (
            new TwoPhaseCompletion2D<T,Descriptor> ( model==freeSurface ),
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        integrateProcessingFunctional (
            new TwoPhaseMacroscopic2D<T,Descriptor> ( rhoDefault, densityRatio, surfaceTension, model ),
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        if ( useSurfaceTension )
        {
            integrateProcessingFunctional (
                new TwoPhaseComputeCurvature2D<T,Descriptor> ( contactAngle, lattice.getBoundingBox() ),
                lattice.getBoundingBox(), twoPhaseArgs, pl );
        }

        /***** New level ******/
        //pl++;

        //integrateProcessingFunctional (
        //        new TwoPhaseInterfaceFilter2D<T,Descriptor>(rhoDefault, densityRatio, model),
        //        lattice.getBoundingBox(), twoPhaseArgs, pl );

        /***** New level ******/
        //pl++;

        // integrateProcessingFunctional (
        //         new TwoPhaseInterfaceFilter2D<T,Descriptor>(rhoDefault, densityRatio, model),
        //         lattice.getBoundingBox(), twoPhaseArgs, pl );

        if ( useSurfaceTension )
        {
            integrateProcessingFunctional (
                new TwoPhaseAddSurfaceTension2D<T,Descriptor> ( surfaceTension ),
                lattice.getBoundingBox(), twoPhaseArgs, pl );
        }

        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new TwoPhaseComputeInterfaceLists2D<T,Descriptor> ( model, rhoDefault, densityRatio, surfaceTension ),
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        //interface->fluid   --     interface->empty
        //interface->empty   --     interface->fluid
        //fluid->interface   --     empty->interface
        integrateProcessingFunctional (
            new TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor> (
                rhoDefault, dynamics->clone(), dynamics2->clone(), force, force2,
                densityRatio, surfaceTension, model ),
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        //empty->interface  --      fluid->interface
        integrateProcessingFunctional (
            new TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor> ( dynamics->clone(), force, densityRatio, model ),
            lattice.getBoundingBox(),
            twoPhaseArgs, pl );

        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new TwoPhaseRemoveFalseInterfaceCells2D<T,Descriptor> ( rhoDefault, model ),
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new TwoPhaseEqualMassExcessReDistribution2D<T,Descriptor> ( model ),
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        integrateProcessingFunctional (
            new TwoPhaseComputeStatistics2D<T,Descriptor>,
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        //integrateProcessingFunctional (
        //    new VerifyTwoPhase2D<T,Descriptor>(),
        //    lattice.getBoundingBox(), twoPhaseArgs, pl);
    }

    Dynamics<T,Descriptor> *dynamics, *dynamics2;
    T rhoDefault;
    T densityRatio;
    T surfaceTension;
    T contactAngle;
    bool useSurfaceTension;
    Array<T,2> force, force2;
    TwoPhaseModel model;
    MultiBlockLattice2D<T, Descriptor> lattice;
    MultiBlockLattice2D<T, Descriptor> *lattice2;
    MultiContainerBlock2D helperLists;
    MultiScalarField2D<T> mass;
    MultiScalarField2D<T>* mass2;
    MultiScalarField2D<int> flag;
    MultiScalarField2D<T> volumeFraction;
    MultiScalarField2D<T> curvature;
    MultiScalarField2D<T> outsideDensity;
    MultiScalarField2D<T> rhoBar;
    MultiScalarField2D<T>* rhoBar2;
    MultiTensorField2D<T,2> j;
    MultiTensorField2D<T,2>* j2;
    MultiTensorField2D<T,2> normal;
    std::vector<MultiBlock2D*> rhoBarJparam, rhoBarJparam2;
    std::vector<MultiBlock2D*> twoPhaseArgs;
};

template<typename T, template<typename U> class Descriptor>
class TwoPhaseOutletMaximumVolumeFraction2D : public BoxProcessingFunctional2D
{
public:
    TwoPhaseOutletMaximumVolumeFraction2D ( T volumeFraction_, TwoPhaseModel model_ )
        : volumeFraction ( volumeFraction_ ),
          model ( model_ )
    { }
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual TwoPhaseOutletMaximumVolumeFraction2D<T,Descriptor>* clone() const
    {
        return new TwoPhaseOutletMaximumVolumeFraction2D<T,Descriptor> ( *this );
    }
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        PLB_ASSERT ( modified.size() >=14 );
        std::fill ( modified.begin(), modified.end(), modif::nothing );
        modified[0] = modif::staticVariables;  // Fluid.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        modified[10] = modif::staticVariables; // Fluid 2.
        modified[11] = modif::staticVariables; // rhoBar2.
        modified[12] = modif::staticVariables; // j2.
        modified[13] = modif::staticVariables; // mass2.
    }
private:
    T volumeFraction;
    TwoPhaseModel model;
};

}  // namespace plb

#endif  // TWO_PHASE_MODEL_2D_H

