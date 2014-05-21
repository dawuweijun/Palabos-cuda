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

#ifndef OFF_LATTICE_BOUNDARY_CONDITION_2D_HH
#define OFF_LATTICE_BOUNDARY_CONDITION_2D_HH

#include "core/globalDefs.h"
#include "offLattice/offLatticeBoundaryCondition2D.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "offLattice/offLatticeBoundaryProfiles2D.h"
#include "triangleToDef.h"

namespace plb {


/* ********** OffLatticeBoundaryCondition2D ********************************** */

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::
    OffLatticeBoundaryCondition2D (
        OffLatticeModel2D<T,BoundaryType>* offLatticeModel_,
        VoxelizedDomain2D<T>& voxelizedDomain_,
        MultiBlockLattice2D<T,Descriptor>& lattice_ )
    : voxelizedDomain(voxelizedDomain_),
      lattice(lattice_),
      boundaryShapeArg(lattice_),
      offLatticeModel(offLatticeModel_),
      offLatticePattern(lattice)
{
    std::vector<MultiBlock2D*> offLatticeIniArg;
    // First argument for compute-off-lattice-pattern.
    offLatticeIniArg.push_back(&offLatticePattern);
    // Remaining arguments for inner-flow-shape.
    offLatticeIniArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeIniArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeIniArg.push_back(&boundaryShapeArg);
    applyProcessingFunctional (
            new OffLatticePatternFunctional2D<T,BoundaryType> (
                offLatticeModel->clone() ),
            offLatticePattern.getBoundingBox(), offLatticeIniArg );
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::
    OffLatticeBoundaryCondition2D (
        OffLatticeModel2D<T,BoundaryType>* offLatticeModel_,
        VoxelizedDomain2D<T>& voxelizedDomain_,
        MultiBlockLattice2D<T,Descriptor>& lattice_,
        MultiParticleField2D<DenseParticleField2D<T,Descriptor> >& particleField_ )
    : offLatticeModel(offLatticeModel_),
      voxelizedDomain(voxelizedDomain_),
      lattice(lattice_),
      boundaryShapeArg(particleField_),
      offLatticePattern(lattice)
{
    std::vector<MultiBlock2D*> offLatticeIniArg;
    // First argument for compute-off-lattice-pattern.
    offLatticeIniArg.push_back(&offLatticePattern);
    // Remaining arguments for inner-flow-shape.
    offLatticeIniArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeIniArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeIniArg.push_back(&boundaryShapeArg);
    applyProcessingFunctional (
            new OffLatticePatternFunctional2D<T,BoundaryType> (
                offLatticeModel->clone() ),
            offLatticePattern.getBoundingBox(), offLatticeIniArg );
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::OffLatticeBoundaryCondition2D (
        OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType> const& rhs )
    : offLatticeModel(rhs.offLatticeModel.clone()),
      voxelizedDomain(rhs.voxelizedDomain),
      lattice(rhs.lattice),
      boundaryShapeArg(rhs.boundaryShapeArg),
      offLatticePattern(rhs.offLatticePattern)
{ }


template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::~OffLatticeBoundaryCondition2D()
{
    delete offLatticeModel;
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
void OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::insert()
{
    std::vector<MultiBlock2D*> offLatticeArg;
    // First three arguments for Guo algorithm.
    offLatticeArg.push_back(&lattice);
    offLatticeArg.push_back(&offLatticePattern);
    // Remaining arguments for inner-flow-shape.
    offLatticeArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeArg.push_back(&boundaryShapeArg);
    plint processorLevel = 1;
    plint numShapeArgs = 3;
    plint numCompletionArgs = 0;
    integrateProcessingFunctional (
            new OffLatticeCompletionFunctional2D<T,Descriptor,BoundaryType> (
                offLatticeModel->clone(), numShapeArgs, numCompletionArgs ),
            boundaryShapeArg.getBoundingBox(), offLatticeArg, processorLevel );
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
void OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::insert (
        std::vector<MultiBlock2D*> const& completionArg )
{
    std::vector<MultiBlock2D*> offLatticeArg;
    // First three arguments for Guo algorithm.
    offLatticeArg.push_back(&lattice);
    offLatticeArg.push_back(&offLatticePattern);
    // Next arguments for inner-flow-shape.
    offLatticeArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeArg.push_back(&boundaryShapeArg);
    // Remaining are optional arguments for completion algorithm.
    plint numCompletionArgs = (plint)completionArg.size();
    for (plint i=0; i<numCompletionArgs; ++i) {
        offLatticeArg.push_back(completionArg[i]);
    }
    plint processorLevel = 1;
    plint numShapeArgs = 3;
    integrateProcessingFunctional (
            new OffLatticeCompletionFunctional2D<T,Descriptor,BoundaryType> (
                offLatticeModel->clone(), numShapeArgs, numCompletionArgs ),
            boundaryShapeArg.getBoundingBox(), offLatticeArg, processorLevel );
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
void OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::apply()
{
    std::vector<MultiBlock2D*> offLatticeArg;
    // First three arguments for Guo algorithm.
    offLatticeArg.push_back(&lattice);
    offLatticeArg.push_back(&offLatticePattern);
    // Remaining arguments for inner-flow-shape.
    offLatticeArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeArg.push_back(&boundaryShapeArg);
    plint numShapeArgs = 3;
    plint numCompletionArgs = 0;
    applyProcessingFunctional (
            new OffLatticeCompletionFunctional2D<T,Descriptor,BoundaryType> (
                offLatticeModel->clone(), numShapeArgs, numCompletionArgs ),
            boundaryShapeArg.getBoundingBox(), offLatticeArg );
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
void OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::apply (
        std::vector<MultiBlock2D*> const& completionArg )
{
    std::vector<MultiBlock2D*> offLatticeArg;
    // First three arguments for Guo algorithm.
    offLatticeArg.push_back(&lattice);
    offLatticeArg.push_back(&offLatticePattern);
    // Next arguments for inner-flow-shape.
    offLatticeArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeArg.push_back(&boundaryShapeArg);
    // Remaining are optional arguments for completion algorithm.
    plint numCompletionArgs = (plint)completionArg.size();
    for (plint i=0; i<numCompletionArgs; ++i) {
        offLatticeArg.push_back(completionArg[i]);
    }
    plint numShapeArgs = 3;
    applyProcessingFunctional (
            new OffLatticeCompletionFunctional2D<T,Descriptor,BoundaryType> (
                offLatticeModel->clone(), numShapeArgs, numCompletionArgs ),
            boundaryShapeArg.getBoundingBox(), offLatticeArg );
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
Array<T,2> OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::getForceOnObject()
{
    std::vector<MultiBlock2D*> arg;
    arg.push_back(&offLatticePattern);
    GetForceOnObjectFunctional2D<T,BoundaryType> functional(offLatticeModel->clone());
    applyProcessingFunctional (
            functional, boundaryShapeArg.getBoundingBox(), arg );
    return functional.getForce();
}

    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiTensorField2D<T,2> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeVelocity(Box2D domain)
{
    std::auto_ptr<MultiTensorField2D<T,2> > velocity(plb::computeVelocity(lattice,domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant<T,2>(*velocity, voxelizedDomain.getVoxelMatrix(), solidFlag,
                       domain, Array<T,2>(T(),T(),T()));
    setToConstant<T,2>(*velocity, voxelizedDomain.getVoxelMatrix(), solidBorderFlag,
                       domain, Array<T,2>(T(),T(),T()));
    return velocity;
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiTensorField2D<T,2> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeVelocity()
{
    return computeVelocity(lattice.getBoundingBox());
}

    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiTensorField2D<T,2> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeVorticity(Box2D domain)
{
    std::auto_ptr<MultiTensorField2D<T,2> > vorticity (
            plb::computeBulkVorticity(*plb::computeVelocity(lattice,domain), domain) );
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant<T,2>(*vorticity, voxelizedDomain.getVoxelMatrix(), solidFlag,
                       domain, Array<T,2>(T(),T(),T()));
    setToConstant<T,2>(*vorticity, voxelizedDomain.getVoxelMatrix(), solidBorderFlag,
                       domain, Array<T,2>(T(),T(),T()));
    return vorticity;
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiTensorField2D<T,2> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeVorticity()
{
    return computeVorticity(lattice.getBoundingBox());
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeVelocityNorm(Box2D domain)
{
    std::auto_ptr<MultiScalarField2D<T> > velocityNorm(plb::computeVelocityNorm(lattice,domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*velocityNorm, voxelizedDomain.getVoxelMatrix(),
                  solidFlag, domain, (T)0);
    setToConstant(*velocityNorm, voxelizedDomain.getVoxelMatrix(),
                  solidBorderFlag, domain, (T)0);
    return velocityNorm;
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeVelocityNorm()
{
    return computeVelocityNorm(lattice.getBoundingBox());
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeVelocityComponent (
            Box2D domain, plint iComp )
{
    std::auto_ptr<MultiScalarField2D<T> > velocityComponent (
            plb::computeVelocityComponent(lattice,domain, iComp));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*velocityComponent, voxelizedDomain.getVoxelMatrix(),
                  solidFlag, domain, (T)0);
    setToConstant(*velocityComponent, voxelizedDomain.getVoxelMatrix(),
                  solidBorderFlag, domain, (T)0);
    return velocityComponent;
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeVelocityComponent(plint iComp)
{
    return computeVelocityComponent(lattice.getBoundingBox(), iComp);
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computePressure(Box2D domain)
{
    std::auto_ptr<MultiScalarField2D<T> > pressure(plb::computeDensity(lattice,domain));
    T averageDensity = computeAverageDensity(domain);
    subtractInPlace(*pressure, averageDensity, domain);
    multiplyInPlace(*pressure, Descriptor<T>::cs2, domain);
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*pressure, voxelizedDomain.getVoxelMatrix(),
                  solidFlag, domain, (T)0);
    setToConstant(*pressure, voxelizedDomain.getVoxelMatrix(),
                  solidBorderFlag, domain, (T)0);
    return pressure;
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computePressure()
{
    return computePressure(lattice.getBoundingBox());
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeDensity(Box2D domain, T solidDensity)
{
    std::auto_ptr<MultiScalarField2D<T> > density(plb::computeDensity(lattice,domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*density, voxelizedDomain.getVoxelMatrix(),
                  solidFlag, domain, solidDensity);
    setToConstant(*density, voxelizedDomain.getVoxelMatrix(),
                  solidBorderFlag, domain, solidDensity);
    return density;
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeDensity(T solidDensity)
{
    return computeDensity(lattice.getBoundingBox(), solidDensity);
}

    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeStrainRateNorm(Box2D domain)
{
    std::auto_ptr<MultiScalarField2D<T> >
         strainRateNorm(computeSymmetricTensorNorm(*plb::computeStrainRateFromStress(lattice,domain)));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*strainRateNorm, voxelizedDomain.getVoxelMatrix(),
                  solidFlag, domain, (T)0);
    setToConstant(*strainRateNorm, voxelizedDomain.getVoxelMatrix(),
                  solidBorderFlag, domain, (T)0);
    return strainRateNorm;
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeStrainRateNorm()
{
    return computeStrainRateNorm(lattice.getBoundingBox());
}


template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeStrainRate(Box2D domain)
{
    std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
         strainRate(plb::computeStrainRateFromStress(lattice,domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    Array<T,SymmetricTensor<T,Descriptor>::n> zeros; zeros.resetToZero();
    setToConstant<T,SymmetricTensor<T,Descriptor>::n>(*strainRate, voxelizedDomain.getVoxelMatrix(),
                  solidFlag, domain, zeros);
    setToConstant<T,SymmetricTensor<T,Descriptor>::n>(*strainRate, voxelizedDomain.getVoxelMatrix(),
                  solidBorderFlag, domain, zeros);
    return strainRate;
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeStrainRate()
{
    return computeStrainRate(lattice.getBoundingBox());
}


template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeShearStressNorm(Box2D domain)
{
    std::auto_ptr<MultiScalarField2D<T> >
         shearStressNorm(computeSymmetricTensorNorm(*plb::computeShearStress(lattice,domain)));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*shearStressNorm, voxelizedDomain.getVoxelMatrix(),
                  solidFlag, domain, (T)0);
    setToConstant(*shearStressNorm, voxelizedDomain.getVoxelMatrix(),
                  solidBorderFlag, domain, (T)0);
    return shearStressNorm;
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiScalarField2D<T> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeShearStressNorm()
{
    return computeShearStressNorm(lattice.getBoundingBox());
}


template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeShearStress(Box2D domain)
{
    std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
         shearStress(plb::computeShearStress(lattice,domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    Array<T,SymmetricTensor<T,Descriptor>::n> zeros; zeros.resetToZero();
    setToConstant<T,SymmetricTensor<T,Descriptor>::n>(*shearStress, voxelizedDomain.getVoxelMatrix(),
                  solidFlag, domain, zeros);
    setToConstant<T,SymmetricTensor<T,Descriptor>::n>(*shearStress, voxelizedDomain.getVoxelMatrix(),
                  solidBorderFlag, domain, zeros);
    return shearStress;
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeShearStress()
{
    return computeShearStress(lattice.getBoundingBox());
}


template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
T OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeAverageVelocityComponent(Box2D domain, plint iComponent)
{
    std::auto_ptr<MultiScalarField2D<T> > density (
            plb::computeVelocityComponent(lattice,domain, iComponent) );
    MultiScalarField2D<int> flagMatrix((MultiBlock2D&)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  fluidBorderFlag, domain, 1);
    return computeAverage(*density, flagMatrix, 1, domain);
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
T OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeAverageDensity()
{
    return computeAverageDensity(voxelizedDomain.getVoxelMatrix().getBoundingBox());
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
T OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeAverageDensity(Box2D domain)
{
    std::auto_ptr<MultiScalarField2D<T> > density (
            plb::computeDensity(lattice,domain) );
    std::auto_ptr<MultiScalarField2D<T> > density2 (
            plb::computeDensity(lattice,domain) );
    MultiScalarField2D<int> flagMatrix((MultiBlock2D&)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  fluidBorderFlag, domain, 1);
    return computeAverage(*density, flagMatrix, 1, domain);
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
T OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeAverageEnergy()
{
    return computeAverageEnergy(voxelizedDomain.getVoxelMatrix().getBoundingBox());
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
T OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeAverageEnergy(Box2D domain)
{
    std::auto_ptr<MultiScalarField2D<T> > energy (
            plb::computeKineticEnergy(lattice,domain) );
    MultiScalarField2D<int> flagMatrix((MultiBlock2D&)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  fluidBorderFlag, domain, 1);
    return computeAverage(*energy, flagMatrix, 1, domain);
}
    
template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
T OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeRMSvorticity()
{
    return computeRMSvorticity(voxelizedDomain.getVoxelMatrix().getBoundingBox());
}

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
T OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType>::computeRMSvorticity(Box2D domain)
{
    std::auto_ptr<MultiScalarField2D<T> > vorticityNormSqr (
            plb::computeNormSqr(*plb::computeBulkVorticity(*plb::computeVelocity(lattice,domain)) ) );
    MultiScalarField2D<int> flagMatrix((MultiBlock2D&)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  fluidBorderFlag, domain, 1);
    return sqrt(computeAverage(*vorticityNormSqr, flagMatrix, 1, domain));
}
    
}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_CONDITION_2D_HH
