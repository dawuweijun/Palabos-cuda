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

#ifndef OFF_LATTICE_MODEL_2D_HH
#define OFF_LATTICE_MODEL_2D_HH

#include "offLattice/offLatticeModel2D.h"
#include "offLattice/voxelizer2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include <algorithm>
#include <cmath>

namespace plb {

template<typename T, class SurfaceData>
OffLatticeModel2D<T,SurfaceData>::OffLatticeModel2D (
        BoundaryShape2D<T,SurfaceData>* shape_, int flowType_ )
    : shape(shape_),
      flowType(flowType_),
      velIsJflag(false),
      partialReplaceFlag(false)
{
    PLB_ASSERT(flowType==voxelFlag::inside || flowType==voxelFlag::outside);
}

template<typename T, class SurfaceData>
OffLatticeModel2D<T,SurfaceData>::OffLatticeModel2D (
        OffLatticeModel2D<T,SurfaceData> const& rhs )
    : shape(rhs.shape->clone()),
      flowType(rhs.flowType),
      velIsJflag(rhs.velIsJflag),
      partialReplaceFlag(rhs.partialReplaceFlag)
{ }

template<typename T, class SurfaceData>
OffLatticeModel2D<T,SurfaceData>& OffLatticeModel2D<T,SurfaceData>::operator= (
        OffLatticeModel2D<T,SurfaceData> const& rhs )
{
    delete shape;
    shape = rhs.shape->clone();
    flowType = rhs.flowType;
    velIsJflag = rhs.velIsJflag;
    partialReplaceFlag = rhs.partialReplaceFlag;
    return *this;
}

template<typename T, class SurfaceData>
OffLatticeModel2D<T,SurfaceData>::~OffLatticeModel2D() {
    delete shape;
}

template<typename T, class SurfaceData>
void OffLatticeModel2D<T,SurfaceData>::provideShapeArguments (
            std::vector<AtomicBlock2D*> args )
{
    BoundaryShape2D<T,SurfaceData>* newShape=shape->clone(args);
    std::swap(shape, newShape);
    delete newShape;
}

template<typename T, class SurfaceData>
plint OffLatticeModel2D<T,SurfaceData>::getTag(plint id) const {
    return shape->getTag(id);
}

template<typename T, class SurfaceData>
bool OffLatticeModel2D<T,SurfaceData>::pointOnSurface (
        Dot2D const& fromPoint, Dot2D const& direction,
        Array<T,2>& locatedPoint, T& distance,
        Array<T,2>& wallNormal, SurfaceData& surfaceData,
        OffBoundary::Type& bdType, plint& id ) const
{
    return shape->gridPointOnSurface (
               fromPoint, direction, locatedPoint, distance, wallNormal, surfaceData, bdType, id );
}

template<typename T, class SurfaceData>
Array<T,2> OffLatticeModel2D<T,SurfaceData>::computeContinuousNormal (
            Array<T,2> const& p, plint id, bool isAreaWeighted ) const
{
    return shape->computeContinuousNormal(p, id, isAreaWeighted);
}

template<typename T, class SurfaceData>
bool OffLatticeModel2D<T,SurfaceData>::intersectsSurface (
        Dot2D const& p1, Dot2D const& p2, plint& id ) const
{
    return shape->intersectsSurface(p1, p2, id);
}

template<typename T, class SurfaceData>
bool OffLatticeModel2D<T,SurfaceData>::isFluid(Dot2D const& location) const
{
    if (flowType==voxelFlag::inside) {
        return shape->isInside(location);
    }
    else {
        return !shape->isInside(location);
    }
}

template<typename T, template<typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::OffLatticeCompletionFunctional2D (
        OffLatticeModel2D<T,SurfaceData>* offLatticeModel_,
        plint numShapeArgs_, plint numCompletionArgs_ )
  : offLatticeModel(offLatticeModel_),
    numShapeArgs(numShapeArgs_),
    numCompletionArgs(numCompletionArgs_)
{ }

template<typename T, template<typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::~OffLatticeCompletionFunctional2D()
{
    delete offLatticeModel;
}

template<typename T, template<typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::
    OffLatticeCompletionFunctional2D (
            OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData> const& rhs)
    : offLatticeModel(rhs.offLatticeModel->clone()),
      numShapeArgs(rhs.numShapeArgs),
      numCompletionArgs(rhs.numCompletionArgs)
{ }

template<typename T, template<typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>&
    OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::operator= (
            OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData> const& rhs )
{
    OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor, class SurfaceData>
void OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::swap(
        OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>& rhs)
{
    std::swap(offLatticeModel, rhs.offLatticeModel);
    std::swap(numShapeArgs, rhs.numShapeArgs);
    std::swap(numCompletionArgs, rhs.numCompletionArgs);
}


template<typename T, template<typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>*
    OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::clone() const
{
    return new OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>(*this);
}

template<typename T, template<typename U> class Descriptor, class SurfaceData>
void OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    PLB_ASSERT( (plint)modified.size() == 2+numShapeArgs+numCompletionArgs );
    modified[0] = modif::staticVariables;  // Lattice.
    modified[1] = modif::nothing; // Container for wet/dry nodes.
    // Possible additional parameters for the shape function and
    //   for the completion algorithm are read-only.
    for (pluint i=2; i<modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template<typename T, template<typename U> class Descriptor, class SurfaceData>
BlockDomain::DomainT OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor, class SurfaceData>
void OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> fields )
{
    PLB_PRECONDITION( (plint)fields.size() == 2+numShapeArgs+numCompletionArgs );
    AtomicBlock2D* lattice = fields[0];
    PLB_ASSERT( lattice );

    AtomicContainerBlock2D* container =  // Container for wet/dry nodes.
        dynamic_cast<AtomicContainerBlock2D*>(fields[1]);
    PLB_ASSERT( container );

    if (numShapeArgs>0) {
        std::vector<AtomicBlock2D*> shapeParameters(numShapeArgs);
        for (plint i=0; i<numShapeArgs; ++i) {
            shapeParameters[i] = fields[i+2];
        }
        offLatticeModel->provideShapeArguments(shapeParameters);
    }
    std::vector<AtomicBlock2D const*> completionParameters(numCompletionArgs);
    for (plint i=0; i<numCompletionArgs; ++i) {
        completionParameters[i] = fields[i+2+numShapeArgs];
    }

    offLatticeModel->boundaryCompletion(*lattice, *container, completionParameters);
}


template<typename T, class SurfaceData>
OffLatticePatternFunctional2D<T,SurfaceData>::
    OffLatticePatternFunctional2D (
            OffLatticeModel2D<T,SurfaceData>* offLatticeModel_ )
  : offLatticeModel(offLatticeModel_)
{ }

template<typename T, class SurfaceData>
OffLatticePatternFunctional2D<T,SurfaceData>::~OffLatticePatternFunctional2D()
{
    delete offLatticeModel;
}

template<typename T, class SurfaceData>
OffLatticePatternFunctional2D<T,SurfaceData>::
    OffLatticePatternFunctional2D (
            OffLatticePatternFunctional2D<T,SurfaceData> const& rhs)
    : offLatticeModel(rhs.offLatticeModel->clone())
{ }

template<typename T, class SurfaceData>
OffLatticePatternFunctional2D<T,SurfaceData>&
    OffLatticePatternFunctional2D<T,SurfaceData>::operator= (
            OffLatticePatternFunctional2D<T,SurfaceData> const& rhs )
{
    OffLatticePatternFunctional2D<T,SurfaceData>(rhs).swap(*this);
    return *this;
}

template<typename T, class SurfaceData>
void OffLatticePatternFunctional2D<T,SurfaceData>::swap(
        OffLatticePatternFunctional2D<T,SurfaceData>& rhs)
{
    std::swap(offLatticeModel, rhs.offLatticeModel);
}


template<typename T, class SurfaceData>
OffLatticePatternFunctional2D<T,SurfaceData>*
    OffLatticePatternFunctional2D<T,SurfaceData>::clone() const
{
    return new OffLatticePatternFunctional2D<T,SurfaceData>(*this);
}

template<typename T, class SurfaceData>
void OffLatticePatternFunctional2D<T,SurfaceData>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;  // Container.
    // Possible additional parameters for the shape function are read-only.
    for (pluint i=1; i<modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template<typename T, class SurfaceData>
BlockDomain::DomainT OffLatticePatternFunctional2D<T,SurfaceData>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, class SurfaceData>
void OffLatticePatternFunctional2D<T,SurfaceData>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> fields )
{
    PLB_PRECONDITION( fields.size() >= 1 );
    AtomicContainerBlock2D* container =
        dynamic_cast<AtomicContainerBlock2D*>(fields[0]);
    PLB_ASSERT( container );
    ContainerBlockData* storeInfo = 
        offLatticeModel->generateOffLatticeInfo();
    container->setData(storeInfo);

    if (fields.size()>1) {
        std::vector<AtomicBlock2D*> shapeParameters(fields.size()-1);
        for (pluint i=0; i<shapeParameters.size(); ++i) {
            shapeParameters[i] = fields[i+1];
        }
        offLatticeModel->provideShapeArguments(shapeParameters);
    }

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                offLatticeModel->prepareCell(Dot2D(iX,iY), *container);
        }
    }
}


template< typename T, class SurfaceData >
GetForceOnObjectFunctional2D<T,SurfaceData>::GetForceOnObjectFunctional2D (
    OffLatticeModel2D<T,SurfaceData>* offLatticeModel_ )
        : offLatticeModel(offLatticeModel_),
          forceId (
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum() )
{ }

template< typename T, class SurfaceData >
GetForceOnObjectFunctional2D<T,SurfaceData>::~GetForceOnObjectFunctional2D()
{
    delete offLatticeModel;
}

template< typename T, class SurfaceData >
GetForceOnObjectFunctional2D<T,SurfaceData>::GetForceOnObjectFunctional2D (
        GetForceOnObjectFunctional2D<T,SurfaceData> const& rhs )
    : PlainReductiveBoxProcessingFunctional2D(rhs),
      offLatticeModel(rhs.offLatticeModel->clone()),
      forceId(rhs.forceId)
{ }

template< typename T, class SurfaceData >
GetForceOnObjectFunctional2D<T,SurfaceData>&
    GetForceOnObjectFunctional2D<T,SurfaceData>::operator= (
            GetForceOnObjectFunctional2D<T,SurfaceData> const& rhs )
{
    delete offLatticeModel; offLatticeModel = rhs.offLatticeModel->clone();
    forceId = rhs.forceId;
    PlainReductiveBoxProcessingFunctional2D::operator=(rhs);
    return *this;
}

template< typename T, class SurfaceData >
void GetForceOnObjectFunctional2D<T,SurfaceData>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> fields )
{
    PLB_PRECONDITION( fields.size() == 1 );
    AtomicContainerBlock2D* offLatticeInfo =
        dynamic_cast<AtomicContainerBlock2D*>(fields[0]);
    PLB_ASSERT( offLatticeInfo );

    Array<T,2> force = offLatticeModel->getLocalForce(*offLatticeInfo);
    this->getStatistics().gatherSum(forceId[0], force[0]);
    this->getStatistics().gatherSum(forceId[1], force[1]);
}

template< typename T, class SurfaceData >
GetForceOnObjectFunctional2D<T,SurfaceData>*
    GetForceOnObjectFunctional2D<T,SurfaceData>::clone() const
{
     return new GetForceOnObjectFunctional2D<T,SurfaceData>(*this);
}

template< typename T, class SurfaceData >
void GetForceOnObjectFunctional2D<T,SurfaceData>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0]=modif::nothing;  // Off-lattice info.
}

template< typename T, class SurfaceData >
BlockDomain::DomainT GetForceOnObjectFunctional2D<T,SurfaceData>::appliesTo() const {
    return BlockDomain::bulk;
}

template< typename T, class SurfaceData >
Array<T,2> GetForceOnObjectFunctional2D<T,SurfaceData>::getForce() const {
    return Array<T,2> (
            this->getStatistics().getSum(forceId[0]),
            this->getStatistics().getSum(forceId[1]) );
}

}  // namespace plb

#endif  // OFF_LATTICE_MODEL_2D_HH
