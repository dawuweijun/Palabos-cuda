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

#ifndef OFF_LATTICE_MODEL_2D_H
#define OFF_LATTICE_MODEL_2D_H

#include "core/globalDefs.h"
#include "offLattice/boundaryShapes2D.h"
namespace plb {

template< typename T, class SurfaceData >
class OffLatticeModel2D {
public:
    OffLatticeModel2D(BoundaryShape2D<T,SurfaceData>* shape_, int flowType_);
    OffLatticeModel2D(OffLatticeModel2D<T,SurfaceData> const& rhs);
    OffLatticeModel2D<T,SurfaceData>& operator= (
            OffLatticeModel2D<T,SurfaceData> const& rhs );
    virtual ~OffLatticeModel2D();
    void provideShapeArguments(std::vector<AtomicBlock2D*> args);
    plint getTag(plint id) const;
    bool pointOnSurface (
            Dot2D const& fromPoint, Dot2D const& direction,
            Array<T,2>& locatedPoint, T& distance,
            Array<T,2>& wallNormal, SurfaceData& surfaceData,
            OffBoundary::Type& bdType, plint& id ) const;
    Array<T,2> computeContinuousNormal (
            Array<T,2> const& p, plint id, bool isAreaWeighted = false ) const;
    bool intersectsSurface (
            Dot2D const& p1, Dot2D const& p2, plint& id ) const;
    bool isFluid(Dot2D const& location) const;
    bool velIsJ() const { return velIsJflag; }
    void setVelIsJ(bool velIsJflag_) { velIsJflag = velIsJflag_; }
    bool getPartialReplace() const { return partialReplaceFlag; }
    void setPartialReplace(bool prFlag) { partialReplaceFlag = prFlag; }
    virtual OffLatticeModel2D<T,SurfaceData>* clone() const =0;
    virtual plint getNumNeighbors() const =0;
    virtual void prepareCell (
            Dot2D const& cellLocation, AtomicContainerBlock2D& container ) =0;
    virtual void boundaryCompletion (
            AtomicBlock2D& lattice,
            AtomicContainerBlock2D& container,
            std::vector<AtomicBlock2D const*> const& args ) =0;
    virtual ContainerBlockData* generateOffLatticeInfo() const =0;
    virtual Array<T,2> getLocalForce(AtomicContainerBlock2D& container) const =0;
private:
    BoundaryShape2D<T,SurfaceData>* shape;
    int flowType;
    bool velIsJflag;
    bool partialReplaceFlag;
};

/// Precompute a list of nodes which are close to the off-lattice boundary
///   (both wet and dry ones).
template<typename T, class SurfaceData>
class OffLatticePatternFunctional2D : public BoxProcessingFunctional2D
{
public:
    /// numNeighbors is the number of neighbors a boundary node must be able
    ///   to access along a given direction.
    OffLatticePatternFunctional2D (
            OffLatticeModel2D<T,SurfaceData>* offLatticeModel_ );
    virtual ~OffLatticePatternFunctional2D();
    OffLatticePatternFunctional2D(OffLatticePatternFunctional2D const& rhs);
    OffLatticePatternFunctional2D& operator= (
            OffLatticePatternFunctional2D const& rhs );
    void swap(OffLatticePatternFunctional2D& rhs);
    virtual OffLatticePatternFunctional2D<T,SurfaceData>* clone() const;

    /// First AtomicBlock: OffLatticeInfo.
    ///   If there are more atomic-blocks then they are forwarded to the
    ///   shape function, to provide additional read-only parameters.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> fields);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    OffLatticeModel2D<T,SurfaceData>* offLatticeModel;
};

template<typename T, template<typename U> class Descriptor, class SurfaceData>
class OffLatticeCompletionFunctional2D : public BoxProcessingFunctional2D
{
public:
    OffLatticeCompletionFunctional2D (
            OffLatticeModel2D<T,SurfaceData>* offLatticeModel_,
            plint numShapeArgs_, plint numCompletionArgs_ );
    virtual ~OffLatticeCompletionFunctional2D();
    OffLatticeCompletionFunctional2D(OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData> const& rhs);
    OffLatticeCompletionFunctional2D& operator= (
            OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData> const& rhs );
    void swap(OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>& rhs);

    /// First AtomicBlock: Lattice; Second AtomicBlock: Container for off-
    ///   lattice info.
    ///   If there are more atomic-blocks then they are forwarded to the
    ///   shape function, to provide additional read-only parameters.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> fields);
    virtual OffLatticeCompletionFunctional2D<T,Descriptor,SurfaceData>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    OffLatticeModel2D<T,SurfaceData>* offLatticeModel;
    plint numShapeArgs, numCompletionArgs;
};

template< typename T, class SurfaceData >
class GetForceOnObjectFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    GetForceOnObjectFunctional2D (
            OffLatticeModel2D<T,SurfaceData>* offLatticeModel_ );
    virtual ~GetForceOnObjectFunctional2D();
    GetForceOnObjectFunctional2D(GetForceOnObjectFunctional2D<T,SurfaceData> const& rhs);
    GetForceOnObjectFunctional2D<T,SurfaceData>& operator= (
            GetForceOnObjectFunctional2D<T,SurfaceData> const& rhs );

    /// First AtomicBlock: Container for off-lattice info.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> fields);
    virtual GetForceOnObjectFunctional2D<T,SurfaceData>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    Array<T,2> getForce() const;
private:
    OffLatticeModel2D<T,Array<T,2> >* offLatticeModel;
    Array<plint,2> forceId;
};

}  // namespace plb

#endif  // OFF_LATTICE_MODEL_2D_H
