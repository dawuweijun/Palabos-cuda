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

/* Main author: Orestis Malaspinas */

#ifndef VOXELIZER_2D_H
#define VOXELIZER_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "offLattice/segmentHash.h"
#include "offLattice/voxelizer.h"
namespace plb {

template<typename T>
std::auto_ptr<MultiScalarField2D<int> > voxelize2D (
        SegmentPolygonMesh2D<T> const& mesh,
        plint symmetricLayer, plint borderWidth );

template<typename T>
std::auto_ptr<MultiScalarField2D<int> > voxelize2D (
        SegmentPolygonMesh2D<T> const& mesh,
        Box2D const& domain, plint borderWidth );

template<typename T>
std::auto_ptr<MultiScalarField2D<int> > voxelize2D (
        SegmentPolygonMesh2D<T> const& mesh,
        Box2D const& domain, plint borderWidth, Box2D seed );

template<typename T>
std::auto_ptr<MultiScalarField2D<int> > revoxelize2D (
        SegmentPolygonMesh2D<T> const& mesh,
        MultiScalarField2D<int>& oldVoxelMatrix,
        MultiContainerBlock2D& hashContainer, plint borderWidth );

template<typename T>
class VoxelizeMeshFunctional2D : public BoxProcessingFunctional2D {
public:
    VoxelizeMeshFunctional2D (
            SegmentPolygonMesh2D<T> const& mesh_ );
    virtual void processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks );
    virtual VoxelizeMeshFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    bool checkIfFacetsCrossed (
            AtomicContainerBlock2D& hashContainer,
            Array<T,2> const& point1, Array<T,2> const& point2,
            T& distance, plint& whichTriangle );
    bool distanceToSurface (
            AtomicContainerBlock2D& hashContainer,
            Array<T,2> const& point, T& distance, bool& isBehind ) const;
    bool createVoxelizationRange (
            Box2D const& domain, ScalarField2D<int>& voxels,
            Array<plint,2>& xRange, Array<plint,2>& yRange );
    bool voxelizeFromNeighbor (
        ScalarField2D<int> const& voxels, AtomicContainerBlock2D& hashContainer,
        Dot2D pos, Dot2D neighbor, int& voxelType );
    void printOffender (
            ScalarField2D<int> const& voxels,
            AtomicContainerBlock2D& hashContainer, Dot2D pos );
private:
    SegmentPolygonMesh2D<T> const& mesh;
};

/// Convert inside flags to innerBoundary, and outside flags to outerBoundary,
///   within a layer of width "borderWidth".
template<typename T>
void detectBorderLine( MultiScalarField2D<T>& voxelMatrix,
                       Box2D const& domain, plint borderWidth );

template<typename T>
class DetectBorderLineFunctional2D : public BoxProcessingFunctional2D_S<T> {
public:
    DetectBorderLineFunctional2D(plint borderWidth_);
    virtual void process(Box2D domain, ScalarField2D<T>& voxels);
    virtual DetectBorderLineFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint borderWidth;
};

} // namespace plb

#endif  // VOXELIZER2D_H
