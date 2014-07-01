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

#ifndef SEGMENT_HASH_H
#define SEGMENT_HASH_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataField2D.h"
#include "particles/particleField2D.h"
#include "offLattice/segmentPolygonMesh2D.h"
namespace plb {

template<typename T>
class SegmentHash {
public:
    SegmentHash(AtomicContainerBlock2D& hashContainer);
    void assignSegments(SegmentPolygonMesh2D<T> const& mesh);
    void bruteReAssignSegments(SegmentPolygonMesh2D<T> const& mesh);
    template<class ParticleFieldT>
    void reAssignSegments (
            SegmentPolygonMesh2D<T> const& mesh, ParticleFieldT& particles,
            std::vector<plint> const& nonParallelVertices );
    void getSegments (
            const Array< T, 2  >& xRange, const Array< T, 2  >& yRange, std::vector< plint >& foundSegments) const;
    void getSegments (
            Box2D const& domain,
            std::vector<plint>& foundSegments ) const;
private:
    ScalarField2D<std::vector<plint> >& segments;
    std::vector<Dot2D>& assignedPositions;
};

template<typename T>
class CreateSegmentHash : public BoxProcessingFunctional2D {
public:
    CreateSegmentHash (
            SegmentPolygonMesh2D<T> const& mesh_ );
    // Field 0: Hash.
    virtual void processGenericBlocks (
                Box2D domain, std::vector<AtomicBlock2D*> fields );
    virtual CreateSegmentHash<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    SegmentPolygonMesh2D<T> const& mesh;
};

template<typename T, class ParticleFieldT>
class ReAssignSegmentHash : public BoxProcessingFunctional2D {
public:
    // The segment hash includes the possibilities for segments
    //   which are too big for proper parallelization, such as, inlet
    //   or outlet segments. The vertices of these segments must be
    //   non-parallel, which means, well defined on each processor. The
    //   argument nonParallelVertices must enumerate at least one vertex
    //   of each non-parallel segment. It is a requirement that all
    //   vertices of all segments connected to the vertices listed in
    //   nonParallelVertices are non-parallel, i.e. define on each
    //   processor.
    ReAssignSegmentHash (
            SegmentPolygonMesh2D<T> const& mesh_,
            std::vector<plint> const& nonParallelVertices_ );
    // Field 0: Hash; Field 1: Particles.
    virtual void processGenericBlocks (
                Box2D domain, std::vector<AtomicBlock2D*> fields );
    virtual ReAssignSegmentHash<T,ParticleFieldT>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    SegmentPolygonMesh2D<T> const& mesh;
    std::vector<plint> nonParallelVertices;
};

template<typename T>
class BruteReAssignSegmentHash : public BoxProcessingFunctional2D {
public:
    BruteReAssignSegmentHash (
            SegmentPolygonMesh2D<T> const& mesh_ );
    // Field 0: Hash.
    virtual void processGenericBlocks (
                Box2D domain, std::vector<AtomicBlock2D*> fields );
    virtual BruteReAssignSegmentHash<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    SegmentPolygonMesh2D<T> const& mesh;
};


}  // namespace plb

#endif  // TRIANGLE_HASH_H

