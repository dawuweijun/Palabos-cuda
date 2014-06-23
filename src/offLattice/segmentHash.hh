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

#ifndef TRIANGLE_HASH_HH
#define TRIANGLE_HASH_HH

#include "core/globalDefs.h"
#include "offLattice/segmentHash.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include <algorithm>

namespace plb {

struct SegmentHashData : public ContainerBlockData {
    SegmentHashData (
            plint nx, plint ny,
            Dot2D const& location )
        : segments(nx,ny)
    {
        segments.setLocation(location);
    }
    virtual SegmentHashData* clone() const {
        return new SegmentHashData(*this);
    }
    ScalarField2D<std::vector<plint> > segments;
    std::vector<Dot2D> assignedPositions;
};


/* ******** class SegmentHash ********************************************* */

template<typename T>
SegmentHash<T>::SegmentHash(AtomicContainerBlock2D& hashContainer)
    : segments (
        dynamic_cast<SegmentHashData*>(hashContainer.getData())->segments ),
      assignedPositions (
        dynamic_cast<SegmentHashData*>(hashContainer.getData())->assignedPositions )
{ }

template<typename T>
void SegmentHash<T>::getSegments (
                Array<T,2> const& xRange,
                Array<T,2> const& yRange,
                std::vector<plint>& foundSegments ) const
{
    // Fit onto the grid by making it bigger, to be sure the segment
    //   is never missed through round-off errors.
    Box2D discreteRange (
            (plint)xRange[0], (plint)xRange[1]+1,
            (plint)yRange[0], (plint)yRange[1]+1 );
    getSegments(discreteRange, foundSegments);
}

template<typename T>
void SegmentHash<T>::getSegments (
        Box2D const& domain,
        std::vector<plint>& foundSegments ) const
{
    Dot2D location(segments.getLocation());
    // Convert to local coordinates.
    Box2D shifted(domain.shift(-location.x,-location.y));
    foundSegments.clear();
    Box2D inters;
    if (intersect(shifted, segments.getBoundingBox(), inters)) {
        for (plint iX=inters.x0; iX<=inters.x1; ++iX) {
            for (plint iY=inters.y0; iY<=inters.y1; ++iY) {
                    std::vector<plint> const& newSegments = segments.get(iX,iY);
                    foundSegments.insert(foundSegments.end(),
                                          newSegments.begin(), newSegments.end());
            }
        }
        std::sort(foundSegments.begin(), foundSegments.end());
        foundSegments.erase( unique(foundSegments.begin(), foundSegments.end()),
                              foundSegments.end() );
    }
}

template<typename T>
void SegmentHash<T>::assignSegments (
        SegmentPolygonMesh2D<T> const& mesh )
{
    Dot2D location(segments.getLocation());
    assignedPositions.clear();
    for (plint iSegment=0; iSegment<mesh.getNumSegments(); ++iSegment)
    {
        Array<T,2> const& vertex0 = mesh.getVertex(iSegment, 0);
        Array<T,2> const& vertex1 = mesh.getVertex(iSegment, 1);

        Array<T,2> xRange (
                     std::min(vertex0[0], vertex1[0]),
                     std::max(vertex0[0], vertex1[0]) );
        Array<T,2> yRange (
                     std::min(vertex0[1], vertex1[1]),
                     std::max(vertex0[1], vertex1[1]) );

        // Fit onto the grid by making it bigger, to be sure the segment
        //   is never missed through round-off errors.
        Box2D discreteRange (
                (plint)xRange[0], (plint)xRange[1]+1,
                (plint)yRange[0], (plint)yRange[1]+1 );
        // Convert to local coordinates.
        discreteRange = discreteRange.shift (
                -location.x, -location.y );
        Box2D inters;
        if (intersect(discreteRange, segments.getBoundingBox(), inters))
        {
            for (plint iX=inters.x0; iX<=inters.x1; ++iX) {
                for (plint iY=inters.y0; iY<=inters.y1; ++iY) {
                        if (segments.get(iX,iY).empty()) {
                            assignedPositions.push_back(Dot2D(iX,iY));
                        }
                        segments.get(iX,iY).push_back(iSegment);
                }
            }
        }
    }
}

template<typename T>
template<class ParticleFieldT>
void SegmentHash<T>::reAssignSegments (
        SegmentPolygonMesh2D<T> const& mesh,
        ParticleFieldT& particles,
        std::vector<plint> const& nonParallelVertices )
{
    // Create domain from which particles are going to be retrieved.
    Box2D domain(segments.getBoundingBox());
    // First of all, remove old segments.
    /*
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                segments.get(iX,iY,iZ).clear();
            }
        }
    }
    */
    for (pluint iAssigned=0; iAssigned<assignedPositions.size(); ++iAssigned) {
        Dot2D pos(assignedPositions[iAssigned]);
        segments.get(pos.x,pos.y).clear();
    }
    assignedPositions.clear();
    // Particles have a bigger envelope than segments (which have
    //   envelope with 2 for Guo for example). The domain must now
    //   be translated into the local coordinates of the particles.
    Dot2D offset = computeRelativeDisplacement(segments, particles);
    domain = domain.shift(offset.x,offset.y);
    // Enlarge by one cell, because segments belong to the hash of
    //   a given AtomicBlock even when one of their vertices is out-
    //   side the AtomicBlock by maximally one cell.
    domain.enlarge(1);
    std::vector<typename ParticleFieldT::ParticleT*> found;
    particles.findParticles(domain, found);
    std::set<plint> segmentIds;
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        plint vertexId = found[iParticle]->getTag();
        std::vector<plint> newSegments (
                mesh.getNeighborSegmentIds(vertexId) );
        segmentIds.insert(newSegments.begin(), newSegments.end());
    }
    for (pluint iVertex=0; iVertex<nonParallelVertices.size(); ++iVertex) {
        plint vertexId = nonParallelVertices[iVertex];
        std::vector<plint> newSegments (
                mesh.getNeighborSegmentIds(vertexId) );
        segmentIds.insert(newSegments.begin(), newSegments.end());
    }

    Dot2D location(segments.getLocation());
    std::set<plint>::const_iterator it = segmentIds.begin();
    for (; it != segmentIds.end(); ++it) {
        plint iSegment = *it;
        Array<T,2> const& vertex0 = mesh.getVertex(iSegment, 0);
        Array<T,2> const& vertex1 = mesh.getVertex(iSegment, 1);

        Array<T,2> xRange (
                     std::min(vertex0[0], vertex1[0]),
                     std::max(vertex0[0], vertex1[0]) );
        Array<T,2> yRange (
                     std::min(vertex0[1], vertex1[1]),
                     std::max(vertex0[1], vertex1[1]) );

        // Fit onto the grid by making it bigger, to be sure the segment
        //   is never missed through round-off errors.
        Box2D discreteRange (
                (plint)xRange[0], (plint)xRange[1]+1,
                (plint)yRange[0], (plint)yRange[1]+1 );
        // Convert to local coordinates.
        discreteRange = discreteRange.shift (
                -location.x, -location.y );
        Box2D inters;
        if (intersect(discreteRange, segments.getBoundingBox(), inters))
        {
            for (plint iX=inters.x0; iX<=inters.x1; ++iX) {
                for (plint iY=inters.y0; iY<=inters.y1; ++iY) {
                    if (segments.get(iX,iY).empty()) {
                        assignedPositions.push_back(Dot2D(iX,iY));
                    }
                    segments.get(iX,iY).push_back(iSegment);
                }
            }
        }
    }
}

template<typename T>
void SegmentHash<T>::bruteReAssignSegments (
        SegmentPolygonMesh2D<T> const& mesh )
{
    PLB_ASSERT(false);
    // First of all, remove old segments.
    Box2D domain(segments.getBoundingBox());
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                segments.get(iX,iY).clear();
        }
    }
    Dot2D location(segments.getLocation());
    for (plint iSegment=0; iSegment<mesh.getNumSegments(); ++iSegment)
    {
        if (mesh.isValidVertex(iSegment,0) &&
            mesh.isValidVertex(iSegment,1) &&
            mesh.isValidVertex(iSegment,2) )
        {
            Array<T,2> const& vertex0 = mesh.getVertex(iSegment, 0);
            Array<T,2> const& vertex1 = mesh.getVertex(iSegment, 1);

            Array<T,2> xRange (
                         std::min(vertex0[0], vertex1[0]),
                         std::max(vertex0[0], vertex1[0]) );
            Array<T,2> yRange (
                         std::min(vertex0[1], vertex1[1]),
                         std::max(vertex0[1], vertex1[1]) );

            // Fit onto the grid by making it bigger, to be sure the segment
            //   is never missed through round-off errors.
            Box2D discreteRange (
                    (plint)xRange[0], (plint)xRange[1]+1,
                    (plint)yRange[0], (plint)yRange[1]+1);
            // Convert to local coordinates.
            discreteRange = discreteRange.shift (
                    -location.x, -location.y);
            Box2D inters;
            if (intersect(discreteRange, segments.getBoundingBox(), inters))
            {
                for (plint iX=inters.x0; iX<=inters.x1; ++iX) {
                    for (plint iY=inters.y0; iY<=inters.y1; ++iY) {
                        segments.get(iX,iY).push_back(iSegment);
                    }
                }
            }
        }
    }
}

/* ******** CreateSegmentHash ************************************ */

template<typename T>
CreateSegmentHash<T>::CreateSegmentHash (
        SegmentPolygonMesh2D<T> const& mesh_ )
    :  mesh(mesh_)
{ }

template<typename T>
void CreateSegmentHash<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock2D* container =
        dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );
    SegmentHashData* hashData
        = new SegmentHashData (
                container->getNx(), container->getNy(),
                container->getLocation() );
    container->setData(hashData);
    SegmentHash<T>(*container).assignSegments(mesh);
}

template<typename T>
CreateSegmentHash<T>* CreateSegmentHash<T>::clone() const {
    return new CreateSegmentHash<T>(*this);
}

template<typename T>
void CreateSegmentHash<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with hash data.
}

template<typename T>
BlockDomain::DomainT CreateSegmentHash<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** ReAssignSegmentHash ************************************ */

template<typename T, class ParticleFieldT>
ReAssignSegmentHash<T,ParticleFieldT>::ReAssignSegmentHash (
        SegmentPolygonMesh2D<T> const& mesh_,
        std::vector<plint> const& nonParallelVertices_ )
    :  mesh(mesh_),
       nonParallelVertices(nonParallelVertices_)
{ }

template<typename T, class ParticleFieldT>
void ReAssignSegmentHash<T,ParticleFieldT>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    AtomicContainerBlock2D* container =
        dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );
    ParticleFieldT* particles =
        dynamic_cast<ParticleFieldT*>(blocks[1]);
    PLB_ASSERT( particles );
    SegmentHash<T>(*container).reAssignSegments (
            mesh,*particles,nonParallelVertices);
}

template<typename T, class ParticleFieldT>
ReAssignSegmentHash<T,ParticleFieldT>* ReAssignSegmentHash<T,ParticleFieldT>::clone() const {
    return new ReAssignSegmentHash<T,ParticleFieldT>(*this);
}

template<typename T, class ParticleFieldT>
void ReAssignSegmentHash<T,ParticleFieldT>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with hash data.
    modified[1] = modif::nothing; // Vertex-Particles.
}

template<typename T, class ParticleFieldT>
BlockDomain::DomainT ReAssignSegmentHash<T,ParticleFieldT>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** BruteReAssignSegmentHash ************************************ */

template<typename T>
BruteReAssignSegmentHash<T>::BruteReAssignSegmentHash (
        SegmentPolygonMesh2D<T> const& mesh_ )
    :  mesh(mesh_)
{ }

template<typename T>
void BruteReAssignSegmentHash<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock2D* container =
        dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );
    SegmentHash<T>(*container).bruteReAssignSegments(mesh);
}

template<typename T>
BruteReAssignSegmentHash<T>* BruteReAssignSegmentHash<T>::clone() const {
    return new BruteReAssignSegmentHash<T>(*this);
}

template<typename T>
void BruteReAssignSegmentHash<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with hash data.
    modified[1] = modif::nothing; // Vertex-Particles.
}

template<typename T>
BlockDomain::DomainT BruteReAssignSegmentHash<T>::appliesTo() const {
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // TRIANGLE_HASH_HH

