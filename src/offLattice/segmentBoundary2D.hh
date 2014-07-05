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

#ifndef INNER_FLOW_BOUNDARY_2D_HH
#define INNER_FLOW_BOUNDARY_2D_HH

#include "core/globalDefs.h"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "offLattice/segmentBoundary2D.h"
#include "offLattice/segmentPolygonMesh2D.h"
#include "offLattice/offLatticeBoundaryProfiles2D.h"
#include "offLattice/segmentToDef.h"
#include "multiBlock/nonLocalTransfer2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "offLattice/makeSparse2D.h"
#include "offLattice/voxelizer2D.h"
#include <cmath>
#include <limits>

namespace plb
{

/******** class BoundaryProfiles2D *****************************************/

template<typename T, class SurfaceData>
BoundaryProfiles2D<T,SurfaceData>::BoundaryProfiles2D()
{
    wallProfile = generateDefaultWallProfile2D<T,SurfaceData>();
    replaceProfile ( 0, wallProfile->clone() );
}

template<typename T, class SurfaceData>
BoundaryProfiles2D<T,SurfaceData>::~BoundaryProfiles2D()
{
    clearProfiles();
    delete wallProfile;
}

template<typename T, class SurfaceData>
BoundaryProfiles2D<T,SurfaceData>::BoundaryProfiles2D ( BoundaryProfiles2D<T,SurfaceData> const& rhs )
{
    wallProfile = rhs.wallProfile->clone();
    typename std::map<plint,BoundaryProfile2D<T,SurfaceData>*>::const_iterator it = rhs.profiles.begin();
    for ( ; it!=rhs.profiles.end(); ++it )
    {
        profiles.insert (
            std::pair<plint,BoundaryProfile2D<T,SurfaceData>*> (
                it->first,it->second->clone() ) );
    }
    inletOutletIds = rhs.inletOutletIds;
    lidNormal = rhs.lidNormal;
    lidCenter = rhs.lidCenter;
    lidRadius = rhs.lidRadius;
}

template<typename T, class SurfaceData>
BoundaryProfiles2D<T,SurfaceData>& BoundaryProfiles2D<T,SurfaceData>::operator= (
    BoundaryProfiles2D<T,SurfaceData> const& rhs )
{
    BoundaryProfiles2D<T,SurfaceData> ( rhs ).swap ( *this );
    return *this;
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::swap ( BoundaryProfiles2D<T,SurfaceData>& rhs )
{
    std::swap ( wallProfile, rhs.wallProfile );
    profiles.swap ( rhs.profiles );
    inletOutletIds.swap ( rhs.inletOutletIds );
    lidNormal.swap ( rhs.lidNormal );
    lidCenter.swap ( rhs.lidCenter );
    lidRadius.swap ( rhs.lidRadius );
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::setWallProfile ( BoundaryProfile2D<T,SurfaceData>* wallProfile_ )
{
    PLB_PRECONDITION ( wallProfile_ );
    delete wallProfile;
    wallProfile = wallProfile_;
    replaceProfile ( 0, wallProfile->clone() );
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::resetProfiles ( std::map<plint,BoundaryProfile2D<T,SurfaceData>*> profiles_ )
{
    PLB_ASSERT ( !profiles_.empty() );
    clearProfiles();
    profiles = profiles_;
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::defineInletOutletTags (
    SegmentBoundary2D<T> const& boundary, plint sortDirection )
{
    inletOutletIds = boundary.getInletOutletIds ( sortDirection );
    boundary.getLidProperties ( sortDirection, lidNormal, lidCenter, lidRadius );
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::adjustInletOutlet (
    SegmentBoundary2D<T> const& boundary, plint sortDirection )
{
    boundary.getLidProperties ( sortDirection, lidNormal, lidCenter, lidRadius );
    for ( pluint iProfile=0; iProfile<inletOutletIds.size(); ++iProfile )
    {
        plint id = inletOutletIds[iProfile];
        typename std::map<plint,BoundaryProfile2D<T,SurfaceData>*>::iterator it = profiles.find ( id );
        PLB_ASSERT ( it != profiles.end() );
        it->second->setNormal ( lidNormal[iProfile] );
        it->second->defineCircularShape ( lidCenter[iProfile], lidRadius[iProfile] );
    }
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::setInletOutlet (
    std::vector<BoundaryProfile2D<T,SurfaceData>*> inletOutlets )
{
    PLB_ASSERT ( inletOutletIds.size() == inletOutlets.size() );
    std::map<plint,BoundaryProfile2D<T,SurfaceData>*> newProfiles;
    newProfiles[0] = wallProfile->clone();
    for ( pluint iProfile=0; iProfile<inletOutlets.size(); ++iProfile )
    {
        plint id = inletOutletIds[iProfile];
#ifdef PLB_DEBUG
        typename std::map<plint,BoundaryProfile2D<T,SurfaceData>*>::const_iterator it = newProfiles.find ( id );
        PLB_ASSERT ( it==newProfiles.end() );
#endif
        newProfiles[id] = inletOutlets[iProfile];
        newProfiles[id]->setNormal ( lidNormal[iProfile] );
        newProfiles[id]->defineCircularShape ( lidCenter[iProfile], lidRadius[iProfile] );

    }
    resetProfiles ( newProfiles );
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::setInletOutlet (
    BoundaryProfile2D<T,SurfaceData>* profile1, BoundaryProfile2D<T,SurfaceData>* profile2 )
{
    std::vector<BoundaryProfile2D<T,SurfaceData>*> inletOutlets ( 2 );
    inletOutlets[0] = profile1;
    inletOutlets[1] = profile2;
    setInletOutlet ( inletOutlets );
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::setInletOutlet (
    BoundaryProfile2D<T,SurfaceData>* profile1, BoundaryProfile2D<T,SurfaceData>* profile2,
    BoundaryProfile2D<T,SurfaceData>* profile3 )
{
    std::vector<BoundaryProfile2D<T,SurfaceData>*> inletOutlets ( 3 );
    inletOutlets[0] = profile1;
    inletOutlets[1] = profile2;
    inletOutlets[2] = profile3;
    setInletOutlet ( inletOutlets );
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::setInletOutlet (
    BoundaryProfile2D<T,SurfaceData>* profile1, BoundaryProfile2D<T,SurfaceData>* profile2,
    BoundaryProfile2D<T,SurfaceData>* profile3, BoundaryProfile2D<T,SurfaceData>* profile4 )
{
    std::vector<BoundaryProfile2D<T,SurfaceData>*> inletOutlets ( 4 );
    inletOutlets[0] = profile1;
    inletOutlets[1] = profile2;
    inletOutlets[2] = profile3;
    inletOutlets[3] = profile4;
    setInletOutlet ( inletOutlets );
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::defineProfile (
    plint tag, BoundaryProfile2D<T,SurfaceData>* profile )
{
    typename std::map<plint,BoundaryProfile2D<T,SurfaceData>*>::iterator it = profiles.find ( tag );
    if ( it!=profiles.end() )
    {
        delete it->second;
    }
    profiles[tag] = profile;
}

template<typename T, class SurfaceData>
BoundaryProfile2D<T,SurfaceData> const& BoundaryProfiles2D<T,SurfaceData>::getProfile (
    SegmentBoundary2D<T> const& boundary, plint iSegment ) const
{
    PLB_ASSERT ( iSegment < ( plint ) boundary.getSegmentTags().size() );
    plint segmentTag = boundary.getSegmentTags() [iSegment];
    typename std::map<plint,BoundaryProfile2D<T,SurfaceData>*>::const_iterator it = profiles.find ( segmentTag );
    if ( it==profiles.end() )
    {
        PLB_ASSERT ( wallProfile );
        return *wallProfile;
    }
    else
    {
        PLB_ASSERT ( it->second );
        return * ( it->second );
    }
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::replaceProfile (
    plint id, BoundaryProfile2D<T,SurfaceData>* newProfile )
{
    typename std::map<plint,BoundaryProfile2D<T,SurfaceData>*>::iterator it = profiles.find ( id );
    if ( it!=profiles.end() )
    {
        delete it->second;
    }
    profiles[id] = newProfile;
}

template<typename T, class SurfaceData>
void BoundaryProfiles2D<T,SurfaceData>::clearProfiles()
{
    typename std::map<plint,BoundaryProfile2D<T,SurfaceData>*>::iterator it = profiles.begin();
    for ( ; it!=profiles.end(); ++it )
    {
        delete it->second;
    }
    profiles.clear();
}

/******** class DEFscaledMesh2D *****************************************/

template<typename T>
DEFscaledMesh2D<T>::DEFscaledMesh2D (
    SegmentSet<T> const& segmentSet_ )
    : margin ( 0 )
{
    initialize ( segmentSet_, 0, 0, Dot2D() );
}

template<typename T>
DEFscaledMesh2D<T>::DEFscaledMesh2D (
    SegmentSet<T> const& segmentSet_,
    plint resolution_, plint referenceDirection_,
    plint margin_, Dot2D location )
    : margin ( margin_ )
{
    initialize ( segmentSet_, resolution_, referenceDirection_, location );
}

template<typename T>
DEFscaledMesh2D<T>::DEFscaledMesh2D (
    SegmentSet<T> const& segmentSet_,
    plint resolution_, plint referenceDirection_,
    plint margin_, plint extraLayer )
    : margin ( margin_ )
{
    plint layer = margin+extraLayer;
    initialize ( segmentSet_, resolution_, referenceDirection_,
                 Dot2D ( layer,layer ) );
}

template<typename T>
DEFscaledMesh2D<T>::~DEFscaledMesh2D()
{
    delete mesh;
}

template<typename T>
void DEFscaledMesh2D<T>::initialize (
    SegmentSet<T> const& segmentSet_, plint resolution_,
    plint referenceDirection_, Dot2D location )
{
    T eps = getEpsilon<T> ( segmentSet_.getPrecision() );

    constructSurfaceMesh<T> (
        segmentSet_.getSegments(),
        vertexList, emanatingEdgeList, edgeList, eps );
    mesh = new SegmentPolygonMesh2D<T> ( vertexList, emanatingEdgeList, edgeList );

    if ( resolution_!=0 )
    {
        // Convert the mesh to lattice units.
        toLatticeUnits2D<T> (
            getMesh(), resolution_, referenceDirection_, physicalLocation, dx );

        Array<T,2> luOffset ( location.x,location.y );
        getMesh().translate ( luOffset );
        physicalLocation -= luOffset*dx;
    }
}

template<typename T>
DEFscaledMesh2D<T>::DEFscaledMesh2D (
    DEFscaledMesh2D<T> const& rhs )
    : vertexList ( rhs.vertexList ),
      emanatingEdgeList ( rhs.emanatingEdgeList ),
      edgeList ( rhs.edgeList ),
      margin ( rhs.margin )
{
    mesh = new SegmentPolygonMesh2D<T> ( vertexList, emanatingEdgeList, edgeList );
}

template<typename T>
DEFscaledMesh2D<T>&
DEFscaledMesh2D<T>::operator= ( DEFscaledMesh2D<T> const& rhs )
{
    DEFscaledMesh2D<T> ( rhs ).swap ( *this );
    return *this;
}

template<typename T>
void DEFscaledMesh2D<T>::swap ( DEFscaledMesh2D<T>& rhs )
{
    vertexList.swap ( rhs.vertexList );
    emanatingEdgeList.swap ( rhs.emanatingEdgeList );
    edgeList.swap ( rhs.edgeList );
    std::swap ( mesh,rhs.mesh );
    std::swap ( margin, rhs.margin );
    std::swap ( physicalLocation, rhs.physicalLocation );
    std::swap ( dx, rhs.dx );
}

template<typename T>
SegmentPolygonMesh2D<T> const& DEFscaledMesh2D<T>::getMesh() const
{
    return *mesh;
}

template<typename T>
SegmentPolygonMesh2D<T>& DEFscaledMesh2D<T>::getMesh()
{
    return *mesh;
}

template<typename T>
plint DEFscaledMesh2D<T>::getMargin() const
{
    return margin;
}


/******** class SegmentBoundary2D *****************************************/

template<typename T>
SegmentBoundary2D<T>::SegmentBoundary2D (
    DEFscaledMesh2D<T> const& defMesh, bool automaticCloseHoles )
    : currentTagNum ( 0 ),
      margin ( defMesh.getMargin() ),
      physicalLocation ( defMesh.getPhysicalLocation() ),
      dx ( defMesh.getDx() )
{
    topology.push ( 1 ); // By default, closed mesh.
    vertexSet.push ( 0 ); // By default, static mesh.

    vertexLists.reserve ( 3 ); // Vertex lists are expensive to copy. Better
    //   pre-allocate a slot for three of them.
    vertexLists.push_back ( defMesh.getVertexList() );
    emanatingEdgeLists[0] = defMesh.getEmanatingEdgeList();
    edgeLists[0] = defMesh.getEdgeList();

    emanatingEdgeLists[1] = emanatingEdgeLists[0];
    edgeLists[1] = edgeLists[0];

    meshes.push_back ( SegmentPolygonMesh2D<T> (vertexLists[0], emanatingEdgeLists[0], edgeLists[0] ) );
    meshes.push_back ( SegmentPolygonMesh2D<T> (vertexLists[0], emanatingEdgeLists[1], edgeLists[1] ) );

    // Prepare the vector "segment type", which later on will inform on
    //   the type of boundary condition implemented by a given segment.
    //   The default, 0, stands for no-slip.
    segmentTagList.resize ( meshes[1].getNumSegments() );
    std::fill ( segmentTagList.begin(), segmentTagList.end(), 0 );
//TODO: fix the function closeHoles;
    if ( automaticCloseHoles )
    {
        closeHoles();
    }
}

template<typename T>
template<typename TMesh>
SegmentBoundary2D<T>::SegmentBoundary2D (
    DEFscaledMesh2D<TMesh> const& defMesh, bool automaticCloseHoles )
    : currentTagNum ( 0 ),
      margin ( defMesh.getMargin() ),
      physicalLocation ( defMesh.getPhysicalLocation() ),
      dx ( defMesh.getDx() )
{
    topology.push ( 1 ); // By default, closed mesh.
    vertexSet.push ( 0 ); // By default, static mesh.

    vertexLists.reserve ( 3 ); // Vertex lists are expensive to copy. Better
    //   pre-allocate a slot for three of them.
    std::vector<Array<TMesh,3> > const& vertexList = defMesh.getVertexList();
    std::vector<Array<T,2> > newVertexList ( vertexList.size() );
    for ( pluint i=0; i<vertexList.size(); ++i )
    {
        newVertexList[i] = Array<T,2> ( vertexList[i] );
    }
    vertexLists.push_back ( newVertexList );
    emanatingEdgeLists[0] = defMesh.getEmanatingEdgeList();
    edgeLists[0] = defMesh.getEdgeList();

    emanatingEdgeLists[1] = emanatingEdgeLists[0];
    edgeLists[1] = edgeLists[0];

    meshes.push_back ( SegmentPolygonMesh2D<T> (
                           vertexLists[0], emanatingEdgeLists[0], edgeLists[0] ) );
    meshes.push_back ( SegmentPolygonMesh2D<T> (
                           vertexLists[0], emanatingEdgeLists[1], edgeLists[1] ) );

    // Prepare the vector "segment type", which later on will inform on
    //   the type of boundary condition implemented by a given segment.
    //   The default, 0, stands for no-slip.
    segmentTagList.resize ( meshes[1].getNumSegments() );
    std::fill ( segmentTagList.begin(), segmentTagList.end(), 0 );

    if ( automaticCloseHoles )
    {
        closeHoles();
    }
}

template<typename T>
void SegmentBoundary2D<T>::closeHoles()
{
    // Close the holes and assign inlets/outlets.
    std::vector<Curve2D> newlids
        = meshes[1].closeHoles();

    // Prepare the vector "segment type", which later on will inform on
    //   the type of boundary condition implemented by a given segment.
    //   The default, 0, stands for no-slip.
    pluint oldNumSegments = segmentTagList.size();
    segmentTagList.resize ( meshes[1].getNumSegments() );
    std::fill ( segmentTagList.begin() +oldNumSegments, segmentTagList.end(), 0 );

    // Assign default functions to inlet/outlets to avoid undefined state.
    tagInletOutlet ( newlids );
    lids.insert ( lids.end(),newlids.begin(),newlids.end() );
    assignLidVertexProperty();
}

template<typename T>
SegmentBoundary2D<T>::~SegmentBoundary2D()
{
    for ( pluint iProp=0; iProp<vertexProperties.size(); ++iProp )
    {
        delete vertexProperties[iProp];
    }
}

template<typename T>
SegmentBoundary2D<T>::SegmentBoundary2D (
    SegmentBoundary2D<T> const& rhs )
    : vertexLists ( rhs.vertexLists ),
      emanatingEdgeLists ( rhs.emanatingEdgeLists ),
      edgeLists ( rhs.edgeLists ),
      segmentTagList ( rhs.segmentTagList ),
      currentTagNum ( rhs.currentTagNum ),
      vertexTagList ( rhs.vertexTagList ),
      vertexProperties ( rhs.vertexProperties.size() ),
      lids ( rhs.lids ),
      margin ( rhs.margin ),
      topology ( rhs.topology ),
      vertexSet ( rhs.vertexSet )
{
    defineMeshes();
    for ( pluint iProp=0; iProp<vertexProperties.size(); ++iProp )
    {
        vertexProperties[iProp] = rhs.vertexProperties[iProp]->clone();
    }
}

template<typename T>
void SegmentBoundary2D<T>::defineMeshes()
{
    meshes.clear();
    for ( pluint iVertices=0; iVertices<vertexLists.size(); ++iVertices )
    {
        meshes.push_back ( SegmentPolygonMesh2D<T> (
                               vertexLists[iVertices], emanatingEdgeLists[0], edgeLists[0] ) );
        meshes.push_back ( SegmentPolygonMesh2D<T> (
                               vertexLists[iVertices], emanatingEdgeLists[1], edgeLists[1] ) );
    }
}

template<typename T>
SegmentBoundary2D<T>&
SegmentBoundary2D<T>::operator= ( SegmentBoundary2D<T> const& rhs )
{
    SegmentBoundary2D<T> ( rhs ).swap ( *this );
    return *this;
}

template<typename T>
void SegmentBoundary2D<T>::swap ( SegmentBoundary2D<T>& rhs )
{
    vertexLists.swap ( rhs.vertexLists );
    emanatingEdgeLists[0].swap ( rhs.emanatingEdgeLists[0] );
    emanatingEdgeLists[1].swap ( rhs.emanatingEdgeLists[1] );
    edgeLists[0].swap ( rhs.edgeLists[0] );
    edgeLists[1].swap ( rhs.edgeLists[1] );
    segmentTagList.swap ( rhs.segmentTagList );
    std::swap ( currentTagNum, rhs.currentTagNum );
    vertexTagList.swap ( rhs.vertexTagList );
    vertexProperties.swap ( rhs.vertexProperties );
    std::swap ( lids, rhs.lids );
    std::swap ( margin, rhs.margin );
    std::swap ( topology, rhs.topology );
    std::swap ( vertexSet, rhs.vertexSet );
    defineMeshes();
}

template<typename T>
SegmentBoundary2D<T> const&
SegmentBoundary2D<T>::select (
    plint whichTopology, plint whichVertices ) const
{
    PLB_PRECONDITION ( whichTopology==0 || whichTopology==1 );
    PLB_PRECONDITION ( whichVertices>=0 && whichVertices < ( plint ) vertexLists.size() );
    topology.top() = whichTopology;
    vertexSet.top() = whichVertices;
    return *this;
}

template<typename T>
SegmentBoundary2D<T> const&
SegmentBoundary2D<T>::pushSelect (
    plint whichTopology, plint whichVertices ) const
{
    PLB_PRECONDITION ( whichTopology==0 || whichTopology==1 );
    PLB_PRECONDITION ( whichVertices>=0 && whichVertices < ( plint ) vertexLists.size() );
    topology.push ( whichTopology );
    vertexSet.push ( whichVertices );
    return *this;
}

template<typename T>
SegmentBoundary2D<T> const&
SegmentBoundary2D<T>::popSelect() const
{
    PLB_PRECONDITION ( topology.size() >= 2 );
    PLB_PRECONDITION ( vertexSet.size() >= 2 );
    topology.pop();
    vertexSet.pop();
    return *this;
}

template<typename T>
void SegmentBoundary2D<T>::getSelection (
    plint& whichTopology, plint& whichVertices ) const
{
    whichTopology = topology.top();
    whichVertices = vertexSet.top();
}

template<typename T>
plint SegmentBoundary2D<T>::currentMesh() const
{
    return 2*vertexSet.top() +topology.top();
}

template<typename T>
SegmentPolygonMesh2D<T> const& SegmentBoundary2D<T>::getMesh() const
{
    return meshes[currentMesh()];
}

template<typename T>
SegmentPolygonMesh2D<T>& SegmentBoundary2D<T>::getMesh()
{
    return meshes[currentMesh()];
}

template<typename T>
plint SegmentBoundary2D<T>::getMargin() const
{
    return margin;
}

template<typename T>
plint SegmentBoundary2D<T>::getTag ( plint iSegment ) const
{
    PLB_ASSERT ( iSegment< ( plint ) segmentTagList.size() );
    return segmentTagList[iSegment];
}

template<typename T>
VertexProperty2D<T> const*
SegmentBoundary2D<T>::getVertexProperty ( plint iVertex ) const
{
    if ( vertexTagList.empty() )
    {
        return 0;
    }
    PLB_ASSERT ( iVertex < ( plint ) vertexTagList.size() );
    return vertexProperties[vertexTagList[iVertex]];
}

template<typename T>
bool SegmentBoundary2D<T>::intersectSegment (
    plint iSegment, AtomicBlock2D* boundaryArg,
    Array<T,2> const& fromPoint, Array<T,2> const& direction,
    Array<T,2>& locatedPoint, T& distance, Array<T,2>& wallNormal ) const
{
    int flag = 0; // Intersection with line segment.
    Array<T,2> point2 ( fromPoint+direction );
    bool doesIntersect =
        getMesh().pointOnSegment ( fromPoint, point2, flag, iSegment,
                                    locatedPoint, wallNormal, distance ) == 1;
    return doesIntersect;
}

template<typename T>
Array<T,2> SegmentBoundary2D<T>::computeContinuousNormal (
    Array<T,2> const& p, plint iSegment, bool isAreaWeighted ) const
{
    return getMesh().computeContinuousNormal ( p, iSegment, isAreaWeighted );
}

template<typename T>
void SegmentBoundary2D<T>::cloneVertexSet ( plint whichVertexSet )
{
    PLB_PRECONDITION ( whichVertexSet < ( plint ) vertexLists.size() &&
                       whichVertexSet >= 0 );
    vertexLists.push_back ( vertexLists[whichVertexSet] );
    plint newVertexSet = ( plint ) vertexLists.size()-1;
    plint numVertexOpen = vertexLists[newVertexSet].size();
    for ( plint iLid=0; ( plint ) iLid<lids.size(); ++iLid )
    {
        numVertexOpen -= lids[iLid].numAddedVertices;
    }
    meshes.push_back ( SegmentPolygonMesh2D<T> (
                           vertexLists[newVertexSet], emanatingEdgeLists[0], edgeLists[0], numVertexOpen ) );
    meshes.push_back ( SegmentPolygonMesh2D<T> (
                           vertexLists[newVertexSet], emanatingEdgeLists[1], edgeLists[1] ) );
}

template<typename T>
std::vector<Curve2D> const&
SegmentBoundary2D<T>::getInletOutlet() const
{
    PLB_PRECONDITION ( topology.top() ==1 );
    return lids;
}

template<typename T>
std::vector<Curve2D>
SegmentBoundary2D<T>::getInletOutlet ( plint sortDirection ) const
{
    PLB_PRECONDITION ( topology.top() ==1 );
    std::vector<Curve2D> lids_copy ( lids );
    std::sort ( lids_copy.begin(), lids_copy.end(), LidLessThan2D<T> ( sortDirection, getMesh() ) );
    return lids_copy;
}

template<typename T>
std::vector<plint> SegmentBoundary2D<T>::getInletOutletIds ( plint sortDirection ) const
{
    std::map<plint,plint> segmentToOriginalLid;
    for ( pluint iLid=0; iLid<lids.size(); ++iLid )
    {
        segmentToOriginalLid[lids[iLid].firstSegment] = iLid;
    }
    std::vector<Curve2D> tmpLids ( lids );
    std::sort ( tmpLids.begin(), tmpLids.end(), LidLessThan2D<T> ( sortDirection, getMesh() ) );
    std::vector<plint> ids ( tmpLids.size() );
    for ( pluint iLid=0; iLid<tmpLids.size(); ++iLid )
    {
        plint originalId = segmentToOriginalLid[tmpLids[iLid].firstSegment];
        ids[iLid] = originalId+1;
    }
    return ids;
}

template<typename T>
void SegmentBoundary2D<T>::getLidProperties (
    plint sortDirection, std::vector<Array<T,2> >& normal,
    std::vector<Array<T,2> >& center, std::vector<T>& radius ) const
{
    // Lid properties can only be computed in a closed mesh, by definition.
    PLB_PRECONDITION ( topology.top() ==1 );
    std::vector<Curve2D> tmpLids ( lids );
    std::sort ( tmpLids.begin(), tmpLids.end(), LidLessThan2D<T> ( sortDirection, getMesh() ) );
    normal.resize ( tmpLids.size() );
    center.resize ( tmpLids.size() );
    radius.resize ( tmpLids.size() );
    for ( pluint iLid=0; iLid<tmpLids.size(); ++iLid )
    {
        normal[iLid] = computeNormal ( getMesh(), tmpLids[iLid] );
        center[iLid] = computeBaryCenter2D ( getMesh(), tmpLids[iLid] );
        radius[iLid] = ( T ) 0.5 * (
                           computeInnerRadius ( getMesh(), tmpLids[iLid] ) +
                           computeOuterRadius ( getMesh(), tmpLids[iLid] ) );
    }
}

template<typename T>
void SegmentBoundary2D<T>::tagInletOutlet (
    std::vector<Curve2D> const& newLids )
{
    // Inlet/Outlet can only be set for a closed mesh, by definition.
    PLB_PRECONDITION ( topology.top() ==1 );

    for ( pluint iLid=0; iLid<newLids.size(); ++iLid )
    {
        ++currentTagNum; // Tag 0 is for default wall portions.
        plint firstSegment = newLids[iLid].firstSegment;
        plint numSegments = newLids[iLid].numSegments;
        for ( plint iSegment = firstSegment; iSegment<firstSegment+numSegments; ++iSegment )
        {
            segmentTagList[iSegment] = currentTagNum;
        }
    }
}

template<typename T>
template<typename DomainFunctional>
plint SegmentBoundary2D<T>::tagDomain ( DomainFunctional functional )
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION ( topology.top() ==1 );
    ++currentTagNum;
    plint newTag = currentTagNum;
    for ( plint iSegment=0; iSegment<getMesh().getNumSegments(); ++iSegment )
    {
        bool isInside = true;
        for ( plint iVertex=0; iVertex<3; ++iVertex )
        {
            Array<T,2> vertex = getMesh().getVertex ( iSegment, iVertex );
            if ( !functional ( vertex ) )
            {
                isInside = false;
                break;
            }
        }
        if ( isInside )
        {
            segmentTagList[iSegment] = newTag;
        }
    }
    return newTag;
}

template<typename T>
template<typename DomainFunctional>
plint SegmentBoundary2D<T>::tagDomain ( DomainFunctional functional, Array<T,2> normal, T angleTolerance, plint previousTag )
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION ( topology.top() ==1 );
    ++currentTagNum;
    plint newTag = currentTagNum;
    for ( plint iSegment=0; iSegment<getMesh().getNumSegments(); ++iSegment )
    {
        if ( previousTag<0 || segmentTagList[iSegment]==previousTag )
        {
            bool isInside = true;
            for ( plint iVertex=0; iVertex<3; ++iVertex )
            {
                Array<T,2> vertex = getMesh().getVertex ( iSegment, iVertex );
                if ( !functional ( vertex ) )
                {
                    isInside = false;
                    break;
                }
            }
            if ( isInside )
            {
                Array<T,2> segmentNormal = getMesh().computeSegmentNormal ( iSegment );
                if ( fabs ( angleBetweenVectors ( normal,segmentNormal ) <angleTolerance ) )
                {
                    segmentTagList[iSegment] = newTag;
                }
            }
        }
    }
    return newTag;
}

template<typename T>
template<typename DomainFunctional>
plint SegmentBoundary2D<T>::setVertexProperty (
    VertexProperty2D<T> const& property, DomainFunctional functional )
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION ( topology.top() ==1 );
    if ( vertexTagList.empty() )
    {
        PLB_ASSERT ( vertexProperties.empty() );
        vertexTagList.resize ( getMesh().getNumVertices() );
        std::fill ( vertexTagList.begin(), vertexTagList.end(), 0 );
        vertexProperties.push_back ( 0 );
    }
    vertexProperties.push_back ( property.clone() );
    plint nextTag = ( plint ) vertexProperties.size()-1;
    for ( plint iVertex=0; iVertex<getMesh().getNumVertices(); ++iVertex )
    {
        Array<T,2> const& vertex = getMesh().getVertex ( iVertex );
        bool isInside = functional ( vertex );
        if ( isInside )
        {
            vertexTagList[iVertex] = nextTag;
        }
    }
    return nextTag;
}

template<typename T>
void SegmentBoundary2D<T>::assignLidVertexProperty()
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION ( topology.top() ==1 );
    std::vector<Curve2D> const& lids = getInletOutlet();
    if ( lids.empty() ) return;

    if ( vertexTagList.empty() )
    {
        PLB_ASSERT ( vertexProperties.empty() );
        vertexTagList.resize ( getMesh().getNumVertices() );
        std::fill ( vertexTagList.begin(), vertexTagList.end(), 0 );
        vertexProperties.push_back ( 0 );
    }

    vertexProperties.push_back ( new InletOutletProperty2D<T> );
    plint lidVertexTag = ( plint ) vertexProperties.size()-1;
    for ( pluint iLid=0; iLid<lids.size(); ++iLid )
    {
        plint centerVertex = lids[iLid].centerVertex;
        vertexTagList[centerVertex] = lidVertexTag;
    }
}


/******** class SegmentFlowShape2D ****************************************/

template< typename T, class SurfaceData >
SegmentFlowShape2D<T,SurfaceData>::SegmentFlowShape2D (
    SegmentBoundary2D<T> const& boundary_,
    BoundaryProfiles2D<T,SurfaceData> const& profiles_ )
    : boundary ( boundary_ ),
      profiles ( profiles_ ),
      voxelFlags ( 0 ), hashContainer ( 0 ), boundaryArg ( 0 )
{ }

template< typename T, class SurfaceData >
bool SegmentFlowShape2D<T,SurfaceData>::isInside (
    Dot2D const& location ) const
{
    PLB_PRECONDITION ( voxelFlags );
    Dot2D localPos = location-voxelFlags->getLocation();
    return voxelFlag::insideFlag ( voxelFlags->get ( localPos.x,localPos.y ) );
}

template< typename T, class SurfaceData >
bool SegmentFlowShape2D<T,SurfaceData>::pointOnSurface (
    Array<T,2> const& fromPoint, Array<T,2> const& direction,
    Array<T,2>& locatedPoint, T& distance,
    Array<T,2>& wallNormal, SurfaceData& surfaceData,
    OffBoundary::Type& bdType, plint& id ) const
{
    PLB_PRECONDITION ( hashContainer ); // Make sure these arguments have
    PLB_PRECONDITION ( boundaryArg );  //   been provided by the user through
    //   the clone function.
    static const T maxDistance = sqrt ( 3 );
    Array<T,2> xRange ( fromPoint[0]-maxDistance, fromPoint[0]+maxDistance );
    Array<T,2> yRange ( fromPoint[1]-maxDistance, fromPoint[1]+maxDistance );
    SegmentHash<T> segmentHash ( *hashContainer );
    std::vector<plint> possibleSegments;
    if ( id>=0 && id<boundary.getMesh().getNumSegments() )
    {
        possibleSegments.push_back ( id );
    }
    else
    {
        segmentHash.getSegments ( xRange, yRange, possibleSegments );
    }

    Array<T,2>  tmpLocatedPoint;
    T           tmpDistance;
    Array<T,2>  tmpNormal;
    T shortestDistance = T();
    plint locatedSegment = -1;

    for ( pluint iPossible=0; iPossible<possibleSegments.size(); ++iPossible )
    {
        plint iSegment = possibleSegments[iPossible];
        if ( boundary.intersectSegment (
                    iSegment, boundaryArg,
                    fromPoint, direction,
                    tmpLocatedPoint, tmpDistance, tmpNormal ) )
        {
            if ( locatedSegment==-1 || tmpDistance<shortestDistance )
            {
                shortestDistance = tmpDistance;
                locatedSegment = iSegment;
                locatedPoint = tmpLocatedPoint;
                distance = tmpDistance;
                wallNormal = tmpNormal;
                profiles.getProfile ( boundary, iSegment ).getData (
                    locatedPoint, iSegment, boundaryArg, surfaceData, bdType );
            }
        }
    }
    if ( locatedSegment != -1 )
    {
        id = locatedSegment;
        return true;
    }
    else
    {
        return false;
    }
}

template<typename T, class SurfaceData>
Array<T,2> SegmentFlowShape2D<T,SurfaceData>::computeContinuousNormal (
    Array<T,2> const& p, plint id, bool isAreaWeighted ) const
{
    return boundary.computeContinuousNormal ( p, id, isAreaWeighted );
}

template<typename T, class SurfaceData>
bool SegmentFlowShape2D<T,SurfaceData>::intersectsSurface (
    Array<T,2> const& p1, Array<T,2> const& p2, plint& id ) const
{
    PLB_PRECONDITION ( hashContainer ); // Make sure these arguments have
    PLB_PRECONDITION ( boundaryArg );  //   been provided by the user through
    //   the clone function.
    static const T maxDistance = sqrt ( 3 );
    Array<T,2> xRange ( p1[0]-maxDistance, p1[0]+maxDistance );
    Array<T,2> yRange ( p1[1]-maxDistance, p1[1]+maxDistance );
    SegmentHash<T> segmentHash ( *hashContainer );
    std::vector<plint> possibleSegments;
    if ( id>=0 && id<boundary.getMesh().getNumSegments() )
    {
        possibleSegments.push_back ( id );
    }
    else
    {
        segmentHash.getSegments ( xRange, yRange, possibleSegments );
    }

    std::vector<plint> selection;
    for ( pluint iPossible=0; iPossible<possibleSegments.size(); ++iPossible )
    {
        plint iSegment = possibleSegments[iPossible];
        int flag = 0;
        Array<T,2> intersection, normal;
        T distance;
        if ( boundary.getMesh().pointOnSegment ( p1, p2, flag, iSegment, intersection, normal, distance ) )
        {
            selection.push_back ( iSegment );
        }
    }
    if ( selection.empty() )
    {
        return false;
    }
    else if ( selection.size() ==1 )
    {
        id = selection[0];
        return true;
    }
    else
    {
        Array<T,2>  locatedPoint;
        T           distance;
        Array<T,2>  wallNormal;
        SurfaceData surfaceData;
        OffBoundary::Type bdType;
        return pointOnSurface ( p1, p2-p1, locatedPoint, distance, wallNormal, surfaceData, bdType, id );
    }
}

template< typename T, class SurfaceData >
plint SegmentFlowShape2D<T,SurfaceData>::getTag ( plint id ) const
{
    return boundary.getTag ( id );
}

template< typename T, class SurfaceData >
bool SegmentFlowShape2D<T,SurfaceData>::distanceToSurface (
    Array<T,2> const& point, T& distance, bool& isBehind ) const
{
    PLB_PRECONDITION ( hashContainer ); // Make sure these arguments have
    PLB_PRECONDITION ( boundaryArg );  //   been provided by the user through
    //   the clone function.
    T maxDistance = sqrt ( 3 );
    Array<T,2> xRange ( point[0]-maxDistance, point[0]+maxDistance );
    Array<T,2> yRange ( point[1]-maxDistance, point[1]+maxDistance );
    SegmentHash<T> segmentHash ( *hashContainer );
    std::vector<plint> possibleSegments;
    segmentHash.getSegments ( xRange, yRange, possibleSegments );

    T    tmpDistance;
    bool tmpIsBehind;
    bool segmentFound = false;

    for ( pluint iPossible=0; iPossible<possibleSegments.size(); ++iPossible )
    {
        plint iSegment = possibleSegments[iPossible];
        boundary.getMesh().distanceToSegment (
            point, iSegment, tmpDistance, tmpIsBehind );
        if ( !segmentFound || tmpDistance<distance )
        {
            distance = tmpDistance;
            isBehind = tmpIsBehind;
            segmentFound = true;
        }
    }
    return segmentFound;
}

template< typename T, class SurfaceData >
SegmentFlowShape2D<T,SurfaceData>*
SegmentFlowShape2D<T,SurfaceData>::clone() const
{
    return new SegmentFlowShape2D<T,SurfaceData> ( *this );
}

template< typename T, class SurfaceData >
SegmentFlowShape2D<T,SurfaceData>*
SegmentFlowShape2D<T,SurfaceData>::clone (
    std::vector<AtomicBlock2D*> args ) const
{
    PLB_PRECONDITION ( args.size() ==3 );
    SegmentFlowShape2D<T,SurfaceData>*
    newShape = new SegmentFlowShape2D<T,SurfaceData> ( *this );
    newShape->voxelFlags = dynamic_cast<ScalarField2D<int>*> ( args[0] );
    newShape->hashContainer = dynamic_cast<AtomicContainerBlock2D*> ( args[1] );
    newShape->boundaryArg = args[2];
    PLB_ASSERT ( newShape->voxelFlags );
    PLB_ASSERT ( newShape->hashContainer );
    PLB_ASSERT ( newShape->boundaryArg );
    return newShape;
}


/******** class VoxelizedDomain2D *****************************************/

template<typename T>
VoxelizedDomain2D<T>::VoxelizedDomain2D (
    SegmentBoundary2D<T> const& boundary_,
    int flowType_, plint extraLayer_, plint borderWidth_,
    plint envelopeWidth_, plint blockSize_, plint gridLevel_, bool dynamicMesh_ )
    : flowType ( flowType_ ),
      borderWidth ( borderWidth_ ),
      boundary ( boundary_ )
{
    PLB_ASSERT ( flowType==voxelFlag::inside || flowType==voxelFlag::outside );
    PLB_ASSERT ( boundary.getMargin() >= borderWidth );
    if ( dynamicMesh_ )
    {
        boundary.pushSelect ( 1,1 ); // Closed, Dynamic.
    }
    else
    {
        boundary.pushSelect ( 1,0 ); // Closed, Static.
    }
    std::auto_ptr<MultiScalarField2D<int> > fullVoxelMatrix (
        voxelize2D ( boundary.getMesh(),
                   boundary.getMargin() +extraLayer_, borderWidth ) );
    fullVoxelMatrix->setRefinementLevel ( gridLevel_ );
    createSparseVoxelMatrix ( *fullVoxelMatrix, blockSize_, envelopeWidth_ );
    createSegmentHash();
    boundary.popSelect();
}

template<typename T>
VoxelizedDomain2D<T>::VoxelizedDomain2D (
    SegmentBoundary2D<T> const& boundary_,
    int flowType_, Box2D const& boundingBox, plint borderWidth_,
    plint envelopeWidth_, plint blockSize_, plint gridLevel_, bool dynamicMesh_ )
    : flowType ( flowType_ ),
      borderWidth ( borderWidth_ ),
      boundary ( boundary_ )
{
    PLB_ASSERT ( flowType==voxelFlag::inside || flowType==voxelFlag::outside );
    PLB_ASSERT ( boundary.getMargin() >= borderWidth );
    if ( dynamicMesh_ )
    {
        boundary.pushSelect ( 1,1 ); // Closed, Dynamic.
    }
    else
    {
        boundary.pushSelect ( 1,0 ); // Closed, Static.
    }
    std::auto_ptr<MultiScalarField2D<int> > fullVoxelMatrix (
        voxelize2D ( boundary.getMesh(), boundingBox, borderWidth ) );
    fullVoxelMatrix->setRefinementLevel ( gridLevel_ );
    createSparseVoxelMatrix ( *fullVoxelMatrix, blockSize_, envelopeWidth_ );
    createSegmentHash();
    boundary.popSelect();
}

template<typename T>
VoxelizedDomain2D<T>::VoxelizedDomain2D (
    SegmentBoundary2D<T> const& boundary_,
    int flowType_, Box2D const& boundingBox, plint borderWidth_,
    plint envelopeWidth_, plint blockSize_, Box2D const& seed, plint gridLevel_, bool dynamicMesh_ )
    : flowType ( flowType_ ),
      borderWidth ( borderWidth_ ),
      boundary ( boundary_ )
{
    PLB_ASSERT ( flowType==voxelFlag::inside || flowType==voxelFlag::outside );
    PLB_ASSERT ( boundary.getMargin() >= borderWidth );
    if ( dynamicMesh_ )
    {
        boundary.pushSelect ( 1,1 ); // Closed, Dynamic.
    }
    else
    {
        boundary.pushSelect ( 1,0 ); // Closed, Static.
    }
    std::auto_ptr<MultiScalarField2D<int> > fullVoxelMatrix (
        voxelize2D ( boundary.getMesh(), boundingBox, borderWidth, seed ) );
    fullVoxelMatrix->setRefinementLevel ( gridLevel_ );
    createSparseVoxelMatrix ( *fullVoxelMatrix, blockSize_, envelopeWidth_ );
    createSegmentHash();
    boundary.popSelect();
}

template<typename T>
void VoxelizedDomain2D<T>::createSparseVoxelMatrix (
    MultiScalarField2D<int>& fullVoxelMatrix,
    plint blockSize_, plint envelopeWidth_ )
{
    if ( blockSize_>0 )
    {
        computeSparseVoxelMatrix ( fullVoxelMatrix, blockSize_, envelopeWidth_ );
    }
    else
    {
        // blockSize=0 means: don't create a sparse representation of the
        //   multi-block. The envelope of the voxel-matrix needs to be extended
        //   for future usage, though.
        extendEnvelopeWidth ( fullVoxelMatrix, envelopeWidth_ );
    }
}

template<typename T>
VoxelizedDomain2D<T>::VoxelizedDomain2D (
    VoxelizedDomain2D<T> const& rhs )
    : boundary ( rhs.boundary ),
      voxelMatrix ( new MultiScalarField2D<int> ( *rhs.voxelMatrix ) ),
      segmentHash ( new MultiContainerBlock2D ( *rhs.segmentHash ) )
{ }

template<typename T>
VoxelizedDomain2D<T>::~VoxelizedDomain2D()
{
    delete voxelMatrix;
    delete segmentHash;
}


template<typename T>
MultiScalarField2D<int>&
VoxelizedDomain2D<T>::getVoxelMatrix()
{
    return *voxelMatrix;
}

template<typename T>
MultiScalarField2D<int> const&
VoxelizedDomain2D<T>::getVoxelMatrix() const
{
    return *voxelMatrix;
}

template<typename T>
MultiContainerBlock2D&
VoxelizedDomain2D<T>::getSegmentHash()
{
    return *segmentHash;
}

template<typename T>
template<class ParticleFieldT>
void VoxelizedDomain2D<T>::adjustVoxelization (
    MultiParticleField2D<ParticleFieldT>& particles, bool dynamicMesh )
{
    if ( dynamicMesh )
    {
        boundary.pushSelect ( 1,1 ); // Closed, Dynamic.
    }
    else
    {
        boundary.pushSelect ( 1,0 ); // Closed, Static.
    }
    reCreateSegmentHash ( particles );
    MultiScalarField2D<int>* newVoxelMatrix =
        revoxelize2D ( boundary.getMesh(), *voxelMatrix, *segmentHash, borderWidth ).release();
    std::swap ( voxelMatrix, newVoxelMatrix );
    delete newVoxelMatrix;
    boundary.popSelect();
}

template<typename T>
void VoxelizedDomain2D<T>::reparallelize ( MultiBlockRedistribute2D const& redistribute )
{
    MultiBlockManagement2D newManagement = redistribute.redistribute ( voxelMatrix->getMultiBlockManagement() );
    MultiScalarField2D<int>* newVoxelMatrix =
        new MultiScalarField2D<int> (
        newManagement,
        voxelMatrix->getBlockCommunicator().clone(),
        voxelMatrix->getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<int>(), 0 );
    copyNonLocal ( *voxelMatrix, *newVoxelMatrix, voxelMatrix->getBoundingBox() );
    std::swap ( voxelMatrix, newVoxelMatrix );
    delete newVoxelMatrix;
    delete segmentHash;
    createSegmentHash();
}

template<typename T>
MultiBlockManagement2D const&
VoxelizedDomain2D<T>::getMultiBlockManagement() const
{
    return voxelMatrix->getMultiBlockManagement();
}

template<typename T>
void VoxelizedDomain2D<T>::computeSparseVoxelMatrix (
    MultiScalarField2D<int>& fullVoxelMatrix,
    plint blockSize, plint envelopeWidth )
{
    // Initialized to zero.
    MultiScalarField2D<int> domainMatrix ( ( MultiBlock2D const& ) fullVoxelMatrix );
    setToConstant ( domainMatrix, fullVoxelMatrix,flowType, domainMatrix.getBoundingBox(), 1 );
    for ( int iLayer=1; iLayer<=boundary.getMargin(); ++iLayer )
    {
        addLayer ( domainMatrix, domainMatrix.getBoundingBox(), iLayer );
    }
    MultiBlockManagement2D sparseBlockManagement =
        computeSparseManagement (
            *plb::reparallelize ( domainMatrix, blockSize,blockSize ),
            envelopeWidth );

    voxelMatrix = new MultiScalarField2D<int> (
        sparseBlockManagement,
        fullVoxelMatrix.getBlockCommunicator().clone(),
        fullVoxelMatrix.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<int>(),
        voxelFlag::undetermined );
    copyNonLocal ( fullVoxelMatrix, *voxelMatrix, voxelMatrix->getBoundingBox() );
}

template<typename T>
void VoxelizedDomain2D<T>::extendEnvelopeWidth (
    MultiScalarField2D<int>& fullVoxelMatrix, plint envelopeWidth )
{
    MultiBlockManagement2D const& oldManagement =
        fullVoxelMatrix.getMultiBlockManagement();
    MultiBlockManagement2D newManagement (
        oldManagement.getSparseBlockStructure(),
        oldManagement.getThreadAttribution().clone(),
        envelopeWidth,
        oldManagement.getRefinementLevel() );
    voxelMatrix = new MultiScalarField2D<int> (
        newManagement,
        fullVoxelMatrix.getBlockCommunicator().clone(),
        fullVoxelMatrix.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<int>(),
        voxelFlag::undetermined );
    plb::copy ( fullVoxelMatrix, *voxelMatrix, fullVoxelMatrix.getBoundingBox() );
}

template<typename T>
void VoxelizedDomain2D<T>::createSegmentHash()
{
    segmentHash = new MultiContainerBlock2D ( *voxelMatrix );
    std::vector<MultiBlock2D*> hashArg;
    hashArg.push_back ( segmentHash );
    applyProcessingFunctional (
        new CreateSegmentHash<T> ( boundary.getMesh() ),
        segmentHash->getBoundingBox(), hashArg );
}

template<typename T>
template<class ParticleFieldT>
void VoxelizedDomain2D<T>::reCreateSegmentHash (
    MultiParticleField2D<ParticleFieldT>& particles )
{
    // The lids are non-parallel, an info which must be provided
    //   to the hash algorithm by means of a list of barycenters.
    //   This is a necessary and sufficient information because
    //   all segments connected to the barycenters are non-parallel,
    //   and all non-parallel segments are connected to a barycenter.
    std::vector<Curve2D> const& lids = boundary.getInletOutlet();
    plint numLids = ( plint ) lids.size();
    std::vector<plint> lidBaryCenters ( numLids );
    for ( plint iLid=0; iLid<numLids; ++iLid )
    {
        lidBaryCenters[iLid] = lids[iLid].centerVertex;
    }

    std::vector<MultiBlock2D*> hashParticleArg;
    hashParticleArg.push_back ( segmentHash );
    hashParticleArg.push_back ( &particles );
    applyProcessingFunctional (
        new ReAssignSegmentHash<T,ParticleFieldT> ( boundary.getMesh(),lidBaryCenters ),
        segmentHash->getBoundingBox(), hashParticleArg );
}


/* ******** DetectBorderLineFunctional2D ************************************* */

template<typename T>
void addLayer ( MultiScalarField2D<T>& matrix,
                Box2D const& domain, T previousLayer )
{
    applyProcessingFunctional ( new AddLayerFunctional2D<T> ( previousLayer ),
                                domain, matrix );
}

template<typename T>
AddLayerFunctional2D<T>::AddLayerFunctional2D ( T previousLayer_ )
    : previousLayer ( previousLayer_ )
{ }

template<typename T>
void AddLayerFunctional2D<T>::process (
    Box2D domain, ScalarField2D<T>& voxels )
{
    for ( plint iX = domain.x0; iX <= domain.x1; ++iX )
    {
        for ( plint iY = domain.y0; iY <= domain.y1; ++iY )
        {
            for ( plint dx=-1; dx<=1; ++dx )
                for ( plint dy=-1; dy<=1; ++dy )
                    if ( ! ( dx==0 && dy==0 ) )
                    {
                        plint nextX = iX + dx;
                        plint nextY = iY + dy;
                        if ( voxels.get ( iX,iY ) ==0 &&
                                voxels.get ( nextX,nextY ) ==previousLayer )
                        {
                            voxels.get ( iX,iY ) = previousLayer+1;
                        }
                    }
        }
    }
}

template<typename T>
AddLayerFunctional2D<T>* AddLayerFunctional2D<T>::clone() const
{
    return new AddLayerFunctional2D<T> ( *this );
}

template<typename T>
void AddLayerFunctional2D<T>::getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT AddLayerFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // INNER_FLOW_BOUNDARY_2D_HH
