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

/* Main author: Dimitrios Kontaxakis */

#ifndef SEGMENT_TO_DEF_HH
#define SEGMENT_TO_DEF_HH

#include "offLattice/segmentToDef.h"
#include "core/util.h"
#include <utility>
#include <limits>
#include <cstdlib>
#include <malloc.h>
#include <queue>

namespace plb
{

template<typename T>
void SegmentToDef<T>::vsAdd (
    Array<T,2> const& coord, plint& index, plint& count )
{
    VsNodeIt it = vertexSet.find ( VertexSetNode ( -1,&coord ) );
    if ( it == vertexSet.end() )
    {
        index = count;
        ++count;
        VertexSetNode newNode ( index,&coord );
        vertexSet.insert ( newNode );
    }
    else
    {
        index = it->i;
    }
}

template<typename T>
void SegmentToDef<T>::vsOrder()
{
    VsNodeConstIt it = vertexSet.begin();
    for ( ; it!=vertexSet.end(); ++it )
    {
        vertexList[it->i] = * ( it->vertex );
    }
}

template<typename T>
inline bool SegmentToDef<T>::VsLessThan::vertexComponentLessThan ( T x, T y )
{
    T tmp = ( T ) 0.5 * ( fabs ( x - y ) + fabs ( y - x ) );
    return ( ( x < y ) && ( tmp > epsilon ) );
}

template<typename T>
inline bool SegmentToDef<T>::VsLessThan::vertexComponentEqual ( T x, T y )
{
    return ( ( !vertexComponentLessThan ( x, y ) ) && ( !vertexComponentLessThan ( y, x ) ) );
}

template<typename T>
inline bool SegmentToDef<T>::VsLessThan::vertexLessThan (
    Array<T,2> const& v1, Array<T,2> const& v2 )
{
    return ( ( vertexComponentLessThan ( v1[0], v2[0] ) ||
               ( vertexComponentEqual ( v1[0], v2[0] ) && vertexComponentLessThan ( v1[1], v2[1] ) ) ||
               ( vertexComponentEqual ( v1[0], v2[0] ) && vertexComponentEqual ( v1[1], v2[1] ) ) ) );
}

template<typename T>
inline bool SegmentToDef<T>::VsLessThan::operator() (
    VertexSetNode const& node1, VertexSetNode const& node2 )
{
    Array<T,2> const& v1 = * ( node1.vertex );
    Array<T,2> const& v2 = * ( node2.vertex );

    return vertexLessThan ( v1, v2 );
}

template<typename T>
plint SegmentToDef<T>::searchEdgeList (
    std::vector<SegmentToDef<T>::EdgeListNode> const& edgeList, plint maxv ) const
{
    for ( pluint iList=0; iList<edgeList.size(); ++iList )
    {
        if ( edgeList[iList].maxv==maxv )
        {
            return iList;
        }
    }
    return -1;
}

template<typename T>
typename SegmentToDef<T>::BvmNodeIt SegmentToDef<T>::bvmAdd ( plint id )
{
    BvmNodeIt it = boundaryVertexMap.find ( id );
    if ( it == boundaryVertexMap.end() )
    {
        it = boundaryVertexMap.insert (
                 std::make_pair ( id,BoundaryVertexMapNode() ) ).first;
    }
    return it;
}

template<typename T>
bool SegmentToDef<T>::bvmCheck() const
{
    BvmNodeConstIt it = boundaryVertexMap.begin();
    for ( ; it !=boundaryVertexMap.end(); ++it )
    {
        BoundaryVertexMapNode const& node = it->second;
        if ( node.v1 == -1 || node.v2 == -1 ||
                node.t1 == -1 || node.t2 == -1 ||
                node.counter != 2 )
        {
            return true;
        }
    }
    return false;
}

template<typename T>
void SegmentToDef<T>::bvmLabel()
{
    plint vertex2, segment2, localEdge2;
    BvmNodeConstIt it = boundaryVertexMap.begin();
    for ( ; it !=boundaryVertexMap.end(); ++it )
    {
        plint vertex = it->first; // Index of the boundary vertex
        BoundaryVertexMapNode const& node = it->second;
        vertex2 = node.v2; // Index of the boundary vertex to which
        // the emanating boundary edge from the
        // boundary vertex `vertex' points to
        segment2 = node.t2; // Index of the adjacent segment to
        // the boundary edge emanating from
        // the boundary vertex `vertex'

        localEdge2 = -1; // Local index (in the segment t2) of the
        // boundary edge that points from `vertex' to
        // `vertex2'
        if ( vertex == edgeList[3*segment2 + 2].pv &&
                vertex2 == edgeList[3*segment2].pv )
        {
            localEdge2 = 0;
        }
        else if ( vertex == edgeList[3*segment2].pv &&
                  vertex2 == edgeList[3*segment2 + 1].pv )
        {
            localEdge2 = 1;
        }
        else if ( vertex == edgeList[3*segment2 + 1].pv &&
                  vertex2 == edgeList[3*segment2 + 2].pv )
        {
            localEdge2 = 2;
        }
        else
        {
            PLB_ASSERT ( false ); // Problem with the boundary of the surface mesh.
        }

        if ( emanatingEdgeList[vertex] == -1 )
            emanatingEdgeList[vertex] = 3*segment2 + localEdge2;
        else
        {
            PLB_ASSERT ( false ); // Problem with the boundary of the surface mesh.
        }
    }
}

template<typename T>
plint& SegmentToDef<T>::globalVertex ( plint segment, plint localVertex )
{
    return edgeList [ 3*segment + ( ( localVertex == 0 ) ? 2 : localVertex-1 ) ].pv;
}

template<typename T>
SegmentToDef<T>::SegmentToDef (
    std::vector<Segment> const& segments, T epsilon )
    : vertexSet ( VsLessThan ( epsilon ) )
{
    numSegments = segments.size();

    numVertices = uniqueVertices ( segments );
    vertexList.resize ( numVertices );
    emanatingEdgeList.resize ( numVertices );

    vsOrder();

    computePointingVertex();

    plint nbe = createEdgeTable();

    if ( fixOrientation() )
    {
#ifdef PLB_DEBUG
        plint nbe_new = createEdgeTable();
#else
        ( void ) createEdgeTable();
#endif
        PLB_ASSERT ( nbe == nbe_new ); // Problem with the surface mesh topology.
    }

    if ( nbe != 0 )
    {
        findBoundaryVertices();
    }

    computeNeighboringEdges();

    computeEmanatingEdges();
}

// Fix the orientation of all segments, so that they are consistently oriented,
//   meaning that all segment normals point either inwards or outwards.
template<typename T>
bool SegmentToDef<T>::fixOrientation()
{
    bool fixedSomething = false;
    std::queue<plint> segmentsToFixNeighbors;
    char *visitedSegments = ( char * ) calloc ( numSegments, sizeof ( char ) );
    for ( plint iSegment = 0; iSegment < numSegments; iSegment++ )
    {
        if ( visitedSegments[iSegment] == 0 )
        {
            segmentsToFixNeighbors.push ( iSegment );
            while ( !segmentsToFixNeighbors.empty() )
            {
                plint segment = segmentsToFixNeighbors.front();
                segmentsToFixNeighbors.pop();
                fixOrientationOfNeighbors ( segment, segmentsToFixNeighbors, visitedSegments,
                                            fixedSomething );
            }
        }
    }

    free ( visitedSegments );
    return fixedSomething;
}

// Fix the orientation of all neighboring segments to a segment of index iSegment,
//   so that ther orientations are consistent to the one of the segment iSegment.
template<typename T>
void SegmentToDef<T>::fixOrientationOfNeighbors ( plint iSegment,
        std::queue<plint>& segmentsToFixNeighbors,
        char* visitedSegments,
        bool& flag )
{
    visitedSegments[iSegment] = 1;

    plint ia = globalVertex ( iSegment, 0 );
    plint ib = globalVertex ( iSegment, 1 );

    plint min_i = std::min ( ia,ib );
    plint max_i = std::max ( ia,ib );
    plint jSegment = -2;
    for ( pluint iEdge = 0; iEdge < edgeTable[min_i].size(); iEdge++ )
        if ( edgeTable[min_i][iEdge].maxv == max_i )
        {
            plint t1 = edgeTable[min_i][iEdge].t1;
            plint t2 = edgeTable[min_i][iEdge].t2;

            if ( iSegment == t1 )
            {
                jSegment = t2;
            }
            else if ( iSegment == t2 )
            {
                jSegment = t1;
            }
            else
            {
                PLB_ASSERT ( false ); // Problem with the surface mesh edges.
            }

            break;
        }

    PLB_ASSERT ( jSegment != -2 ); // Problem with the surface mesh.

    if ( jSegment != -1 )
    {
        if ( visitedSegments[jSegment] == 0 )
        {
            plint j0 = globalVertex ( jSegment, 0 );
            plint j1 = globalVertex ( jSegment, 1 );

            if ( ia == j0 && ib == j1 )
            {
                std::swap ( segmentIndices[jSegment][0], segmentIndices[jSegment][1] );
                globalVertex ( jSegment,0 ) = segmentIndices[jSegment][0];
                globalVertex ( jSegment,1 ) = segmentIndices[jSegment][1];
                flag = true;
            }
#if 0
            else if ( ia == j1 && ib == j2 )
            {
                std::swap ( segmentIndices[jSegment][1], segmentIndices[jSegment][2] );
                globalVertex ( jSegment,1 ) = segmentIndices[jSegment][1];
                globalVertex ( jSegment,2 ) = segmentIndices[jSegment][2];
                flag = true;
            }
            else if ( ia == j2 && ib == j0 )
            {
                std::swap ( segmentIndices[jSegment][2], segmentIndices[jSegment][0] );
                globalVertex ( jSegment,2 ) = segmentIndices[jSegment][2];
                globalVertex ( jSegment,0 ) = segmentIndices[jSegment][0];
                flag = true;
            }
            else if ( ( ia == j1 && ib == j0 ) ||
                      ( ia == j2 && ib == j1 ) ||
                      ( ia == j0 && ib == j2 ) )
            {
                flag = false;
            }
#endif
            else
            {
                PLB_ASSERT ( false ); // Problem with the surface mesh.
            }

            visitedSegments[jSegment] = 1;
            segmentsToFixNeighbors.push ( jSegment );
        }
    }
}


// Make vertices unique through a unique labelling scheme.
//   Return the number of unique vertices.
template<typename T>
plint SegmentToDef<T>::uniqueVertices ( std::vector<Segment> const& segments )
{

    edgeList.resize ( 3*numSegments );

    segmentIndices.resize ( numSegments );
    plint index=0;
    plint count=0;
    for ( plint iSegment=0; iSegment<numSegments; ++iSegment )
    {
        vsAdd ( segments[iSegment][0], index, count );
        segmentIndices[iSegment][0] = index;

        vsAdd ( segments[iSegment][1], index, count );
        segmentIndices[iSegment][1] = index;

        vsAdd ( segments[iSegment][2], index, count );
        segmentIndices[iSegment][2] = index;
    }
    return count;
}

template<typename T>
void SegmentToDef<T>::computePointingVertex()
{
    for ( plint iSegment=0; iSegment<numSegments; ++iSegment )
    {
        globalVertex ( iSegment,0 ) = segmentIndices[iSegment][0];
        globalVertex ( iSegment,1 ) = segmentIndices[iSegment][1];
        globalVertex ( iSegment,2 ) = segmentIndices[iSegment][2];
    }
}

/// Create hash-table of edges, distinguish boundary edges
/// and count boundary edges.
/// The purpose of this table is to identify the edges and
/// to separate between boundary and interior edges.
template<typename T>
plint SegmentToDef<T>::createEdgeTable()
{

    /* Create the edge hash table */
    edgeTable.resize ( 0 );
    edgeTable.resize ( vertexList.size() );
    /* Local indices of edge endpoints */
    int epoint[3][2] = { {0, 1}, {1, 2}, {2,0} };

    plint nbe = 0; /* Number of boundary edges */

    for ( plint iSegment = 0; iSegment < numSegments; ++iSegment )
    {
        for ( plint localEdge = 0; localEdge < 3; ++localEdge )
        {
            plint i0 = globalVertex ( iSegment, epoint[localEdge][0] );
            plint i1 = globalVertex ( iSegment, epoint[localEdge][1] );
            plint min_i = std::min ( i0,i1 );
            plint max_i = std::max ( i0,i1 );

            plint iList = searchEdgeList ( edgeTable[min_i], max_i );
            if ( iList==-1 )
            {
                EdgeListNode newNode;
                newNode.maxv = max_i;
                newNode.t1   = iSegment;
                newNode.t2   = -1;
                edgeTable[min_i].push_back ( newNode );
                ++nbe;
            }
            else
            {
                EdgeListNode& node = edgeTable[min_i][iList];
                if ( node.t2 == -1 )
                {
                    node.t2 = iSegment;
                    --nbe;
                }
                else
                {
                    PLB_ASSERT ( false ); // The surface mesh contains non-manifold edges.
                }
            }
        }
    }
    return nbe;
}


// Create the map of boundary vertices with additional information on the
//   adjacent boundary edges and their orientation.
template<typename T>
void SegmentToDef<T>::findBoundaryVertices()
{
    BoundaryVertexMap boundaryVertexMap;
    BvmNodeIt node_i = boundaryVertexMap.end();
    BvmNodeIt node_j = boundaryVertexMap.end();
    plint jVertex, segment;
    plint iVertex;
    int lock;
    for ( iVertex=0, lock=1;
            iVertex< ( plint ) edgeTable.size(); ++iVertex, lock=1 )
    {
        for ( pluint iEdge=0; iEdge<edgeTable[iVertex].size(); ++iEdge )
        {
            if ( edgeTable[iVertex][iEdge].t2 == -1 )
            {
                if ( lock )
                {
                    node_i = bvmAdd ( iVertex );
                    lock = 0;
                }
                jVertex = edgeTable[iVertex][iEdge].maxv;
                segment = edgeTable[iVertex][iEdge].t1;
                node_j = bvmAdd ( jVertex );
                node_i->second.counter++;
                node_j->second.counter++;

                plint localVertex_i = -1;
                if ( iVertex == globalVertex ( segment, 0 ) )
                    localVertex_i = 0;
                else if ( iVertex == globalVertex ( segment, 1 ) )
                    localVertex_i = 1;
                else if ( iVertex == globalVertex ( segment, 2 ) )
                    localVertex_i = 2;
                else
                    PLB_ASSERT ( false ); // Problem with the boundary of the surface mesh.

                plint localVertex_j = -1;
                if ( jVertex == globalVertex ( segment, 0 ) )
                    localVertex_j = 0;
                else if ( jVertex == globalVertex ( segment, 1 ) )
                    localVertex_j = 1;
                else if ( jVertex == globalVertex ( segment, 2 ) )
                    localVertex_j = 2;
                else
                    PLB_ASSERT ( false ); // Problem with the boundary of the surface mesh.

                if ( ( localVertex_i == 0 && localVertex_j == 1 ) ||
                        ( localVertex_i == 1 && localVertex_j == 2 ) ||
                        ( localVertex_i == 2 && localVertex_j == 0 ) )
                {
                    node_i->second.v2 = jVertex;
                    node_i->second.t2 = segment;
                    node_j->second.v1 = iVertex;
                    node_j->second.t1 = segment;
                }
                else if ( ( localVertex_i == 1 && localVertex_j == 0 ) ||
                          ( localVertex_i == 2 && localVertex_j == 1 ) ||
                          ( localVertex_i == 0 && localVertex_j == 2 ) )
                {
                    node_i->second.v1 = jVertex;
                    node_i->second.t1 = segment;
                    node_j->second.v2 = iVertex;
                    node_j->second.t2 = segment;
                }
            }
        }
    }
#ifdef PLB_DEBUG
    bool failure = bvmCheck();
#else
    ( void ) bvmCheck();
#endif
    PLB_ASSERT ( !failure ); // Problem with the boundary of the surface mesh.
}

/// Use all previous information to fill in the missing fields in
///   the directed edge data.
template<typename T>
void SegmentToDef<T>::computeNeighboringEdges()
{
    for ( plint iVertex=0; iVertex< ( plint ) edgeTable.size(); ++iVertex )
    {
        for ( plint iEdge=0; iEdge< ( plint ) edgeTable[iVertex].size(); ++iEdge )
        {
            EdgeListNode& node = edgeTable[iVertex][iEdge];
            plint jVertex      = node.maxv;
            plint segment1    = node.t1;
            plint segment2    = node.t2;
            plint localEdge1   = -1;

            plint id0 = globalVertex ( segment1, 0 );
            plint id1 = globalVertex ( segment1, 1 );
            plint id2 = globalVertex ( segment1, 2 );

            if ( ( iVertex == id0 && jVertex == id1 ) ||
                    ( iVertex == id1 && jVertex == id0 ) )
                localEdge1 = 0;
            else if ( ( iVertex == id1 && jVertex == id2 ) ||
                      ( iVertex == id2 && jVertex == id1 ) )
                localEdge1 = 1;
            else if ( ( iVertex == id2 && jVertex == id0 ) ||
                      ( iVertex == id0 && jVertex == id2 ) )
                localEdge1 = 2;
            else
                PLB_ASSERT ( false ); // Problem with the surface mesh connectivity.

            if ( segment2 != -1 ) /* Internal edge */
            {
                plint va1 = globalVertex ( segment1, localEdge1 );
                plint vb1 = globalVertex ( segment1, ( ( localEdge1 != 2 ) ? localEdge1+1 : 0 ) );

                computeInternalNeighboringEdges (
                    iVertex, jVertex, segment1, segment2,
                    localEdge1, va1, vb1 );
            }
            else   /* Boundary edge */
            {
                // By definition, on boundary edges, the neighboring
                //   edge has negative id equal to -1.
                edgeList[3*segment1 + localEdge1].ne = -1;
            }
        }
    }
}

template<typename T>
void SegmentToDef<T>::computeInternalNeighboringEdges (
    plint iVertex, plint jVertex, plint segment1, plint segment2,
    plint localEdge1, plint va1, plint vb1 )
{
    plint localEdge2 = -1;

    plint id0 = globalVertex ( segment2, 0 );
    plint id1 = globalVertex ( segment2, 1 );
    plint id2 = globalVertex ( segment2, 2 );

    if ( ( iVertex == id0 && jVertex == id1 ) ||
            ( iVertex == id1 && jVertex == id0 ) )
    {
        localEdge2 = 0;
    }
    else if ( ( iVertex == id1 && jVertex == id2 ) ||
              ( iVertex == id2 && jVertex == id1 ) )
    {
        localEdge2 = 1;
    }
    else if ( ( iVertex == id2 && jVertex == id0 ) ||
              ( iVertex == id0 && jVertex == id2 ) )
    {
        localEdge2 = 2;
    }
    else
    {
        PLB_ASSERT ( false ); // Problem with the surface mesh connectivity.
    }

    plint va2 = globalVertex ( segment2, localEdge2 );
    plint vb2 = globalVertex ( segment2, ( ( localEdge2 != 2 ) ? localEdge2 + 1 : 0 ) );

    if ( va1 != vb2 || vb1 != va2 )
    {
        PLB_ASSERT ( false ); // Problem with the surface mesh orientation.
    }

    edgeList[3*segment1 + localEdge1].ne = 3*segment2 + localEdge2;
    edgeList[3*segment2 + localEdge2].ne = 3*segment1 + localEdge1;
}

template<typename T>
void SegmentToDef<T>::computeEmanatingEdges()
{
    // Handle emanating edges.
    for ( pluint iVertex=0; iVertex<vertexList.size(); ++iVertex )
    {
        emanatingEdgeList[iVertex] = -1;
    }
    // For every boundary vertex choose as an emanating edge
    // the one that belongs to the boundary to facilitate
    // boundary operations.
    bvmLabel();

    for ( plint iSegment = 0; iSegment < numSegments; ++iSegment )
    {
        for ( plint localEdge = 0; localEdge < 3; ++localEdge )
        {
            plint& emanatingEdge =
                emanatingEdgeList[globalVertex ( iSegment, localEdge )];
            if ( emanatingEdge == -1 )
            {
                emanatingEdge = 3*iSegment + localEdge;
            }
        }
    }
}

template<typename T>
void SegmentToDef<T>::generateOnce (
    std::vector<Array<T,2> >& vertexList_,
    std::vector<plint>& emanatingEdgeList_,
    std::vector<Edge2D>& edgeList_ )
{
    vertexList.swap ( vertexList_ );
    emanatingEdgeList.swap ( emanatingEdgeList_ );
    edgeList.swap ( edgeList_ );
}

} // namespace plb

#endif  // SEGMENT_TO_DEF_HH
