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

/* Main author: Huijie Zhang */

#ifndef SEGMENT_TO_DEF_H
#define SEGMENT_TO_DEF_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "offLattice/segmentPolygonMesh2D.h"
#include "offLattice/segmentSet.h"
#include <vector>
#include <map>
#include <set>
#include <queue>

namespace plb {

template<typename T>
class SegmentToDef {
public:
    typedef typename SegmentSet<T>::Segment Segment;
private:  // This class should only be used by the function \sa constructSurfaceMesh().
    SegmentToDef(std::vector<Segment> const& segments, T epsilon=std::numeric_limits<float>::epsilon());
    void generateOnce (
        std::vector<Array<T,2> >& vertexList_,
        std::vector<plint>& emanatingEdgeList_,
        std::vector<Edge2D>& edgeList_ );
private:
    struct EdgeListNode {
        plint maxv;   /* Integer key value representing the maximum of the global
                         indices of the two vertices that constitute a mesh edge */
        plint t1, t2; /* Indices of the adjacent segments to the specific
                         edge. The value -1 in t2 indicates that the edge has no
                         segment neighbor and therefore is a boundary edge */
    };
    struct VertexSetNode {
        VertexSetNode(plint i_, Array<T,2> const* vertex_)
            : i(i_), vertex(vertex_)
        { }
        plint i;                  // Global index of the vertex
        Array<T,2> const* vertex; // Pointer to vertex coordinates
    };
    struct BoundaryVertexMapNode {
        BoundaryVertexMapNode()
            : v1(-1), t1(-1), v2(-1), t2(-1), counter(0)
        { }
        plint v1; /* Global index of the vertex that the boundary edge
                     that points to the boundary vertex of index equal to
                     the map key starts from */
        plint t1; /* Index of the adjacent segment of the boundary edge
                     that points to the boundary vertex of index equal to
                     the map key */

        plint v2; /* Global index of the vertex that the boundary edge
                     that starts from the boundary vertex of index equal
                     to the map key points to */
        plint t2; /* Index of the adjacent segment of the boundary edge
                     that starts from the boundary vertex of index equal
                     to the map key */

        plint counter; /* Counter of boundary edges attached to the boundary
                          vertex of index equal to the map key. If `counter'
                          is not equal to two then there is a problem with
                          the mesh */
    };

    class VsLessThan {
    public:
        VsLessThan(T epsilon_=std::numeric_limits<float>::epsilon())
            : epsilon(epsilon_)
        { }
        bool operator()(VertexSetNode const& node1, VertexSetNode const& node2);
    private:
        bool vertexComponentLessThan(T x, T y);
        bool vertexComponentEqual(T x, T y);
        bool vertexLessThan(Array<T,2> const& v1, Array<T,2> const& v2);
    private:
        T epsilon;
    };
    typedef std::set<VertexSetNode,VsLessThan> VertexSet;
    typedef typename VertexSet::iterator VsNodeIt;
    typedef typename VertexSet::const_iterator VsNodeConstIt;

    typedef std::map<plint, BoundaryVertexMapNode> BoundaryVertexMap;
    typedef typename BoundaryVertexMap::iterator BvmNodeIt;
    typedef typename BoundaryVertexMap::const_iterator BvmNodeConstIt;
private:
    void vsAdd( Array<T,2> const& coord, plint& index, plint& count );
    void vsOrder();
    plint searchEdgeList (
        std::vector<EdgeListNode> const& edgeList, plint maxv ) const;
    BvmNodeIt bvmAdd(plint id);
    bool bvmCheck() const;
    void bvmLabel();
    plint& globalVertex(plint segment, plint localVertex);
    plint uniqueVertices(std::vector<Segment> const& segments);
    void computePointingVertex();
    plint createEdgeTable();
    void findBoundaryVertices();
    void computeNeighboringEdges();
    void computeInternalNeighboringEdges (
        plint iVertex, plint jVertex, plint segment1, plint segment2,
        plint localEdge1, plint va1, plint vb1 );
    void computeEmanatingEdges();
    bool fixOrientation();
    void fixOrientationOfNeighbors(plint iSegment, std::queue<plint>& segmentsToFixNeighbors,
                                   char* visitedSegments, bool& flag);
private:
    std::vector<Array<plint,2> > segmentIndices;
    VertexSet vertexSet;
    std::vector<std::vector<EdgeListNode> > edgeTable;
    BoundaryVertexMap boundaryVertexMap;
private:
    std::vector<Array<T,2> > vertexList;
    std::vector<plint> emanatingEdgeList;
    std::vector<Edge2D> edgeList;
    plint numSegments, numVertices;
template<typename U> friend 
void constructSurfaceMesh (
        std::vector<typename SegmentToDef<U>::Segment> const& segments,
        std::vector<Array<U,2> >& vertexList,
        std::vector<plint>& emanatingEdgeList,
        std::vector<Edge2D>& edgeList, U epsilon );
};

template<typename T>
void constructSurfaceMesh (
        std::vector<typename SegmentToDef<T>::Segment> const& segments,
        std::vector<Array<T,2> >& vertexList,
        std::vector<plint>& emanatingEdgeList,
        std::vector<Edge2D>& edgeList,
        T epsilon=std::numeric_limits<float>::epsilon() )
{
    SegmentToDef<T>(segments, epsilon).generateOnce (
            vertexList, emanatingEdgeList, edgeList );
}

} // namespace plb

#endif  // SEGMENT_TO_DEF_H
