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

#ifndef SEGMENT_POLYGON_MESH_H
#define SEGMENT_POLYGON_MESH_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "core/array.h"
#include <string>
#include <vector>

namespace plb {

/// Holds one of the three edges of a triangle.
struct Edge2D {
    plint pv;  // Pointing Vertex.
    plint ne;  // Neighboring Edge.
};

/// Refers to a lid, i.e. a collection of segments, constructed
///   by filling a hole in the mesh.
struct Curve2D{
    // Comment: if the close-holes algorithm is to be made more
    //   general one day, you will need to (1) assign the proper
    //   value to numAddedVertices, and (2) replace centerVertex
    //   by something else and change all codes which access
    //   centerVertex.
    Curve2D() : numAddedVertices(1) { }
    plint firstSegment;
    plint numSegments;
    std::vector<plint> boundaryVertices;
    plint centerVertex;
    plint numAddedVertices;
};

/// Represent a polygon mesh, made up of adjacent segments,
///   in a directed-edge format.
/** This format conveniently provides access to the surface topology
 * through segments and vertices, but is less convenient for accessing
 * it directly through edges.
 *
 * The segments are numbered continuously from 0 to getNumSegments()-1,
 * and the vertices from 0 to getNumVertices()-1.
 *
 * Once the mesh is constructed, you can change the position of its
 * vertices, but you should not change its topology, which is static.
 * One exception: the method closeHoles changes the topology.
 *
 * For more information on the directed-edges format, check the paper
 * "Directed edges - A scalable representation for triangle meshes",
 * Journal of Graphics Tools (3), 1998, pp 1 - 11  
 **/
template<typename T>
class SegmentPolygonMesh2D {
public:
    /// Typedefs for backward compatibility.
    typedef plb::Curve2D Curve2D;
    typedef plb::Edge2D Edge2D;
public:
    /// The ownership over the parameters vertexList, emanatingEdgeList, and
    ///   edgeList is not taken by the class SegmentPolygonMesh2D, and they
    ///   are not being copied. They must be managed from outside and kept
    ///   alive up to the usage time of the SegmentPolygonMesh2D instance.
    ///   The user may choose to disregard vertices at the end of vertexList,
    ///   and to only use a selected amount at the beginning: do this through
    ///   the parameter numVertices. A negative value for numVertices means:
    ///   use all vertices in vertexList.
    SegmentPolygonMesh2D (
            std::vector<Array<T,2> >& vertexList_,
            std::vector<plint>& emanatingEdgeList_,
            std::vector<Edge2D>& edgeList_,
            plint numVertices_ = -1 );

    /// Return number of segments on the surface (they are numbered 
    ///   continuously from 0 to getNumSegments() -1.
    plint getNumSegments() const {
        return numSegments;
    }
    /// Return number of vertices on the surface (they are numbered 
    ///   continuously from 0 to getNumVertices() -1.
    plint getNumVertices() const {
        return numVertices;
    }
    /// Reset the position of a vertex to a new value;
    void replaceVertex(plint iVertex, Array<T,2> const& newPosition);

    /// Reset all vertices of a mesh to a default value.
    /** This function is used when debugging programs in which the access to the
     *  mesh is parallelized. Default values can be used to point out undefined
     *  state. In debug mode, vertices are identified to be invalid if all three
     *  coordinates are distinctly less than -1.
     */
    void resetVertices(Array<T,2> const& defaultVertex);

    /// Get one of the three vertices of a given triangle.
    ///   "localVertex" must be 0, 1 or 2.
    Array<T,2> const& getVertex(plint iSegment, int localVertex) const;
    /// Get the coordinates of one of the vertices by providing its global index.
    Array<T,2> const& getVertex(plint iVertex) const;
    bool isValidVertex(plint iSegment, int localVertex) const;
    bool isValidVertex(plint iVertex) const;

    /// Compute the minimum and maximum vertex positions in every direction.
    void computeBoundingBox (
            Array<T,2>& xRange,
            Array<T,2>& yRange ) const;

    /// Translate the surface mesh.
    void translate(Array<T,2> const& vector);
    /// Scale the surface mesh.
    void scale(T alpha);
    /// Rotate the surface mesh.
    ///   The arguments of this function are the Euler angles in radians.
    ///   The so-called "x-convention" is used, in which the rotation is
    ///   given by the three angles (phi, theta, psi), where:
    ///   1.] The first rotation is by an angle phi about the z-axis,
    ///   2.] The second rotation is by an angle theta in [0, pi] about
    ///       the new x-axis,
    ///   3.] The third rotation is by an angle psi about the new z-axis.
    void rotate(T theta);

    /// Smooth the surface mesh.
    ///   The triangular surface mesh is smoothed by using a spatial
    ///   averaging algorithm. Interior vertices are treated differently
    ///   than boundary ones. The mesh is smoothed as many times as
    ///   the integer argument "maxiter" indicates. If "isMeasureWeighted"
    ///   is true, then triangle areas and edge lengths are used in the
    ///   spatial averaging procedure. The smoothing method uses a relaxation
    ///   algorithm with a relaxation parameter 0 <= relax <= 1.
    void smooth(plint maxiter = 1, T relax = 1.0, bool isMeasureWeighted = false);

    /// Get the global vertex index of a specific vertex local to a triangle.
    plint getVertexId(plint iSegment, plint localVertex) const;
    /// Get a list of the neighboring vertices of a given vertex.
    std::vector<plint> getNeighborVertexIds(plint iVertex) const;
    /// Get a list of the neighboring vertices of a given edge.
    std::vector<plint> getNeighborVertexIds(plint iVertex, plint jVertex) const;
    /// Get a list of the neighboring segments of a given vertex.
    std::vector<plint> getNeighborSegmentIds(plint iVertex) const;
    /// Get a list of the adjacent segments of a given triangle.
    ///   The adjacent segments to a specific triangle are defined
    ///   as those that have common edges with the given triangle.
    std::vector<plint> getAdjacentSegmentIds(plint iSegment) const;
    /// Get a list of the adjacent segments of a given edge.
    ///   The adjacent segments to a specific edge are defined
    ///   as those that have the specific edge common.
    std::vector<plint> getAdjacentSegmentIds(plint iVertex, plint jVertex) const;

    /// Compute the normal vector for a given triangle. If "isAreaWeighted" is false,
    ///   then the normal has length equal to one. If "isAreaWeighted" is true, then
    ///   the normal has length equal to twice the area of the triangle.
    Array<T,2> computeSegmentNormal(plint iSegment, bool isAreaWeighted = false) const;
    Array<T,2> computeSegmentNormal(
            plint iVertex, plint jVertex, plint kVertex, bool isAreaWeighted = false) const;

    /// Compute the normal vector at a given edge. If "isAreaWeighted" is true, then
    ///   the triangle normals used in the computation are area weighted.
    Array<T,2> computeEdgeNormal(
            plint iVertex, plint jVertex, bool isAreaWeighted = false) const;
    /// Compute the normal vector at a given vertex. If "isAreaWeighted" is false, then
    ///   the vertex normal is calculated as the simple normalized sum of the individual
    ///   unit normals of the segments which share the specific vertex. If
    ///   "isAreaWeighted" is true, then the vertex normal is calculated as the
    ///   normalized sum of the individual area weighted normals of the segments which
    ///   share the vertex under consideration.
    Array<T,2> computeVertexNormal(plint iVertex, bool isAreaWeighted = false) const;

    /// Compute a "continuous" normal at a point "p" which belongs to the triangle
    ///   with index "iSegment". The computed normal vector is normalized (has length
    ///   equal to one) and is continuous as its base point moves across the segments.
    ///   The point "p" must belong to the interior or to the boundary of the triangle
    ///   with index "iSegment". The currently implemented algorithm uses a barycentric
    ///   coordinate representation of the normal with respect to the three vertex
    ///   normals of the given triangle. If "isAreaWeighted" is true, then the vertex
    ///   normals used in the barycentric representation are area weighted.
    Array<T,2> computeContinuousNormal(
            Array<T,2> const& p, plint iSegment, bool isAreaWeighted = false) const;

    /// Compute the area of a given triangle.
    T computeSegmentArea(plint iSegment) const;
    T computeSegmentArea(plint iVertex, plint jVertex, plint kVertex) const;
    /// Compute the one third of the sum of the areas of the segments that share
    ///   the specific edge.
    T computeEdgeArea(plint iVertex, plint jVertex) const;
    /// Compute the one third of the sum of the areas of all segments that
    ///   share the given vertex.
    T computeVertexArea(plint iVertex) const;

    /// Compute the length of a given edge.
    T computeEdgeLength(plint iVertex, plint jVertex) const;

    /// Compute the dihedral angle of an edge (in radians).
    ///   By convention, if the edge is a boundary edge, the
    ///   dihedral angle returned by the function is 0.
    T computeDihedralAngle(plint iVertex, plint jVertex) const;

    /// Compute the one sixth of the sum of the heights of the
    ///   two segments incident to the specific edge.
    T computeEdgeTileSpan(plint iVertex, plint jVertex) const;

    /// Export the surface mesh as an ASCII STL file.
    void writeAsciiSTL(std::string fname) const;
    /// Export the surface mesh as an binary STL file.
    void writeBinarySTL(std::string fname) const;

    /// Return true if the vertex belongs to the boundary
    ///   or false if the vertex belongs to the interior
    ///   of the mesh.
    bool isBoundaryVertex(plint iVertex) const;
    /// Return true if the vertex belongs to the interior
    ///   or false if the vertex belongs to the boundary
    ///   of the mesh.
    bool isInteriorVertex(plint iVertex) const;

    /// Return true if the edge belongs to the boundary
    ///   or false if the edge belongs to the interior
    ///   of the mesh.
    bool isBoundaryEdge(plint iVertex, plint jVertex) const;
    /// Return true if the edge belongs to the interior
    ///   or false if the edge belongs to the boundary
    ///   of the mesh.
    bool isInteriorEdge(plint iVertex, plint jVertex) const;

    /// Function to compute the intersection between a triangle and a line segment
    ///   between points "point1" and "point2" (if flag = 0), or between a triangle
    ///   and a half-line starting at "point1" and containing "point2" (if flag = 1),
    ///   or between a triangle and a whole line containing both "point1" and
    ///   "point2" (if flag = 2). "intersection", "normal" and "distance" are
    ///   objects whose states are changed by this function. These states are undefined,
    ///   and cannot be used by the caller function, when the return value of this
    ///   function is not 1.
    ///   Returns 1 if an intersection is found and the intersection is inside
    ///   the triangle, 0 if there is no intersection or the intersection is
    ///   outside the triangle, and -1 is the line is parallel to the triangle
    ///   and intersects with an infinity of points.
    int pointOnSegment (
            Array<T,2> const& point1, Array<T,2> const& point2, int flag, plint iSegment,
            Array<T,2>& intersection, Array<T,2>& normal, T& distance ) const;

    /// The following function is a specialized version of the previous one. It simply checks
    ///   if a line segment between points "point1" and "point2" intersects a triangle. This
    ///   function is written for optimization purposes, and returns "true" if an intersection
    ///   is found and the intersection is inside the triangle, or "false" if there is no
    ///   intersection or the intersection is outside the triangle or if the line segment
    ///   belongs to the triangle and intersects with an infinity of points.
    bool segmentIntersectsSegment (
            Array<T,2> const& point1, Array<T,2> const& point2, plint iSegment ) const;

    /// Compute the distance between a point and a line which goes through a
    ///   given edge. Returns also a flag indicating whether the location on
    ///   the line which is closest to the point is inside the edge or not.
    void distanceToEdgeLine (
            Array<T,2> const& point, plint iSegment, plint whichEdge,
            T& distance, bool& intersectionIsInside ) const;

    /// Compute the distance between a point and a plane which goes through a
    ///   given triangle. Returns two flags. The first indicates whether the
    ///   location on the plane which is closest to the point is inside the
    ///   triangle or not. The second indicates whether the point is "behind"
    ///   the plane, i.e. on the side which is opposed to the triangle normal.
    void distanceToSegmentPlane (
            Array<T,2> const& point, plint iSegment,
            T& distance, bool& intersectionIsInside, bool& pointIsBehind ) const;

    /// Compute the distance between the point and a triangle. The location
    ///   on the triangle the distance to which is being computed is chosen
    ///   as follows. The point is first projected on the plane defined by
    ///   the triangle. If the obtained location is inside the triangle it is
    ///   being selected. Otherwise, it is projected on the three lines defined
    ///   by the edges, and the nearest one is picked out. If it is inside the
    ///   edge it is being selected. Otherwise, the answer consists of the
    ///   shortest distance between the point and the three vertices. Note
    ///   that none of the three cases (plane, edge, vertex) is degenerate:
    ///   there exists an infinity of points in space for which the nearest
    ///   distance to the triangle is on a vertex.
    void distanceToSegment (
            Array<T,2> const& point, plint iSegment,
            T& distance, bool& pointIsBehind ) const;

    /// Function to reverse the orientation of the given surface mesh.
    void reverseOrientation();

    /// Detect holes in the mesh and fill each of them by creating a lid,
    ///   i.e. a collection of segments created by adding a single new
    ///   vertex at the barycenter of the hole, and connecting it radially
    ///   with each boundary vertex of the hole. Returns a Lid2Dstructure
    ///   for each created lid.
    std::vector<Curve2D> closeHoles();
    /// Returns, for each detected hole, a list of vertices.
    std::vector<std::vector<plint> > detectHoles();
    void avoidIntegerPositions();
    void avoidIntegerPosition(plint iVertex);
    void inflate(T amount=1.e-3);
public:
    /// Get a handle to the vertices.
    std::vector<Array<T,2> >const& vertices() const {
        PLB_PRECONDITION(vertexList);
        return *vertexList;
    }
    /// Get a handle to the emanating edges.
    std::vector<plint> const& emanatingEdges() const {
        PLB_PRECONDITION(emanatingEdgeList);
        return *emanatingEdgeList;
    }
    /// Get a handle to the edges.
    std::vector<Edge2D> const& edges() const {
        PLB_PRECONDITION(edgeList);
        return *edgeList;
    }
    /// Get one of the three vertices of a given triangle (non-const).
    ///   "localVertex" must be 0, 1 or 2.
    Array<T,2>& getVertex(plint iSegment, int localVertex);
    /// Get the coordinates of one of the vertices by providing its global index (non-const).
    Array<T,2>& getVertex(plint iVertex);
public:
    /// Export the surface mesh as an HTML file.
    void writeHTML(std::string fname);
    void writeHTML(std::string fname, std::string title, T phys_dx, Array<T,2> phys_location);
private:
    /// Ensure the fact that a vertex is in a valid range by an assertion,
    ///   if STL vertices are bounded by IOpolicy().
    void assertVertex(Array<T,2> const& vertex) const;
    /// Get a non-const handle to the vertices.
    std::vector<Array<T,2> >& vertices() {
        PLB_PRECONDITION(vertexList);
        return *vertexList;
    }
    /// Get a non-const handle to the emanating edges.
    std::vector<plint>& emanatingEdges() {
        PLB_PRECONDITION(emanatingEdgeList);
        return *emanatingEdgeList;
    }
    /// Get a non-const handle to the edges.
    std::vector<Edge2D>& edges() {
        PLB_PRECONDITION(edgeList);
        return *edgeList;
    }
private:
    /// Close a single hole by creating a lid.
    Curve2D closeHole(std::vector<plint> const& hole);
private:
    plint prev(plint iEdge) const;
    plint next(plint iEdge) const;
    plint changeEdgeId(plint iEdge) const;
private:
    std::vector<Array<T,2> >* vertexList;
    std::vector<plint>* emanatingEdgeList;
    std::vector<Edge2D>* edgeList;
    plint numSegments, numVertices;
public:
    static const T eps0, eps1;
};

template<typename T>
class LidLessThan2D {
public:
    LidLessThan2D(plint mainDirection_, SegmentPolygonMesh2D<T> const& mesh_)
        : mainDirection(mainDirection_),
          mesh(mesh_)
    { }
    bool operator()(Curve2D const& lid1, Curve2D const& lid2) const
    {
        plint dim1 = mainDirection;
        plint dim2 = (dim1+1)%3;
        plint dim3 = (dim2+1)%3;
        // Compare on bary-center for main coordinate, and, in the unlikely
        //   case that these two floats are equal, on the two subsequent
        //   coordinates.
        Array<T,2> baryCenter1 = computeBaryCenter(mesh, lid1);
        Array<T,2> baryCenter2 = computeBaryCenter(mesh, lid2);
        return
          (baryCenter1[dim1] < baryCenter2[dim1]) || (
                  equiv(baryCenter1[dim1],baryCenter2[dim1]) && (
                      (baryCenter1[dim2] < baryCenter2[dim2]) || (
                          equiv(baryCenter1[dim2],baryCenter2[dim2]) &&
                          (baryCenter1[dim3] < baryCenter2[dim3]) ) ) ) ;
    }
    static bool equiv(T a, T b) {
        return !(a<b || b<a);
    }
private:
    plint mainDirection;
    SegmentPolygonMesh2D<T> const& mesh;
};

/// Returns the scaling factor.
template<typename T>
T toLatticeUnits (
        SegmentPolygonMesh2D<T>& mesh,
        plint resolution, plint referenceDirection );
template<typename T>
void toLatticeUnits2D (
        SegmentPolygonMesh2D<T>& mesh,
        plint resolution, plint referenceDirection,
        Array<T,2>& location, T& dx);


/* ******* Lid2Doperations ************************************************** */

template<typename T>
Array<T,2> computeBaryCenter2D (
        SegmentPolygonMesh2D<T> const& mesh, Curve2D const& lid );

template<typename T>
void computeBoundingBox (
        SegmentPolygonMesh2D<T> const& mesh, Curve2D const& lid,
        Array<T,2>& xLim, Array<T,2>& yLim);

template<typename T>
T computeInnerRadius (
        SegmentPolygonMesh2D<T> const& mesh, Curve2D const& lid );

template<typename T>
T computeOuterRadius (
        SegmentPolygonMesh2D<T> const& mesh, Curve2D const& lid );
//TODO in 2D condition, we should compute the length of the curve
/*
template<typename T>
T computeArea (
        SegmentPolygonMesh2D<T> const& mesh, Curve2D const& lid );
*/
template<typename T>
Array<T,2> computeNormal (
        SegmentPolygonMesh2D<T> const& mesh, Curve2D const& lid );

template<typename T>
void reCenter (
        SegmentPolygonMesh2D<T>& mesh, Curve2D const& lid );

} // namespace plb

#endif  // SEGMENT_POLYGON_MESH_H
