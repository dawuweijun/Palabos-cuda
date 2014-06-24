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

#ifndef SEGMENT_SET_H
#define SEGMENT_SET_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "core/array.h"
#include <string>
#include <vector>
#include <cstdio>

namespace plb {

template<typename T>
class SegmentSet {
public:
    typedef Array<Array<T, 2>,2> Segment;
public:
    SegmentSet(Precision precision_ = FLT);
    SegmentSet(std::vector<Segment> const& segments_, Precision precision_ = FLT);
    SegmentSet(std::string fname, Precision precision_ = FLT, SurfaceGeometryFileFormat fformat = STL);
    std::vector<Segment> const& getSegments() const;
    Precision getPrecision() const { return precision; }
    void setPrecision(Precision precision_);

    /// Translate the Segment set polygon mesh.
    void translate(Array<T,2> const& vector);
    /// Scale the segment set polygon mesh.
    void scale(T alpha);
    /// Rotate the segment set polygon mesh.
    ///   The arguments of this function are the Euler angles in radians.
    ///   The so-called "x-convention" is used, in which the rotation is
    ///   given by the three angles (phi, theta, psi), where:
    ///   1.] The first rotation is by an angle phi about the z-axis,
    ///   2.] The second rotation is by an angle theta in [0, pi] about
    ///       the new x-axis,
    ///   3.] The third rotation is by an angle psi about the new z-axis.
    void rotate(T theta);
    
    /// Erase the current segment set surface mesh, and merge into it the new meshes.
    ///   This function currently does not check for duplicate
    ///   segments in the new meshes, and does not handle
    ///   cases when the resulting merged surface mesh is
    ///   non-manifold. It is in the user's jurisdiction to
    ///   make sure that the resulting mesh is well defined.
    ///   Practically, this can be achieved if the new segment set
    ///   meshes given as arguments are mutually disjoint.
    void merge(std::vector<SegmentSet<T>*> meshes);

    /// Append to the current segment mesh, the mesh that is passed as an argument.
    ///   This function currently does not check for duplicate
    ///   segments in the two meshes, and does not handle
    ///   cases when the resulting merged surface mesh is
    ///   non-manifold. It is in the user's jurisdiction to
    ///   make sure that the resulting mesh is well defined.
    ///   Practically, this can be achieved if the two segment set
    ///   meshes are mutually disjoint.
    void append(SegmentSet<T> const& mesh);

    /// Refine the current segment set surface mesh by splitting each
    ///   original segment into four smaller segments constructed by
    ///   joining the middle points of the original edges. The old mesh
    ///   is deleted.
    void refine();

    /// A very simple orientation reversing function.
    void reverseOrientation();

    /// Export the mesh as an ASCII STL file.
    void writeAsciiSTL(std::string fname) const;
    /// Export the mesh as an binary STL file.
    void writeBinarySTL(std::string fname) const;

    /// Cut the current segment set mesh by a plane "plane" which is
    ///   defined by a point and a normal. This cutting operation will
    ///   produce a new segment set "newSegmentSet" which is
    ///   the one that belongs to the half-space that the normal of the
    ///   cutting plane points outside of. The function returns 1 if
    ///   the cut was successful, 0 if there was no intersection between
    ///   the original segment set and the plane provided, and -1 if an
    ///   error occured.
    int cutWithLine(Line2D<T> const& plane, SegmentSet<T>& newSegmentSet) const;

    /// Cut the current segment set mesh by a plane "plane" which is
    ///   defined by a point and a normal and a cuboid "cuboid" which is
    ///   defined by a lower left corner and a right upper corner. This
    ///   cutting operation will cut off all segments, or parts of segments
    ///   that are fully contained in the cuboid and are positioned in
    ///   the half-space that the normal of the cutting plane points into.
    ///   Caution is needed, because this function does not use the cuboid
    ///   to cut, but only to select parts of the original segment set,
    ///   to be then cut by the plane. Obviously, at least a part of the
    ///   cutting plane must be contained in the cuboid, for the cutting
    ///   to make any sense. If not used wisely, this function can lead to
    ///   broken STL files. The function returns 1 if the cut was successful,
    ///   0 if there was no intersection between the original segment set
    ///   the plane and the cuboid provided, and -1 if an error occured.
    int cutWithLine(Line2D<T> const& plane, Cuboid2D<T> const& cuboid,
            SegmentSet<T>& newSegmentSet) const;

    T getMinEdgeLength() const { return minSegLength; }
    T getMaxEdgeLength() const { return maxSegLength; }

    Cuboid2D<T> getBoundingCuboid() const { return boundingCuboid; }

private:
    void readSTL(std::string fname);
    void readAsciiSTL(FILE* fp);
    void readBinarySTL(FILE* fp);
    void check(Segment& segment, Array<T,2> const& n);
    bool checkNoAbort(Segment& segment, Array<T,2> const& n);
    void computeSegmentLength(pluint iSegment, T& segLength) const;
    void computeMinMaxSegments();
    void computeBoundingCuboid();
    /// Cut the current segment by a plane "plane" which is
    ///   defined by a point and a normal. This cutting operation will
    ///   add or not one or more segments to the segment set "newSegmentSet".
    ///   The segments added belong to the half-space that the normal of
    ///   the plane given points outside of. The function returns 1 if the
    ///   cut was successful, 0 if there was no intersection between the
    ///   original segment set and the plane and -1 if an error occured.
    int cutSegmentWithLine(Line2D<T> const& plane, Segment const& segment,
            SegmentSet<T>& newSegmentSet) const;

private:
    std::vector<Segment> segments;
    T minSegLength, maxSegLength;
    Cuboid2D<T> boundingCuboid;
    /// Enumeration constant that sets the precision for segment set checks.
    /// Single precision (float): FLT
    /// Double precision (double): DBL
    /// Extended precision (long double): LDBL
    Precision precision;
};

} // namespace plb

#endif  // SEGMENT_SET_H
