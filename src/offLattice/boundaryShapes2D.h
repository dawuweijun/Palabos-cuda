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

#ifndef BOUNDARY_SHAPES_2D_H
#define BOUNDARY_SHAPES_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/atomicBlock2D.h"
#include "core/geometry2D.h"
#include "core/array.h"

namespace plb {
/*
struct OffBoundary {
    enum Type {dirichlet, neumann, freeSlip, constRhoInlet, densityNeumann, flux, isolation};
};
*/
/// Description of a shape as the boundary of a given volume.
template<typename T, class SurfaceData>
struct BoundaryShape2D {
    virtual ~BoundaryShape2D() { }
    /// Decide whether a given discrete lattice node is inside the solid shape.
    /** The reason why there is an isInside function instead of a isOutside
     *  one is that while there exists only one voxelFlag::inside flag,
     *  there are several flags which count as outside: undetermined,
     *  outside, and borderline. This is particularly important for the
     *  undetermined flag, because the outer envelopes in a multi-block
     *  structure have no special meaning and are default-initialized to
     *  undetermined.
     **/
    virtual bool isInside(Dot2D const& location) const =0;
    /// Get the distance to the wall, and the velocity value on the wall,
    ///   from a given real-valued position (in lattice units), and along
    ///   a given direction. Returns true if there exists an intersection
    ///   along the indicated direction.
    ///   ATTENTION: The id is an in-and-out value. If you know the right
    ///   id of the surface intersection (e.g. the triangle ID in case of
    ///   a triangular mesh) you can provide it to speed up the calculation.
    ///   However, if you don't know it, you MUST provide the value -1, 
    ///   because otherwise the result might be wrong.
    virtual bool pointOnSurface (
            Array<T,2> const& fromPoint, Array<T,2> const& direction,
            Array<T,2>& locatedPoint, T& distance,
            Array<T,2>& wallNormal, SurfaceData& surfaceData,
            OffBoundary::Type& bdType, plint& id ) const =0;
    /// Get the distance to the wall, and the data on the wall,
    ///   from a given discrete lattice node, and along a given direction.
    ///   Returns true if there exists an intersection along the indicated
    ///   direction.
    ///   ATTENTION: The id is an in-and-out value. If you know the right
    ///   id of the surface intersection (e.g. the triangle ID in case of
    ///   a triangular mesh) you can provide it to speed up the calculation.
    ///   However, if you don't know it, you MUST provide the value -1, 
    ///   because otherwise the result might be wrong.
    virtual bool gridPointOnSurface (
            Dot2D const& fromPoint, Dot2D const& direction,
            Array<T,2>& locatedPoint, T& distance,
            Array<T,2>& wallNormal, SurfaceData& surfaceData,
            OffBoundary::Type& bdType, plint& id ) const
    {
        return pointOnSurface(Array<T,2>((T)fromPoint.x, (T)fromPoint.y),
                              Array<T,2>((T)direction.x, (T)direction.y),
                              locatedPoint, distance, wallNormal, surfaceData, bdType, id);
    }
    /// Given a point p on the surface of the shape, determine its "continuous normal".
    ///   If the shape is for example piecewise linear, the normal is adjusted to vary
    ///   continuously over the surface.
    virtual Array<T,2> computeContinuousNormal (
            Array<T,2> const& p, plint id, bool isAreaWeighted = false ) const =0;
    /// Say if a given segment intersects the surface.
    virtual bool intersectsSurface (
            Array<T,2> const& p1, Array<T,2> const& p2, plint& id ) const =0;
    /// Say if a given segment with integer-valued endpoints intersects the surface.
    virtual bool gridIntersectsSurface (
            Dot2D const& p1, Dot2D const& p2, plint& id ) const
    {
        return intersectsSurface( Array<T,2>((T)p1.x, (T)p1.y),
                                  Array<T,2>((T)p2.x, (T)p2.y), id );
    }
    /// Get the shortest distance to the wall. Returns true in case of success.
    ///   The flag isBehind indicates if the point is "behind" the wall, i.e.
    ///   in the direction opposite to the wall normal.
    virtual bool distanceToSurface( Array<T,2> const& point,
                                    T& distance, bool& isBehind ) const =0;
    /// Get the shortest distance to the wall. Returns true in case of success.
    ///   The flag isBehind indicates if the point is "behind" the wall, i.e.
    ///   in the direction opposite to the wall normal.
    virtual bool gridDistanceToSurface( Dot2D const& point,
                                        T& distance, bool& isBehind ) const
    {
        return distanceToSurface( Array<T,2>((T)point.x, (T)point.y),
                                  distance, isBehind );
    }
    /// Return the tag (id of boundary-portion with specific boundary,
    ///   condition); the id is the one returned by pointOnSurface.
    virtual plint getTag(plint id) const =0;
    /// Plain clone function.
    virtual BoundaryShape2D<T,SurfaceData>* clone() const =0;
    /// In case the shape class needs additional meshed data, this clone function
    ///   offers the possibility to provide the data.
    virtual BoundaryShape2D<T,SurfaceData>* clone(std::vector<AtomicBlock2D*> args) const {
        return clone();
    }
};

template<typename T, class SurfaceData>
class BoundaryShapeIntersection2D : public BoundaryShape2D<T,SurfaceData> {
public:
    BoundaryShapeIntersection2D();
    virtual ~BoundaryShapeIntersection2D();
    BoundaryShapeIntersection2D(BoundaryShapeIntersection2D<T,SurfaceData> const& rhs);
    BoundaryShapeIntersection2D<T,SurfaceData>& operator=(BoundaryShapeIntersection2D<T,SurfaceData> const& rhs);
    void swap(BoundaryShapeIntersection2D<T,SurfaceData>& rhs);
    void addComponent(BoundaryShape2D<T,SurfaceData>* component);

    virtual bool isInside(Array<T,2> const& location) const;
    virtual bool isInside(Array<T,2> const& location, T epsilon) const;
    virtual bool pointOnSurface (
            Array<T,2> const& fromPoint, Array<T,2> const& direction,
            Array<T,2>& locatedPoint, T& distance,
            Array<T,2>& wallNormal, SurfaceData& surfaceData,
            OffBoundary::Type& bdType, plint& id ) const;
    virtual bool intersectsSurface (
            Array<T,2> const& p1, Array<T,2> const& p2, plint& id ) const;
    virtual bool distanceToSurface( Array<T,2> const& point,
                                    T& distance, bool& isBehind ) const;
    virtual plint getTag(plint id) const;
    virtual BoundaryShapeIntersection2D<T,SurfaceData>* clone() const;
private:
    std::vector<BoundaryShape2D<T,SurfaceData>*> components;
};

}  // namespace plb

#endif  // BOUNDARY_SHAPES_2D_H

