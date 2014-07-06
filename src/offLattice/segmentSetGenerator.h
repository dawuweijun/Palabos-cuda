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

#ifndef SEGMENT_SET_GENERATOR_H
#define SEGMENT_SET_GENERATOR_H

#include "core/globalDefs.h"
#include "offLattice/segmentSet.h"
#include <vector>
#include <memory>

namespace plb {

/// Create and return a sphere as a set of segments. The center and radius of the sphere
///   must be given as arguments. A minimum number of segments for the sphere triangulation
///   must be provided as well. This number is suggestive for the resolution. The
///   actual number of segments can be greater than the one provided.
template<typename T>
SegmentSet<T>* constructCircle(Array<T,2> const& center, T radius, plint minNumOfSegments);


template<typename T>
SegmentSet<T> patchTubes(SegmentSet<T> const& geometryWithOpenings, plint sortDirection, std::vector<T> patchLengths);

/// Create and return a rectangle. The rectangle is on the x-y plane, and its lower left
///   corner is at the origin of the axes. It's sides have length "lx" and "ly", while
///   the number of points for the segment are "nx" and "ny" on the x and y axis,
///   respectively. This means that the total number of segments is 2*(nx-1)*(ny-1).
template<typename T>
SegmentSet<T> constructRectangle(T lx, T ly);

} // namespace plb

#endif  // SEGMENT_SET_GENERATOR_H

