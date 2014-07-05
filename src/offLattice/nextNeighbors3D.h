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

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef NEXT_NEIGHBORS_3D_H
#define NEXT_NEIGHBORS_3D_H

#include "core/globalDefs.h"

namespace plb {

template <typename T>
struct NextNeighbor3D {
    static const int numNeighbors=26;
    static const int c[numNeighbors][3];
    static const T d1;
    static const T d2;
    static const T d3;
    static const T d[numNeighbors];
    static const T id1;
    static const T id2;
    static const T id3;
    static const T invD[numNeighbors];
};

template<typename T, template<typename U> class Descriptor>
struct NextNeighborPop3D {
    NextNeighborPop3D();
    int ids[NextNeighbor3D<T>::numNeighbors];
};

template<typename T, template<typename U> class Descriptor>
inline plint nextNeighborPop3D(plint iNeighbor) {
    static NextNeighborPop3D<T,Descriptor> instance;
    return instance.ids[iNeighbor];
}

}  // namespace plb

#endif  // NEXT_NEIGHBORS_3D_H
