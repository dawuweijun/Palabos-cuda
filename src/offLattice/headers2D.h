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
 * Groups all the include files for the 2D off-lattice directory.
 */

#include "offLattice/marchingCube.h"
#include "offLattice/segmentHash.h"
#include "offLattice/nextNeighbors2D.h"
#include "offLattice/segmentToDef.h"
#include "offLattice/segmentPolygonMesh2D.h"
#include "offLattice/voxelizer2D.h"
#include "offLattice/makeSparse2D.h"
#include "offLattice/segmentHash.h"
#include "offLattice/offLatticeBoundaryProcessor2D.h"
#include "offLattice/offLatticeBoundaryProfiles2D.h"
#include "offLattice/offLatticeBoundaryCondition2D.h"
#include "offLattice/boundaryShapes2D.h"
#include "offLattice/segmentBoundary2D.h"
#include "offLattice/offLatticeModel2D.h"
#include "offLattice/guoOffLatticeModel2D.h"
#include "offLattice/bouzidiOffLatticeModel2D.h"
#include "offLattice/guoAdvDiffOffLatticeModel2D.h"
#include "offLattice/segmentSetGenerator.h"
#include "offLattice/immersedWalls2D.h"
#include "offLattice/filippovaHaenel2D.h"

#ifndef PLB_BGP
#include "offLattice/generalizedOffLatticeModel2D.h"
#endif

