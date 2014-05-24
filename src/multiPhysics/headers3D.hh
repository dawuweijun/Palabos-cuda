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
 * Groups all the generic 3D implementation files in the directory multiPhysics.
 */

#include "multiPhysics/boussinesqThermalProcessor3D.hh"
#include "multiPhysics/advectionDiffusion3D.hh"
#include "multiPhysics/interparticlePotential.hh"
#include "multiPhysics/shanChenProcessor3D.hh"
#include "multiPhysics/thermalDataAnalysis3D.hh"
#include "multiPhysics/heLeeProcessor3D.hh"
#include "multiPhysics/freeSurfaceModel3D.hh"
#include "multiPhysics/freeSurfaceBoundaryCondition3D.hh"
#include "multiPhysics/freeSurfaceInitializer3D.hh"
#include "multiPhysics/freeSurfaceAnalysis3D.hh"
#include "multiPhysics/multiFreeSurfaceModel3D.hh"
#include "multiPhysics/createBubbles3D.hh"
#include "multiPhysics/bubbleHistory3D.hh"
#include "multiPhysics/bubbleMatch3D.hh"
#include "multiPhysics/twoPhaseModel3D.hh"
#include "multiPhysics/REVProcessor3D.hh"
#include "multiPhysics/REVHelperFunctional3D.hh"
#include "multiPhysics/REVHelperWrapper3D.hh"
