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
 * #include "core/globalDefs.h"
 * Groups all the 2D include files in the boundaryConditions directory.
 */

#include "boundaryCondition/boundaryDynamics.h"
#include "boundaryCondition/regularizedBoundaryDynamics.h"
#include "boundaryCondition/regularizedBoundaryDynamics2D.h"
#include "boundaryCondition/equilibriumBoundaryDynamics.h"
#include "boundaryCondition/inamuroAnalyticalDynamics.h"
#include "boundaryCondition/inamuroBoundary2D.h"
#include "boundaryCondition/zouHeBoundary2D.h"
#include "boundaryCondition/zouHeDynamics.h"
#include "boundaryCondition/boundaryCondition2D.h"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor2D.h"
#include "boundaryCondition/wrappedLocalBoundaryProcessor2D.h"
#include "boundaryCondition/neumannCondition2D.h"
#include "boundaryCondition/bounceBackModels.h"
#include "boundaryCondition/bounceBackModels2D.h"
#include "boundaryCondition/spongeZones2D.h"

#ifndef PLB_BGP
#include "boundaryCondition/generalizedBoundaryDynamics.h"
#include "boundaryCondition/generalizedBoundaryCondition2D.h"
#endif

