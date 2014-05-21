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
 * Utilities for 2D multi data distributions -- header file.
 */

#ifndef REDISTRIBUTION_2D_H
#define REDISTRIBUTION_2D_H

#include "parallelism/mpiManager.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockManagement2D.h"

namespace plb {

struct MultiBlockRedistribute2D {
    virtual ~MultiBlockRedistribute2D() { }
    virtual MultiBlockManagement2D redistribute(MultiBlockManagement2D const& original) const=0;
};

class RandomRedistribute2D : public MultiBlockRedistribute2D {
public:
    RandomRedistribute2D(pluint rseed_=10);
    virtual MultiBlockManagement2D redistribute(MultiBlockManagement2D const& original) const;
private:
    pluint rseed;
};

}  // namespace plb

#endif  // REDISTRIBUTION_2D_H
