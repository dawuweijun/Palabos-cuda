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

#ifndef PARTICLE_NON_LOCAL_TRANSFER_2D_H
#define PARTICLE_NON_LOCAL_TRANSFER_2D_H

#include "core/globalDefs.h"
#include "multiBlock/localMultiBlockInfo2D.h"
#include "particles/multiParticleField2D.h"


namespace plb {

template<class ParticleFieldT>
void copy (                                                                   
        MultiParticleField2D<ParticleFieldT> const& from, Box2D const& fromDomain,
        MultiParticleField2D<ParticleFieldT>& to, Box2D const& toDomain );

} // namespace plb

#endif  // PARTICLE_NON_LOCAL_TRANSFER_2D_H