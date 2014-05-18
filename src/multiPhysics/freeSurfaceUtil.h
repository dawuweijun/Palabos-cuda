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

#ifndef FREE_SURFACE_UTIL_H
#define FREE_SURFACE_UTIL_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiBlockLattice2D.h"

#include <vector>
#include <set>
#include <string>

namespace plb {

	/// Constants used in a free surface flag matrix for cell tagging.
	namespace twoPhaseFlag {
		enum Flag {empty=0, interface=1, fluid=2, wall=4, protect=5, protectEmpty=6};
		inline std::string flagToString(int flag) {
			switch(flag) {
				case empty:     return "empty";
				case interface: return "interface";
				case fluid:     return "fluid";
				case wall:      return "wall";
				case protect:   return "protect";
				case protectEmpty:   return "protectEmpty";
				default: PLB_ASSERT( false );
			}
		}
		inline Flag invert(int flag) {
			switch(flag) {
				case empty:     return fluid;
				case interface: return interface;
				case fluid:     return empty;
				case wall:      return wall;
				case protect:   return protect;
				case protectEmpty: return protectEmpty;
				default: PLB_ASSERT( false );
			}
		}
		inline bool isWet(int flag) {
			return flag==interface || flag==fluid || flag==protect;
		}
		inline bool isFullWet(int flag) {
			return flag==fluid || flag==protect;
		}
		inline bool isEmpty(int flag) {
			return flag==empty || flag==protectEmpty;
		}
	}

	/// Data structure for holding lists of cells along the free surface in an AtomicContainerBlock.
	template< typename T,template<typename U> class Descriptor>
	struct InterfaceLists : public ContainerBlockData {
		typedef Array<plint,Descriptor<T>::d> Node;
		/// Holds all nodes which have excess mass.
		std::map<Node,T> massExcess;
		/// Holds all nodes that need to change status from interface to fluid.
		std::set<Node>   interfaceToFluid;
		/// Holds all nodes that need to change status from interface to empty.
		std::set<Node>   interfaceToEmpty;
		/// Holds all nodes that need to change status from empty to interface.
		std::set<Node>   emptyToInterface;

		virtual InterfaceLists<T,Descriptor>* clone() const {
			return new InterfaceLists<T,Descriptor>(*this);
		}
	};

}  // namespace plb

#endif  // FREE_SURFACE_UTIL_H

