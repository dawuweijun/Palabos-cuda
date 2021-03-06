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

#ifndef REV_PROCESSORS_3D_H
#define REV_PROCESSORS_3D_H

#include "core/globalDefs.h"
#include "core/block3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
namespace plb{
template<typename T1, template<typename U> class Descriptor, typename T2=T1>
class BrinkmanProcessor3D: public BoxProcessingFunctional3D_LT<T1,Descriptor,T2,9>
{
public:
    virtual void process ( Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                           TensorField3D<T2,9>& negNiuInvsK );
    virtual BrinkmanProcessor3D<T1,Descriptor,T2> *clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::staticVariables;
        modified[1] = modif::staticVariables;
    }
};

template<typename T, template<typename U> class Descriptor>
class BrinkmanProcessor3DL: public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    BrinkmanProcessor3DL (const Array<T,9> &negNiuInvsK_ ):negNiuInvsK(negNiuInvsK_){};
    virtual void process ( Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual BrinkmanProcessor3DL<T,Descriptor> *clone() const{
      return new BrinkmanProcessor3DL<T,Descriptor>(negNiuInvsK);
    };
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::staticVariables;
    }
private:
    Array<T,9> negNiuInvsK;
};
}// namespace plb
#endif
