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

#ifndef REV_HELPER_FUNCTIONAL_3D_H
#define REV_HELPER_FUNCTIONAL_3D_H
#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
namespace plb
{
template <typename T1,typename T2=T1>
class BoxTensorRotationFunctional3D:public BoxProcessingFunctional3D_TT<T1,3,T2,9>
{

public:

    BoxTensorRotationFunctional3D ( const Array<T1,9>& rawTensor_,const Array<T1,3>& dir_ ) :rawTensor ( rawTensor_ ),rawKDir ( dir_ )
    {
        PLB_PRECONDITION ( rawKDir[0]*rawKDir[0]+rawKDir[1]*rawKDir[1]+rawKDir[2]*rawKDir[2]==1. );
    };
    virtual void process ( Box3D domain, TensorField3D<T1,3>& orient,TensorField3D<T2,9>& invKField );
    virtual BoxTensorRotationFunctional3D* clone() const
    {
        return new BoxTensorRotationFunctional3D ( rawTensor,rawKDir );
    };

    void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    };

private:
    Array<T1,9> rawTensor;
    Array<T1,3> rawKDir;
};
}//namespace plb
#endif

