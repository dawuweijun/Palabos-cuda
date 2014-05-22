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

#ifndef REV_HELPER_FUNCTIONAL_2D_H
#define REV_HELPER_FUNCTIONAL_2D_H
#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
namespace plb
{
template <typename T1 ,typename T2>
class ComputeLocalInvK2D:public BoxProcessingFunctional2D_TT<T1,2,T2,4>
{
public:
    ComputeLocalInvK2D ( T1 K_[2][2],T1 dir_[2] ) :K ( K_ ), rawKDir ( dir_ )
    {
        PLB_PRECONDITION ( rawKDir[0]*rawKDir[0]+rawKDir[1]*rawKDir[1]==1. );
        //此处计算invK
        T1 det=K_[0][0]*K_[1][1]-K_[0][1]*K_[1][0];

        PLB_ASSERT ( abs ( det ) >0. );

        T1 inv_det=1./abs ( det );

        invK[0] = inv_det *K_[1][1];
        invK[1] = -inv_det *K_[0][1];
        invK[2] = -inv_det *K_[1][0];
        invK[3] = inv_det *K_[0][0];

    };
    void process ( Box2D domain, TensorField2D<T1,2>& orient,TensorField2D<T2,4>& invK );
    virtual ComputeLocalInvK2D* clone() const
    {
        return new ComputeLocalInvK2D ( K );
    };
private:
    T1 K[2][2],invK[2][2] ,rawKDir[2];
};
}//namespace plb
#endif

