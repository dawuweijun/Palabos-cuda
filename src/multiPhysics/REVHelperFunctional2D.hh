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

#ifndef REV_HELPER_FUNCTIONAL_2D_HH
#define REV_HELPER_FUNCTIONAL_2D_HH
#include "multiPhysics/REVHelperFunctional2D.h"
#include "core/array.h"
namespace plb
{
template <typename T1,typename T2>
void ComputeLocalInvK2D<T1,T2>::process ( Box2D domain, TensorField2D<T1,2>& orient,TensorField2D<T2,4>& invK )
{
    T1 rmat[2][2];

    for ( int iX=domain.x0; iX<=domain.x1; iX++ )
    {
        for ( int iY=domain.y0; iY<=domain.y1; iY++ )
        {
            const Array< T1, 2> &to_B=orient.get ( iX,iY );
            double sqrNorm=to_B[0]*to_B[0]+to_B[1]*to_B[1];
            Array< T2, 4 > &tmpinvk=invK.get ( iX,iY );
            tmpinvk.resetToZero();
            PLB_PRECONDITION ( sqrNorm==0.||sqrNorm==1. );
            if ( sqrNorm!=0. )
            {
                //计算角度
                T1 sinTheta = rawKDir[0]*to_B[1]-rawKDir[1]*to_B[0];
                T1 cosTheta = rawKDir[0]*to_B[0]+rawKDir[1]*to_B[1];
                //计算顺时针旋转矩阵
                rmat[0][0]=cosTheta;
                rmat[0][1]=sinTheta;
                rmat[1][0]=-sinTheta;
                rmat[1][1]=cosTheta;
                //计算渗透率矩阵的逆矩阵
                tmpinvk[0]=invK[0][0]*rmat[0][0]+invK[0][1]*rmat[1][0];
                tmpinvk[1]=invK[0][0]*rmat[0][1]+invK[0][1]*rmat[1][1];
                tmpinvk[2]=invK[1][0]*rmat[0][0]+invK[1][1]*rmat[1][0];
                tmpinvk[3]=invK[1][0]*rmat[0][1]+invK[1][1]*rmat[1][1];
            }
        }
    }
}

}//namespace plb
#endif



