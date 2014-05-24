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

#ifndef REV_HELPER_FUNCTIONAL_3D_HH
#define REV_HELPER_FUNCTIONAL_3D_HH
#include "multiPhysics/REVHelperFunctional3D.h"
#include "core/array.h"
namespace plb
{
template <typename T1,typename T2>
void BoxLocalInvKFunctional3D<T1,T2>::process ( Box3D domain, TensorField3D<T1,3>& orient,TensorField3D<T2,9>& invKField )
{
    T1 rmat[3][3];
    T2 Axis[3],normAxis[3];
    for ( int iX=domain.x0; iX<=domain.x1; iX++ )
    {
        for ( int iY=domain.y0; iY<=domain.y1; iY++ )
        {
            for ( int iZ=domain.z0; iZ<=domain.y1; iZ++ )
            {
                const Array< T1, 3> &to_B=orient.get ( iX,iY,iZ );
                double sqrNorm=to_B[0]*to_B[0]+to_B[1]*to_B[1]+to_B[2]*to_B[2];
                Array< T2, 9 > &tmpinvk=invKField.get ( iX,iY ,iZ );
                tmpinvk.resetToZero();

                PLB_PRECONDITION ( std::abs ( sqrNorm-0.5 ) <0.500001 );

                if ( sqrNorm!=0. )
                {
                    //计算出from rawDir to B 的单位旋转轴axis：
                    Axis[0]=rawKDir[1]*to_B[2]-rawKDir[2]*to_B[1];
                    Axis[1]=rawKDir[2]*to_B[0]-rawKDir[0]*to_B[2];
                    Axis[2]=rawKDir[0]*to_B[1]-rawKDir[1]*to_B[0];
                    //计算角度,sinTheta与符号无关，由Axis向量保证
                    T2 sinTheta = Axis[0]*Axis[0]+Axis[1]*Axis[1]+Axis[2]*Axis[2];
                    T2 cosTheta = rawKDir[0]*to_B[0]+rawKDir[1]*to_B[1]+rawKDir[2]*to_B[2];
                    T2 invNorm =1.0/sqrt ( sinTheta );
                    //归一化
                    normAxis[0]=invNorm*Axis[0];
                    normAxis[1]=invNorm*Axis[1];
                    normAxis[2]=invNorm*Axis[2];

                    //计算顺时针旋转矩阵
                    sinTheta=-sinTheta;
                    T2 negCos=1. -cosTheta;

                    rmat[0][0]=cosTheta+normAxis[0]*normAxis[0]* negCos;
                    rmat[0][1]=normAxis[0]*normAxis[1]* negCos - normAxis[2]*sinTheta;
                    rmat[0][2]=normAxis[0]*normAxis[2]* negCos + normAxis[1]*sinTheta;

                    rmat[1][0]=normAxis[1]*normAxis[0]* negCos + normAxis[2]*sinTheta;
                    rmat[1][1]=cosTheta + normAxis[1]*normAxis[1]* negCos;
                    rmat[1][2]=normAxis[1]*normAxis[2]* negCos - normAxis[0]*sinTheta;

                    rmat[2][0]=normAxis[2]*normAxis[0]* negCos - normAxis[1]*sinTheta;
                    rmat[2][1]=normAxis[2]*normAxis[1]* negCos + normAxis[0]*sinTheta;
                    rmat[2][2]=cosTheta + normAxis[2]*normAxis[2]* negCos;

                    //计算渗透率矩阵的逆矩阵
                    tmpinvk[0] = invK[0][0]*rmat[0][0]+invK[0][1]*rmat[1][0]+invK[0][2]*rmat[2][0];
                    tmpinvk[1] = invK[0][0]*rmat[0][1]+invK[0][1]*rmat[1][1]+invK[0][2]*rmat[2][1];
                    tmpinvk[2] = invK[0][0]*rmat[0][2]+invK[0][1]*rmat[1][2]+invK[0][2]*rmat[2][2];
                    tmpinvk[3] = invK[1][0]*rmat[0][0]+invK[1][1]*rmat[1][0]+invK[1][2]*rmat[2][0];
                    tmpinvk[4] = invK[1][0]*rmat[0][1]+invK[1][1]*rmat[1][1]+invK[1][2]*rmat[2][1];
                    tmpinvk[5] = invK[1][0]*rmat[0][2]+invK[1][1]*rmat[1][2]+invK[1][2]*rmat[2][2];
                    tmpinvk[6] = invK[2][0]*rmat[0][0]+invK[2][1]*rmat[1][0]+invK[2][2]*rmat[2][0];
                    tmpinvk[7] = invK[2][0]*rmat[0][1]+invK[2][1]*rmat[1][1]+invK[2][2]*rmat[2][1];
                    tmpinvk[8] = invK[2][0]*rmat[0][2]+invK[2][1]*rmat[1][2]+invK[2][2]*rmat[2][2];
                }
            }
        }
    }
}
}//namespace plb
#endif //REV_HELPER_FUNCTIONAL_3D_HH

