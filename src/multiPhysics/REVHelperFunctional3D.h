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
class BoxLocalInvKFunctional3D:public BoxProcessingFunctional3D_TT<T1,3,T2,9>
{

public:
	typedef Array< T1, 3  > Vector3D;
	typedef Array< BoxLocalInvKFunctional3D::Vector3D,3 > Matrix3D;

    BoxLocalInvKFunctional3D ( const Matrix3D& K_,const Vector3D& dir_ ):K(K_),rawKDir(dir_)
    {
        PLB_PRECONDITION ( rawKDir[0]*rawKDir[0]+rawKDir[1]*rawKDir[1]+rawKDir[2]*rawKDir[2]==1. );
        //此处计算invK
        T1 det=K[0][0] * ( K[1][1] * K[2][2] - K[1][2] * K[2][1] ) -
               K[0][1] * ( K[1][0] * K[2][2] - K[1][2] * K[2][0] ) +
               K[0][2] * ( K[1][0] * K[2][1] - K[1][1] * K[2][0] );

        PLB_ASSERT ( abs ( det ) >0. );

        T1 inv_det =1./abs ( det );


        invK[0][0] = inv_det * ( K[1][1] * K[2][2] - K[1][2] * K[2][1] );
        invK[0][1] = inv_det * ( K[0][2] * K[2][1] - K[0][1] * K[2][2] );
        invK[0][2] = inv_det * ( K[0][1] * K[1][2] - K[0][2] * K[1][1] );
        invK[1][0] = inv_det * ( K[1][2] * K[2][0] - K[1][0] * K[2][2] );
        invK[1][1] = inv_det * ( K[0][0] * K[2][2] - K[0][2] * K[2][0] );
        invK[1][2] = inv_det * ( K[0][2] * K[1][0] - K[0][0] * K[1][2] );
        invK[2][0] = inv_det * ( K[1][0] * K[2][1] - K[1][1] * K[2][0] );
        invK[2][1] = inv_det * ( K[0][1] * K[2][0] - K[0][0] * K[2][1] );
        invK[2][2] = inv_det * ( K[0][0] * K[1][1] - K[0][1] * K[1][0] );

    };
    virtual void process ( Box3D domain, TensorField3D<T1,3>& orient,TensorField3D<T2,9>& invK );
    virtual BoxLocalInvKFunctional3D* clone() const
    {
        return new BoxLocalInvKFunctional3D ( K,rawKDir );
    };

    void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {

    };

private:
    Matrix3D K ;
    T1 invK[3][3];
    Vector3D rawKDir;
};
}//namespace plb
#endif

