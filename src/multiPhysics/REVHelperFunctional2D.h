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
#include <cmath>
namespace plb
{
template <typename T1 ,typename T2>
class BoxLocalInvKFunctional2D:public BoxProcessingFunctional2D_TT<T1,2,T2,4>
{
public:

    typedef Array<T1,2> Vector2D;
    typedef Array<BoxLocalInvKFunctional2D::Vector2D,2> Matrix2D;

    BoxLocalInvKFunctional2D ( const Matrix2D &K_, const Vector2D &dir ) :K ( K_ ),rawKDir ( dir )
    {
        PLB_PRECONDITION ( std::abs ( rawKDir[0]*rawKDir[0]+rawKDir[1]*rawKDir[1]-1.0 ) <1.e-8 );
        std::cout<<rawKDir[0]<<rawKDir[1]<<std::endl;
        //此处计算invK
        T1 det=K_[0][0]*K_[1][1]-K_[0][1]*K_[1][0];
        PLB_PRECONDITION ( std::abs ( det ) >0. );

        T1 inv_det=1./std::abs ( det );
        invK[0][0] = inv_det *K_[1][1];
        invK[0][1] = -inv_det *K_[0][1];
        invK[1][0] = -inv_det *K_[1][0];
        invK[1][1] = inv_det *K_[0][0];
    };
    void process ( Box2D domain, TensorField2D<T1,2>& orient,TensorField2D<T2,4>& invKField );
    virtual BoxLocalInvKFunctional2D* clone() const
    {
        return new BoxLocalInvKFunctional2D<T1,T2> ( K, rawKDir );
    };
    void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    };
private:
    Matrix2D K;
    T1 invK[2][2];
    Vector2D rawKDir;
};
}//namespace plb
#endif


