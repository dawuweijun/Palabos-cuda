#ifndef REV_HELPER_WRAPPER_2D_HH
#define REV_HELPER_WRAPPER_2D_HH
#include "multiPhysics/REVHelperWrapper2D.h"
namespace plb
{
/* *************** Local InvK ******************************************* */
template <typename T1,typename T2>
void computeLocalInvK2D ( MultiTensorField2D<T1,2>& orient,
                          MultiTensorField2D<T2,4>& invK ,
                          const typename BoxLocalInvKFunctional2D<T1,T2>::Matrix2D &K,
                          const typename BoxLocalInvKFunctional2D<T1,T2>::Vector2D &dir )
{
    applyProcessingFunctional (
        new BoxLocalInvKFunctional2D<T1,T2> ( K,dir ), orient.getBoundingBox(), orient, invK );
};
template <typename T1,typename T2>
std::auto_ptr<MultiTensorField2D<T2,4> >computeLocalInvK2D (
    MultiTensorField2D< T1, 2  >& orient,
    const typename BoxLocalInvKFunctional2D<T1,T2>::Matrix2D &K,
    const typename BoxLocalInvKFunctional2D<T1,T2>::Vector2D &dir )
{
    MultiTensorField2D<T2,4>* invK = new MultiTensorField2D<T2,4> ( orient.getNx(), orient.getNy() );
    computeLocalInvK2D ( orient,invK,K,dir );
    return std::auto_ptr<MultiTensorField2D<T2,4> > ( invK );
};

}//namespace plb
#endif //REV_HELPER_WRAPPER_2D_HH
