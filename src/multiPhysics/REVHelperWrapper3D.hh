#ifndef REV_HELPER_WRAPPER_3D_HH
#define REV_HELPER_WRAPPER_3D_HH
#include "multiPhysics/REVHelperWrapper3D.h"
namespace plb
{
/* *************** Local InvK ******************************************* */
template <typename T1,typename T2>
void computeLocalInvK3D ( MultiTensorField3D<T1,3>& orient,
                          MultiTensorField3D<T2,9>& invK ,
                          const typename BoxLocalInvKFunctional3D<T1,T2>::Matrix3D &K,
                          const typename BoxLocalInvKFunctional3D<T1,T2>::Vector3D &dir )
{
    applyProcessingFunctional (
        new BoxLocalInvKFunctional3D<T1,T2> ( K,dir ), orient.getBoundingBox(), orient, invK );
};
template <typename T1,typename T2>
std::auto_ptr<TensorField3D<T2,9> >computeLocalInvK3D (
    MultiTensorField3D<T1,3>& orient,
    const typename BoxLocalInvKFunctional3D<T1,T2>::Matrix3D &K,
    const typename BoxLocalInvKFunctional3D<T1,T2>::Vector3D &dir )
{
    TensorField3D<T2,9>* invK = new TensorField3D<T2,9> ( orient.getNx(), orient.getNy() );
    computeLocalInvK3D ( orient,invK,K,dir );
    return std::auto_ptr<TensorField3D<T2,9> > ( invK );
};

}//namespace plb
#endif //REV_HELPER_WRAPPER_3D_HH
