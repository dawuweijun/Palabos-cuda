#ifndef REV_HELPER_WRAPPER_3D_HH
#define REV_HELPER_WRAPPER_3D_HH
#include "multiPhysics/REVHelperWrapper3D.h"
namespace plb
{
/* *************** Local InvK ******************************************* */
template <typename T1,typename T2>
void computeLocalInvK3D ( MultiTensorField3D<T1,3>& orient,
                          MultiTensorField3D<T2,9>& localInvK ,
                          const Array<T1,9> &rawInvK,
                          const Array<T1,3> &rawDir )
{
    applyProcessingFunctional (
        new BoxTensorRotationFunctional3D<T1,T2> ( rawInvK,rawDir ), orient.getBoundingBox(), orient, localInvK );
};
template <typename T1,typename T2>
std::auto_ptr<TensorField3D<T2,9> >computeLocalInvK3D (
    MultiTensorField3D<T1,3>& orient,
    const Array<T1,9> &rawInvK,
    const Array<T1,3> &rawDir )
{
    TensorField3D<T2,9>* localInvK = new TensorField3D<T2,9> ( orient.getNx(), orient.getNy() );
    computeLocalInvK3D ( orient,localInvK,rawInvK,rawDir );
    return std::auto_ptr<TensorField3D<T2,9> > ( localInvK );
};

}//namespace plb
#endif //REV_HELPER_WRAPPER_3D_HH
