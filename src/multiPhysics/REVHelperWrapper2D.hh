#ifndef REV_HELPER_WRAPPER_2D_HH
#define REV_HELPER_WRAPPER_2D_HH
#include "multiPhysics/REVHelperWrapper2D.h"
namespace plb
{
/* *************** Local InvK ******************************************* */
template <typename T1,typename T2>
void computeLocalInvK2D ( MultiTensorField2D<T1,2>& orient,
                          MultiTensorField2D<T2,4>& localInvK ,
                          const Array<T1,4> &rawInvK,
                          const Array<T1,2> &rawDir )
{
    applyProcessingFunctional (
        new BoxTensorRotationFunctional2D<T1,T2> ( rawInvK,rawDir ), orient.getBoundingBox(), orient, localInvK );
};
template <typename T1,typename T2>
std::auto_ptr<MultiTensorField2D<T2,4> >computeLocalInvK2D (
    MultiTensorField2D< T1, 2  >& orient,
    const Array<T1,4> &rawInvK,
    const Array<T1,2> &rawDir )
{
    MultiTensorField2D<T2,4>* localInvK = new MultiTensorField2D<T2,4> ( orient.getNx(), orient.getNy() );
    computeLocalInvK2D ( orient,*localInvK,rawInvK,rawDir );
    return std::auto_ptr<MultiTensorField2D<T2,4> > ( localInvK );
};

}//namespace plb
#endif //REV_HELPER_WRAPPER_2D_HH
