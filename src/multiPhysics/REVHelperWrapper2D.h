#ifndef REV_HELPER_WRAPPER_2D_H
#define REV_HELPER_WRAPPER_2D_H
#include "multiPhysics/REVHelperFunctional2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "auto_ptr.h"
namespace plb
{
template <typename T1,typename T2>
void computeLocalInvK2D ( MultiTensorField2D< T1, 2  >& orient,
                          MultiTensorField2D< T2, 4  >& localInvK,
                          const typename BoxTensorRotationFunctional2D< T1, T2 >::Matrix2D& rawInvK,
                          const typename BoxTensorRotationFunctional2D< T1, T2 >::Vector2D& rawDir );

template <typename T1,typename T2>
std::auto_ptr<MultiTensorField2D<T2,4> > computeLocalInvK2D (
    MultiTensorField2D<T1,2>& orient,
    const typename BoxTensorRotationFunctional2D<T1,T2>::Matrix2D &rawInvK,
    const typename BoxTensorRotationFunctional2D<T1,T2>::Vector2D &rawDir );

}//namespace plb
#endif //REV_HELPER_WRAPPER_2D_H
