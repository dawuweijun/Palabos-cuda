#ifndef CUDAMANAGED_H
#define CUDAMANAGED_H

#include"plbDebug.h"
#include<cuda.h>
#include <cuda_runtime_api.h>

namespace plb
{
class  cudaManaged
{
    /*
     * 重载 new 操作，自动分配 cudaManaged Memory
     */
public:
    void *operator new ( size_t len )
    {
        void *ptr;
        PLB_ASSERT(cudaMallocManaged ( &ptr, len,cudaMemAttachGlobal)!=cudaErrorNotSupported);
        return ptr;
    }
    /*
     * 重载 delete 操作，自动删除 cudaManaged Memory;
     */
    void operator delete ( void *ptr )
    {
        cudaFree ( ptr );
    }
};
}
#endif//CUDAMANAGED_H
