#ifndef CUDAMANAGED_H
#define CUDAMANAGED_H

#include<cuda.h>
#include <cuda_runtime_api.h>
extern int cudaMallocManaged ( void * ptr,size_t len );
// extern int cudaFree ( void * ptr );
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
        cudaMallocManaged ( &ptr, len );
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
