///////////////////////////////////////////////////////////////////////////////
//
// Palabos CUDA Portable C/C++ Utilities
//
// Copyright 2012 NVIDIA Corporation
//
// License: Apache License, v2.0 http://www.apache.org/licenses/LICENSE-2.0.html
//
// The home for Malabos is https://github.com/harrism/hemi
//
///////////////////////////////////////////////////////////////////////////////
// Please see the file README.md (https://github.com/harrism/hemi/README.md)
// for full documentation and discussion.
///////////////////////////////////////////////////////////////////////////////
#ifndef PLB_CUDA_H
#define PLB_CUDA_H

#include <stdio.h>
#include <assert.h>
#include "host_defines.h"
#include "cuda_runtime_api.h"
#include "driver_types.h"
namespace plb
{

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda ( cudaError_t result )
{
#if defined(DEBUG) || defined(_DEBUG)
    if ( result != cudaSuccess )
    {
        fprintf ( stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString ( result ) );
        assert ( result == cudaSuccess );
    }
#endif
    return result;
}

// Convenience function for checking CUDA error state including
// errors caused by asynchronous calls (like kernel launches). Note that
// this causes device synchronization, but is a no-op in release builds.
inline
cudaError_t checkCudaErrors()
{
    cudaError_t result = cudaSuccess;
    checkCuda ( result = cudaGetLastError() ); // runtime API errors
#if defined(DEBUG) || defined(_DEBUG)
    result = cudaDeviceSynchronize(); // async kernel launch errors
    if ( result != cudaSuccess )
        fprintf ( stderr, "CUDA Launch Error: %s\n", cudaGetErrorString ( result ) );
#endif
    return result;
}

#ifdef __CUDACC__ // CUDA compiler
	#define PLB_CUDA_COMPILER              // to detect CUDACC compilation
	#define PLB_LOC_STRING "Device"

	#ifdef __CUDA_ARCH__
		#define PLB_DEV_CODE                 // to detect device compilation
	#endif

	#define PLB_KERNEL(name)               __global__ void name ## _kernel
	#define PLB_KERNEL_NAME(name)          name ## _kernel

	#if defined(DEBUG) || defined(_DEBUG)
		#define PLB_KERNEL_LAUNCH(name, gridDim, blockDim, sharedBytes, streamId, ...) \
			do {                                                                     \
				name ## _kernel<<< (gridDim), (blockDim), (sharedBytes), (streamId) >>>\
				(__VA_ARGS__);                                                     \
				checkCudaErrors();                                                     \
			} while(0)
	#else
		#define PLB_KERNEL_LAUNCH(name, gridDim, blockDim, sharedBytes, streamId, ...) \
			name ## _kernel<<< (gridDim) , (blockDim), (sharedBytes), (streamId) >>>(__VA_ARGS__)
	#endif

	#define PLB_HOST_ONLY \
	#ifdef __CUDACC__ \
		#undef __CUDACC__ \
		#error "This function only can be called by host" \
	#endif

	#define PLB_DEV_CALLABLE               __host__ __device__
	#define PLB_DEV_CALLABLE_INLINE        __host__ __device__ inline
	#define PLB_DEV_CALLABLE_MEMBER        __host__ __device__
	#define PLB_DEV_CALLABLE_INLINE_MEMBER __host__ __device__ inline

// Constants: declares both a device and a host copy of this constant
// static and extern flavors can be used to declare static and extern
// linkage as required.
	#define PLB_DEFINE_CONSTANT(def, value) \
		__constant__ def ## _devconst = value; \
		def ## _hostconst = value
	#define PLB_DEFINE_STATIC_CONSTANT(def, value) \
		static __constant__ def ## _devconst = value; \
		static def ## _hostconst = value
	#define PLB_DEFINE_EXTERN_CONSTANT(def) \
		extern __constant__ def ## _devconst; \
		extern def ## _hostconst

// use to access device constant explicitly
	#define PLB_DEV_CONSTANT(name) name ## _devconst

// use to access a constant defined with PLB_DEFINE_*_CONSTANT
// automatically chooses either host or device depending on compilation
	#ifdef PLB_DEV_CODE
		#define PLB_CONSTANT(name) name ## _devconst
	#else
		#define PLB_CONSTANT(name) name ## _hostconst
	#endif

	#if !defined(PLB_ALIGN)
		#define PLB_ALIGN(n) __align__(n)
	#endif

#else             // host compiler
	#define PLB_HOST_COMPILER              // to detect non-CUDACC compilation
	#define PLB_LOC_STRING "Host"

	#define PLB_KERNEL(name)               void name
	#define PLB_KERNEL_NAME(name)          name
	#define PLB_KERNEL_LAUNCH(name, gridDim, blockDim, sharedBytes, streamId, ...) name(__VA_ARGS__)


	#define PLB_HOST_ONLY

	#define PLB_DEV_CALLABLE
	#define PLB_DEV_CALLABLE_INLINE        inline
	#define PLB_DEV_CALLABLE_MEMBER
	#define PLB_DEV_CALLABLE_INLINE_MEMBER inline

	#define PLB_DEFINE_CONSTANT(def, value) def ## _hostconst = value
	#define PLB_DEFINE_STATIC_CONSTANT(def, value) static def ## _hostconst = value
	#define PLB_DEFINE_EXTERN_CONSTANT(def) extern def ## _hostconst

	#undef PLB_DEV_CONSTANT // requires NVCC, so undefined here!
	#define PLB_CONSTANT(name) name ## _hostconst

	#if !defined(PLB_ALIGN)
		#if defined(__GNUC__)
			#define PLB_ALIGN(n) __attribute__((aligned(n)))
		#elif defined(_MSC_VER)
			#define PLB_ALIGN(n) __declspec(align(n))
		#else
			#error "Please provide a definition of PLB_ALIGN for your host compiler!"
		#endif
	#endif

#endif

/*
 * TODO:如何保证BlockLattice中Cell的偏移量在Host与Device上保持一致?
 * ATTENTION:
 * 	1,BlockLattice相对全局有一个偏移量，Cell相对单个BlockLattice有一个偏移量;
 * 	2,单个BlockLattice上的数据可能在不同的cuda block中执行;
 */

// Returns the offset of the current thread's element within the grid for device
// code compiled with NVCC, or 0 for sequential host code.
PLB_DEV_CALLABLE_INLINE
int hemiGetElementOffset()
{
#ifdef PLB_DEV_CODE
    return blockIdx.x * blockDim.x + threadIdx.x;
#else
    return 0;
#endif
}

// Returns the offset of the current thread's X element within the grid for device
// code compiled with NVCC, or 0 for sequential host code.
PLB_DEV_CALLABLE_INLINE
int hemiGetElementXOffset()
{
    return hemiGetElementOffset();
}

// Returns the offset of the current thread's Y element within the grid for device
// code compiled with NVCC, or 0 for sequential host code.
PLB_DEV_CALLABLE_INLINE
int hemiGetElementYOffset()
{
#ifdef PLB_DEV_CODE
    return blockIdx.y * blockDim.y + threadIdx.y;
#else
    return 0;
#endif
}

// Returns the stride of the current grid (blockDim.x * gridDim.x) for device
// code compiled with NVCC, or 1 for sequential host code.
PLB_DEV_CALLABLE_INLINE
int hemiGetElementStride()
{
#ifdef PLB_DEV_CODE
    return blockDim.x * gridDim.x;
#else
    return 1;
#endif
}

// Returns the stride of the current grid (blockDim.x * gridDim.x) for device
// code compiled with NVCC, or 1 for sequential host code.
PLB_DEV_CALLABLE_INLINE
int hemiGetElementXStride()
{
    return hemiGetElementStride();
}

// Returns the stride of the current grid (blockDim.x * gridDim.x) for device
// code compiled with NVCC, or 1 for sequential host code.
PLB_DEV_CALLABLE_INLINE
int hemiGetElementYStride()
{
#ifdef PLB_DEV_CODE
    return blockDim.y * gridDim.y;
#else
    return 1;
#endif
}

}//END NAMESPACE PLB

#endif // PLB_CUDA_H


