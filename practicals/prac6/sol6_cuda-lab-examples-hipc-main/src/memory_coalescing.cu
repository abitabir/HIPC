#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
    if (result != cudaSuccess) {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
#endif
    return result;
}

__global__ void stride(double* a, int s)
{
    int i = (blockDim.x * blockIdx.x + threadIdx.x) * s;
    a[i] = a[i] + 1;
}

void runTest(int deviceId, int nMB)
{
    int blockSize = 256;
    float ms;

    double *d_a;
    cudaEvent_t startEvent, stopEvent;
      
    int n = nMB*1024*1024/sizeof(double);

    // NB:  d_a(33*nMB) for stride case
    checkCuda( cudaMalloc(&d_a, n * 33 * sizeof(double)) );

    checkCuda( cudaEventCreate(&startEvent) );
    checkCuda( cudaEventCreate(&stopEvent) );

    printf("Stride, Bandwidth (GB/s):\n");

    stride<<<n/blockSize, blockSize>>>(d_a, 1); // warm up
    
    for (int i = 1; i <= 32; i++) {
        checkCuda( cudaMemset(d_a, 0, n * sizeof(double)) );

        checkCuda( cudaEventRecord(startEvent,0) );
        stride<<<n/blockSize, blockSize>>>(d_a, i);
        checkCuda( cudaEventRecord(stopEvent,0) );
        checkCuda( cudaEventSynchronize(stopEvent) );

        checkCuda( cudaEventElapsedTime(&ms, startEvent, stopEvent) );
        printf("%d\t%f\n", i, 2*nMB/ms);
    }

    checkCuda( cudaEventDestroy(startEvent) );
    checkCuda( cudaEventDestroy(stopEvent) );
    cudaFree(d_a);
}

int main(int argc, char **argv)
{
    int nMB = 4;
    int deviceId = 0;

    cudaDeviceProp prop;

    checkCuda( cudaSetDevice(deviceId) );
    checkCuda( cudaGetDeviceProperties(&prop, deviceId) );
    printf("Device: %s\n", prop.name);
    printf("Transfer size (MB): %d\n", nMB);

    runTest(deviceId, nMB);
}
