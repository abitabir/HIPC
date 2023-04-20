#include <stdio.h>
#include <cuda.h>

__global__ void helloWorld() {
    printf("I am block (%d), thread (%d)\n", blockIdx.x, threadIdx.x);
}

int main() {
    int grid_size = 3;
    int block_size = 3;

    helloWorld<<<grid_size, block_size>>>();

    cudaDeviceSynchronize();

    return 0;
}
