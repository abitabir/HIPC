#include <stdio.h>
#include <cuda.h>

void queryDevice() {
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    printf("Number of GPUs: %d\n", deviceCount);

    int drvVer;
    int runVer;
    cudaDriverGetVersion(&drvVer);
    cudaRuntimeGetVersion(&runVer);
    cudaDeviceProp deviceProperties;
    printf("CUDA Driver Version / Runtime Version: %d.%d / %d.%d\n", drvVer/1000, (drvVer%100)/10, runVer/1000, (runVer%100)/10);

    int i;
    for (i = 0; i < deviceCount; i++) {
        cudaGetDeviceProperties(&deviceProperties, i);
        printf("Name: %s\n", deviceProperties.name);
    }
}

int main() {
    queryDevice();
}
